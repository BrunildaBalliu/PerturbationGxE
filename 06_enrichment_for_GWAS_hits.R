#%%%%%%%%%%%%%%% Libraries and functions
x=c("cowplot","data.table", "dplyr", "ggplot2",  "ggrepel", "magrittr" , "reshape2", "corrplot","stringr","gridBase","grid","qvalue")
garbage=suppressMessages(lapply(x, require, character.only = TRUE))

source('../Functions.R')
source('scripts/00_MoMeIR_Functions.R')

#%%%%%%%%%%%%%%% Perturbation and cell order and colors
if(1){
  perturb_order=c("GLUC", "INSU", "IGF1", "SB20", "SP60", "U012", "WORT", "IL-6", "TGFB1", "TNFa", "DEXA", "ADIP", "LEPT", "ATOR", "ROSI", "METF", "DECA", "LAUR", "RETA", "IBMX", "ISOP")
  cell_order=c("HMCL-7304", "SGBS", "HEPG2")
  cell_labels=c("Muscle", "Fat", "Liver")
  cell_col=c("#E69F00","#999999","#56B4E9")
  color_scheme = c("white", "#FCD9DD","#F04257", "#9DE1FB", "#26BCF7", "#D690DF", "#BD52CB")
}


#%%%%%%%%%%%%%%% DE and matched non-DE results
if(1){
 MatchedRes=NULL
  
  for(Cell in cell_labels){
    for (treat in perturb_order){
      tmp2 = read.csv(file = paste0('results/CoexpressionNetworks/GeneMatching/',Cell,"_",treat,".csv"))  
      MatchedRes=rbind(MatchedRes,tmp2)
    }
  }
  
  DE_results=MatchedRes %>% select(cell, Treatment, geneID, sig5FDR_hFDR3)
}

#%%%%%%%%%%%%%%% GWAS catalog results
# Clean GWAS catalogue
if(0){
gwas_catalog=fread(input = 'results/GWAS_catalog/gwas_catalog_v1.0.2-associations_e100_r2020-06-17.tsv', header = T, sep = '\t',stringsAsFactors = F) 
gwas_catalog_trait_annotation=readxl::read_xlsx(path = 'results/GWAS_catalog/gwas_catalog_trait-mappings_r2020-06-17.xlsx',col_names = T) 
gwas_catalog_trait_annotation %<>% filter(!duplicated(`Disease trait`))
gwas_catalog_annotated = merge(x = gwas_catalog, y = gwas_catalog_trait_annotation %>% select(`Disease trait` , `Parent term` ),by.x = "DISEASE/TRAIT",by.y = "Disease trait")

# Filter out parent terms not to be included in analysis
gwas_catalog_annotated %<>% 
  select(`Parent term`,MAPPED_TRAIT,MAPPED_GENE, CONTEXT, INTERGENIC) %>% 
  rename(PARENT_TRAIT=`Parent term`, TRAIT=MAPPED_TRAIT,GENE=MAPPED_GENE)  %>%
  filter(!PARENT_TRAIT %in% c("Other measurement","Other disease","Other trait")) 

# Clean up mapped gene names
# MAPPED GENE(S)*: Gene(s) mapped to the strongest SNP. If the SNP is located within a gene, that gene is listed. If the SNP is located within multiple genes, these genes are listed separated by commas. If the SNP is intergenic, the upstream and downstream genes are listed, separated by a hyphen.
entries_with_commas=grep(pattern = ",", x = gwas_catalog_annotated$GENE)
gwas_catalog_entries_with_commas=NULL
for(i in entries_with_commas){
  genes_i=unique(unlist(strsplit(x = gwas_catalog_annotated$GENE[i],split = ', ')))
  tmp=do.call("rbind", replicate(length(genes_i), data.frame(gwas_catalog_annotated[i,]), simplify = FALSE))
  tmp[,"GENE"]=genes_i
  gwas_catalog_entries_with_commas=rbind(gwas_catalog_entries_with_commas,tmp)
}

entries_with_hyphen=grep(pattern = "-", x = gwas_catalog_annotated$GENE)
gwas_catalog_entries_with_hyphen=NULL
for(i in entries_with_hyphen){
  genes_i=unique(unlist(strsplit(x = gwas_catalog_annotated$GENE[i],split = ' - ')))
  tmp=do.call("rbind", replicate(length(genes_i), data.frame(gwas_catalog_annotated[i,]), simplify = FALSE))
  tmp[,"GENE"]=genes_i
  gwas_catalog_entries_with_hyphen=rbind(gwas_catalog_entries_with_hyphen,tmp)
}

gwas_catalog_clean=rbind(gwas_catalog_annotated[-c(entries_with_commas,entries_with_hyphen),], gwas_catalog_entries_with_commas, gwas_catalog_entries_with_hyphen)
gwas_catalog_clean %<>% filter(!duplicated(paste(PARENT_TRAIT,TRAIT,GENE)))

write.table(x = gwas_catalog_clean, file = 'results/GWAS_catalog/gwas_catalog_v1.0.2-associations_e100_r2020-06-17_clean.tsv',sep = '\t', quote = F,row.names = F,col.names = T) 

}

gwas_catalog_clean=fread(input = 'results/GWAS_catalog/gwas_catalog_v1.0.2-associations_e100_r2020-06-17_clean.tsv', header = T, sep = '\t',stringsAsFactors = F) 
mapped_trait=gwas_catalog_clean %>% group_by(TRAIT) %>% summarise(n=length(unique(GENE))) %>% filter(n>=100)  %>% arrange(desc(n))
parent_term=gwas_catalog_clean %>% group_by(PARENT_TRAIT) %>% summarise(n=length(unique(GENE))) %>% filter(n>=100)  %>% arrange(desc(n))

#%%%%%%%%%%%%%%% Test for enrichment by trait
if(0){
  if(0){
    GWAS_enrich_each_trait=NULL
    for(trait in mapped_trait$TRAIT){
      gwas_res_trait_i=gwas_catalog_clean %>% filter(TRAIT == trait)
      for(cells in cell_labels){
        for(perturb in perturb_order){
          print(c(trait,cells,perturb))
          tmp=droplevels(DE_results %>% filter(cell==cells & Treatment==perturb))  %>% 
            mutate(GWASgene=1*(geneID %in% unique(gwas_res_trait_i$GENE))) 
          
          tab=table(DE=tmp$sig5FDR_hFDR3,GWAS=tmp$GWASgene)
          if(nrow(tab) ==2 & ncol(tab) ==2) GWAS_enrich_each_trait=rbind(GWAS_enrich_each_trait,c(gwas_res_trait_i$PARENT_TRAIT[1],trait,cells,perturb,unlist(fisher.test(x = tab)[c("estimate","conf.int","p.value")]),c(tab)))
        }}}
    colnames(GWAS_enrich_each_trait) = c("Parent_Trait","Trait","Cell", "Treatment", "OR", "conf.int1","conf.int2","p.value","n_noGWASnoDE","n_noGWASDE","n_GWASnoDE","n_GWASDE")
    write.table(x = GWAS_enrich_each_trait,file = 'results/GWAS_catalog/GWAS_enrich_each_trait_FisherExactTest.tsv',append = F,quote = F,sep = '\t',row.names = F, col.names = T)
  }
  
  
  GWAS_enrich_each_trait=fread(input = 'results/GWAS_catalog/GWAS_enrich_each_trait_FisherExactTest.tsv',sep = '\t', header = T, stringsAsFactors = F)
  
  GWAS_enrich_each_trait %<>%
    mutate(cell = factor(x = Cell,levels = cell_labels), 
           perturbation = factor(x = Treatment,levels = perturb_order),
           bh_p=p.adjust(p.value,method = "BH"),
           # Trait=factor(x = Trait,levels = trait_levels, labels = trait_labels),
           sign=factor(x = cut(x = bh_p,breaks = c(0,.05,.10,1),include.lowest = T), 
                       levels = c("[0,0.05]", "(0.05,0.1]", "(0.1,1]"))) 
  
  GWAS_enrich_each_trait$Trait[ GWAS_enrich_each_trait$Trait=="granulocyte percentage of myeloid white cells"] = "granulocyte % of myeloid white cells"
  GWAS_enrich_each_trait$Trait[ GWAS_enrich_each_trait$Trait=="monocyte percentage of leukocytes"] = "monocyte % of leukocytes"
  GWAS_enrich_each_trait$Trait[ GWAS_enrich_each_trait$Trait=="birth weight, parental genotype effect measurement"] = "birth weight, parental genotype effect"

  traits_2_keep=as.character(unique((GWAS_enrich_each_trait %>% filter(bh_p <=.1))$Trait))
  
  p_GWAS_trait=ggplot(data = GWAS_enrich_each_trait   %>%  
           mutate(OR=ifelse(test = bh_p<=.1,yes = OR,no = NA))%>%
           filter(Trait %in% traits_2_keep) ,
         mapping = aes(x = perturbation, y = OR, col=cell)) + 
    geom_point(size=4,position=position_dodge(.2)) + 
    facet_wrap(~Parent_Trait+Trait, ncol = 4,scales = "free_x") + 
    theme_bw()  + 
    theme(axis.title.y = element_text(size = 20, hjust = .5),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12, colour = "black", angle=90, vjust=.5, hjust = 1),
          axis.text.y = element_text(size = 15, colour = "black"),
          plot.margin = unit(c(5.5, 5.5, 5.5, 20), "points"),
          legend.position = "top",
          legend.background = element_blank(),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          strip.text = element_text(size=12))  + 
    scale_color_manual(values = cell_col) +
    geom_hline(yintercept = 1) + 
    ylab("OR for GWAS enrichment")  

  if(saveFig) ggsave(filename = 'manuscript/Figures/Fig04_GWASenrich_bytrait.pdf', plot = p_GWAS_trait, width = 15, height = 15)
}

#%%%%%%%%%%%%%%% Test for enrichment by trait type
if(1){
  GWAS_enrich_each_trait_type=data.frame(matrix(data = NA,ncol =  11, nrow = 1), stringsAsFactors = F)
  colnames(GWAS_enrich_each_trait_type) = c("Parent_Trait","Cell", "Treatment", "OR", "conf.int1","conf.int2","p.value","n_noGWASnoDE","n_noGWASDE","n_GWASnoDE","n_GWASDE")
  
  for(trait in parent_term$PARENT_TRAIT){
    gwas_res_trait_i=gwas_catalog_clean %>% filter(PARENT_TRAIT == trait)
    for(cells in cell_labels){
      for(perturb in perturb_order){
        print(c(trait,cells,perturb))
        tmp=droplevels(DE_results %>% filter(cell==cells & Treatment==perturb))  %>% 
          mutate(GWASgene=1*(geneID %in% unique(gwas_res_trait_i$GENE))) 
        
        tab=table(DE=tmp$sig5FDR_hFDR3,GWAS=tmp$GWASgene)
        if(nrow(tab) ==2 & ncol(tab) ==2) GWAS_enrich_each_trait_type=rbind(GWAS_enrich_each_trait_type,c(trait,cells,perturb,unlist(fisher.test(x = tab)[c("estimate","conf.int","p.value")]),c(tab)))
      }}}
  
  GWAS_enrich_each_trait_type=GWAS_enrich_each_trait_type[-1,]
  GWAS_enrich_each_trait_type %<>% 
    mutate(OR=as.numeric(OR), conf.int1=as.numeric(conf.int1),        
           conf.int2=as.numeric(conf.int2), p.value=as.numeric(p.value),            
           BHPvalue=p.adjust(p = as.numeric(p.value),method = "BH"))
  
  write.table(x = GWAS_enrich_each_trait_type,file = 'results/GWAS_catalog/GWAS_enrich_each_trait_type_FisherExactTest.tsv',append = F,quote = F,sep = '\t',row.names = F, col.names = T)
  # WriteXLS::WriteXLS(x = GWAS_enrich_each_trait_type, ExcelFileName = 'manuscript/Tables/TableS09_GWAS_enrich_each_trait_type.xls',row.names = F,col.names = T,BoldHeaderRow = T,FreezeRow = 1)
  
}

#%%%%%%%%%%%%%%% Test for enrichment across all traits
if(1){
  GWAS_enrich_all_traits=data.frame(matrix(data = NA,nrow = length(cell_labels) * length(perturb_order), ncol = 10,
                                           dimnames = list(NULL,c("cell", "Treatment", "OR", "conf.int1","conf.int2","p.value","n_noGWASnoDE","n_noGWASDE","n_GWASnoDE","n_GWASDE"))))
  GWAS_enrich_all_traits$cell=rep(cell_labels,each=length(perturb_order))
  GWAS_enrich_all_traits$Treatment=rep(perturb_order,times=length(cell_labels))
  
  for(cells in cell_labels){
    for(perturb in perturb_order){
      print(c(cells,perturb))
      tmp=droplevels(DE_results %>% filter(cell==cells & Treatment==perturb))  %>% 
        mutate(GWASgene=1*(geneID %in% unique(gwas_catalog_annotated$`REPORTED GENE(S)`))) 
      
      tab=table(DE=tmp$sig5FDR_hFDR3,GWAS=tmp$GWASgene)
      
      if(nrow(tab) ==2 & ncol(tab) ==2) GWAS_enrich_all_traits[GWAS_enrich_all_traits$cell==cells & GWAS_enrich_all_traits$Treatment==perturb,-c(1:2)] = c(unlist(fisher.test(x = tab)[c("estimate","conf.int","p.value")]),c(tab))
      
    }
  }
  
  GWAS_enrich_all_traits %<>% na.omit() %>% mutate(Treatment = factor(Treatment,levels = perturb_order),
                                                   cell=factor(x = cell,levels = cell_labels),
                                                   Sig=factor(x = p.adjust(p = p.value,method = "BH")<=.05,
                                                              levels = c(FALSE,TRUE), 
                                                              labels = c(">.05", "<=.05")))
  
  traits2keep=as.character((GWAS_enrich_all_traits %>% group_by(Treatment) %>% summarise(Sig=any(Sig=="<=.05")) %>% filter(Sig))$Treatment)
  
  d=.5
  p_GWAS = ggplot(data = GWAS_enrich_all_traits %>% filter(Treatment %in% traits2keep ),
         mapping = aes(x = Treatment, y = OR, col=cell)) + 
    geom_point(mapping = aes(size=Sig),position=position_dodge(d)) + 
    # geom_errorbar(aes(ymin= conf.int1, ymax=conf.int2), width=.1, position=position_dodge(d)) +
    theme_bw()  +
    theme(axis.title.y = element_text(size = 20, hjust = .5),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 15, colour = "black", angle=90, vjust=.5, hjust = 1),
          axis.text.y = element_text(size = 15, colour = "black"),
          plot.margin = unit(c(5.5, 5.5, 5.5, 20), "points"),
          legend.position = "top",
          legend.background = element_blank(),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20))  + 
    scale_color_manual(values = cell_col) +
    geom_hline(yintercept = 1) + 
    ylab("OR for GWAS enrichment") + 
    guides(color=FALSE, size=guide_legend(title = "Adjusted P-value"))
}





