library(Seurat)
library(SeuratDisk)
library(magrittr)
library(tidyverse)
library(org.Hs.eg.db)
setwd("~/DATA/BRAIN/AD")




seu <- LoadH5Seurat("./dlpfc.h5seurat")

meta <- read.csv("./dlpfc_meta.csv",row.names = 1)
seu@meta.data <- meta[rownames(seu@meta.data ),]






seu$Lewy.body.disease.pathology %<>% fct_collapse(reference = c("Reference","Not Identified (olfactory bulb assessed)","Not Identified (olfactory bulb not assessed)"))

seu$Lewy.body.disease.pathology %<>% fct_relevel(c( "reference", "Neocortical (Diffuse)","Olfactory bulb only",
                                                   "Limbic (Transitional)", "Brainstem-predominant", "Amygdala-predominant"
))


library(rtracklayer)
gene <- rtracklayer::import("./gencode.v32.primary_assembly.annotation.gtf") %>% as.data.frame()
gene %<>% select(gene_name,gene_id,seqnames) %>% distinct() %>% mutate(gene_id = str_remove(gene_id,"\\..*"))
gene %<>% filter(str_detect(seqnames,"chr\\d+"))

geneID_to_sym <- gene$gene_name %>% setNames(gene$gene_id)

geneSymToID <-  gene$gene_id %>% setNames(gene$gene_name)
qs::qsave(seu, file = "dfpfc.qs",nthreads = 20)
seu_org <- qs::qread("./dfpfc.qs")

count <- GetAssayData(seu_org,slot = "count")
rownames(count) <- geneID_to_sym[rownames(count)]
count <- count[rownames(count) %in% gene$gene_name,]

seu <- CreateSeuratObject(counts = count,meta.data = seu_org@meta.data)

seu@reductions <- seu_org@reductions

qs::qsave(seu, file = "dfpfc.qs")

seu <- qs::qread("./dfpfc.qs")

seu %<>% NormalizeData() %>% ScaleData(features = c("NOSTRIN", "FLT1", "CLDN5", "SLC1A3", "MUSTN1", "AQP4", "FGFR3","GFAP", "OPALIN", "MOG", "FYB1", "APBB1IP", "PTPRC", "TYROBP","CSPG4", "PDGFRA", "MYT1", "SLC17A7", "RORB", "LINC00507", "THEMIS","FEZF2", "CCN2", "FAM84B","POU3F1","NPSR1","GAD1", "GAD2", "PVALB", "LHX6", "SST", "LAMP5", "PAX6", "VIP"))

qs::qsave(seu, file = "dfpfc.qs",nthreads = 20)

Lewy.body.disease.pathology <- c( "reference", "Neocortical (Diffuse)","Olfactory bulb only","Limbic (Transitional)", "Brainstem-predominant", "Amygdala-predominant"
)

Lewy.body.disease.pathology_color = MetBrewer::met.brewer(name = "Tiepolo",n = length(unique(seu$Lewy.body.disease.pathology)),direction = -1) %>% setNames(c( "reference", "Neocortical (Diffuse)","Olfactory bulb only","Limbic (Transitional)", "Brainstem-predominant", "Amygdala-predominant"
))


subclass_color <- MetBrewer::met.brewer(name = "Juarez",n = length(unique(seu$Subclass))) %>% setNames(seu@meta.data %>% select(Class,Subclass) %>% distinct() %>% arrange(Class,Subclass) %$% Subclass )

subclass_color <- c(Chandelier = "#A82203", Lamp5 = "#8A392C", `Lamp5 Lhx6` = "#6C5055",
                    Pax6 = "#4F677E", Pvalb = "#317EA7", Sncg = "#328FB4", Sst = "#5F9697",
                    `Sst Chodl` = "#8D9E7A", Vip = "#BAA55C", `L2/3 IT` = "#E7AD3F",
                    `L4 IT` = "#EBA03D", `L5 ET` = "#E38F41", `L5 IT` = "#DC7D46",
                    `L5/6 NP` = "#D46C4A", `L6 CT` = "#CA5F4C", `L6 IT` = "#B26546",
                    `L6 IT Car3` = "#9B6B40", L6b = "#837239", Astrocyte = "#6C7833",
                    Endothelial = "#567238", `Microglia-PVM` = "#406443", OPC = "#2B554F",
                    Oligodendrocyte = "#15475B", VLMC = "#003967")

#



{DimPlot(seu,group.by = "Lewy.body.disease.pathology",split.by = "Lewy.body.disease.pathology",shuffle = T,label.box = T,label = T,repel = T,raster = F,label.size = 8) + theme(legend.position = "none")  +  scale_color_manual(values = Lewy.body.disease.pathology_color)  + scale_fill_manual(values = Lewy.body.disease.pathology_color)+ cowplot::theme_nothing()}%T>% ggsave(filename = "dlpfc_global_Lewy.body.disease.pathology_umap.png",height = 7,width = 30,dpi = 300)







{ggplot(seu@meta.data %>% group_by(Lewy.body.disease.pathology,Subclass,Specimen.ID) %>% summarise(count = n()) %>% group_by(Specimen.ID) %>% mutate(per = count/sum(count)),mapping = aes(Lewy.body.disease.pathology,y = per)) + geom_boxplot(aes(fill = Subclass)) + ggsignif::geom_signif(comparisons = map(Lewy.body.disease.pathology[-c(1)],~c("reference",.)),map_signif_level = T,step_increase = 0.05,margin_top = -0.3 ) + facet_wrap(~Subclass,scales = "free_y") + cowplot::theme_cowplot() + theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))} %T>% ggsave(filename = "dlpfc_global_subclass_lewy_ratio.pdf",height = 12,width = 15,dpi = 300)

{ggplot(seu@meta.data %>% group_by(Lewy.body.disease.pathology,Class,Specimen.ID) %>% summarise(count = n()) %>% group_by(Specimen.ID) %>% mutate(per = count/sum(count)),mapping = aes(Lewy.body.disease.pathology,y = per)) + geom_boxplot(aes(fill = Class)) + facet_wrap(~Class,scales = "free_y") + ggsignif::geom_signif(comparisons = map(Lewy.body.disease.pathology[-c(1)],~c("reference",.)),step_increase = 0.05,margin_top = -0.3,map_signif_level = T) + cowplot::theme_cowplot()+ theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))  } %T>% ggsave(filename = "dlpfc_global_class_lewy_ratio.pdf",height = 10,width = 15,dpi = 300)


seu_list <- seu %>% SplitObject("Subclass")
library(parallel)
x = seu_list[["L6 CT"]]
degTable <- parallel::mclapply(seu_list,mc.cores = length(seu_list),function(x) {
  l = Lewy.body.disease.pathology[[1]]
  res <- parallel::mclapply(Lewy.body.disease.pathology[-c(1)],mc.cores = length(Lewy.body.disease.pathology)-1,function(l){
    res <- Seurat::FindMarkers(x,group.by = "Lewy.body.disease.pathology",ident.1 = l,ident.2 = NULL)
    res$cluster <- l
    res$gene <- rownames(res)
    res$subclass <- x$Subclass[[1]]
    res
  }) %>% purrr::reduce(rbind)
}
) %>% purrr::reduce(rbind)

degTable %>% write.csv("dlpfc_Lewy.body.disease.pathology_deg.csv")



{ggplot( degTable  %>% filter(p_val_adj < 0.2,avg_log2FC >0)  %>%
           dplyr::select(subclass,cluster) %>% table %>% as.data.frame(),aes(x = subclass,y = Freq, group = cluster, fill = cluster)) + geom_bar(stat = "identity",position = "dodge") + scale_y_log10() + cowplot::theme_cowplot() +  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + MetBrewer::scale_fill_met_d(name = "Navajo")} %>%
  ggsave(filename = "dlpfc_lewy_num_of_deg_up.pdf",height = 6,width = 15,dpi = 300)

{ggplot( degTable  %>% filter(p_val_adj < 0.2,avg_log2FC <0)  %>%
           dplyr::select(subclass,cluster) %>% table %>% as.data.frame(),aes(x = subclass,y = Freq, group = cluster, fill = cluster)) + geom_bar(stat = "identity",position = "dodge") + scale_y_log10() +  cowplot::theme_cowplot() +   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + MetBrewer::scale_fill_met_d(name = "Navajo")} %>%
  ggsave(filename = "dlpfc_lewy_num_of_deg_down.pdf",height = 6,width = 15,dpi = 300)

library(ClusterGVis)


celltype = "L5 ET"
for (celltype in c("Astrocyte", "Chandelier", "Endothelial", "L2/3 IT", "L4 IT",
  "L5 ET", "L5 IT",
  "L5/6 NP", "L6 CT", "L6 IT", "L6 IT Car3",
  "L6b", "Lamp5", "Lamp5 Lhx6", "Microglia-PVM", "Oligodendrocyte",
  "OPC", "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip", "VLMC"
)) {
  
  tmpSeu <- seu_list[[celltype]]
  Idents(tmpSeu) <- "Lewy.body.disease.pathology"
  print(timestamp())
  print(celltype)
  
  print("up")
  markers <- degTable %>% filter(avg_log2FC >0, subclass == celltype,p_val < 0.05) %>%
    dplyr::group_by(gene) %>%
    dplyr::top_n(n = 1,wt = avg_log2FC) %>%
    ungroup %>%
    dplyr::top_n(n = 15, wt = avg_log2FC) %>%
    mutate(unluckGene = str_detect(gene,"\\.[0-9]$"))

  markers %<>% mutate(cluster = factor(cluster,levels = Lewy.body.disease.pathology)) %>% arrange(cluster,unluckGene)

  st.dataPlot <- prepareDataFromscRNA(object = tmpSeu,
                                      diffData = markers,
                                      showAverage = T)
  
  


  enrich <- parallel::mclapply(unique(st.dataPlot$wide.res$cluster) %>% sort,mc.cores = 4,function(i){

    res <- clusterProfiler::enrichGO(gene =geneSymToID[st.dataPlot$wide.res %>% filter(cluster == i) %$% gene] ,ont = "BP",
                                     pAdjustMethod = "BH",
                                     keyType  = "ENSEMBL", pvalueCutoff = 1,
                                     qvalueCutoff = 0.2,
                                     OrgDb = org.Hs.eg.db)
    if (is.null(res)) {
      return(res)
    }
    res@result$group <- str_c("C",i)

    res@result$ratio <- sapply(res@result$GeneRatio, function(x) eval(parse(text = x)))*100

    res@result %>% arrange(pvalue) %>% slice_head(n = 5) %>% dplyr::select(group,Description,pvalue,ratio)

  }) %>% keep(~!is.null(.)) %>% purrr::reduce(rbind)



  

  markGenes <- markers %>% filter(avg_log2FC >0)  %>%
    dplyr::group_by(cluster) %>%
    arrange(unluckGene,desc(abs(avg_log2FC))) %>%
    slice_head(n = 10) %$% gene

  


  pdf(str_c("dlpfc_",str_replace(celltype,"/","_"),"_lewy_up.pdf"),height = 10/nrow(markers)*10,width = 2*(ncol(st.dataPlot$wide.res)-2),onefile = F)
  visCluster(object = st.dataPlot,
             plot.type = "both",
             column_names_rot = 45,
             markGenes = markGenes,
             sample.col =  Lewy.body.disease.pathology_color,
             markGenes.side = "left",
             annoTerm.data = enrich,
             line.side = "left",ctAnno.col =  Lewy.body.disease.pathology_color[markers$cluster %>% unique],
             cluster.order = c(1:length(unique(markers))),
             go.col =  Lewy.body.disease.pathology_color[unique(markers$cluster)[enrich$group %>% str_extract("\\d+") %>% as.numeric()]]
             
             
  )
  dev.off()

  
  print("down")
  markers <- degTable %>% filter(avg_log2FC <0, subclass == celltype,p_val < 0.05) %>%
    dplyr::group_by(gene) %>%
    dplyr::top_n(n = 1,wt = avg_log2FC) %>%
    ungroup %>%
    dplyr::top_n(n = 15, wt = avg_log2FC) %>%
    mutate(unluckGene = str_detect(gene,"\\.[0-9]$"))

  markers %<>% mutate(cluster = factor(cluster,levels = Lewy.body.disease.pathology)) %>% arrange(cluster,unluckGene)

  st.dataPlot <- prepareDataFromscRNA(object = tmpSeu,
                                      diffData = markers,
                                      showAverage = T)
  tmp <- st.dataPlot
  tmp$wide.res <- st.dataPlot$wide.res %>% filter(str_detect(gene,"\\.[0-9]$",negate = T))


  enrich <- parallel::mclapply(unique(st.dataPlot$wide.res$cluster) %>% sort,mc.cores = 4,function(i){

    res <- clusterProfiler::enrichGO(gene =geneSymToID[st.dataPlot$wide.res %>% filter(cluster == i) %$% gene] ,ont = "BP",
                                     pAdjustMethod = "BH",
                                     keyType  = "ENSEMBL", pvalueCutoff = 1,
                                     qvalueCutoff = 0.2,
                                     OrgDb = org.Hs.eg.db)
    if (is.null(res)) {
      return(res)
    }
    res@result$group <- str_c("C",i)

    res@result$ratio <- sapply(res@result$GeneRatio, function(x) eval(parse(text = x)))*100

    res@result %>% arrange(pvalue) %>% slice_head(n = 5) %>% dplyr::select(group,Description,pvalue,ratio)

  }) %>% keep(~!is.null(.)) %>% purrr::reduce(rbind)



  

  markGenes <- markers %>% filter(avg_log2FC <0)  %>%
    dplyr::group_by(cluster) %>%
    arrange(cluster ,unluckGene,desc(abs(avg_log2FC))) %>%
    slice_head(n = 3) %$% gene

  


  pdf(str_c("dlpfc_",str_replace(celltype,"/","_"),"_lewy_down.pdf"),height = 10/nrow(markers)*10,width = 2*(ncol(st.dataPlot$wide.res)-2),onefile = F)
  visCluster(object = st.dataPlot,
             plot.type = "both",
             column_names_rot = 45,
             markGenes = markGenes,
             sample.col =  Lewy.body.disease.pathology_color,
             markGenes.side = "left",
             annoTerm.data = enrich,
             line.side = "left",
             ctAnno.col =  Lewy.body.disease.pathology_color[markers$cluster %>% unique],
             cluster.order = c(1:length(unique(markers))),
             go.col = Lewy.body.disease.pathology_color[unique(markers$cluster)[enrich$group %>% str_extract("\\d+") %>% as.numeric()]]
             
             
  )
  dev.off()
}




