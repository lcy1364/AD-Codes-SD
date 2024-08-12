library(Seurat)
library(SeuratDisk)
library(magrittr)
library(tidyverse)
library(org.Hs.eg.db)
setwd("~/DATA/BRAIN/AD")

meta <- read.csv("./dlpfc_meta.csv",row.names = 1)


meta$ADNC %<>% fct_collapse(N = c("Reference","Not AD"))
meta$ADNC <- str_c("ADNC_",meta$ADNC)


meta$APOE4.status[meta$APOE4.status == "Reference" | meta$APOE4.status == "N"] = "Apoe4_N"
meta$APOE4.status[meta$APOE4.status == "Y"] = "Apoe4_Y"


meta$Lewy.body.disease.pathology %<>% fct_collapse(N = c("Reference","Not Identified (olfactory bulb assessed)","Not Identified (olfactory bulb not assessed)"))
meta$Lewy.body.disease.pathology <- str_c("Lewy_",meta$Lewy.body.disease.pathology)


meta <- meta %>% unite(col = "ADNC.APOE4.LEWY",,ADNC,APOE4.status,Lewy.body.disease.pathology,remove = F,sep = "/")
meta$ADNC.APOE4.LEWY %<>% factor( levels= c(
  "ADNC_N/Apoe4_N/Lewy_N",
  "ADNC_N/Apoe4_N/Lewy_Neocortical (Diffuse)",
  "ADNC_N/Apoe4_N/Lewy_Limbic (Transitional)",
  "ADNC_Low/Apoe4_N/Lewy_N",
  "ADNC_Low/Apoe4_N/Lewy_Neocortical (Diffuse)",
  "ADNC_Low/Apoe4_N/Lewy_Limbic (Transitional)",
  "ADNC_Intermediate/Apoe4_N/Lewy_N",
  "ADNC_Intermediate/Apoe4_N/Lewy_Neocortical (Diffuse)",
  "ADNC_Intermediate/Apoe4_N/Lewy_Limbic (Transitional)",
  "ADNC_Intermediate/Apoe4_Y/Lewy_N",
  "ADNC_Intermediate/Apoe4_Y/Lewy_Neocortical (Diffuse)",
  "ADNC_High/Apoe4_N/Lewy_Neocortical (Diffuse)",
  "ADNC_High/Apoe4_N/Lewy_N",
  "ADNC_High/Apoe4_N/Lewy_Limbic (Transitional)",
  "ADNC_High/Apoe4_N/Lewy_Amygdala-predominant",
  "ADNC_High/Apoe4_Y/Lewy_N",
  "ADNC_High/Apoe4_Y/Lewy_Neocortical (Diffuse)",
  "ADNC_High/Apoe4_Y/Lewy_Limbic (Transitional)",
  "ADNC_High/Apoe4_Y/Lewy_Amygdala-predominant"))


library(rtracklayer)
gene <- rtracklayer::import("./gencode.v32.primary_assembly.annotation.gtf") %>% as.data.frame()
gene %<>% dplyr::select(gene_name,gene_id,seqnames) %>% distinct() %>% mutate(gene_id = str_remove(gene_id,"\\..*"))
gene %<>% filter(str_detect(seqnames,"chr\\d+"))

geneID_to_sym <- gene$gene_name %>% setNames(gene$gene_id)
geneSymToID <-  gene$gene_id %>% setNames(gene$gene_name)





seu_org <- LoadH5Seurat("./dlpfc.h5seurat")
count <- GetAssayData(seu_org,slot = "count")
rownames(count) <- geneID_to_sym[rownames(count)]
count <- count[rownames(count) %in% gene$gene_name,]
seu <- CreateSeuratObject(counts = count,meta.data = meta)
seu@reductions <- seu_org@reductions
rm(seu_org)
rm(count)
gc()

seu %<>% NormalizeData() %>% ScaleData(features = c("NOSTRIN", "FLT1", "CLDN5", "SLC1A3", "MUSTN1", "AQP4", "FGFR3","GFAP", "OPALIN", "MOG", "FYB1", "APBB1IP", "PTPRC", "TYROBP","CSPG4", "PDGFRA", "MYT1", "SLC17A7", "RORB", "LINC00507", "THEMIS","FEZF2", "CCN2", "FAM84B","POU3F1","NPSR1","GAD1", "GAD2", "PVALB", "LHX6", "SST", "LAMP5", "PAX6", "VIP"))

qs::qsave(seu, file = "dlpfc.qs",nthreads = 20)


subclassToClass <- seu@meta.data %>% dplyr::select(Class,Subclass) %>% distinct() %>% {setNames(.$Class,nm = .$Subclass)}


seu <- seu[,seu$Lewy.body.disease.pathology %>% str_detect("Amy|Limbic|Neocortical|Lewy_N")]

seu$ADNC.APOE4.LEWY %<>% factor(levels = ADNC.APOE4.LEWY_order)

Idents(seu) <- seu$ADNC.APOE4.LEWY


ADNC.APOE4.LEWY_order_color = MetBrewer::met.brewer(name = "Tiepolo",n = length(ADNC.APOE4.LEWY_order),direction = -1) %>% setNames(ADNC.APOE4.LEWY_order)

subclass_color <- c(Chandelier = "#A82203", Lamp5 = "#8A392C", `Lamp5 Lhx6` = "#6C5055",
                    Pax6 = "#4F677E", Pvalb = "#317EA7", Sncg = "#328FB4", Sst = "#5F9697",
                    `Sst Chodl` = "#8D9E7A", Vip = "#BAA55C", `L2/3 IT` = "#E7AD3F",
                    `L4 IT` = "#EBA03D", `L5 ET` = "#E38F41", `L5 IT` = "#DC7D46",
                    `L5/6 NP` = "#D46C4A", `L6 CT` = "#CA5F4C", `L6 IT` = "#B26546",
                    `L6 IT Car3` = "#9B6B40", L6b = "#837239", Astrocyte = "#6C7833",
                    Endothelial = "#567238", `Microglia-PVM` = "#406443", OPC = "#2B554F",
                    Oligodendrocyte = "#15475B", VLMC = "#003967")







{ggplot(seu@meta.data %>% group_by(ADNC.APOE4.LEWY,Class,Specimen.ID) %>% summarise(count = n()) %>% group_by(Specimen.ID) %>% mutate(per = count/sum(count)),mapping = aes(factor(ADNC.APOE4.LEWY,levels = rev(ADNC.APOE4.LEWY_order)),y = per)) + geom_boxplot(aes(fill = Class)) + facet_wrap(~Class,scales = "free_y",ncol = 1) + scale_x_discrete(position = "top") + ggsignif::geom_signif(comparisons = map(ADNC.APOE4.LEWY_order[-c(1)],~c("ADNC_N/Apoe4_N/Lewy_N",.)),step_increase = 0.05,margin_top = -0.3,map_signif_level = T) + cowplot::theme_cowplot()+ coord_flip() + theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1))  } %T>% ggsave(filename = "dlpfc_class_class_ADNC.APOE4.LEWY_ratio.pdf",height = 15,width = 15,dpi = 300)


{
  ggplot(
    seu@meta.data %>% group_by(ADNC.APOE4.LEWY, Subclass, Specimen.ID) %>% summarise(count = n()) %>% group_by(Specimen.ID) %>% mutate(per = count /sum(count))  %>% filter(Subclass == "L4 IT"),
    mapping = aes(ADNC.APOE4.LEWY,y = per)
  ) + geom_boxplot(aes(fill = Subclass)) + ggsignif::geom_signif(
    comparisons = map(ADNC.APOE4.LEWY_order[-c(1)],  ~ c("ADNC_N/Apoe4_N/Lewy_N", .)),
    step_increase = 0.05,
    margin_top = -0.3,
    map_signif_level = T
  ) + cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(
    angle = 270,
    hjust = 0,
    vjust = 1
  ))
} %T>% ggsave(filename = "dlpfc_class_L4_IT_ADNC.APOE4.LEWY_ratio.pdf",height = 15,width = 15,dpi = 300)


{
  ggplot(
    seu@meta.data %>% group_by(ADNC.APOE4.LEWY, Subclass, Specimen.ID) %>% summarise(count = n()) %>% group_by(Specimen.ID) %>% mutate(per = count /sum(count))  %>% filter(str_detect(Subclass,"Micro|Astro")),
    mapping = aes(factor(
      ADNC.APOE4.LEWY, levels = ADNC.APOE4.LEWY_order
    ), y = per)
  ) + geom_boxplot(aes(fill = Subclass)) + ggsignif::geom_signif(
    comparisons = map(ADNC.APOE4.LEWY_order[-c(1)],  ~ c("ADNC_N/Apoe4_N/Lewy_N", .)),
    step_increase = 0.05,
    margin_top = -0.3,
    map_signif_level = T
  ) + cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(
      angle = 270,
      hjust = 0,
      vjust = 1
    ))
} %T>% ggsave(filename = "dlpfc_class_Micro_Astro_ADNC.APOE4.LEWY_ratio.pdf",height = 15,width = 15,dpi = 300)

 
seu_list <- seu %>% SplitObject("Subclass")

x = seu_list[["L6 CT"]]

degTable <-
  parallel::mclapply(seu_list, mc.cores = length(seu_list), function(x) {
    res <-
      parallel::mclapply(ADNC.APOE4.LEWY_order[-c(1)], mc.cores = length(ADNC.APOE4.LEWY_order) -
                           1, function(l) {
                             res <- Seurat::FindMarkers(x, ident.1 = l, ident.2 = NULL)
                             res$cluster <- l
                             res$gene <- rownames(res)
                             res$subclass <- x$Subclass[[1]]
                             res
                           }) %>% purrr::reduce(rbind)

  }) %>% purrr::reduce(rbind)

degTable %>% write.csv("dlpfc_ADNC.APOE4.LEWY_deg.csv")


{
    ggplot(mapping = aes(
      x =  factor(cluster, levels = ADNC.APOE4.LEWY_order),
      y = Freq
    )) + geom_bar(
      data = degTable  %>% mutate(Up = subclassToClass[subclass]) %>% filter(p_val_adj < 0.2, avg_log2FC > 0)  %>% dplyr::select(Up, cluster) %>% table %>% as.data.frame() %>% filter(str_detect(cluster, "ADNC_N", negate = T)) %>%  mutate(Freq = log10(Freq+1)),
      mapping = aes(fill = Up,group = Up),
      stat = "identity",
      position = "dodge"
    ) + MetBrewer::scale_fill_met_d(name = "Navajo") +
    ggnewscale::new_scale_fill() +
    geom_bar(
      data = degTable  %>% mutate(Down = subclassToClass[subclass]) %>% filter(p_val_adj < 0.2, avg_log2FC <  0)  %>% dplyr::select(Down, cluster) %>% table %>% as.data.frame() %>% filter(str_detect(cluster, "ADNC_N", negate = T)) %>%  mutate(Freq = -log10(Freq+1)),
      mapping = aes(fill = Down,group = Down),
      stat = "identity",
      position = "dodge"
    )  +  cowplot::theme_cowplot() +   theme(axis.text.x = element_text(
      angle = 270,
      hjust = 0,
      vjust = 1
    ))  +
    MetBrewer::scale_fill_met_d(name = "Hokusai2") +
    theme(axis.line.x = element_blank())+ geom_hline(yintercept = 0) + coord_flip() + scale_y_continuous(breaks = -3:3,labels = str_c("10^",as.character(abs(-3:3))))
  } %>%
    ggsave(
      filename = "dlpfc_ADNC.APOE4.LEWY_num_of_deg.pdf",
      height = 10,
      width = 13,
      dpi = 300
    )



{
  ggplot(mapping = aes(
    x =  factor(cluster, levels = ADNC.APOE4.LEWY_order),
    y = Freq
  )) + geom_bar(
    data = degTable  %>% mutate(Up = subclass) %>% filter(p_val_adj < 0.2, avg_log2FC > 0)  %>% filter( subclass == "L4 IT") %>% dplyr::select(Up, cluster) %>% table %>% as.data.frame() %>% filter(str_detect(cluster, "ADNC_N", negate = T)) %>%  mutate(Freq = log10(Freq+1)),
    mapping = aes(fill = Up,group = Up),
    stat = "identity",
    position = "dodge"
  ) + MetBrewer::scale_fill_met_d(name = "Navajo") +
    ggnewscale::new_scale_fill() +
    geom_bar(
      data = degTable  %>% mutate(Down = subclass) %>% filter(p_val_adj < 0.2, avg_log2FC <  0) %>% filter( subclass == "L4 IT") %>% dplyr::select(Down, cluster) %>% table %>% as.data.frame() %>% filter(str_detect(cluster, "ADNC_N", negate = T)) %>%  mutate(Freq = -log10(Freq+1)),
      mapping = aes(fill = Down,group = Down),
      stat = "identity",
      position = "dodge"
    )  +  cowplot::theme_cowplot() +   theme(axis.text.x = element_text(
      angle = 270,
      hjust = 0,
      vjust = 1
    ))  +
    MetBrewer::scale_fill_met_d(name = "Hokusai2") +
    theme(axis.line.x = element_blank())+ geom_hline(yintercept = 0) + coord_flip() + scale_y_continuous(breaks = -3:3,labels = str_c("10^",as.character(abs(-3:3))))
    } %>%
  ggsave(
    filename = "dlpfc_ADNC.APOE4.LEWY_L4_IT_num_of_deg.pdf",
    height = 10,
    width = 13,
    dpi = 300
  )



{
  ggplot(mapping = aes(
    x =  factor(cluster, levels = ADNC.APOE4.LEWY_order),
    y = Freq
  )) + geom_bar(
    data = degTable  %>% mutate(Up = subclass) %>% filter(p_val_adj < 0.2, avg_log2FC > 0)  %>% filter(str_detect(subclass,"Micro|Astro")) %>% dplyr::select(Up, cluster) %>% table %>% as.data.frame() %>% filter(str_detect(cluster, "ADNC_N", negate = T)) %>%  mutate(Freq = log10(Freq+1)),
    mapping = aes(fill = Up,group = Up),
    stat = "identity",
    position = "dodge"
  ) + MetBrewer::scale_fill_met_d(name = "Navajo") +
    ggnewscale::new_scale_fill() +
    geom_bar(
      data = degTable  %>% mutate(Down = subclass) %>% filter(p_val_adj < 0.2, avg_log2FC <  0) %>% filter(str_detect(subclass,"Micro|Astro")) %>% dplyr::select(Down, cluster) %>% table %>% as.data.frame() %>% filter(str_detect(cluster, "ADNC_N", negate = T)) %>%  mutate(Freq = -log10(Freq+1)),
      mapping = aes(fill = Down,group = Down),
      stat = "identity",
      position = "dodge"
    )  +  cowplot::theme_cowplot() +   theme(axis.text.x = element_text(
      angle = 270,
      hjust = 0,
      vjust = 1
    ))  +
    MetBrewer::scale_fill_met_d(name = "Hokusai2") +
    theme(axis.line.x = element_blank())+ geom_hline(yintercept = 0) + coord_flip() + scale_y_continuous(breaks = -3:3,labels = str_c("10^",as.character(abs(-3:3))))
  } %>%
  ggsave(
    filename = "dlpfc_ADNC.APOE4.LEWY_Micro_Astro_num_of_deg.pdf",
    height = 10,
    width = 13,
    dpi = 300
  )



diseaseRNA <- read.delim("./brain_disease") %>% filter(str_detect(Disease,"Alzheimer"), Type == "mRNA") %$% Gene.symbol %>% unique
length(diseaseRNA)

{
  ggplot(mapping = aes(
    x =  factor(cluster, levels = ADNC.APOE4.LEWY_order),
    y = Freq
  )) + geom_bar(
    data = degTable  %>% mutate(Up = subclassToClass[subclass]) %>% filter(p_val_adj < 0.2, avg_log2FC > 0, gene %in% diseaseRNA)  %>% dplyr::select(Up, cluster) %>% table %>% as.data.frame() %>% filter(str_detect(cluster, "ADNC_N", negate = T)) %>%  mutate(Freq = Freq/length(diseaseRNA)),
    mapping = aes(fill = Up,group = Up),
    stat = "identity",
    position = "dodge"
  ) + MetBrewer::scale_fill_met_d(name = "Navajo") +
    ggnewscale::new_scale_fill() +
    geom_bar(
      data = degTable  %>% mutate(Down = subclassToClass[subclass]) %>% filter(p_val_adj < 0.2, avg_log2FC <0, gene %in% diseaseRNA)  %>% dplyr::select(Down, cluster) %>% table %>% as.data.frame() %>% filter(str_detect(cluster, "ADNC_N", negate = T)) %>% mutate(Freq = -Freq/length(diseaseRNA)),
      mapping = aes(fill = Down,group = Down),
      stat = "identity",
      position = "dodge"
    )  +  cowplot::theme_cowplot() +   theme(axis.text.x = element_text(
      angle = 270,
      hjust = 0,
      vjust = 1
    ))  +
    MetBrewer::scale_fill_met_d(name = "Hokusai2") +
    theme(axis.line.x = element_blank())+ geom_hline(yintercept = 0) + coord_flip()
} %T>%
  ggsave(
    filename = "dlpfc_ADNC.APOE4.LEWY_AD_RiskGenePer.pdf",
    height = 10,
    width = 13,
    dpi = 300
  )



{
  ggplot(mapping = aes(
    x =  factor(cluster, levels = ADNC.APOE4.LEWY_order),
    y = Freq
  )) + geom_bar(
    data = degTable  %>% mutate(Up = subclass) %>% filter(p_val_adj < 0.2, avg_log2FC > 0, gene %in% diseaseRNA)  %>% filter( subclass == "L4 IT") %>% dplyr::select(Up, cluster) %>% table %>% as.data.frame() %>% filter(str_detect(cluster, "ADNC_N", negate = T)) %>%   mutate(Freq = Freq/length(diseaseRNA)),
    mapping = aes(fill = Up,group = Up),
    stat = "identity",
    position = "dodge"
  ) + MetBrewer::scale_fill_met_d(name = "Navajo") +
    ggnewscale::new_scale_fill() +
    geom_bar(
      data = degTable  %>% mutate(Down = subclass) %>% filter(p_val_adj < 0.2, avg_log2FC <  0, gene %in% diseaseRNA) %>% filter( subclass == "L4 IT") %>% dplyr::select(Down, cluster) %>% table %>% as.data.frame() %>% filter(str_detect(cluster, "ADNC_N", negate = T)) %>%   mutate(Freq = -Freq/length(diseaseRNA)),
      mapping = aes(fill = Down,group = Down),
      stat = "identity",
      position = "dodge"
    )  +  cowplot::theme_cowplot() +   theme(axis.text.x = element_text(
      angle = 270,
      hjust = 0,
      vjust = 1
    ))  +
    MetBrewer::scale_fill_met_d(name = "Hokusai2") +
    theme(axis.line.x = element_blank())+ geom_hline(yintercept = 0) + coord_flip()
  } %T>%
  ggsave(
    filename = "dlpfc_ADNC.APOE4.LEWY_L4_IT_AD_RiskGenePer.pdf",
    height = 7,
    width = 10,
    dpi = 300
  )



{
  ggplot(mapping = aes(
    x =  factor(cluster, levels = ADNC.APOE4.LEWY_order),
    y = Freq
  )) + geom_bar(
    data = degTable  %>% mutate(Up = subclass) %>% filter(p_val_adj < 0.2, avg_log2FC > 0, gene %in% diseaseRNA)  %>% filter(str_detect(subclass,"Micro|Astro")) %>% dplyr::select(Up, cluster) %>% table %>% as.data.frame() %>% filter(str_detect(cluster, "ADNC_N", negate = T))  %>% mutate(Freq = Freq/length(diseaseRNA)),
    mapping = aes(fill = Up,group = Up),
    stat = "identity",
    position = "dodge"
  ) + MetBrewer::scale_fill_met_d(name = "Navajo") +
    ggnewscale::new_scale_fill() +
    geom_bar(
      data = degTable  %>% mutate(Down = subclass) %>% filter(p_val_adj < 0.2, avg_log2FC <  0, gene %in% diseaseRNA) %>% filter(str_detect(subclass,"Micro|Astro")) %>% dplyr::select(Down, cluster) %>% table %>% as.data.frame() %>% filter(str_detect(cluster, "ADNC_N", negate = T)) %>%  mutate(Freq = -Freq/length(diseaseRNA)),
      mapping = aes(fill = Down,group = Down),
      stat = "identity",
      position = "dodge"
    )  +  cowplot::theme_cowplot() +   theme(axis.text.x = element_text(
      angle = 270,
      hjust = 0,
      vjust = 1
    ))  +
    MetBrewer::scale_fill_met_d(name = "Hokusai2") +
    theme(axis.line.x = element_blank())+ geom_hline(yintercept = 0) + coord_flip()
  } %T>%
  ggsave(
    filename = "dlpfc_ADNC.APOE4.LEWY_Micro_Astro_AD_RiskGenePer.pdf",
    height = 10,
    width = 13,
    dpi = 300
  )
library(ClusterGVis)


celltype = "Astrocyte"
for (celltype in c("Astrocyte", "Chandelier", "Endothelial", "L2/3 IT", "L4 IT",
                   "L5 ET", "L5 IT",
                   "L5/6 NP", "L6 CT", "L6 IT", "L6 IT Car3",
                   "L6b", "Lamp5", "Lamp5 Lhx6", "Microglia-PVM", "Oligodendrocyte",
                   "OPC", "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip", "VLMC"
)) {
  
  tmpSeu <- seu_list[[celltype]]
  Idents(tmpSeu) <- "APOE4.ADNC"
  print(timestamp())
  print(celltype)
  
  print("up")
  markers <- degTable %>% filter(avg_log2FC >0, subclass == celltype,p_val < 0.05) %>%
    dplyr::group_by(gene) %>%
    dplyr::top_n(n = 1,wt = avg_log2FC) %>%
    ungroup %>%
    dplyr::top_n(n = 15, wt = avg_log2FC) %>%
    mutate(unluckGene = str_detect(gene,"\\.[0-9]$"))

  markers %<>% mutate(cluster = factor(cluster,levels = apoe4)) %>% arrange(cluster,unluckGene)

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

  


  pdf(str_c("dlpfc_",str_replace(celltype,"/","_"),"_apoe_up.pdf"),height = 10/nrow(markers)*10,width = 2*(ncol(st.dataPlot$wide.res)-2),onefile = F)
  visCluster(object = st.dataPlot,
             plot.type = "both",
             column_names_rot = 45,
             markGenes = markGenes,
             sample.col =  apoe4.body.disease.pathology_color,
             markGenes.side = "left",
             annoTerm.data = enrich,
             line.side = "left",
             ctAnno.col =  apoe4.body.disease.pathology_color[markers$cluster %>% unique],
             cluster.order = c(1:length(unique(markers))),
             go.col =  apoe4.body.disease.pathology_color[unique(markers$cluster)[enrich$group %>% str_extract("\\d+") %>% as.numeric()]]
             
             
  )
  dev.off()

  
  print("down")
  markers <- degTable %>% filter(avg_log2FC <0, subclass == celltype,p_val < 0.05) %>%
    dplyr::group_by(gene) %>%
    dplyr::top_n(n = 1,wt = avg_log2FC) %>%
    ungroup %>%
    dplyr::top_n(n = 15, wt = avg_log2FC) %>%
    mutate(unluckGene = str_detect(gene,"\\.[0-9]$"))

  markers %<>% mutate(cluster = factor(cluster,levels = apoe4)) %>% arrange(cluster,unluckGene)

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

  


  pdf(str_c("dlpfc_",str_replace(celltype,"/","_"),"_apoe_down.pdf"),height = 10/nrow(markers)*10,width = 2*(ncol(st.dataPlot$wide.res)-2),onefile = F)
  visCluster(object = st.dataPlot,
             plot.type = "both",
             column_names_rot = 45,
             markGenes = markGenes,
             sample.col =  apoe4.body.disease.pathology_color,
             markGenes.side = "left",
             annoTerm.data = enrich,
             line.side = "left",
             ctAnno.col =  apoe4.body.disease.pathology_color[markers$cluster %>% unique],
             cluster.order = c(1:length(unique(markers))),
             go.col = apoe4.body.disease.pathology_color[unique(markers$cluster)[enrich$group %>% str_extract("\\d+") %>% as.numeric()]],
             
             
  )
  dev.off()
}




data <- meta %>% filter(!is.na(ADNC.APOE4.LEWY)) %>% group_by(ADNC,APOE4.status,Lewy.body.disease.pathology,ADNC.APOE4.LEWY,Class,Specimen.ID) %>% summarise(count = n()) %>% group_by(Specimen.ID) %>% mutate(per = count/sum(count))



model <- lm(per ~ ADNC * APOE4.status *Lewy.body.disease.pathology* Class, data=data)

summary_model <- summary(model)
coefficients <- summary_model$coefficients
r_squared <- summary_model$r.squared



plot <- ggplot(data, aes(x=ADNC.APOE4.LEWY, y=per, group = Class,color = Class)) +
  geom_point() +  
  
  geom_smooth(method = "loess", se=FALSE) +  
  labs(title='ADNC', x='ADNC * APOE', y='percentage', color='Category') +
  theme_minimal()+ labs(subtitle = str_c("R squared ",round(r_squared,digits = 2))) + theme(axis.text.x = element_text(angle = -45,hjust = 0))

plot %>% ggsave(filename = "dlpfc_regression_adnc_apoe_lewy_global.png",height = 5,width = 7)



ct = unique(data$Class)[[1]]
for (ct in unique(data$Class)) {
  tmp = data %>% filter(Class == ct)

  
  model <- lm(per ~ ADNC * APOE4.status *Lewy.body.disease.pathology, data=tmp)

  summary_model <- summary(model)
  coefficients <- summary_model$coefficients
  r_squared <- summary_model$r.squared
  

  plot <- ggplot(tmp, aes(x= ADNC.APOE4.LEWY, y=per, group = Class,color = Class)) +
    geom_point() +  
    
    geom_smooth(method = "loess", se=FALSE) +  
    labs(title='ADNC', x='ADNC * APOE', y='percentage', color='Category') +
    theme_minimal()+ labs(subtitle = str_c("R squared ",round(r_squared,digits = 2))) + theme(axis.text.x = element_text(angle = -45,hjust = 0))

  plot %>% ggsave(filename = str_c("dlpfc_regression_adnc_apoe_lewy ",ct,"_.png"),height = 5,width = 7)
}



