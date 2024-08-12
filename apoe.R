library(Seurat)
library(SeuratDisk)
library(magrittr)
library(tidyverse)
library(org.Hs.eg.db)
setwd("~/DATA/BRAIN/AD")



meta <- read.csv("./dlpfc_meta.csv",row.names = 1)


meta$APOE4.status[meta$APOE4.status == "Reference" | meta$APOE4.status == "N"] = "Apoe4_N"
meta$APOE4.status[meta$APOE4.status == "Y"] = "Apoe4_Y"


meta$ADNC %<>% fct_collapse(N = c("Reference","Not AD"))
meta$ADNC <- str_c("ADNC_",meta$ADNC)


seu@meta.data <- meta[rownames(seu@meta.data ),]


seu@meta.data <- seu@meta.data %>% unite(col = "APOE4.ADNC",APOE4.status,ADNC,remove = F,sep = "/")

apoe4_order <-c( "Apoe4_N/ADNC_N", "Apoe4_N/ADNC_Low","Apoe4_N/ADNC_Intermediate","Apoe4_N/ADNC_High", "Apoe4_Y/ADNC_Intermediate","Apoe4_Y/ADNC_High")

seu$APOE4.ADNC %<>% fct_relevel(apoe4_order)
                                                   
                                                   


subclassToClass <- seu@meta.data %>% dplyr::select(Class,Subclass) %>% distinct() %>% {setNames(.$Class,nm = .$Subclass)}



apoe4.body.disease.pathology_color = MetBrewer::met.brewer(name = "Tiepolo",n = length(unique(seu$APOE4.ADNC)),direction = -1) %>% setNames(apoe4)


subclass_color <- c(Chandelier = "#A82203", Lamp5 = "#8A392C", `Lamp5 Lhx6` = "#6C5055",
                    Pax6 = "#4F677E", Pvalb = "#317EA7", Sncg = "#328FB4", Sst = "#5F9697",
                    `Sst Chodl` = "#8D9E7A", Vip = "#BAA55C", `L2/3 IT` = "#E7AD3F",
                    `L4 IT` = "#EBA03D", `L5 ET` = "#E38F41", `L5 IT` = "#DC7D46",
                    `L5/6 NP` = "#D46C4A", `L6 CT` = "#CA5F4C", `L6 IT` = "#B26546",
                    `L6 IT Car3` = "#9B6B40", L6b = "#837239", Astrocyte = "#6C7833",
                    Endothelial = "#567238", `Microglia-PVM` = "#406443", OPC = "#2B554F",
                    Oligodendrocyte = "#15475B", VLMC = "#003967")


{ggplot(seu@meta.data %>% group_by(APOE4.ADNC,Class,Specimen.ID) %>% summarise(count = n()) %>% group_by(Specimen.ID) %>% mutate(per = count/sum(count)),mapping = aes(APOE4.ADNC,y = per)) + geom_boxplot(aes(fill = Class)) + facet_wrap(~Class,scales = "free_y") + ggsignif::geom_signif(comparisons = map(apoe4_order[-c(1)],~c("Apoe4_N/ADNC_N",.)),step_increase = 0.05,margin_top = -0.3,map_signif_level = T) + cowplot::theme_cowplot()+ theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))  } %T>% ggsave(filename = "dlpfc_global_class_ADNC.APOE4_ratio.pdf",height = 10,width = 15,dpi = 300)



{
  ggplot(
    seu@meta.data %>% group_by(APOE4.ADNC, Subclass, Specimen.ID) %>% summarise(count = n()) %>% group_by(Specimen.ID) %>% mutate(per = count /sum(count))  %>% filter(Subclass == "L4 IT"),
    mapping = aes(APOE4.ADNC,y = per)
  ) + geom_boxplot(aes(fill = Subclass)) + ggsignif::geom_signif(
    comparisons = map(apoe4_order[-c(1)],  ~ c("Apoe4_N/ADNC_N", .)),
    step_increase = 0.05,
    margin_top = -0.3,
    map_signif_level = T
  ) + cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(
      angle = 270,
      hjust = 0,
      vjust = 1
    ))
} %T>% ggsave(filename = "dlpfc_class_L4_IT_ADNC.APOE4_ratio.pdf",height = 15,width = 10,dpi = 300)


{
  ggplot(
    seu@meta.data %>% group_by(APOE4.ADNC, Subclass, Specimen.ID) %>% summarise(count = n()) %>% group_by(Specimen.ID) %>% mutate(per = count /sum(count))  %>% filter(str_detect(Subclass,"Micro|Astro")),
    mapping = aes(factor(
      APOE4.ADNC, levels = apoe4_order
    ), y = per)
  ) + geom_boxplot(aes(fill = Subclass)) + ggsignif::geom_signif(
    comparisons = map(apoe4_order[-c(1)],  ~ c("Apoe4_N/ADNC_N", .)),
    step_increase = 0.05,
    margin_top = -0.3,
    map_signif_level = T
  ) + cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(
      angle = 270,
      hjust = 0,
      vjust = 1
    ))
} %T>% ggsave(filename = "dlpfc_class_Micro_Astro_ADNC.APOE4_ratio.pdf",height = 15,width = 10,dpi = 300)

seu_list <- seu %>% SplitObject("Subclass")
x = seu_list[["L6 CT"]]

degTable <- parallel::mclapply(seu_list,mc.cores = length(seu_list),function(x) {
  
  res <- parallel::mclapply(apoe4_order[-c(1)],mc.cores = length(apoe4_order)-1,function(l){
    res <- Seurat::FindMarkers(x,group.by = "APOE4.ADNC",ident.1 = l,ident.2 = NULL)
    res$cluster <- l
    res$gene <- rownames(res)
    res$subclass <- x$Subclass[[1]]
    res
  }) %>% purrr::reduce(rbind)
}
) %>% purrr::reduce(rbind)

degTable %>% write.csv("dlpfc_apoe4_deg.csv")


{
  ggplot(mapping = aes(
    x =  factor(cluster, levels = apoe4_order),
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
    filename = "dlpfc_ADNC.APOE4_num_of_deg.pdf",
    height = 10,
    width = 13,
    dpi = 300
  )



{
  ggplot(mapping = aes(
    x =  factor(cluster, levels = apoe4_order),
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
    filename = "dlpfc_ADNC.APOE4_L4_IT_num_of_deg.pdf",
    height = 7,
    width = 10,
    dpi = 300
  )



{
  ggplot(mapping = aes(
    x =  factor(cluster, levels = apoe4_order),
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
    filename = "dlpfc_ADNC.APOE4_Micro_Astro_num_of_deg.pdf",
    height = 10,
    width = 13,
    dpi = 300
  )



diseaseRNA <- read.delim("./brain_disease") %>% filter(str_detect(Disease,"Alzheimer"), Type == "mRNA") %$% Gene.symbol %>% unique
length(diseaseRNA)

{
  ggplot(mapping = aes(
    x =  factor(cluster, levels = apoe4_order),
    y = Freq
  )) + geom_bar(
    data = degTable  %>% mutate(Up = subclassToClass[subclass]) %>% filter(p_val_adj < 0.2, avg_log2FC > 0, gene %in% diseasRNA)  %>% dplyr::select(Up, cluster) %>% table %>% as.data.frame() %>% filter(str_detect(cluster, "ADNC_N", negate = T)) %>%  mutate(Freq = Freq/length(diseaseRNA)),
    mapping = aes(fill = Up,group = Up),
    stat = "identity",
    position = "dodge"
  ) + MetBrewer::scale_fill_met_d(name = "Navajo") +
    ggnewscale::new_scale_fill() +
    geom_bar(
    data = degTable  %>% mutate(Down = subclassToClass[subclass]) %>% filter(p_val_adj < 0.2, avg_log2FC <0, gene %in% diseasRNA)  %>% dplyr::select(Down, cluster) %>% table %>% as.data.frame() %>% filter(str_detect(cluster, "ADNC_N", negate = T)) %>% mutate(Freq = -Freq/length(diseaseRNA)),
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
    filename = "dlpfc_ADNC.APOE4_AD_RiskGenePer.pdf",
    height = 10,
    width = 13,
    dpi = 300
  )



{
  ggplot(mapping = aes(
    x =  factor(cluster, levels = apoe4_order),
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
    filename = "dlpfc_ADNC.APOE4_L4_IT_num_of_deg.pdf",
    height = 7,
    width = 10,
    dpi = 300
  )



{
  ggplot(mapping = aes(
    x =  factor(cluster, levels = apoe4_order),
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
    filename = "dlpfc_ADNC.APOE4_Micro_Astro_num_of_deg.pdf",
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


data <- meta %>% group_by(ADNC,APOE4.status,Class,Specimen.ID) %>% summarise(count = n()) %>% group_by(Specimen.ID) %>% mutate(per = count/sum(count))

data$ADNC %>% unique %>% sort %>% rev %>% dput
data$ADNC %<>% factor(levels = c("ADNC_N", "ADNC_Low", "ADNC_Intermediate", "ADNC_High"))

data$APOE4.status  %<>% factor(levels = c("Apoe4_N", "Apoe4_Y"))


model <- lm(per ~ ADNC * APOE4.status * Class, data=data)

summary_model <- summary(model)
coefficients <- summary_model$coefficients
r_squared <- summary_model$r.squared



plot <- ggplot(data, aes(x=interaction(ADNC,APOE4.status), y=per, group = Class,color = Class)) +
  geom_point() +  
  
  geom_smooth(method = "loess", se=FALSE) +  
  labs(title='ADNC', x='ADNC * APOE', y='percentage', color='Category') +
  theme_minimal()+ labs(subtitle = str_c("R squared ",round(r_squared,digits = 2))) + theme(axis.text.x = element_text(angle = 45,hjust = 1))

plot %>% ggsave(filename = "dlpfc_adnc_regression_apoe_global.png",height = 5,width = 7)



ct = unique(data$Class)[[1]]
for (ct in unique(data$Class)) {
  tmp = data %>% filter(Class == ct)

  
  model <- lm(per ~ ADNC * APOE4.status, data=tmp)

  summary_model <- summary(model)
  coefficients <- summary_model$coefficients
  r_squared <- summary_model$r.squared
  

  plot <- ggplot(tmp, aes(x=interaction(ADNC,APOE4.status), y=per, group = Class,color = Class)) +
    geom_point() +  
    
    geom_smooth(method = "loess", se=FALSE) +  
    labs(title='ADNC', x='ADNC * APOE', y='percentage', color='Category') +
    theme_minimal()+ labs(subtitle = str_c("R squared ",round(r_squared,digits = 2))) + theme(axis.text.x = element_text(angle = 45,hjust = 1))

  plot %>% ggsave(filename = str_c("dlpfc_adnc_regression_apoe ",ct,"_.png"),height = 5,width = 7)

}



