library(tidyverse)
setwd("~/DATA/BRAIN/AD")
options(bitmapType = "cairo")
degFiles <- list.files(".","deg.csv")

region = degFiles[[1]] %>% str_extract("^[a-z]{1,}")
degFiles <- list.files(".","deg.csv")

region = degFiles[[1]] %>% str_extract("^[a-z]{1,}")

diseaseGenes <- read.delim("./Brain disease-gene associations.txt") %>% filter(str_detect(Disease,"Alzheimer")) %$% Gene.symbol %>% unique

ct = "Astrocyte"
for (ct in c("Astrocyte","Microglia-PVM")) {

  
  file1 <- degFiles %>% keep(~str_detect(.,"adnc_deg.csv")) %>% read.csv() %>% filter(subclass %in% ct, p_val_adj < 0.05)
  file1$comparision <- "ADNC"

  file2 <- degFiles %>% keep(~str_detect(.,"ADNC.APOE4.LEWY_deg.csv")) %>% read.csv() %>% filter(subclass %in% ct, p_val_adj < 0.05)
  file2$comparision <- "ADNC.APOE4.LEWY"

  file3 <- degFiles %>% keep(~str_detect(.,"apoe4_deg.csv")) %>% read.csv() %>% filter(subclass %in% ct, p_val_adj < 0.05)
  file3$comparision <- "APOE4"

  data <- list(ADNC = file1$gene %>% unique, APOE4 = file3$gene %>% unique, ADNC.APOE4.LEWY = file3$gene %>% unique, AD = diseaseGenes)


  
  library(ggvenn)
  ggsave(str_c(region,"_",ct,"_","AD_genes_Venn.pdf"),ggvenn(data, c("ADNC", "APOE4","ADNC.APOE4.LEWY","AD")),height = 6,width = 6)

  library(UpSetR)
  
  cairo_pdf(str_c(region,"_",ct,"_","AD_genes_Upset.pdf"),height = 6,width = 6)
  print(upset(fromList(data),order.by = "freq"))
  dev.off()

  sharedGenes <- data %>% purrr::reduce(intersect)

  
  
}
