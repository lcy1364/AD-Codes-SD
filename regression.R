setwd("~/DATA/BRAIN/AD")

meta <- read.csv("./dlpfc_meta.csv",row.names = 1)


meta$ADNC %<>% fct_collapse(N = c("Reference","Not AD"))
meta$ADNC <- str_c("ADNC_",meta$ADNC)
meta$ADNC %<>% factor(levels = c("ADNC_N", "ADNC_Low", "ADNC_Intermediate", "ADNC_High"))


meta$APOE4.status[meta$APOE4.status == "Reference" | meta$APOE4.status == "N"] = "Apoe4_N"
meta$APOE4.status[meta$APOE4.status == "Y"] = "Apoe4_Y"
meta$APOE4.status  %<>% factor(levels = c("Apoe4_N", "Apoe4_Y"))


meta$Lewy.body.disease.pathology %<>% fct_collapse(N = c("Reference","Not Identified (olfactory bulb assessed)","Not Identified (olfactory bulb not assessed)"))
meta$Lewy.body.disease.pathology <- str_c("Lewy_",meta$Lewy.body.disease.pathology)


meta <- meta %>% unite(col = "ADNC.APOE4.LEWY",,ADNC,APOE4.status,Lewy.body.disease.pathology,remove = F,sep = "/")
meta$ADNC.APOE4.LEWY %<>% factor(levels = c(
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

meta %<>% filter(!is.na(ADNC.APOE4.LEWY))



data <- meta %>% group_by(ADNC,Class,Specimen.ID) %>% summarise(count = n()) %>% group_by(Specimen.ID) %>% mutate(per = count/sum(count))


model <- lm(per ~ ADNC * Class, data=data)

summary_model <- summary(model)
coefficients <- summary_model$coefficients
r_squared <- summary_model$r.squared



plot <- ggplot(data, aes(x = ADNC, y = per, group = Class,color = Class)) +
  geom_point() +  
  
  geom_smooth(method='loess', se=FALSE) +  
  labs(title='ADNC', x='ADNC', y='percentage', color='Category') +
  theme_minimal()+ labs(subtitle = str_c("R squared ",round(r_squared,digits = 2))) +
  theme(axis.text.x = element_text(angle = -45,hjust = 0))
plot %>% ggsave(filename = "dlpfc_adnc_regression_global.png",height = 5,width = 7)



ct = unique(data$Class)[[1]]
for (ct in unique(data$Class)) {
  tmp = data %>% filter(Class == ct)

  
  model <- lm(per ~  ADNC, data=tmp)
  summary_model <- summary(model)
  coefficients <- summary_model$coefficients
  r_squared <- summary_model$r.squared
  

  plot1 <- ggplot(tmp, aes(x=ADNC, y=per, group = Class,color = Class)) +
    geom_point() +  
    
    geom_smooth(method='lm', formula = y ~ x, se=FALSE) +  
    labs(title='ADNC', x='ADNC', y='percentage', color='Category') + facet_wrap(~Class) +
    theme_minimal()+ labs(subtitle = str_c("R squared ",round(r_squared,digits = 2))) +
    theme(axis.text.x = element_text(angle = -45,hjust = 0))

  plot1 %>% ggsave(filename = str_c("dlpfc_adnc_regression ",ct,"_.png"),height = 5,width = 7)

}


data <- meta  %>% group_by(ADNC,Subclass,Specimen.ID) %>% summarise(count = n()) %>% group_by(Specimen.ID) %>% mutate(per = count/sum(count)) %>%  filter(str_detect(Subclass,"Micro|Astro"))

ct = unique(data$Subclass)[[1]]
for (ct in unique(data$Subclass)) {
  tmp = data %>% filter(Subclass == ct)

  
  model <- lm(per ~  ADNC, data=tmp)
  summary_model <- summary(model)
  coefficients <- summary_model$coefficients
  r_squared <- summary_model$r.squared
  

  plot1 <- ggplot(tmp, aes(x=ADNC, y=per, group = Subclass,color = Subclass)) +
    geom_point() +  
    
    geom_smooth(method='lm', formula = y ~ x, se=FALSE) +  
    labs(title='ADNC', x='ADNC', y='percentage', color='Category') +
    theme_minimal()+ labs(subtitle = str_c("R squared ",round(r_squared,digits = 2))) +
    theme(axis.text.x = element_text(angle = -45,hjust = 0))

  plot1 %>% ggsave(filename = str_c("dlpfc_adnc_regression ",ct,"_.png"),height = 5,width = 7)

}


data <- meta %>% group_by(ADNC,APOE4.status,Class,Specimen.ID) %>% summarise(count = n()) %>% group_by(Specimen.ID) %>% mutate(per = count/sum(count))


model <- lm(per ~ ADNC * APOE4.status * Class, data=data)

summary_model <- summary(model)
coefficients <- summary_model$coefficients
r_squared <- summary_model$r.squared



plot <- ggplot(data, aes(x=interaction(ADNC,APOE4.status), y=per, group = Class,color = Class)) +
  geom_point() +  
  
  geom_smooth(method = "loess", se=FALSE) +  
  labs(title='ADNC * APOE', x='ADNC * APOE', y='percentage', color='Category') +
  theme_minimal()+ labs(subtitle = str_c("R squared ",round(r_squared,digits = 2))) + theme(axis.text.x = element_text(angle = -45,hjust = 0))

plot %>% ggsave(filename = "dlpfc_adnc_apoe_regression_global.png",height = 5,width = 7)



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
    labs(title='ADNC * APOE', x='ADNC * APOE', y='percentage', color='Category') +
    theme_minimal()+ labs(subtitle = str_c("R squared ",round(r_squared,digits = 2))) + theme(axis.text.x = element_text(angle = -45,hjust = 0))

  plot %>% ggsave(filename = str_c("dlpfc_adnc_apoe_regression ",ct,"_.png"),height = 5,width = 7)
}



data <- meta  %>% group_by(ADNC,APOE4.status,Subclass,Specimen.ID) %>% summarise(count = n()) %>% group_by(Specimen.ID) %>% mutate(per = count/sum(count)) %>%  filter(str_detect(Subclass,"Micro|Astro"))

ct = unique(data$Subclass)[[1]]
for (ct in unique(data$Subclass)) {
  tmp = data %>% filter(Subclass == ct)

  
  model <- lm(per ~  ADNC * APOE4.status, data=tmp)
  summary_model <- summary(model)
  coefficients <- summary_model$coefficients
  r_squared <- summary_model$r.squared
  

  plot1 <- ggplot(tmp, aes(x=interaction(ADNC,APOE4.status), y=per, group = Subclass,color = Subclass)) +
    geom_point() +  
    
    geom_smooth(method='lm', formula = y ~ x, se=FALSE) +  
    labs(title='ADNC * APOE', x='ADNC * APOE', y='percentage', color='Category') +
    theme_minimal()+ labs(subtitle = str_c("R squared ",round(r_squared,digits = 2))) +
    theme(axis.text.x = element_text(angle = -45,hjust = 0))

  plot1 %>% ggsave(filename = str_c("dlpfc_adnc_apoe_regression ",ct,"_.png"),height = 5,width = 7)

}


data <- meta %>% filter(!is.na(ADNC.APOE4.LEWY)) %>% group_by(ADNC,APOE4.status,Lewy.body.disease.pathology,ADNC.APOE4.LEWY,Class,Specimen.ID) %>% summarise(count = n()) %>% group_by(Specimen.ID) %>% mutate(per = count/sum(count))



model <- lm(per ~ ADNC * APOE4.status *Lewy.body.disease.pathology* Class, data=data)

summary_model <- summary(model)
coefficients <- summary_model$coefficients
r_squared <- summary_model$r.squared



plot <- ggplot(data, aes(x=ADNC.APOE4.LEWY, y=per, group = Class,color = Class)) +
  geom_point() +  
  
  geom_smooth(method = "loess", se=FALSE) +  
  labs(title='ADNC * APOE * Lewy', x='ADNC * APOE * Lewy', y='percentage', color='Category') +
  theme_minimal()+ labs(subtitle = str_c("R squared ",round(r_squared,digits = 2))) + theme(axis.text.x = element_text(angle = -45,hjust = 0))

plot %>% ggsave(filename = "dlpfc_adnc_apoe_lewy_regression_global.png",height = 5,width = 7)



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
    labs(title='ADNC * APOE * Lewy', x='ADNC * APOE * Lewy', y='percentage', color='Category') +
    theme_minimal()+ labs(subtitle = str_c("R squared ",round(r_squared,digits = 2))) + theme(axis.text.x = element_text(angle = -45,hjust = 0))

  plot %>% ggsave(filename = str_c("dlpfc_adnc_apoe_lewy_regression_adnc_apoe_lewy ",ct,"_.png"),height = 5,width = 7)

}



data <- meta %>% filter(!is.na(ADNC.APOE4.LEWY)) %>% group_by(ADNC,APOE4.status,Lewy.body.disease.pathology,ADNC.APOE4.LEWY,Subclass,Specimen.ID) %>% summarise(count = n()) %>% group_by(Specimen.ID) %>% mutate(per = count/sum(count)) %>%  filter(str_detect(Subclass,"Micro|Astro"))

ct = unique(data$Subclass)[[1]]
for (ct in unique(data$Subclass)) {
  tmp = data %>% filter(Subclass == ct)

  
  model <- lm(per ~  ADNC * APOE4.status * Lewy.body.disease.pathology, data=tmp)
  summary_model <- summary(model)
  coefficients <- summary_model$coefficients
  r_squared <- summary_model$r.squared
  

  plot1 <- ggplot(tmp, aes(x=ADNC.APOE4.LEWY, y=per, group = Subclass,color = Subclass)) +
    geom_point() +  
    
    geom_smooth(method='lm', formula = y ~ x, se=FALSE) +  
    labs(title='ADNC * APOE * Lewy', x='ADNC * APOE * Lewy', y='percentage', color='Category') +
    theme_minimal()+ labs(subtitle = str_c("R squared ",round(r_squared,digits = 2))) +
    theme(axis.text.x = element_text(angle = -45,hjust = 0)) 

  plot1 %>% ggsave(filename = str_c("dlpfc_adnc_apoe_lewy_regression_adnc_apoe_lewy ",ct,"_.pdf"),height = 5,width = 8)

}

