set_wd()
libztw()

library("survminer")
library("survival")
library(broom)


# load data ---------------------------------------------------------------
filelists <- list.files(pattern = "sherlock.*RData",)
filelists <- filelists[!str_detect(filelists,'landscape')]
lapply(filelists,load,.GlobalEnv)


load('sherlock_landscape_v2.RData')
load('HRD_LOH.RData')
load('quiet_tumor.RData')
data_hrd <- HRD_LOH %>% mutate(Type="HRD_LOH")

data_quiet <- 
  sherlock_samples_unique %>% 
  select(Subject,Tumor_Barcode) %>% 
  left_join(quiet_tumor) %>% 
  filter(group=="Quiet tumors (N=33)") %>% 
  mutate(Gene="Quiet Tumor",Alteration="",Type="@") %>% 
  select(colnames(data_hrd))



data_rtk_p53 <- 
  sherlock_samples_unique %>% 
  select(Subject,Tumor_Barcode) %>%
  left_join(sherlock_rtk) %>% 
  left_join(sherlock_p53 %>% select(Tumor_Barcode,P53=MDM2_TP53)) %>% 
  mutate(Group=paste0(RTK_Altered_Status,"\nP53 ",P53)) %>% mutate(Group=factor(Group,levels = c("RTK-RAS-\nP53 N", "RTK-RAS+\nP53 N", "RTK-RAS-\nP53 Y", "RTK-RAS+\nP53 Y")))


# load sherlock mt content
load('sherlock_mtcontent.RData')


data_all <- bind_rows(
  data_top %>% mutate(Alteration=as.character(Alteration))%>% as_tibble(),
  data_sftb%>% as_tibble(),
  data_scna%>% as_tibble(),
  data_arm%>% as_tibble(),
  data_fusion%>% as_tibble(),
  data_feature1 %>% mutate(Type="@") %>% filter(!(Gene %in% c( "SCNA_Cluster" ,"MSI")))%>% as_tibble(),
  data_hrd%>% as_tibble(),
  data_hrd %>% mutate(Gene="HRD") %>% unique()%>% as_tibble(),
  data_quiet%>% as_tibble(),
  data_arm %>% filter(Gene=="Chr19q"|Gene=="Chr19p") %>% mutate(Gene="19p/19q") %>% unique()%>% as_tibble(),
  sherlock_mtcontent %>% arrange(desc(MT_Ratio)) %>% slice(1:116) %>% select(Subject,Tumor_Barcode) %>% mutate(Gene="mtDNA-high",Alteration="mtDNA-high",Type="@") %>% as_tibble()
)%>% mutate(Key=paste0(Gene,"|",Type)) 

#'ERBB2','MET','KRAS','ALK','EGFR'
# data_all2 <- bind_rows(
#   data_top %>% mutate(Alteration=as.character(Alteration))  %>% as_tibble(),
#   data_sftb%>% as_tibble(),
#   data_scna%>% as_tibble(),
#   data_arm%>% as_tibble(),
#   data_fusion %>%  as_tibble(),
#   data_feature1 %>% mutate(Type="@") %>% filter(!(Gene %in% c( "SCNA_Cluster" ,"MSI","Kataegis","WGD","Histology","p53_deficiency", "RTK-RAS+")))%>% as_tibble(),
#   data_hrd%>% as_tibble()
# )%>% mutate(Key=paste0(Gene,"|",Type))

data_all2 <- bind_rows(
  data_top %>% mutate(Alteration=as.character(Alteration)) %>% filter(Gene!="TP53") %>% as_tibble(),
  data_sftb%>% as_tibble(),
  data_scna%>% filter(Gene!="MDM2") %>% as_tibble(),
  data_arm%>% as_tibble(),
  data_fusion %>%  as_tibble(),
  data_feature1 %>% mutate(Type="@") %>% filter(!(Gene %in% c( "SCNA_Cluster" ,"MSI","Kataegis","WGD","Histology", "RTK-RAS+")))%>% as_tibble(),
  data_hrd%>% as_tibble()
)%>% mutate(Key=paste0(Gene,"|",Type))

tmp <- data_all2 %>% select(Subject,Key) %>% unique() %>% count(Key) %>% filter(n>0.05*232) %>% pull(Key)

data_all2 <- data_all2 %>% filter(Key %in% tmp)

## remove chr19 HD
# data_all %>% filter(Gene=="19p/19q")
# chr19hd <- data_all %>% filter(Gene=="Chr19_HD") %>% pull(Tumor_Barcode)
# 
# data_all <- bind_rows(
#   data_all %>% filter(Gene!="19p/19q"),
#   data_all %>% filter(Gene=="19p/19q") %>% filter(!Tumor_Barcode %in% chr19hd)
# )


suvdata <- sherlock_samples_unique %>% 
  select(Subject,Tumor_Barcode,Study) %>% 
  left_join(sherlock_clinical) %>% 
  mutate(passive_smoking=factor(passive_smoking,levels = c("N","Y"))) %>% 
  rename(age=age_at_diagnosis) %>% 
  mutate(Study=if_else(Study=="Laval - Quebec",Study, 'Other study'))

# mutate(stage = fct_collapse(stage,
#                             IA="IA",
#                             IB="IB",
#                             II="II",
#                             III=c('III','IIIA','IIIB'))
# ) %>% 


suvdata %>% select(Subject:gender)%>% filter_all(any_vars(is.na(.)))
#sherlock_suvdata <- suvdata
#save(sherlock_suvdata,file='sherlock_suvdata.RData')
load('sherlock_histology.RData')
sherlock_histology <- sherlock_histology %>% mutate(Histology=if_else(Histology=="Carcinoids","Others",Histology)) %>% mutate(Histology=factor(Histology,levels = c("Adenocarcinomas", "Others")))
tmp <- suvdata %>% left_join(sherlock_histology)
# blank survival ----------------------------------------------------------
fit <- coxph(Surv(survival_months, death) ~ age+gender+stage+passive_smoking+grade+Study+Histology, data = tmp) ### overall 
ggforest(fit)
ggsave('Figures/blank_suv.pdf',width = 8,height = 6)


## calculated all p-value # 

#suvdata <- suvdata %>% select(Subject,Tumor_Barcode,Study,gender,age,stage,grade,passive_smoking,death,survival_months,stage_info)

allvariates <-  unique(data_all2$Key) 
resultall <- tibble(Key=character(),p=numeric(),lab=character())
for(key in allvariates){
  tmp <- data_all2 %>%
    mutate(Key=paste0(Gene,"|",Type)) %>% 
    filter(Key==key) %>% 
    select(Tumor_Barcode) %>% 
    mutate(Key="Y")
  
  suvdata_tmp <- suvdata %>% 
    left_join(tmp) %>%
    mutate(Key=if_else(is.na(Key),"N",Key)) %>% 
    mutate(Key=factor(Key,levels = c("N","Y")))
  
  
  resdata <- SurvZTW(suvdata = suvdata_tmp,plot = FALSE)
  resultall <-  tibble(Key=key,p=as.numeric(resdata[1]),lab=resdata[2]) %>% 
    bind_rows(resultall)
  #SurvZTW(suvdata = suvdata_tmp,plot = TRUE,keyname="RTK")
}

resultall <- resultall %>% mutate(FDR=p.adjust(p,method = 'BH')) %>% arrange(p) 
resultall %>% View()

resultall %>% filter(str_detect(Key,"SCNA_Arm")) %>% mutate(FDR=p.adjust(p,method = 'BH')) %>% arrange(p)
resultall %>% filter(str_detect(Key,"LOH")) %>% mutate(FDR=p.adjust(p,method = 'BH')) %>% arrange(p)


save(resultall,file='resultall.RData',version = 2)

tmp <- data_all %>% mutate(Key=paste0(Gene,"|",Type)) %>% select(Subject,Key) %>% mutate(value=1) %>% unique() %>% pivot_wider(names_from = Key,values_from=value,values_fill=0 )

suvdata %>% select(Subject,age,gender,stage,death,survival_months)

ff <- resultall$Key[1:5]
data_all3 <- data_all2 %>% filter(Key %in% ff) %>% select(Subject,Key) %>% unique() %>% mutate(Score=1) %>% group_by(Subject) %>% summarise(Key=sum(Score)) %>% ungroup() %>% 
  select(Subject,Key)
tmp <- left_join(suvdata,data_all3) %>% mutate(Key=if_else(is.na(Key),0,Key))
fit <- coxph(Surv(survival_months, death) ~ Key+age+gender+stage, data = tmp) ### overall 
broom::tidy(fit,conf.int=T) %>% filter(term=="Key") %>% mutate(p.value=scientific(p.value))
SurvZTWm(suvdata = (tmp %>% mutate(Key=if_else(Key>3,3,Key))%>%  mutate(Key=as.factor(Key))),plot = TRUE,keyname="Score",filename = 'Figures/Survival/Risk_survival.pdf')


# plot signficant group # 

dataprep <- function(key){
  tmp <- data_all %>%
    mutate(Key=paste0(Gene,"|",Type)) %>% 
    filter(Key==key) %>% 
    select(Tumor_Barcode) %>% 
    mutate(Key="Y")
  
  suvdata_tmp <- suvdata %>% 
    left_join(tmp) %>%
    mutate(Key=if_else(is.na(Key),"N",Key)) %>% 
    mutate(Key=factor(Key,levels = c("N","Y")))
  
  return(suvdata_tmp)
}


SurvZTW(suvdata = dataprep('p53_deficiency|@'),plot = TRUE,keyname="p53_deficiency",filename = 'Figures/Survival/p53_deficiency_survival.pdf')
SurvZTW(suvdata = dataprep('CHEK2|HRD_LOH'),plot = TRUE,keyname="CHEK2 LOH",filename = 'Figures/Survival/CHEK2_LOH_survival.pdf')
SurvZTW(suvdata = dataprep('ATM|HRD_LOH'),plot = TRUE,keyname="ATM LOH",filename = 'Figures/Survival/ATM_LOH_survival.pdf')
SurvZTW(suvdata = dataprep('Chr22q|SCNA_Arm'),plot = TRUE,keyname="Chr22q Loss",filename = 'Figures/Survival/Chr22q_loss_survival.pdf')
SurvZTW(suvdata = dataprep('Chr15q|SCNA_Arm'),plot = TRUE,keyname="Chr15q Loss",filename = 'Figures/Survival/Chr15q_loss_survival.pdf')
SurvZTW(suvdata = dataprep('Chr19p|SCNA_Arm'),plot = TRUE,keyname="Chr19p Loss",filename = 'Figures/Survival/Chr19p_loss_survival.pdf')
SurvZTW(suvdata = dataprep('Chr19q|SCNA_Arm'),plot = TRUE,keyname="Chr19q Loss",filename = 'Figures/Survival/Chr19q_loss_survival.pdf')

SurvZTW(suvdata = dataprep('19p/19q|SCNA_Arm'),plot = TRUE,keyname="Chr19p/q Loss",filename = 'Figures/Survival/19pq_loss_survival.pdf')
SurvZTW(suvdata = dataprep('SETD2|Mutation_Driver'),plot = TRUE,keyname="SETD2 Mutation",filename = 'Figures/Survival/SETD2_Mutation_Driver_survival.pdf')
SurvZTW(suvdata = dataprep('HLA_LOH|@'),plot = TRUE,keyname="HLA LOH",filename = 'Figures/Survival/HLA_LOH_survival.pdf')
SurvZTW(suvdata = dataprep('Chr19_HD|@'),plot = TRUE,keyname="Chr19 HD ",filename = 'Figures/Survival/Chr19_HD_survival.pdf')
SurvZTW(suvdata = dataprep('Quiet Tumor|@'),plot = TRUE,keyname="Quiet Tumor",filename = 'Figures/Survival/Quiet_Tumor_survival.pdf')
SurvZTW(suvdata = dataprep('WGD|@'),plot = TRUE,keyname="WGD ",filename = 'Figures/Survival/WGD_survival.pdf')


SurvZTW(suvdata = dataprep('RTK-RAS+|@'),plot = TRUE,keyname="RTK-RAS",filename = 'Figures/Survival/RTK_survival.pdf')


SurvZTW(suvdata = dataprep('mtDNA-high|@'),plot = TRUE,keyname="mtDNA-high",filename = 'Figures/Survival/mtDNA_survival.pdf')

SurvZTW(suvdata = dataprep('CTNNB1|Mutation_Driver'),plot = TRUE,keyname="CTNNB1",filename = 'Figures/Survival/CTNNB1_survival.pdf')

SurvZTW(suvdata = dataprep('Risk Score|High'),plot = TRUE,keyname="Risk_High",filename = 'Figures/Survival/Risk_score_survival.pdf')
save(data_all,suvdata,file="sherlock_survival.RData")


## carcinoid survival
tmp <- suvdata %>% left_join(sherlock_histology) %>% filter(Histology!="Others") %>% mutate(Key=if_else(Histology=="Adenocarcinomas","N","Y")) %>% mutate(Key=factor(Key,levels = c("N","Y")))
SurvZTW(suvdata = tmp,plot = TRUE,keyname="Carcinoids",filename = 'Figures/Survival/Carcinoids_survival.pdf')



# P53 deficiency ----------------------------------------------------------
suvdata_tmp <- sherlock_p53 %>% 
  mutate(Gene=if_else(is.na(Gene),"WT",Gene)) %>% 
  mutate(Gene=factor(Gene,levels = c('WT','MDM2','TP53'))) %>% 
  select(Tumor_Barcode,Key=Gene) %>% 
  right_join(suvdata)

SurvZTWm(suvdata = suvdata_tmp,plot = TRUE,keyname="P53 deficiency",pvalsize = 3,filename = 'Figures/Survival/P53_deficiency_group_survival.pdf')

## RTK group ##
tmp <- sherlock_rtk_info %>% 
  filter(Hugo_Symbol %in% c('EGFR','KRAS','ALK',"ERBB2","MET"),Variant_Type!="CNV") %>% 
  select(Subject,Hugo_Symbol) %>% 
  unique() 
#  group_by(Subject) %>% 
#  summarise(Hugo_Symbol=paste0(Hugo_Symbol,collapse = ',')) %>%
# filter(!str_detect(Hugo_Symbol,","))
# tmp$Hugo_Symbol[tmp$Subject=="NSLC-0021"] <- "ALK"
# tmp$Hugo_Symbol[tmp$Subject=="NSLC-0186"] <- "EGFR"
# tmp$Hugo_Symbol[tmp$Subject=="NSLC-0156"] <- "EGFR"
# tmp$Hugo_Symbol[tmp$Subject=="NSLC-0005"] <- "EGFR"
# tmp <- sherlock_rtk_info %>% 
#   filter(Variant_Type == "Fusion") %>% 
#   select(Subject,Hugo_Symbol) %>% 
#   unique() %>% 
#   mutate(Hugo_Symbol="Fusion") %>% 
#   bind_rows(tmp)

suvdata_tmp <- suvdata %>% 
  left_join(tmp) %>% 
  mutate(Key=if_else(is.na(Hugo_Symbol),"WT",Hugo_Symbol)) %>% 
  mutate(Key=factor(Key,levels = c('WT','ERBB2','MET','ALK','KRAS','EGFR')))

SurvZTWm(suvdata = suvdata_tmp,plot = TRUE,keyname="RTK-RAS",pvalsize = 2.5,filename = 'Figures/Survival/rtk_group_survival.pdf')


## combined RTK and TP53 ## 
tmp <- sherlock_samples_unique %>% 
  select(Tumor_Barcode) %>% 
  left_join(sherlock_rtk) %>% 
  left_join(sherlock_p53 %>% select(Tumor_Barcode,P53='MDM2_TP53')) %>% 
  mutate(Key=paste0(RTK_Altered_Status," & P53 ",P53)) %>% mutate(Key=factor(Key,levels = c("RTK-RAS- & P53 N", "RTK-RAS+ & P53 N", "RTK-RAS- & P53 Y", "RTK-RAS+ & P53 Y")))

suvdata_tmp <- suvdata %>% left_join(tmp)
SurvZTWm(suvdata = suvdata_tmp,plot = TRUE,keyname="RTK-RAS & P53",pvalsize = 2.5,filename = 'Figures/Survival/rtk_p53_group_survival.pdf')



# CNV clsuter 
load('cnvclust.RData')
tmp <- sherlock_samples_unique %>% 
  select(Tumor_Barcode) %>% 
  left_join(cnvclust) %>% 
  mutate(Key=CNV_Clust) %>% mutate(Key=factor(Key,levels = c("C3","C2","C1")))

suvdata_tmp <- suvdata %>% left_join(tmp)
SurvZTWm(suvdata = suvdata_tmp,plot = TRUE,keyname="CNV_Clust",pvalsize = 2.5,filename = 'Figures/Survival/cnv_cluster_survival.pdf')

load('cnvclust.RData')
tmp <- sherlock_samples_unique %>% 
  select(Tumor_Barcode) %>% 
  left_join(cnvclust) %>% 
  mutate(Key=if_else(CNV_Clust=="C3","Piano","Non-Piano")) %>%
  mutate(Key=factor(Key,levels = c("Non-Piano","Piano")))

suvdata_tmp <- suvdata %>% left_join(tmp)
SurvZTWm(suvdata = suvdata_tmp,plot = TRUE,keyname="CNV_Clust",pvalsize = 2.5,filename = 'Figures/Survival/cnv_cluster_survival2.pdf')

load('sherlock_histology.RData')
adtmp <- sherlock_histology %>% filter(Histology=="Adenocarcinomas") %>% pull(Tumor_Barcode)
tmp <- sherlock_samples_unique %>% 
  select(Tumor_Barcode) %>% 
  left_join(cnvclust) %>% 
  filter(Tumor_Barcode %in% adtmp) %>% 
  mutate(Key=if_else(CNV_Clust=="C3","Piano","Non-Piano")) %>%
  mutate(Key=factor(Key,levels = c("Piano","Non-Piano")))

suvdata_tmp <- suvdata %>% left_join(tmp) %>% filter(Tumor_Barcode %in% adtmp)
SurvZTWm(suvdata = suvdata_tmp,plot = TRUE,keyname="CNV_Clust",pvalsize = 2.5,filename = 'Figures/Survival/cnv_cluster_survival2_AD.pdf')







## TL length
tmp <- sherlock_samples_unique %>% 
  left_join(
    sherlock_tl %>% mutate(TLratio=Tumor_TL/Normal_TL) 
  ) %>%
  arrange(TLratio) %>% 
  mutate(TLgroup="")

tmp$TLgroup[1:77] <- 'TL-Low'
tmp$TLgroup[78:154] <- 'TL-Middle'
tmp$TLgroup[155:232] <- 'TL-High'


suvdata_tmp <- suvdata %>% 
  left_join(tmp %>% select(Tumor_Barcode,Key=TLgroup)) %>% 
  mutate(Key=factor(Key,levels = c('TL-Low','TL-Middle','TL-High')))

SurvZTWm(suvdata = suvdata_tmp,plot = TRUE,keyname="Telomere length",pvalsize = 2.5,filename = 'Figures/Survival/TL_group_survival.pdf')





# Function ----------------------------------------------------------------
SurvZTW <- function(suvdata,plot=TRUE,legend.labs=NULL,keyname="Strata",filename=NULL,width = 7,height = 5.5,keyterm="KeyY"){
  fit <- coxph(Surv(survival_months, death) ~ Key+age+gender+stage, data = suvdata) ### overall 
  suvpvalue <- broom::tidy(fit,conf.int=T) %>% filter(term==keyterm) %>% mutate(lab=paste0("P = ",format(p.value,digits = 2,scientific = FALSE),"\nHR = ",format(exp(estimate),digits = 2,scientific = FALSE)," (95% CI: ",format(exp(conf.low),digits = 2,scientific = FALSE),"-",format(exp(conf.high),digits = 2,scientific = FALSE),")")) %>% select(term,p.value,lab)
  if(plot){
    fit2 <- coxph(Surv(survival_months, death) ~ strata(Key)+age+gender+stage, data = suvdata) ## strata curve ##
    ##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
    suvfit <- survfit(fit2)
    #plot(suvfit)
    suvfit$call$formula <- as.formula("Surv(survival_months, death) ~ strata(Key)+age+gender+stage")
    
    if(is.null(legend.labs)) 
    {
      legend.labs <- c(paste0("N (",suvfit$n[1],")"),paste0("Y (",suvfit$n[2],")"))
    }
    keyname=keyname
    ggsurv <- ggsurvplot(
      suvfit,                     # survfit object with calculated statistics.
      data = suvdata,             # data used to fit survival curves.
      risk.table = FALSE,       # show risk table.
      pval = suvpvalue$lab,             # show p-value of log-rank test.
      pval.size=4,
      #pval = TRUE,
      #pval.method
      #conf.int = TRUE,         # show confidence intervals for 
      # point estimates of survival curves.
      palette = pal_jama()(2),
      #xlim = c(0,210),         # present narrower X axis, but not affect
      # survival estimates.
      xlab = "Time in months",   # customize X axis label.
      break.time.by = 50,     # break X axis in time intervals by 500.
      ggtheme = theme_ipsum_rc(axis_title_size = 14,base_size = 12,axis_title_just = "m",axis = "XY",axis_col = "black"), # customize plot and risk table with a theme.
      font.family="Roboto Condensed" ,
      risk.table.y.text.col = T,# colour risk table text annotations.
      risk.table.height = 0.25, # the height of the risk table
      risk.table.y.text = FALSE,# show bars instead of names in text annotations
      # in legend of risk table.
      censor=TRUE,
      ncensor.plot = FALSE,      # plot the number of censored subjects at time t
      ncensor.plot.height = 0.25,
      conf.int.style = "step",  # customize style of confidence intervals
      surv.median.line = "none",
      legend.labs = legend.labs
    )
    # ggsurv$plot <- ggsurv$plot + labs(
    #   title    = "Survival curve",                     
    #   subtitle = "Based on Kaplan-Meier estimates"
    # )
    ggsurv <- ggpar(
      ggsurv,
      # font.title    = c(15, "bold", "darkblue"),         
      # font.subtitle = c(14, "bold.italic", "purple"), 
      # font.caption  = c(14, "plain", "orange"),        
      # font.x        = c(14, "bold.italic", "red"),          
      # font.y        = c(14, "bold.italic", "darkred"),      
      # font.xtickslab = c(12, "plain", "darkgreen"),
      #legend = c(0.8,1.1)
      legend = "top"
    )
    ggsurv$plot <- ggsurv$plot+
      scale_x_continuous(breaks = pretty_breaks(6),expand = c(0,0))+
      scale_y_continuous(breaks = pretty_breaks(6),expand = c(0,0))+
      theme(legend.direction = "horizontal")+
      labs(color=keyname)+
      coord_cartesian(clip = 'off')
    ggsurv$plot <- flush_ticks(ggsurv$plot)
    ggsurv
    
    if(is.null(filename)){
      filename=paste0(keyname,"-Survival.pdf")
    }
    
    ggsave(filename,width = width,height = height,device = cairo_pdf)
    
    
  }
  return(c(suvpvalue$p.value,suvpvalue$lab))
}

SurvZTWm <- function(suvdata,plot=TRUE,legend.labs=NULL,keyname="Strata",filename=NULL,width = 7,height = 5.5,pvalsize=4,pvalab=NULL){
  
  fit <- coxph(Surv(survival_months, death) ~ Key+age+gender+stage, data = suvdata) ### overall 
  suvpvalue <- tidy(fit,conf.int = T) %>% filter(str_detect(term,"Key")) %>% mutate(lab=paste0("P = ",format(p.value,digits = 2,scientific = FALSE),"\nHR = ",format(exp(estimate),digits = 2,scientific = FALSE)," (95% CI: ",format(exp(conf.low),digits = 2,scientific = FALSE),"-",format(exp(conf.high),digits = 2,scientific = FALSE),")")) %>% select(term,p.value,lab) %>% 
    mutate(term=str_remove(term,"Key")) %>% mutate(lab=paste0(term,":\n",lab))
  if(is.null(pvalab)){
    pvalab <- paste0(suvpvalue$lab,collapse = "\n")
  }
  
  fit2 <- coxph(Surv(survival_months, death) ~ strata(Key)+age+gender+stage, data = suvdata) ## strata curve ##
  ##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
  suvfit <- survfit(fit2)
  #plot(suvfit)
  suvfit$call$formula <- as.formula("Surv(survival_months, death) ~ strata(Key)+age+gender+stage")
  
  if(is.null(legend.labs)) 
  {
    legend.labs <- paste0(levels(suvdata$Key)," (", suvfit$n,")")
  }
  keyname=keyname
  ggsurv <- ggsurvplot(
    suvfit,                     # survfit object with calculated statistics.
    data = suvdata,             # data used to fit survival curves.
    risk.table = FALSE,       # show risk table.
    pval = pvalab,             # show p-value of log-rank test.
    pval.size=pvalsize,
    #pval = TRUE,
    #pval.method
    #conf.int = TRUE,         # show confidence intervals for 
    # point estimates of survival curves.
    palette = pal_jama()(length(levels(suvdata$Key))),
    #xlim = c(0,210),         # present narrower X axis, but not affect
    # survival estimates.
    xlab = "Time in months",   # customize X axis label.
    break.time.by = 50,     # break X axis in time intervals by 500.
    ggtheme = theme_ipsum_rc(axis_title_size = 14,base_size = 12,axis_title_just = "m",axis = "XY",axis_col = "black"), # customize plot and risk table with a theme.
    font.family="Roboto Condensed" ,
    risk.table.y.text.col = T,# colour risk table text annotations.
    risk.table.height = 0.25, # the height of the risk table
    risk.table.y.text = FALSE,# show bars instead of names in text annotations
    # in legend of risk table.
    censor=TRUE,
    ncensor.plot = FALSE,      # plot the number of censored subjects at time t
    ncensor.plot.height = 0.25,
    conf.int.style = "step",  # customize style of confidence intervals
    surv.median.line = "none",
    legend.labs = legend.labs
  )
  # ggsurv$plot <- ggsurv$plot + labs(
  #   title    = "Survival curve",                     
  #   subtitle = "Based on Kaplan-Meier estimates"
  # )
  ggsurv <- ggpar(
    ggsurv,
    # font.title    = c(15, "bold", "darkblue"),         
    # font.subtitle = c(14, "bold.italic", "purple"), 
    # font.caption  = c(14, "plain", "orange"),        
    # font.x        = c(14, "bold.italic", "red"),          
    # font.y        = c(14, "bold.italic", "darkred"),      
    # font.xtickslab = c(12, "plain", "darkgreen"),
    #legend = c(0.8,1.1)
    legend = "top"
  )
  ggsurv$plot <- ggsurv$plot+
    scale_x_continuous(breaks = pretty_breaks(6),expand = c(0,0))+
    scale_y_continuous(breaks = pretty_breaks(6),expand = c(0,0))+
    theme(legend.direction = "horizontal")+
    labs(color=keyname)+
    coord_cartesian(clip = 'off')
  ggsurv$plot <- flush_ticks(ggsurv$plot)
  ggsurv
  
  if(is.null(filename)){
    filename=paste0(keyname,"-Survival.pdf")
  }
  
  ggsave(filename,width = width,height = height,device = cairo_pdf)
  
  return(suvpvalue)
}




















# COXNET ------------------------------------------------------------------
#rm(list=ls())
set_wd()
libztw()
load('sherlock_survival.RData')

suvdata_tmp <- suvdata %>% filter(!is.na(death),!is.na(survival_months),!is.na(age),!is.na(gender), !is.na(stage))
#data_all_tmp <- data_all %>% filter(Key %in% c("MDM2|SCNA_Focal","TP53|Mutation_Driver","ERBB2|Mutation_Driver","MET|Mutation_Driver","ALK|Fusion","KRAS|Mutation_Driver","EGFR|Mutation_Driver","CHEK2|HRD_LOH","ATM|HRD_LOH","HLA_LOH|@","Chr22q|SCNA_Arm","Chr15q|SCNA_Arm", "SETD2|Mutation_Driver"))

data_all_tmp <- data_all %>% filter(!(Key %in% c("HRD|HRD_LOH","p53_deficiency|@" ,"RTK-RAS+|@","Histology|@","Quiet Tumor|@")))
tmpkey <- data_all_tmp %>% select(Subject,Key) %>% unique() %>% count(Key) %>% arrange(n) %>% filter(n>(0.03*232)) %>% pull(Key)
data_all_tmp <- data_all_tmp %>% filter(Key %in% tmpkey) %>% unique()

tmp <- data_all_tmp %>% filter(Tumor_Barcode %in% suvdata_tmp$Tumor_Barcode) %>% 
  select(Tumor_Barcode,Key) %>% mutate(Value=1) %>% 
  unique()

tmpdata <- suvdata_tmp %>% select(Tumor_Barcode,age,gender,stage) %>% 
  mutate(gender=if_else(gender=="Female",1L,0L)) %>% 
  mutate(stage=as.integer(stage)) %>% 
  left_join(tmp) %>% unique() %>% 
  pivot_wider(names_from = Key,values_from=Value,values_fill=0) %>% select(-`NA`)

suvdata_tmp$Tumor_Barcode==tmpdata$Tumor_Barcode
patient.data <- NULL
patient.data <- list(x=as.matrix(tmpdata[,-1]),time= suvdata_tmp$survival_months, status=suvdata_tmp$death)
dataname <- colnames(tmpdata)[-1]

library("glmnet")
library("survival")
betainfo_all <- tibble(name=character(),beta=numeric(),lambda=numeric(),Seq=integer(),p.value=numeric(),HR=numeric(),HR_low=numeric(), HR_high=numeric())
for(i in 1:1000){
  #print(i)
  betainfo <-NULL
  cv.fit <- cv.glmnet(patient.data$x, Surv(patient.data$time, patient.data$status), family="cox", maxit = 10000)
  fit <- glmnet(patient.data$x, Surv(patient.data$time,patient.data$status), family =  "cox", maxit = 10000)
  #plot(cv.fit)
  cv.fit$lambda.min
  Coefficients <- coef(fit, s = cv.fit$lambda.min)
  Active.Index <- which(Coefficients != 0)
  Active.Coefficients  <- Coefficients[Active.Index]
  Active.Index
  betainfo <- tibble(name=dataname[Active.Index],beta=Active.Coefficients) %>% mutate(lambda=cv.fit$lambda.min,Seq=i)
  seqdata <- NULL
  suppressMessages(seqdata <- tmpdata %>% select(-age,-gender,-stage) %>% select(Tumor_Barcode,one_of(betainfo$name)) %>% pivot_longer(cols = -Tumor_Barcode) %>% left_join(betainfo %>% select(name,beta)) %>% mutate(value=value*beta) %>% group_by(Tumor_Barcode) %>% summarise(score=sum(value)) %>% mutate(Key=if_else(score>median(score),"Y","N")) %>% left_join(suvdata_tmp))
  
  fit <- coxph(Surv(survival_months, death) ~ Key+age+gender+stage, data = seqdata) ### overall 
  suvpvalue <- tidy(fit,conf.int = T) %>% filter(str_detect(term,"Key")) %>% mutate(Seq=i,HR=exp(estimate),HR_low=exp(conf.low),HR_high=exp(conf.high)) %>% select(Seq,p.value,contains('HR'))
  
  betainfo <- left_join(betainfo,suvpvalue)
  
  betainfo_all <- bind_rows(betainfo_all,betainfo)
  
}
save(betainfo_all,file='betainfo_all.RData')


betainfo_all %>% left_join(
  betainfo_all %>% count(Seq)
) %>% 
  group_by(name,beta,lambda,p.value,HR,HR_low,HR_high) %>% slice(1) %>% View()



# plot
besti <- 24
betainfo <- betainfo_all %>% filter(Seq==besti) %>% select(name:Seq)
betainfo
seqdata <- tmpdata %>% select(-age,-gender,-stage) %>% select(Tumor_Barcode,one_of(betainfo$name)) %>% pivot_longer(cols = -Tumor_Barcode) %>% left_join(betainfo %>% select(name,beta)) %>% mutate(value=value*beta) %>% group_by(Tumor_Barcode) %>% summarise(score=sum(value)) %>% mutate(Key=if_else(score>median(score),"Y","N")) %>% left_join(suvdata_tmp)
SurvZTW(suvdata = seqdata,plot = TRUE,keyname="Risk Score",filename = 'Figures/Survival/test_survival.pdf')


## polar plots

betainfo %>% 
  filter(name!="stage") %>% 
  mutate(Key=if_else(beta>0,"+","-")) %>% 
  mutate(name=fct_reorder(name,beta)) %>% 
  ggplot( aes(name, abs(beta), fill = Key)) +
  geom_bar(width = 0.97, stat = "identity", color = "white") +
  scale_y_continuous(limits = c(0,0.55)) +
  theme_minimal() +
  theme(axis.line.y=element_line(colour =  "#cccccc"),axis.ticks.y = element_line())+
  scale_fill_igv()+
  coord_polar(start = 4.5) 

ggsave(filename = 'Figures/betainfo_polar_plot.pdf',width = 5,height = 5,device = cairo_pdf)

tmpdata %>% select(-age,-gender,-stage) %>% select(Tumor_Barcode,one_of(betainfo$name)) %>% select(-Tumor_Barcode) %>% 
  count(`MDM2|SCNA_Focal`, `TP53|Mutation_Driver`, `ATM|HRD_LOH`, `CHEK2|HRD_LOH`, `Chr15q|SCNA_Arm`, `Chr22q|SCNA_Arm`, `SETD2|Mutation_Driver`) %>% 
  arrange(desc(n)) %>% mutate(Seq=seq_along(n)) %>% 
  pivot_longer(cols = -c(Seq,n)) %>% left_join(betainfo %>% select(name,beta)) %>% mutate(value=value*beta) %>% group_by(Seq,n) %>% summarise(score=sum(value)) %>% 
  ungroup() %>% 
  arrange(score) %>% view()


seqdata %>% filter(score!=0) %>% ggplot(aes(score))+geom_density()+geom_vline(xintercept = c(betainfo$beta))



# Old analyses ------------------------------------------------------------
#
# RTK ---------------------------------------------------------------------

fit <- coxph(Surv(survival_months, death) ~ RTK_Altered_Status+age+gender+stage, data = suvdata) ### overall 
ggforest(fit)
summary(fit)
suvpvalue <- tidy(fit) %>% filter(term=="RTK_Altered_StatusRTK/RAS +") %>% pull(p.value)
ggsave('RTK-RAS-Survival-forest.pdf',width = 7,height = 4,device = cairo_pdf)

fit2 <- coxph(Surv(survival_months, death) ~ strata(RTK_Altered_Status)+age+gender+stage, data = suvdata) ## strata curve ##
##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
suvfit <- survfit(fit2)
plot(suvfit)
suvfit$call$formula <- as.formula("Surv(survival_months, death) ~ strata(RTK_Altered_Status)+age+gender+stage")

ggsurv <- ggsurvplot(
  suvfit,                     # survfit object with calculated statistics.
  data = suvdata,             # data used to fit survival curves.
  risk.table = FALSE,       # show risk table.
  pval = round(suvpvalue,digits = 3),             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,210),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv"
)
ggsurv$plot <- ggsurv$plot + labs(
  title    = "Survival curves",                     
  subtitle = "Based on Kaplan-Meier estimates"
)
ggsurv <- ggpar(
  ggsurv,
  font.title    = c(16, "bold", "darkblue"),         
  font.subtitle = c(15, "bold.italic", "purple"), 
  font.caption  = c(14, "plain", "orange"),        
  font.x        = c(14, "bold.italic", "red"),          
  font.y        = c(14, "bold.italic", "darkred"),      
  font.xtickslab = c(12, "plain", "darkgreen"),
  legend = "top"
)

ggsave('RTK-RAS-Survival.pdf',width = 10,height = 6,device = cairo_pdf)

### kataegis survival #####
fit <- coxph(Surv(survival_months, death) ~ Kataegis+age+gender+stage, data = suvdata) ### overall 
ggforest(fit)
summary(fit)
suvpvalue <- tidy(fit) %>% filter(term=="KataegisY") %>% pull(p.value)
ggsave('Kataegis-Survival-forest.pdf',width = 7,height = 4,device = cairo_pdf)

fit2 <- coxph(Surv(survival_months, death) ~ strata(Kataegis)+age+gender+stage, data = suvdata) ## strata curve ##

suvfit <- survfit(fit2)
plot(suvfit)
suvfit$call$formula <- as.formula("Surv(survival_months, death) ~ strata(Kataegis)+age+gender+stage")

ggsurv <- ggsurvplot(
  suvfit,                     # survfit object with calculated statistics.
  data = suvdata,             # data used to fit survival curves.
  risk.table = FALSE,       # show risk table.
  pval = round(suvpvalue,3),             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#7570b3", "#e7298a"),
  xlim = c(0,210),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("Kataegis=N", "Kataegis=Y")    # change legend labels.
)
ggsurv$plot <- ggsurv$plot + labs(
  title    = "Survival curves",                     
  subtitle = "Based on Kaplan-Meier estimates"
)
ggsurv <- ggpar(
  ggsurv,
  font.title    = c(16, "bold", "darkblue"),         
  font.subtitle = c(15, "bold.italic", "purple"), 
  font.caption  = c(14, "plain", "orange"),        
  font.x        = c(14, "bold.italic", "red"),          
  font.y        = c(14, "bold.italic", "darkred"),      
  font.xtickslab = c(12, "plain", "darkgreen"),
  legend = "top"
)

ggsurv

ggsave('Kataegis-Survival.pdf',width = 10,height = 6,device = cairo_pdf)


#### MSI suvrvial ####
fit <- coxph(Surv(survival_months, death) ~ MSI_Status+age+gender+stage, data = suvdata) ### overall 
ggforest(fit)
summary(fit)
suvpvalue <- tidy(fit) %>% filter(term=="MSI_StatusMSI") %>% pull(p.value)
ggsave('MSI-Survival-forest.pdf',width = 7,height = 4,device = cairo_pdf)

fit2 <- coxph(Surv(survival_months, death) ~ strata(MSI_Status)+age+gender+stage, data = suvdata) ## strata curve ##

suvfit <- survfit(fit2)
plot(suvfit)
suvfit$call$formula <- as.formula("Surv(survival_months, death) ~ strata(MSI_Status)+age+gender+stage")

ggsurv <- ggsurvplot(
  suvfit,                     # survfit object with calculated statistics.
  data = suvdata,             # data used to fit survival curves.
  risk.table = FALSE,       # show risk table.
  pval = paste0("p=",as.character(scientific(suvpvalue))),            # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#ca0020", "#0571b0"),
  xlim = c(0,210),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv"  # add the median survival pointer.
  # legend.labs = 
  #   c("MSS", "MSI")    # change legend labels.
)
ggsurv$plot <- ggsurv$plot + labs(
  title    = "Survival curves",                     
  subtitle = "Based on Kaplan-Meier estimates"
)
ggsurv <- ggpar(
  ggsurv,
  font.title    = c(16, "bold", "darkblue"),         
  font.subtitle = c(15, "bold.italic", "purple"), 
  font.caption  = c(14, "plain", "orange"),        
  font.x        = c(14, "bold.italic", "red"),          
  font.y        = c(14, "bold.italic", "darkred"),      
  font.xtickslab = c(12, "plain", "darkgreen"),
  legend = "top"
)

ggsurv

ggsave('MSI-Survival.pdf',width = 10,height = 6,device = cairo_pdf)


# WGD ---------------------------------------------------------------------

fit <- coxph(Surv(survival_months, death) ~ WGD_Status+age+gender+stage, data = suvdata) ### overall 
ggforest(fit)
summary(fit)
suvpvalue <- tidy(fit) %>% filter(term=="WGD_StatusY") %>% pull(p.value)
ggsave('WGD-Survival-forest.pdf',width = 7,height = 4,device = cairo_pdf)

fit2 <- coxph(Surv(survival_months, death) ~ strata(WGD_Status)+age+gender+stage, data = suvdata) ## strata curve ##
##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
suvfit <- survfit(fit2)
plot(suvfit)
suvfit$call$formula <- as.formula("Surv(survival_months, death) ~ strata(WGD_Status)+age+gender+stage")

ggsurv <- ggsurvplot(
  suvfit,                     # survfit object with calculated statistics.
  data = suvdata,             # data used to fit survival curves.
  risk.table = FALSE,       # show risk table.
  pval = round(suvpvalue,digits = 3),             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,210),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("WGD N", "WGD Y")    # change legend labels.
)
ggsurv$plot <- ggsurv$plot + labs(
  title    = "Survival curves",                     
  subtitle = "Based on Kaplan-Meier estimates"
)
ggsurv <- ggpar(
  ggsurv,
  font.title    = c(16, "bold", "darkblue"),         
  font.subtitle = c(15, "bold.italic", "purple"), 
  font.caption  = c(14, "plain", "orange"),        
  font.x        = c(14, "bold.italic", "red"),          
  font.y        = c(14, "bold.italic", "darkred"),      
  font.xtickslab = c(12, "plain", "darkgreen"),
  legend = "top"
)

ggsave('WGD-Survival.pdf',width = 10,height = 6,device = cairo_pdf)


# Chr19loss ---------------------------------------------------------------------

fit <- coxph(Surv(survival_months, death) ~ Chr19loss_Clust+age+gender+stage, data = suvdata) ### overall 
ggforest(fit)
summary(fit)
suvpvalue <- tidy(fit) %>% filter(term=="Chr19loss_ClustY") %>% pull(p.value)
ggsave('chr19loss-Survival-forest.pdf',width = 7,height = 4,device = cairo_pdf)

fit2 <- coxph(Surv(survival_months, death) ~ strata(Chr19loss_Clust)+age+gender+stage, data = suvdata) ## strata curve ##
##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
suvfit <- survfit(fit2)
plot(suvfit)
suvfit$call$formula <- as.formula("Surv(survival_months, death) ~ strata(Chr19loss_Clust)+age+gender+stage")

ggsurv <- ggsurvplot(
  suvfit,                     # survfit object with calculated statistics.
  data = suvdata,             # data used to fit survival curves.
  risk.table = FALSE,       # show risk table.
  pval = round(suvpvalue,digits = 3),             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,210),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("Chr19loss N", "Chr19loss Y")    # change legend labels.
)
ggsurv$plot <- ggsurv$plot + labs(
  title    = "Survival curves",                     
  subtitle = "Based on Kaplan-Meier estimates"
)
ggsurv <- ggpar(
  ggsurv,
  font.title    = c(16, "bold", "darkblue"),         
  font.subtitle = c(15, "bold.italic", "purple"), 
  font.caption  = c(14, "plain", "orange"),        
  font.x        = c(14, "bold.italic", "red"),          
  font.y        = c(14, "bold.italic", "darkred"),      
  font.xtickslab = c(12, "plain", "darkgreen"),
  legend = "top"
)

ggsave('Chr19loss-Survival.pdf',width = 10,height = 6,device = cairo_pdf)


# CNV_Cluster ---------------------------------------------------------------------

fit <- coxph(Surv(survival_months, death) ~ CNV_Clust+age+gender+stage, data = suvdata) ### overall 
ggforest(fit)
summary(fit)
suvpvalue <- tidy(fit) %>% filter(term=="CNV_ClustC2") %>% pull(p.value)
ggsave('CNV_Clust-Survival-forest.pdf',width = 7,height = 4,device = cairo_pdf)

fit2 <- coxph(Surv(survival_months, death) ~ strata(CNV_Clust)+age+gender+stage, data = suvdata) ## strata curve ##
##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
suvfit <- survfit(fit2)
plot(suvfit)
suvfit$call$formula <- as.formula("Surv(survival_months, death) ~ strata(CNV_Clust)+age+gender+stage")

ggsurv <- ggsurvplot(
  suvfit,                     # survfit object with calculated statistics.
  data = suvdata,             # data used to fit survival curves.
  risk.table = FALSE,       # show risk table.
  pval = round(suvpvalue,digits = 3),             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#E7B800", "#2E9FDF",'darkblue'),
  xlim = c(0,210),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv"  # add the median survival pointer.
  # legend.labs = 
  #   c("WGD N", "WGD Y")    # change legend labels.
)
ggsurv$plot <- ggsurv$plot + labs(
  title    = "Survival curves",                     
  subtitle = "Based on Kaplan-Meier estimates"
)
ggsurv <- ggpar(
  ggsurv,
  font.title    = c(16, "bold", "darkblue"),         
  font.subtitle = c(15, "bold.italic", "purple"), 
  font.caption  = c(14, "plain", "orange"),        
  font.x        = c(14, "bold.italic", "red"),          
  font.y        = c(14, "bold.italic", "darkred"),      
  font.xtickslab = c(12, "plain", "darkgreen"),
  legend = "top"
)

ggsave('CNV_Clust-Survival.pdf',width = 10,height = 6,device = cairo_pdf)

# SFTP ---------------------------------------------------------------------

fit <- coxph(Surv(survival_months, death) ~ SFTP+age+gender+stage, data = suvdata) ### overall 
ggforest(fit)
summary(fit)
suvpvalue <- tidy(fit) %>% filter(term=="SFTPY") %>% pull(p.value)
ggsave('SFTP-Survival-forest.pdf',width = 7,height = 4,device = cairo_pdf)

fit2 <- coxph(Surv(survival_months, death) ~ strata(SFTP)+age+gender+stage, data = suvdata) ## strata curve ##
##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
suvfit <- survfit(fit2)
plot(suvfit)
suvfit$call$formula <- as.formula("Surv(survival_months, death) ~ strata(SFTP)+age+gender+stage")

ggsurv <- ggsurvplot(
  suvfit,                     # survfit object with calculated statistics.
  data = suvdata,             # data used to fit survival curves.
  risk.table = FALSE,       # show risk table.
  pval = round(suvpvalue,digits = 3),             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,210),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("SFTP N", "SFTP Y")    # change legend labels.
)
ggsurv$plot <- ggsurv$plot + labs(
  title    = "Survival curves",                     
  subtitle = "Based on Kaplan-Meier estimates"
)
ggsurv <- ggpar(
  ggsurv,
  font.title    = c(16, "bold", "darkblue"),         
  font.subtitle = c(15, "bold.italic", "purple"), 
  font.caption  = c(14, "plain", "orange"),        
  font.x        = c(14, "bold.italic", "red"),          
  font.y        = c(14, "bold.italic", "darkred"),      
  font.xtickslab = c(12, "plain", "darkgreen"),
  legend = "top"
)

ggsave('SFTP-Survival.pdf',width = 10,height = 6,device = cairo_pdf)

# TP53 ---------------------------------------------------------------------

fit <- coxph(Surv(survival_months, death) ~ TP53+age+gender+stage, data = suvdata) ### overall 
ggforest(fit)
summary(fit)
suvpvalue <- tidy(fit) %>% filter(term=="TP53Y") %>% pull(p.value)
ggsave('TP53-Survival-forest.pdf',width = 7,height = 4,device = cairo_pdf)

fit2 <- coxph(Surv(survival_months, death) ~ strata(TP53)+age+gender+stage, data = suvdata) ## strata curve ##
##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
suvfit <- survfit(fit2)
plot(suvfit)
suvfit$call$formula <- as.formula("Surv(survival_months, death) ~ strata(TP53)+age+gender+stage")

ggsurv <- ggsurvplot(
  suvfit,                     # survfit object with calculated statistics.
  data = suvdata,             # data used to fit survival curves.
  risk.table = FALSE,       # show risk table.
  pval = round(suvpvalue,digits = 3),             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,210),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("TP53 N", "TP53 Y")    # change legend labels.
)
ggsurv$plot <- ggsurv$plot + labs(
  title    = "Survival curves",                     
  subtitle = "Based on Kaplan-Meier estimates"
)
ggsurv <- ggpar(
  ggsurv,
  font.title    = c(16, "bold", "darkblue"),         
  font.subtitle = c(15, "bold.italic", "purple"), 
  font.caption  = c(14, "plain", "orange"),        
  font.x        = c(14, "bold.italic", "red"),          
  font.y        = c(14, "bold.italic", "darkred"),      
  font.xtickslab = c(12, "plain", "darkgreen"),
  legend = "top"
)

ggsave('TP53-Survival.pdf',width = 10,height = 6,device = cairo_pdf)


# TERT ---------------------------------------------------------------------

fit <- coxph(Surv(survival_months, death) ~ TERT+age+gender+stage, data = suvdata) ### overall 
ggforest(fit)
summary(fit)
suvpvalue <- tidy(fit) %>% filter(term=="TERTY") %>% pull(p.value)
ggsave('TERT-Survival-forest.pdf',width = 7,height = 4,device = cairo_pdf)

fit2 <- coxph(Surv(survival_months, death) ~ strata(TERT)+age+gender+stage, data = suvdata) ## strata curve ##
##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
suvfit <- survfit(fit2)
plot(suvfit)
suvfit$call$formula <- as.formula("Surv(survival_months, death) ~ strata(TERT)+age+gender+stage")

ggsurv <- ggsurvplot(
  suvfit,                     # survfit object with calculated statistics.
  data = suvdata,             # data used to fit survival curves.
  risk.table = FALSE,       # show risk table.
  pval = round(suvpvalue,digits = 3),             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,210),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("TERT N", "TERT Y")    # change legend labels.
)
ggsurv$plot <- ggsurv$plot + labs(
  title    = "Survival curves",                     
  subtitle = "Based on Kaplan-Meier estimates"
)
ggsurv <- ggpar(
  ggsurv,
  font.title    = c(16, "bold", "darkblue"),         
  font.subtitle = c(15, "bold.italic", "purple"), 
  font.caption  = c(14, "plain", "orange"),        
  font.x        = c(14, "bold.italic", "red"),          
  font.y        = c(14, "bold.italic", "darkred"),      
  font.xtickslab = c(12, "plain", "darkgreen"),
  legend = "top"
)

ggsave('TERT-Survival.pdf',width = 10,height = 6,device = cairo_pdf)


# MDM2 ---------------------------------------------------------------------

fit <- coxph(Surv(survival_months, death) ~ MDM2+age+gender+stage, data = suvdata) ### overall 
ggforest(fit)
summary(fit)
suvpvalue <- tidy(fit) %>% filter(term=="MDM2Y") %>% pull(p.value)
ggsave('MDM2-Survival-forest.pdf',width = 7,height = 4,device = cairo_pdf)

fit2 <- coxph(Surv(survival_months, death) ~ strata(MDM2)+age+gender+stage, data = suvdata) ## strata curve ##
##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
suvfit <- survfit(fit2)
plot(suvfit)
suvfit$call$formula <- as.formula("Surv(survival_months, death) ~ strata(MDM2)+age+gender+stage")

ggsurv <- ggsurvplot(
  suvfit,                     # survfit object with calculated statistics.
  data = suvdata,             # data used to fit survival curves.
  risk.table = FALSE,       # show risk table.
  pval = round(suvpvalue,digits = 3),             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,210),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("MDM2 N", "MDM2 Y")    # change legend labels.
)
ggsurv$plot <- ggsurv$plot + labs(
  title    = "Survival curves",                     
  subtitle = "Based on Kaplan-Meier estimates"
)
ggsurv <- ggpar(
  ggsurv,
  font.title    = c(16, "bold", "darkblue"),         
  font.subtitle = c(15, "bold.italic", "purple"), 
  font.caption  = c(14, "plain", "orange"),        
  font.x        = c(14, "bold.italic", "red"),          
  font.y        = c(14, "bold.italic", "darkred"),      
  font.xtickslab = c(12, "plain", "darkgreen"),
  legend = "top"
)

ggsave('MDM2-Survival.pdf',width = 10,height = 6,device = cairo_pdf)

# CDKN2A ---------------------------------------------------------------------

fit <- coxph(Surv(survival_months, death) ~ CDKN2A+age+gender+stage, data = suvdata) ### overall 
ggforest(fit)
summary(fit)
suvpvalue <- tidy(fit) %>% filter(term=="CDKN2AY") %>% pull(p.value)
ggsave('CDKN2A-Survival-forest.pdf',width = 7,height = 4,device = cairo_pdf)

fit2 <- coxph(Surv(survival_months, death) ~ strata(CDKN2A)+age+gender+stage, data = suvdata) ## strata curve ##
##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
suvfit <- survfit(fit2)
plot(suvfit)
suvfit$call$formula <- as.formula("Surv(survival_months, death) ~ strata(CDKN2A)+age+gender+stage")

ggsurv <- ggsurvplot(
  suvfit,                     # survfit object with calculated statistics.
  data = suvdata,             # data used to fit survival curves.
  risk.table = FALSE,       # show risk table.
  pval = round(suvpvalue,digits = 3),             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,210),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("CDKN2A N", "CDKN2A Y")    # change legend labels.
)
ggsurv$plot <- ggsurv$plot + labs(
  title    = "Survival curves",                     
  subtitle = "Based on Kaplan-Meier estimates"
)
ggsurv <- ggpar(
  ggsurv,
  font.title    = c(16, "bold", "darkblue"),         
  font.subtitle = c(15, "bold.italic", "purple"), 
  font.caption  = c(14, "plain", "orange"),        
  font.x        = c(14, "bold.italic", "red"),          
  font.y        = c(14, "bold.italic", "darkred"),      
  font.xtickslab = c(12, "plain", "darkgreen"),
  legend = "top"
)

ggsave('CDKN2A-Survival.pdf',width = 10,height = 6,device = cairo_pdf)


# HLA ---------------------------------------------------------------------

fit <- coxph(Surv(survival_months, death) ~ HLA+age+gender+stage, data = suvdata) ### overall 
ggforest(fit)
summary(fit)
suvpvalue <- tidy(fit) %>% filter(term=="HLAY") %>% pull(p.value)
ggsave('HLA-Survival-forest.pdf',width = 7,height = 4,device = cairo_pdf)

fit2 <- coxph(Surv(survival_months, death) ~ strata(HLA)+age+gender+stage, data = suvdata) ## strata curve ##
##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
suvfit <- survfit(fit2)
plot(suvfit)
suvfit$call$formula <- as.formula("Surv(survival_months, death) ~ strata(HLA)+age+gender+stage")

ggsurv <- ggsurvplot(
  suvfit,                     # survfit object with calculated statistics.
  data = suvdata,             # data used to fit survival curves.
  risk.table = FALSE,       # show risk table.
  pval = round(suvpvalue,digits = 3),             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,210),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("HLA N", "HLA Y")    # change legend labels.
)
ggsurv$plot <- ggsurv$plot + labs(
  title    = "Survival curves",                     
  subtitle = "Based on Kaplan-Meier estimates"
)
ggsurv <- ggpar(
  ggsurv,
  font.title    = c(16, "bold", "darkblue"),         
  font.subtitle = c(15, "bold.italic", "purple"), 
  font.caption  = c(14, "plain", "orange"),        
  font.x        = c(14, "bold.italic", "red"),          
  font.y        = c(14, "bold.italic", "darkred"),      
  font.xtickslab = c(12, "plain", "darkgreen"),
  legend = "top"
)

ggsave('HLA-Survival.pdf',width = 10,height = 6,device = cairo_pdf)








# MDM2+ TP53 ---------------------------------------------------------------------
suvdata <- sherlock_suvdata
load('sherlock_p53.RData')
suvdata <- suvdata %>% left_join(sherlock_p53 %>% select(Tumor_Barcode,MDM2_TP53) %>% mutate(MDM2_TP53=factor(MDM2_TP53,levels=c('N','Y'))))

fit <- coxph(Surv(survival_months, death) ~ MDM2_TP53+age+gender+stage, data = suvdata) ### overall 
ggforest(fit)
summary(fit)
suvpvalue <- tidy(fit) %>% filter(term=="MDM2_TP53Y") %>% pull(p.value)
ggsave('p53-Survival-forest.pdf',width = 7,height = 4,device = cairo_pdf)

fit2 <- coxph(Surv(survival_months, death) ~ strata(MDM2_TP53)+age+gender+stage, data = suvdata) ## strata curve ##
##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
suvfit <- survfit(fit2)
plot(suvfit)
suvfit$call$formula <- as.formula("Surv(survival_months, death) ~ strata(MDM2_TP53)+age+gender+stage")

ggsurv <- ggsurvplot(
  suvfit,                     # survfit object with calculated statistics.
  data = suvdata,             # data used to fit survival curves.
  risk.table = FALSE,       # show risk table.
  pval = scientific(suvpvalue,digits = 3),             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,210),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("MDM2_TP53 N", "MDM2_TP53 Y")    # change legend labels.
)
ggsurv$plot <- ggsurv$plot + labs(
  title    = "Survival curves",                     
  subtitle = "Based on Kaplan-Meier estimates"
)
ggsurv <- ggpar(
  ggsurv,
  font.title    = c(16, "bold", "darkblue"),         
  font.subtitle = c(15, "bold.italic", "purple"), 
  font.caption  = c(14, "plain", "orange"),        
  font.x        = c(14, "bold.italic", "red"),          
  font.y        = c(14, "bold.italic", "darkred"),      
  font.xtickslab = c(12, "plain", "darkgreen"),
  legend = "top"
)

ggsave('p53-Survival.pdf',width = 10,height = 6,device = cairo_pdf)






# chr19p ---------------------------------------------------------------------
suvdata <- sherlock_suvdata
load('sherlock_19pq.RData')
suvdata <- suvdata %>% left_join(sherlock_19pq %>% select(Tumor_Barcode,Chr19p))

fit <- coxph(Surv(survival_months, death) ~ Chr19p+age+gender+stage, data = suvdata) ### overall 
ggforest(fit)
summary(fit)
suvpvalue <- tidy(fit) %>% filter(term=="Chr19pY") %>% pull(p.value)
ggsave('Chr19p-Survival-forest.pdf',width = 7,height = 4,device = cairo_pdf)

fit2 <- coxph(Surv(survival_months, death) ~ strata(Chr19p)+age+gender+stage, data = suvdata) ## strata curve ##
##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
suvfit <- survfit(fit2)
plot(suvfit)
suvfit$call$formula <- as.formula("Surv(survival_months, death) ~ strata(Chr19p)+age+gender+stage")

ggsurv <- ggsurvplot(
  suvfit,                     # survfit object with calculated statistics.
  data = suvdata,             # data used to fit survival curves.
  risk.table = FALSE,       # show risk table.
  pval = round(suvpvalue,digits = 3),             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,210),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("Chr19p N", "Chr19p Y")    # change legend labels.
)
ggsurv$plot <- ggsurv$plot + labs(
  title    = "Survival curves",                     
  subtitle = "Based on Kaplan-Meier estimates"
)
ggsurv <- ggpar(
  ggsurv,
  font.title    = c(16, "bold", "darkblue"),         
  font.subtitle = c(15, "bold.italic", "purple"), 
  font.caption  = c(14, "plain", "orange"),        
  font.x        = c(14, "bold.italic", "red"),          
  font.y        = c(14, "bold.italic", "darkred"),      
  font.xtickslab = c(12, "plain", "darkgreen"),
  legend = "top"
)

ggsave('Chr19p-Survival.pdf',width = 10,height = 6,device = cairo_pdf)


