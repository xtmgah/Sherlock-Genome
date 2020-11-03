set_wd()
libztw()

filelists <- list.files(pattern = "sherlock.*RData",)
filelists <- filelists[!str_detect(filelists,'landscape')]
lapply(filelists,load,.GlobalEnv)

load('HRD_LOH.RData')
HRD_LOH <- sherlock_samples_unique %>% 
  select(Tumor_Barcode) %>% 
  left_join(HRD_LOH) %>% 
  select(-Type,-Subject) %>% 
  pivot_wider(names_from = 'Gene',values_from = 'Alteration') %>% 
  select(-`NA`) %>% 
  mutate_all(~replace_na(., 'nLOH')) 

HRD_LOH <- HRD_LOH %>% mutate(HRD_LOH=if_else(ATM=="nLOH" & BRCA2=="nLOH" & RAD51D=="nLOH" & BRCA1=="nLOH" & CHEK2=="nLOH" & PALB2=="nLOH", "nLOH","LOH"))
HRD_LOH[HRD_LOH=="LOH"] <- "Y"
HRD_LOH[HRD_LOH=="nLOH"] <- "N"


load('quiet_tumor.RData')


data_all <- sherlock_samples_unique %>% 
  select(Subject,Tumor_Barcode) %>% 
  left_join(sherlock_rtk) %>% 
  left_join(
    sherlock_p53 %>% select(Tumor_Barcode,P53=MDM2_TP53)
  ) %>% 
  left_join(
    quiet_tumor %>% mutate(Quiet=if_else(group=="Other tumors (N=199)",'N','Y')) %>% select(Subject,Quiet)
  ) %>% 
  left_join(
    sherlock_tmb %>% select(Subject,TMB)
  ) %>% 
  left_join(
    sherlock_sis_numbers %>% select(Tumor_Barcode,SV_Num)
  ) %>% 
  left_join(
    sherlock_cnvinfo %>% select(Tumor_Barcode,CNV_Coverage)
  ) %>% left_join(
    sherlock_tl %>% mutate(TLratio=Tumor_TL/Normal_TL) %>% select(Tumor_Barcode,TLratio)
  ) %>% left_join(
    sherlock_kataegis %>% select(-Loci)
  ) %>% left_join(
    sherlock_wgd %>% select(Tumor_Barcode,WGD_Status)
  ) %>% 
  mutate(HLA_LOH=if_else(Subject %in% sherlock_lohhla$Subject,"Y","N")) %>% 
  left_join(
    sherlock_cnvinfo %>% select(Tumor_Barcode,Chr19_HD=Chr19loss_Clust)
  ) %>% 
  left_join(HRD_LOH)


data_all <- data_all %>% replace_na(list(TMB=0.01,SV_Num=0,CNV_Coverage=0,Kataegis="N"))
data_all %>% filter_all(any_vars(is.na(.)))

save(data_all,file="sherlock_group_analysis2.RData")

# RTK group ---------------------------------------------------------------
data1 <- data_all %>%
  select(Tumor_Barcode,RTK_Altered_Status,`TMB (log2)`=TMB,`SCNA (%)`=CNV_Coverage,`SV (log2)`=SV_Num,`T/N TL ratio`=TLratio) %>% 
  pivot_longer(cols = -c(Tumor_Barcode,RTK_Altered_Status)) %>% 
  mutate(name=fct_inorder(name))
data1_labe <- data1 %>% group_by(name) %>% 
  do(tidy(wilcox.test(value~RTK_Altered_Status,data=.))) %>% 
  mutate(label=paste0(name,'\n','P=',format(p.value,digits = 2,scientific = T))) %>% 
  select(name,label)
  
data1 <- data_all %>%
  mutate(CNV_Coverage=100*CNV_Coverage,TMB=log2(TMB),SV_Num=log2(SV_Num)) %>% 
  select(Tumor_Barcode,RTK_Altered_Status,`TMB (log2)`=TMB,`SCNA (%)`=CNV_Coverage,`SV (log2)`=SV_Num,`T/N TL ratio`=TLratio) %>% 
  pivot_longer(cols = -c(Tumor_Barcode,RTK_Altered_Status)) %>% 
  left_join(data1_labe) %>% 
  mutate(name=fct_inorder(label))

data2_labe <- data_all %>% 
  select(RTK_Altered_Status,`Kataegis (%)`=Kataegis,`WGD status (%)`=WGD_Status,HLA_LOH,Chr19_HD,HRD_LOH,`BRCA2 LOH`=BRCA2,BRCA1,PALB2,ATM) %>% 
  pivot_longer(cols = -RTK_Altered_Status) %>% 
  group_by(name) %>% 
  nest() %>% 
  mutate(matrix=map(data,table)) %>% 
  mutate(test=map(matrix,fisher.test)) %>% 
  mutate(result=map(test,tidy)) %>% 
  select(name,result) %>% 
  unnest() %>% 
  ungroup() %>% 
  mutate(label=paste0(name,'\n','P=',format(p.value,digits = 2,scientific = T),'\nOR=',format(estimate,digits = 2))) %>% 
  select(name,label,p.value)

data2 <- data_all %>% 
  select(RTK_Altered_Status,`Kataegis (%)`=Kataegis,`WGD status (%)`=WGD_Status,`BRCA2 LOH`=BRCA2) %>% 
  pivot_longer(cols = -RTK_Altered_Status) %>% 
  group_by(RTK_Altered_Status,name) %>% 
  count(value) %>% 
  mutate(freq=round(n/sum(n),2)) %>% 
  ungroup() %>% 
  left_join(data2_labe) %>% 
  mutate(name=fct_reorder(label,p.value))


data_all %>% 
  select(Tumor_Barcode,RTK_Altered_Status) %>% 
  left_join(sherlock_SBS %>% mutate(APOBEC=SBS2+SBS13)) %>% 
  pivot_longer(cols = c(contains('SBS'),'APOBEC')) %>% 
  group_by(name) %>% 
  do(tidy(wilcox.test(value~RTK_Altered_Status,data=.))) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  mutate(FDR=p.adjust(p.value,method = 'BH'))


data3 <- data_all %>% 
  select(Tumor_Barcode,RTK_Altered_Status) %>% 
  left_join(sherlock_SBS) %>% 
  group_by(RTK_Altered_Status) %>% 
  summarise_at(vars(contains("SBS")),sum) %>% 
  gather(Signature,Weight,contains("SBS")) %>% 
  mutate(name="Mutation Signature (%)")

p1 <- groupplot1(data1)
p2 <- groupplot2(data2)+theme(plot.margin = margin(l=-1.3,unit = 'cm'))
p3 <- groupplot3(data3)+theme(plot.margin = margin(l=-0.2,unit = 'cm'))
pall <- plot_grid(p1,p2,p3,align = 'h',axis = 'tb',nrow = 1,rel_widths = c(4,1.8,1.3))
ggsave('Figures/RTK_group.pdf',plot = pall,width = 12,height = 8,device = cairo_pdf)



## enriched Cluster and unclustered SVs
load('/Volumes/data/NSLC/Clonality/Palimpsest_old/SV_signature/SV_signature.RData')

tmp <- data_all %>% dplyr::select(Tumor_Barcode,RTK_Altered_Status) %>% 
  left_join(SVsignatures_exp$sig_nums %>% rownames_to_column(var="Tumor_Barcode")) %>% 
  rename(`Clustered SVs`=Signature.2,`Non-clustered SVs`=Signature.1) %>% 
  pivot_longer(cols = -c(Tumor_Barcode,RTK_Altered_Status))

tmplab <- tmp %>% group_by(name) %>% 
  do(tidy(wilcox.test(value~RTK_Altered_Status,data=.))) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(label=paste0(name,'<br>*P* = ',format(p.value,digits  = 2))) %>% 
  ungroup()
  
library(ggtext)
tmp %>% 
  left_join(tmplab) %>% 
  ggplot(aes(RTK_Altered_Status,value))+
  geom_boxplot(fill = gray(0.85),outlier.shape=NA,width=0.6)+
  geom_jitter(position=position_jitter(w=0.18,h=0.1),size=2,shape=21,fill="gray60")+  
  labs(x="",y="Number of SVs")+
  facet_wrap(~label,scales = "free_y",nrow = 1)+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(base_size = 14,axis_title_size = 14,axis_title_just = 'm',grid = "XYxy",ticks = TRUE,axis = "")+
  theme(panel.spacing.x=unit(0.5, "lines"))+
  panel_border(size = 0.5,color = "black")+
  theme(strip.background = element_blank(),strip.text = element_markdown(hjust = 0.5))
  
ggsave(file="Figures/RTK_group_SVs_Cluster.pdf",width = 6,height =7,device = cairo_pdf)





# P53 group ---------------------------------------------------------------
data1 <- data_all %>%
  select(Tumor_Barcode,P53,`TMB (log2)`=TMB,`SCNA (%)`=CNV_Coverage,`SV (log2)`=SV_Num,`T/N TL ratio`=TLratio) %>% 
  pivot_longer(cols = -c(Tumor_Barcode,P53)) %>% 
  mutate(name=fct_inorder(name))
data1_labe <- data1 %>% group_by(name) %>% 
  do(tidy(wilcox.test(value~P53,data=.))) %>% 
  mutate(label=paste0(name,'\n','P=',format(p.value,digits = 2,scientific = T))) %>% 
  select(name,label)

data1 <- data_all %>%
  mutate(CNV_Coverage=100*CNV_Coverage,TMB=log2(TMB),SV_Num=log2(SV_Num)) %>% 
  select(Tumor_Barcode,P53,`TMB (log2)`=TMB,`SCNA (%)`=CNV_Coverage,`SV (log2)`=SV_Num) %>% 
  pivot_longer(cols = -c(Tumor_Barcode,P53)) %>% 
  left_join(data1_labe) %>% 
  mutate(name=fct_inorder(label))

data2_labe <- data_all %>% 
  select(P53,`Kataegis (%)`=Kataegis,`WGD status (%)`=WGD_Status,`HLA LOH`=HLA_LOH,Chr19_HD,HRD_LOH,`BRCA1 LOH`=BRCA1,BRCA2,PALB2,ATM) %>% 
  pivot_longer(cols = -P53) %>% 
  group_by(name) %>% 
  nest() %>% 
  mutate(matrix=map(data,table)) %>% 
  mutate(test=map(matrix,fisher.test)) %>% 
  mutate(result=map(test,tidy)) %>% 
  select(name,result) %>% 
  unnest() %>% 
  ungroup() %>% 
  mutate(label=paste0(name,'\n','P=',format(p.value,digits = 2,scientific = T),'\nOR=',format(estimate,digits = 2))) %>% 
  select(name,label,estimate ,p.value)

data2 <- data_all %>% 
  select(P53,`Kataegis (%)`=Kataegis,`WGD status (%)`=WGD_Status,`HLA LOH`=HLA_LOH,`BRCA1 LOH`=BRCA1) %>% 
  pivot_longer(cols = -P53) %>% 
  group_by(P53,name) %>% 
  count(value) %>% 
  mutate(freq=round(n/sum(n),2)) %>% 
  ungroup() %>% 
  left_join(data2_labe) %>% 
  mutate(name=fct_reorder(label,p.value))


data_all %>% 
  select(Tumor_Barcode,P53) %>% 
  left_join(sherlock_SBS %>% mutate(APOBEC=SBS2+SBS13)) %>% 
  pivot_longer(cols = c(contains('SBS'),APOBEC)) %>% 
  group_by(name) %>% 
  do(tidy(wilcox.test(value~P53,data=.))) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  mutate(FDR=p.adjust(p.value,method = 'BH'))


data3 <- data_all %>% 
  select(Tumor_Barcode,P53) %>% 
  left_join(sherlock_SBS) %>% 
  group_by(P53) %>% 
  summarise_at(vars(contains("SBS")),sum) %>% 
  gather(Signature,Weight,contains("SBS")) %>% 
  mutate(name="Mutation Signature (%)")

p1 <- groupplot1(data1,var = P53)
p2 <- groupplot2(data2,var = P53)+theme(plot.margin = margin(l=-1.3,unit = 'cm'))
p3 <- groupplot3(data3,var = P53)+theme(plot.margin = margin(l=-0.2,unit = 'cm'))
pall <- plot_grid(p1,p2,p3,align = 'h',axis = 'tb',nrow = 1,rel_widths = c(4,1.2,1.3))
ggsave('Figures/p53_group.pdf',plot = pall,width = 10,height = 8,device = cairo_pdf)


## enriched Cluster and unclustered SVs
load('/Volumes/data/NSLC/Clonality/Palimpsest_old/SV_signature/SV_signature.RData')

tmp <- data_all %>% dplyr::select(Tumor_Barcode,P53) %>% 
  left_join(SVsignatures_exp$sig_nums %>% rownames_to_column(var="Tumor_Barcode")) %>% 
  rename(`Clustered SVs`=Signature.2,`Non-clustered SVs`=Signature.1) %>% 
  pivot_longer(cols = -c(Tumor_Barcode,P53))

tmplab <- tmp %>% group_by(name) %>% 
  do(tidy(wilcox.test(value~P53,data=.))) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(label=paste0(name,'<br>*P* = ',format(p.value,digits  = 2))) %>% 
  ungroup()

library(ggtext)
tmp %>% 
  left_join(tmplab) %>% 
  ggplot(aes(P53,value))+
  geom_boxplot(fill = gray(0.85),outlier.shape=NA,width=0.6)+
  geom_jitter(position=position_jitter(w=0.18,h=0.1),size=2,shape=21,fill="gray60")+  
  labs(x="",y="Number of SVs")+
  facet_wrap(~label,scales = "free_y",nrow = 1)+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(base_size = 14,axis_title_size = 14,axis_title_just = 'm',grid = "XYxy",ticks = TRUE,axis = "")+
  theme(panel.spacing.x=unit(0.5, "lines"))+
  panel_border(size = 0.5,color = "black")+
  theme(strip.background = element_blank(),strip.text = element_markdown(hjust = 0.5))

ggsave(file="Figures/P53_group_SVs_Cluster.pdf",width = 6,height =7,device = cairo_pdf)



# Quiet Tumor -------------------------------------------------------------
data1 <- data_all %>%
  select(Tumor_Barcode,Quiet,`TMB (log2)`=TMB,`SCNA (%)`=CNV_Coverage,`SV (log2)`=SV_Num,`T/N TL ratio`=TLratio) %>% 
  pivot_longer(cols = -c(Tumor_Barcode,Quiet)) %>% 
  mutate(name=fct_inorder(name))
data1_labe <- data1 %>% group_by(name) %>% 
  do(tidy(wilcox.test(value~Quiet,data=.))) %>% 
  mutate(label=paste0(name,'\n','P=',format(p.value,digits = 2,scientific = T))) %>% 
  select(name,label)

data1 <- data_all %>%
  mutate(CNV_Coverage=100*CNV_Coverage,TMB=log2(TMB),SV_Num=log2(SV_Num)) %>% 
  select(Tumor_Barcode,Quiet,`TMB (log2)`=TMB,`SCNA (%)`=CNV_Coverage,`SV (log2)`=SV_Num,`T/N TL ratio`=TLratio) %>% 
  pivot_longer(cols = -c(Tumor_Barcode,Quiet)) %>% 
  left_join(data1_labe) %>% 
  mutate(name=fct_inorder(label))

data2_labe <- data_all %>% 
  select(Quiet,`Kataegis (%)`=Kataegis,`WGD status (%)`=WGD_Status,`HLA LOH`=HLA_LOH,Chr19_HD,`HRD LOH`=HRD_LOH,`BRCA1 LOH`=BRCA1,BRCA2,PALB2,ATM) %>% 
  pivot_longer(cols = -Quiet) %>% 
  group_by(name) %>% 
  nest() %>% 
  mutate(matrix=map(data,table)) %>% 
  mutate(test=map(matrix,fisher.test)) %>% 
  mutate(result=map(test,tidy)) %>% 
  select(name,result) %>% 
  unnest() %>% 
  ungroup() %>% 
  mutate(label=paste0(name,'\n','P=',format(p.value,digits = 2,scientific = T),'\nOR=',format(estimate,digits = 2))) %>% 
  select(name,label,estimate ,p.value)

data2 <- data_all %>% 
  select(Quiet,`Kataegis (%)`=Kataegis,`WGD status (%)`=WGD_Status,`HRD LOH`=HRD_LOH) %>% 
  pivot_longer(cols = -Quiet) %>% 
  group_by(Quiet,name) %>% 
  count(value) %>% 
  mutate(freq=round(n/sum(n),2)) %>% 
  ungroup() %>% 
  left_join(data2_labe) %>% 
  mutate(name=fct_reorder(label,p.value))


data_all %>% 
  select(Tumor_Barcode,Quiet) %>% 
  left_join(sherlock_SBS) %>% 
  pivot_longer(cols = contains('SBS')) %>% 
  group_by(name) %>% 
  do(tidy(wilcox.test(value~Quiet,data=.))) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  mutate(FDR=p.adjust(p.value,method = 'BH'))


data3 <- data_all %>% 
  select(Tumor_Barcode,Quiet) %>% 
  left_join(sherlock_SBS) %>% 
  group_by(Quiet) %>% 
  summarise_at(vars(contains("SBS")),sum) %>% 
  gather(Signature,Weight,contains("SBS")) %>% 
  mutate(name="Mutation Signature (%)")

p1 <- groupplot1(data1,var = Quiet)
p2 <- groupplot2(data2,var = Quiet)+theme(plot.margin = margin(l=-1.3,unit = 'cm'))
p3 <- groupplot3(data3,var = Quiet)+theme(plot.margin = margin(l=-0.2,unit = 'cm'))
pall <- plot_grid(p1,p2,p3,align = 'h',axis = 'tb',nrow = 1,rel_widths = c(4,1.8,1.3))
ggsave('Figures/Quiet_Tumor_group.pdf',plot = pall,width = 12,height = 8,device = cairo_pdf)



## TL, Subclone ratio, NPRCC, purity

load('sherlock_samples_unique.RData')
tmp <- data_all %>% 
  select(Tumor_Barcode,Quiet,`T/N TL ratio`=TLratio) %>% 
  left_join(
    sherlock_samples_unique %>% select(Tumor_Barcode,`Tumor purity`=Tumor_Purity,NRPCC,`Subclonal mutation ratio`=Subclonal_Mutation_Ratio)
  ) %>% 
  pivot_longer(cols = -c(Tumor_Barcode,Quiet))

tmplab <- tmp %>% group_by(name) %>% 
  do(tidy(wilcox.test(value~Quiet,data=.))) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(label=paste0(name,'<br>*P* = ',format(p.value,digits  = 2))) %>% 
  ungroup()

tmplab <- tmplab %>% mutate(name=factor(name,levels = unique(tmp$name)[c(2,4,3,1)])) %>% arrange(name) %>% mutate(label=fct_inorder(label))
load('sherlock_histology.RData')

histcolor <- c('#a6611a','#f1b6da','#d01c8b')
names(histcolor) <- unique(sherlock_histology$Histology)

library(ggtext)
tmp %>% left_join(sherlock_histology) %>% 
  left_join(tmplab) %>% 
  ggplot(aes(Quiet,value))+
  geom_boxplot(fill = gray(0.85),outlier.shape=NA,width=0.6)+
  geom_jitter(aes(fill=Histology),position=position_jitter(w=0.18,h=0.1),size=3,shape=21)+  
  labs(x="Pianissimo tumor",y="")+
  facet_wrap(~label,scales = "free_y",nrow = 1)+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(base_size = 14,axis_title_size = 14,axis_title_just = 'm',grid = "XYxy",ticks = TRUE,axis = "")+
  theme(panel.spacing.x=unit(0.5, "lines"))+
  panel_border(size = 0.5,color = "black")+
  scale_fill_manual(values = histcolor)+
  theme(strip.background = element_blank(),strip.text = element_markdown(hjust = 0.5))

ggsave(file="Figures/Quiet_Tumor_features.pdf",width = 9,height =8,device = cairo_pdf)






# RTK and P53 combined ----------------------------------------------------

data_all <- data_all %>% mutate(Group=paste0(RTK_Altered_Status,"\nP53 ",P53)) %>% mutate(Group=factor(Group,levels = c("RTK-RAS-\nP53 N", "RTK-RAS+\nP53 N", "RTK-RAS-\nP53 Y", "RTK-RAS+\nP53 Y")))
data1 <- data_all %>% filter(Group %in% c("RTK-RAS+\nP53 Y","RTK-RAS-\nP53 N")) %>% 
  select(Tumor_Barcode,Group,`TMB (log2)`=TMB,`SCNA (%)`=CNV_Coverage,`SV (log2)`=SV_Num,`T/N TL ratio`=TLratio) %>% 
  pivot_longer(cols = -c(Tumor_Barcode,Group)) %>% 
  mutate(name=fct_inorder(name)) 
data1_labe <- data1 %>% group_by(name) %>% 
  do(tidy(wilcox.test(value~Group,data=.))) %>% 
  mutate(label=paste0(name,'\n','P=',format(p.value,digits = 2,scientific = T))) %>% 
  select(name,label)

data1 <- data_all %>%
  mutate(CNV_Coverage=100*CNV_Coverage,TMB=log2(TMB),SV_Num=log2(SV_Num)) %>% 
  select(Tumor_Barcode,Group,`TMB (log2)`=TMB,`SCNA (%)`=CNV_Coverage,`SV (log2)`=SV_Num,`T/N TL ratio`=TLratio) %>% 
  pivot_longer(cols = -c(Tumor_Barcode,Group)) %>% 
  left_join(data1_labe) %>% 
  mutate(name=fct_inorder(label))

data2_labe <- data_all %>% filter(Group %in% c("RTK-RAS+\nP53 Y","RTK-RAS-\nP53 N")) %>% 
  select(Group,`Kataegis (%)`=Kataegis,`WGD status (%)`=WGD_Status,`HLA LOH`=HLA_LOH,Chr19_HD,`HRD LOH`=HRD_LOH,`BRCA2 LOH`=BRCA2,`BRCA1 LOH`=BRCA1,PALB2,ATM) %>% 
  pivot_longer(cols = -Group) %>% 
  group_by(name) %>% 
  nest() %>% 
  mutate(matrix=map(data,table)) %>% 
  mutate(test=map(matrix,fisher.test)) %>% 
  mutate(result=map(test,tidy)) %>% 
  select(name,result) %>% 
  unnest() %>% 
  ungroup() %>% 
  mutate(label=paste0(name,'\n','P=',format(p.value,digits = 2,scientific = T),'\nOR=',format(estimate,digits = 2))) %>% 
  select(name,label,p.value)

data2 <- data_all %>% 
  select(Group,`Kataegis (%)`=Kataegis,`WGD status (%)`=WGD_Status,`BRCA2 LOH`=BRCA2,`BRCA1 LOH`=BRCA1,`HLA LOH`=HLA_LOH,`HRD LOH`=HRD_LOH) %>% 
  pivot_longer(cols = -Group) %>% 
  group_by(Group,name) %>% 
  count(value) %>% 
  mutate(freq=round(n/sum(n),2)) %>% 
  ungroup() %>% 
  left_join(data2_labe) %>% 
  mutate(name=fct_reorder(label,p.value))


data_all %>% filter(Group %in% c("RTK-RAS+\nP53 Y","RTK-RAS-\nP53 N")) %>% 
  select(Tumor_Barcode,Group) %>% 
  left_join(sherlock_SBS) %>% 
  pivot_longer(cols = contains('SBS')) %>% 
  group_by(name) %>% 
  do(tidy(wilcox.test(value~Group,data=.))) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  mutate(FDR=p.adjust(p.value,method = 'BH'))


data3 <- data_all %>% 
  select(Tumor_Barcode,Group) %>% 
  left_join(sherlock_SBS) %>% 
  group_by(Group) %>% 
  summarise_at(vars(contains("SBS")),sum) %>% 
  gather(Signature,Weight,contains("SBS")) %>% 
  mutate(name="Mutation Signature (%)")

p1 <- groupplot1(data1,var=Group)
p2 <- groupplot2(data2,var=Group)+theme(plot.margin = margin(l=-1.3,unit = 'cm'))
p3 <- groupplot3(data3,var=Group)+theme(plot.margin = margin(l=-0.2,unit = 'cm'))
pall <- plot_grid(p1,p2,p3,align = 'h',axis = 'tb',nrow = 1,rel_widths = c(4,2,1.1))
ggsave('Figures/rtk_p53_group.pdf',plot = pall,width = 18,height = 8,device = cairo_pdf)



# C1/C2/C3 ----------------------------------------------------------------
load('sherlock_group_analysis2.RData')
load('sherlock_cnvinfo.RData')
load('sherlock_histology.RData')
load('sherlock_samples_unique.RData')
load('sherlock_signature.RData')
tmp <- sherlock_cnvinfo %>% left_join(sherlock_histology) %>% mutate(Group=if_else(CNV_Clust %in% c("C1","C2"),CNV_Clust,if_else(Histology=="Adenocarcinomas","C3-LUAD",if_else(Histology=="Carcinoids","C3-Carcinoids",NA_character_)))) %>% select(Tumor_Barcode,Group) %>% filter(!is.na(Group)) %>% 
  mutate(Group=factor(Group,levels = c("C2","C1","C3-LUAD","C3-Carcinoids")))

data_all <- tmp %>% left_join(data_all)
data_all <- data_all %>% left_join(sherlock_samples_unique %>% select(Tumor_Barcode,Subclonal_Mutation_Ratio))

data1 <- data_all %>% filter(Group %in% c("C3-LUAD","C2")) %>% mutate(Group=factor(Group,levels = c("C2","C3-LUAD"))) %>% 
  select(Tumor_Barcode,Group,`TMB (log2)`=TMB,`SCNA (%)`=CNV_Coverage,`SV (log2)`=SV_Num,`T/N TL ratio`=TLratio,`Subclonal Mutation Ratio`=Subclonal_Mutation_Ratio) %>% 
  pivot_longer(cols = -c(Tumor_Barcode,Group)) %>% 
  mutate(name=fct_inorder(name)) 

data1_labe <- data1 %>% group_by(name) %>% 
  do(tidy(wilcox.test(value~Group,data=.))) %>% 
  mutate(label=paste0(name,'\n','P=',format(p.value,digits = 2,scientific = T))) %>% 
  select(name,label)

data1 <- data_all %>%
  mutate(CNV_Coverage=100*CNV_Coverage,TMB=log2(TMB),SV_Num=log2(SV_Num)) %>% 
  select(Tumor_Barcode,Group,`TMB (log2)`=TMB,`SCNA (%)`=CNV_Coverage,`SV (log2)`=SV_Num,`T/N TL ratio`=TLratio,`Subclonal Mutation Ratio`=Subclonal_Mutation_Ratio) %>% 
  pivot_longer(cols = -c(Tumor_Barcode,Group)) %>% 
  left_join(data1_labe) %>% 
  mutate(name=fct_inorder(label))

data2_labe <- data_all %>% filter(Group %in% c("C3-LUAD","C2")) %>% mutate(Group=factor(Group,levels = c("C2","C3-LUAD"))) %>% 
  select(Group,`Kataegis (%)`=Kataegis,`WGD status (%)`=WGD_Status,`HLA LOH`=HLA_LOH,Chr19_HD,`HRD LOH`=HRD_LOH,`BRCA2 LOH`=BRCA2,`BRCA1 LOH`=BRCA1,PALB2,ATM) %>% 
  pivot_longer(cols = -Group) %>% 
  group_by(name) %>% 
  nest() %>% 
  mutate(matrix=map(data,table)) %>% 
  mutate(test=map(matrix,fisher.test)) %>% 
  mutate(result=map(test,tidy)) %>% 
  select(name,result) %>% 
  unnest(cols = c(result)) %>% 
  ungroup() %>% 
  mutate(label=paste0(name,'\n','P=',format(p.value,digits = 2,scientific = T),'\nOR=',format(estimate,digits = 2))) %>% 
  select(name,label,p.value)

data2 <- data_all %>% 
  select(Group,`Kataegis (%)`=Kataegis,`WGD status (%)`=WGD_Status,`BRCA2 LOH`=BRCA2,`BRCA1 LOH`=BRCA1,`HLA LOH`=HLA_LOH,`HRD LOH`=HRD_LOH) %>% 
  pivot_longer(cols = -Group) %>% 
  group_by(Group,name) %>% 
  count(value) %>% 
  mutate(freq=round(n/sum(n),2)) %>% 
  ungroup() %>% 
  left_join(data2_labe) %>% 
  mutate(name=fct_reorder(label,p.value))


data_all %>% filter(Group %in% c("C3-LUAD","C2")) %>% mutate(Group=factor(Group,levels = c("C2","C3-LUAD"))) %>% 
  select(Tumor_Barcode,Group) %>% 
  left_join(sherlock_SBS) %>% 
  pivot_longer(cols = contains('SBS')) %>% 
  group_by(name) %>% 
  do(tidy(wilcox.test(value~Group,data=.))) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  mutate(FDR=p.adjust(p.value,method = 'BH'))


data3 <- data_all %>% 
  select(Tumor_Barcode,Group) %>% 
  left_join(sherlock_SBS) %>% 
  group_by(Group) %>% 
  summarise_at(vars(contains("SBS")),sum) %>% 
  gather(Signature,Weight,contains("SBS")) %>% 
  mutate(name="Mutation Signature (%)")


load('/Volumes/data/NSLC/Mutation/Signautres_final2/Analysis/com_sbs1536.RData')
data4 <- data_all %>% 
  select(Tumor_Barcode,Group) %>%
  left_join(sbs1536data %>% rename(Tumor_Barcode=Samples)) %>% 
  group_by(Group) %>% 
  summarise_at(vars(one_of(c('SBS','M6','M18','M52'))),sum) %>% 
  gather(Signature,Weight,one_of(c('SBS','M6','M18','M52'))) %>% 
  mutate(name="Mutation Signature (%)")

source('~/Box/Copied/cache/Biowulf/Sigvisualfunc.R')
msigcolor <- c('#bababa','#66c2a5','#fc8d62','#8da0cb')
names(msigcolor) <- colnames(sbs1536data)[-1]
SBScolor <- c(SBScolor,msigcolor)



p1 <- groupplot1(data1,var=Group)
p2 <- groupplot2(data2,var=Group)+theme(plot.margin = margin(l=-1.3,unit = 'cm'))
p3 <- groupplot3(data3,var=Group,SBScolor = SBScolor)+theme(plot.margin = margin(l=-0.2,unit = 'cm'))
p4 <- groupplot3(data4,var=Group,SBScolor = SBScolor)+theme(plot.margin = margin(l=-0.2,unit = 'cm'))
pall <- plot_grid(p1,p2,p3,p4,align = 'h',axis = 'tb',nrow = 1,rel_widths = c(5,2,1.1,1))
ggsave('Figures/SCNA_group.pdf',plot = pall,width = 18,height = 8,device = cairo_pdf)












## enriched Cluster and unclustered SVs
load('/Volumes/data/NSLC/Clonality/Palimpsest_old/SV_signature/SV_signature.RData')

tmp <- data_all %>% dplyr::select(Tumor_Barcode,Group) %>% 
  left_join(SVsignatures_exp$sig_nums %>% rownames_to_column(var="Tumor_Barcode")) %>% 
  rename(`Clustered SVs`=Signature.2,`Non-clustered SVs`=Signature.1) %>% 
  pivot_longer(cols = -c(Tumor_Barcode,Group))

tmplab <- tmp %>% filter(Group %in% c("RTK-RAS+\nP53 Y","RTK-RAS-\nP53 N")) %>% 
  group_by(name) %>% 
  do(tidy(wilcox.test(value~Group,data=.))) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(label=paste0(name,'<br>*P* = ',format(p.value,digits  = 2))) %>% 
  ungroup()

library(ggtext)
tmp %>% 
  left_join(tmplab) %>% 
  ggplot(aes(Group,value))+
  geom_boxplot(fill = gray(0.85),outlier.shape=NA,width=0.6)+
  geom_jitter(position=position_jitter(w=0.18,h=0.1),size=2,shape=21,fill="gray60")+  
  labs(x="",y="Number of SVs")+
  facet_wrap(~label,scales = "free_y",nrow = 1)+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(base_size = 14,axis_title_size = 14,axis_title_just = 'm',grid = "XYxy",ticks = TRUE,axis = "")+
  theme(panel.spacing.x=unit(0.5, "lines"))+
  panel_border(size = 0.5,color = "black")+
  theme(strip.background = element_blank(),strip.text = element_markdown(hjust = 0.5))

ggsave(file="Figures/rtk_p53_group_SVs_Cluster.pdf",width = 10,height =7,device = cairo_pdf)


# Function ---------------------------------------------------------------------

groupplot1 <- function(data,var=RTK_Altered_Status){ 
  library(drlib)
  rmeidan <- data %>% group_by({{var}},name) %>% summarise(Median=median(value,na.rm = TRUE)) %>% arrange(Median)
  npoint=length(levels(rmeidan$Group))
  p1 <- ggplot(data = data, aes(x = {{var}}, y = (value), group = reorder_within(Tumor_Barcode,value,name)))+
    geom_boxplot(size=2,fatten=NULL,col="grey40")+
    geom_point(data = rmeidan, aes(y = (Median), x = {{var}}), shape = 95, inherit.aes = FALSE, color = "red", size = 10)+
    facet_wrap(~name,scales = 'free_y',nrow = 1)+
    scale_y_continuous(breaks = pretty_breaks(),expand = expand_scale(mult = c(0.02,0.02)))+
    theme_ipsum_rc(base_size = 10,axis_title_size = 12,axis = '',grid = FALSE,axis_title_just = "m",ticks = TRUE,strip_text_size = 10)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.y = element_text(margin = margin(r = 5)),panel.spacing = unit(0.4, "lines"),axis.ticks.x = element_blank())+
    geom_vline(xintercept = 1:(npoint-1)+0.5,colour="#cccccc",linetype=2)+
    labs(x="",y="")+
    #coord_cartesian(clip = 'off')+
    panel_border(color = 'black',size = 0.3)
  #ggsave('Figures/tmp.pdf',plot = p1,width = 8,height = 8,device = cairo_pdf)
  return(p1)
}



groupplot2 <- function(data,var=RTK_Altered_Status){ 
  p2 <- data %>% 
    ggplot(aes({{var}},freq,fill=value))+
    geom_col()+
    geom_text(aes(label=paste0(freq*100,"%")),position = position_stack(vjust = 0.5),size=3)+ 
    facet_wrap(~name)+
    theme_ipsum_rc(base_size = 10,axis_title_size = 12,axis = 'y',axis_title_just = "m",ticks = FALSE,axis_col = "black",strip_text_size = 10)+
    theme(axis.title.x = element_blank(),legend.position = "top",legend.box.background = element_blank(),legend.box.spacing = unit(0.6,"cm"),legend.key = element_rect(size =5,colour = "white"),legend.key.height = unit(0.8,"cm"),axis.ticks.y = element_line(colour = "black"),legend.key.size = unit(0.5, "cm"),legend.text = element_text(size = 10),legend.box.margin=margin(l=-8),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),panel.grid.major.y = element_line(),axis.title.y = element_text(margin = margin(r = 5)),panel.spacing = unit(0.8, "lines"),axis.ticks.x = element_blank())+
    labs(y="")+
    scale_fill_d3()+
    labs(fill='')+
    scale_y_continuous(breaks = pretty_breaks(),label=percent,expand = expand_scale(mult = c(0,0)))+
    scale_x_discrete(expand = c(0.2,0))
  return(p2)
}

# Sigantures --------------------------------------------------------------
groupplot3 <- function(data,var=RTK_Altered_Status,SBScolor=SBScolor){ 
  #source('~/Box/Copied/cache/Biowulf/Sigvisualfunc.R')
  p3 <- data %>% 
    ggplot(aes({{var}},Weight,fill=factor(Signature,levels = names(SBScolor))))+
    geom_bar(stat="identity",position="fill",col="gray95",size=0)+
    facet_wrap(~name)+
    theme_ipsum_rc(base_size = 10,axis_title_just = "m",axis_title_size = 12,axis = "y",axis_col = "black",strip_text_size = 10)+
    theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major=element_line(),legend.position = "right",legend.box.background = element_blank(),legend.box.spacing = unit(0.6,"cm"),legend.key = element_rect(size =0),legend.key.height = unit(0.8,"cm"),axis.ticks.y = element_line(colour = "black"),legend.key.size = unit(0.5, "cm"),legend.text = element_text(size = 10),legend.box.margin=margin(l=-10),axis.ticks.x = element_blank())+
    scale_y_continuous(breaks = pretty_breaks(),label=percent,expand = expand_scale(mult = c(0,0)))+
    labs(x="",y="")+
    scale_x_discrete(expand = c(0.2,0))+
    scale_fill_manual("",values = SBScolor)
  return(p3)
}











