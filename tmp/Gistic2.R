setwd("/Volumes/data/NSLC/SCNA/Battenberg_Gistic/")
#rm(list=ls())
library(tidyverse)
library(maftools)
library(hrbrthemes)
library(scales)
library(cowplot)
library(ggrepel)
library(valr)

#load('../purity_ploid/sampleinfo3.RData')
load('~/NIH-Work/EAGLE_NSLC/FirstPaper/Biowulf/sherlock_samples.RData')
load('/Volumes/data/Squamous_Cell_Carcinoma/DriveGene/lung_drivegene.RData')
load('cytoband_hg19.RData')
hg19centro <- cyto %>% 
  filter(Centro!=0,Chromosome %in% c(1:22)) %>% 
  select(Chromosome,Start_Position=Start,End_Position=End) %>% 
  group_by(Chromosome) %>% 
  arrange(End_Position) %>% 
  slice(1) %>% 
  select(Chromosome,Centro=End_Position) %>% 
  ungroup()



##### CNV data #### 
am.gistic = readGistic(gisticAllLesionsFile = './result/all_lesions.conf_90.txt',gisticAmpGenesFile = 'result/amp_genes.conf_90.txt',gisticDelGenesFile = 'result/del_genes.conf_90.txt',gisticScoresFile = 'result/scores.gistic', isTCGA = FALSE,cnLevel = "all")
#gcp = gisticChromPlot(gistic = am.gistic, markBands = "all",ref.build = 'hg19',fdr=0.1)

gcp = GistChrom_ztw(gistic = am.gistic, markBands = "all",ref.build = 'hg19',fdr=1)


# pdf('somaticInteractions_cnv.pdf', width = 14,height = 10)
# gisticOncoPlot(gistic = am.gistic, clinicalData = getClinicalData(x = allw.maf.orignal), clinicalFeatures =c('Data','Ulceration','Gender') , sortByAnnotation = TRUE, top = 30)
# dev.off()




# link drive gene  --------------------------------------------------------------

amp <- read_delim('result/amp_genes.conf_90.txt',delim = "\t",col_names = T) %>% select(-X26)
amp1 <- amp %>% rename(variant=cytoband) %>% slice(1:3) %>% gather(cytoband,value,-variant) %>% spread(variant,value)
amp2 <- amp %>% rename(variant=cytoband) %>% slice(4:dim(amp)[1]) %>% gather(cytoband,value,-variant) %>% select(-variant) %>% filter(!is.na(value)) %>% rename(gene=value)
amp <- amp1 %>% left_join(amp2)
amp <- amp %>% mutate(gene=str_replace_all(gene,"[\\[\\]]","")) %>% left_join(amp %>% count(cytoband)) %>% mutate(drive_gene=if_else(gene %in% lung_drivegene$gene,gene,"")) %>% mutate(`q value`=as.numeric(`q value`)) %>% arrange(`q value`,n,cytoband) 
amp <- amp %>% group_by(cytoband,`q value`,`wide peak boundaries`,n) %>% summarise(drive_gene=paste0(drive_gene,collapse = ","),gene=paste0(gene,collapse = ",")) %>% ungroup() %>% mutate(drive_gene=str_remove(drive_gene,"^\\,"),drive_gene=str_remove(drive_gene,"\\,$"),drive_gene=str_replace_all(drive_gene,",,*",","))%>% arrange(`q value`)
amp %>% write_excel_csv("amp_drivegene.csv",col_names = T)


del <- read_delim('result/del_genes.conf_90.txt',delim = "\t",col_names = T) %>% select(-X20)
del1 <- del %>% rename(variant=cytoband) %>% slice(1:3) %>% gather(cytoband,value,-variant) %>% spread(variant,value)
del2 <- del %>% rename(variant=cytoband) %>% slice(4:dim(del)[1]) %>% gather(cytoband,value,-variant) %>% select(-variant) %>% filter(!is.na(value)) %>% rename(gene=value)
del <- del1 %>% left_join(del2)
del <- del %>% mutate(gene=str_replace_all(gene,"[\\[\\]]","")) %>% left_join(del %>% count(cytoband))%>% mutate(drive_gene=if_else(gene %in% lung_drivegene$gene,gene,"")) %>% mutate(`q value`=as.numeric(`q value`)) %>% arrange(`q value`,n,cytoband) 
del <- del %>% group_by(cytoband,`q value`,`wide peak boundaries`,n) %>% summarise(drive_gene=paste0(drive_gene,collapse = ","),gene=paste0(gene,collapse = ",")) %>% ungroup() %>% mutate(drive_gene=str_remove(drive_gene,"^\\,"),drive_gene=str_remove(drive_gene,"\\,$"),drive_gene=str_replace_all(drive_gene,",,*",",")) %>% arrange(`q value`)
del %>% write_excel_csv("del_drivegene.csv",col_names = T)




# combind cnv data and remove small peaks ------------------------------------------------------
focal_cnv <- left_join(
  gcp$g.lin %>% mutate(Wide_Peak=paste0('chr',Chromosome,':',Start_Position,'-',End_Position)),
  bind_rows(amp,del) %>% mutate(drive_gene=str_remove(drive_gene,",$")) %>% select(Cytoband=cytoband,Wide_Peak=`wide peak boundaries`,TotalGene=n,Drive_Gene=drive_gene)
) %>% mutate(PeakWidth=End_Position-Start_Position+1) 

flanklen <- 500000
focal_remove <- focal_cnv %>% 
  filter(PeakWidth<50000) %>% 
  select(chrom=Chromosome,start=Start_Position,end=End_Position) %>% 
  mutate(start=start-flanklen,end=end+flanklen) %>% 
  bind_rows(
    cyto %>% 
      filter(Centro!=0,Chromosome %in% c(1:22)) %>% 
      select(chrom=Chromosome,start=Start,end=End)
  ) 

scorenames <- colnames(gcp$gis.scores )
gcp$gis.scores <- gcp$gis.scores %>% 
  select(chrom=Chromosome,start=Start_Position,end=End_Position,everything()) %>% 
  mutate(Variant_Classification=as.character(Variant_Classification),chrom=as.character(chrom)) %>% 
  bed_intersect(focal_remove,invert = T) %>% 
  select(Chromosome=chrom,Start_Position=start,End_Position=end,everything()) %>% 
  select(one_of(scorenames))

scorenames <- colnames(focal_cnv)
focal_cnv <- focal_cnv %>% 
  select(chrom=Chromosome,start=Start_Position,end=End_Position,everything()) %>% 
  mutate(Variant_Classification=as.character(Variant_Classification),chrom=as.character(chrom)) %>% 
  bed_intersect(focal_remove,invert = T) %>% 
  select(Chromosome=chrom,Start_Position=start,End_Position=end,everything()) %>% 
  select(one_of(scorenames))




# Gistic plot -------------------------------------------------------------

# hg19centro <- hg19centro %>% 
#   left_join(as_tibble(gcp$chr.tbl) %>% select(Chromosome=chr,start)) %>% 
#   mutate(Centro=Pos+start-1) %>% 
#   select(Chromosome,Centro)

chrlevels <- as.character(c(1:22))

hg19backgroup <- 
  cyto %>% group_by(Chromosome) %>% summarise(Start=min(Start)+1,End=max(End)) %>% 
  filter(Chromosome %in% c(1:22)) %>% mutate(Chromosome=factor(Chromosome,levels = chrlevels))

# Amplificaiton -----------------------------------------------------------

ampdata <- bind_rows(
  gcp$gis.scores %>% filter(Variant_Classification=="Amp") %>% select(Chromosome,pos=Start_Position,G_Score,fdr),
  gcp$gis.scores %>% filter(Variant_Classification=="Amp") %>% select(Chromosome,pos=End_Position ,G_Score,fdr),
  hg19centro %>% mutate(fdr=0,G_Score=0) %>% select(Chromosome,pos=Centro,fdr,G_Score),
  hg19backgroup %>% mutate(fdr=0,G_Score=0) %>%  select(Chromosome,pos=Start,fdr,G_Score),
  hg19backgroup %>% mutate(fdr=0,G_Score=0) %>%  select(Chromosome,pos=End,fdr,G_Score)
) %>% arrange(pos) %>% 
  filter(Chromosome %in% c(1:22)) 

fdrmax <- max(ampdata$fdr)
ampdata_limit <- ampdata %>% group_by(Chromosome) %>% summarise(min=min(pos),max=max(pos)) %>% mutate(fdr_max=25,pos=NA_integer_) %>% mutate(Chromosome=factor(Chromosome,levels = chrlevels)) %>% arrange(Chromosome)

ampdata_limit$col <- rep(c('white','gray90'),11)

secbreacks <- NULL
secbreacks[1] <- gcp$gis.scores %>% filter(Variant_Classification=="Amp") %>% arrange(G_Score) %>% filter(G_Score>=0.9) %>% slice(1) %>% pull(fdr)
secbreacks[2] <- gcp$gis.scores %>% filter(Variant_Classification=="Amp") %>% arrange(G_Score) %>% filter(G_Score>=0.8) %>% slice(1) %>% pull(fdr)
secbreacks[3] <-gcp$gis.scores %>% filter(Variant_Classification=="Amp") %>% arrange(G_Score) %>% filter(G_Score>=0.6) %>% slice(1) %>% pull(fdr)
secbreacks[4] <-gcp$gis.scores %>% filter(Variant_Classification=="Amp") %>% arrange(G_Score) %>% filter(G_Score>=0.4) %>% slice(1) %>% pull(fdr)
secbreacks[5] <-gcp$gis.scores %>% filter(Variant_Classification=="Amp") %>% arrange(G_Score) %>% filter(G_Score>=0.2) %>% slice(1) %>% pull(fdr)
secbreacks[6] <- -log10(0.25)
names(secbreacks) <- c(0.9,0.8,0.6,0.4,0.2,0)
names(secbreacks)[6] <- gcp$gis.scores %>% filter(Variant_Classification=="Amp") %>% arrange(fdr) %>% filter(fdr>= -log10(0.25)) %>% slice(1) %>% pull(G_Score) %>% round(digits = 2)

p_amp <- ampdata %>% 
  ggplot()+
  geom_rect(data=ampdata_limit,aes(xmin=min,xmax=max,ymin=0,ymax=64),fill=ampdata_limit$col)+
  geom_step(aes(pos,fdr),size=0.5,colour="red",na.rm = TRUE,direction = 'vh')+
  facet_grid(~factor(Chromosome,levels = chrlevels),scales = "free_x",space = 'free')+
  geom_vline(data=hg19centro,aes(xintercept=Centro),size=0.2,colour="gray20",linetype=4)+
  geom_hline(data=hg19centro %>% mutate(fdrcutoff=2),aes(yintercept=fdrcutoff),size=0.2,colour="green4",linetype=1)+
  scale_x_continuous(breaks = pretty_breaks(),expand = c(0,0),sec.axis = dup_axis())+
  scale_y_continuous(breaks = c(1,2,4,8,16,32),expand = expand_scale(add = c(0, 0)), sec.axis = dup_axis(breaks = secbreacks,labels = names(secbreacks),name = 'G-Score'),labels = math_format(10^-.x),trans='log2')+
  labs(x="",y="q-value")+
  theme(text = element_text(family = "Roboto Condensed"),panel.spacing = unit(0, "lines"),panel.grid.major.y = element_line(colour = 'gray50',linetype = 5),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 14),strip.placement = "inside",strip.switch.pad.wrap = unit(0, "cm"),strip.text = element_blank(),strip.background = element_blank(),axis.line= element_line(colour = 'black'),axis.title.y.right = element_text(size=14),axis.text.y.right = element_text(size = 12),axis.title.x.top = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank() )+
  coord_cartesian(ylim = c(0.35,64))+panel_border(colour = 'gray80',size = 0.2)


focal_cnv_amp <- focal_cnv %>% filter(Variant_Classification=="Amp") %>% mutate(pos=(Start_Position+End_Position)/2,fdr=-log10(qvalues))%>% select(Chromosome,pos,fdr,Cytoband,Drive_Gene,TotalGene) %>% mutate(gene=Cytoband) %>% mutate(Chromosome=factor(Chromosome,levels = chrlevels))

focal_cnv_amp$gene[focal_cnv_amp$Cytoband=="12q15"]="\n12q15/MDM2"
focal_cnv_amp$gene[focal_cnv_amp$Cytoband=="5p15.33"]="5p15.33/TERT"
focal_cnv_amp$gene[focal_cnv_amp$Cytoband=="8q24.21"]="8q24.21/MYC\n"
focal_cnv_amp$gene[focal_cnv_amp$Cytoband=="14q13.3"]="14q13.3/NKX2-1"
focal_cnv_amp$gene[focal_cnv_amp$Cytoband=="16p13.2"]="16p13.2/GRIN2A"
focal_cnv_amp$gene[focal_cnv_amp$Cytoband=="1q21.3"]="1q21.3/MCL1"
focal_cnv_amp$gene[focal_cnv_amp$Cytoband=="7q36.1"]="7q36.1/BRAF"
#focal_cnv_amp$gene[focal_cnv_amp$Cytoband=="16q12.1"]="16q12.1/N4BP1"

focal_cnv_amp$gene[focal_cnv_amp$Cytoband=="20q11.23"]="20q11.23"
focal_cnv_amp$gene[focal_cnv_amp$Cytoband=="20q13.33"]="20q13.33  " 

focal_cnv_amp$gene[focal_cnv_amp$Cytoband=="1p34.2"]=""
focal_cnv_amp$gene[focal_cnv_amp$Cytoband=="3q26.2"]=""
focal_cnv_amp$gene[focal_cnv_amp$Cytoband=="13q33.3"]=""
focal_cnv_amp$gene[focal_cnv_amp$Cytoband=="16q12.1"]=""
focal_cnv_amp$gene[focal_cnv_amp$Cytoband=="17q24.3"]=""
focal_cnv_amp <- bind_rows(
  focal_cnv_amp,
  tibble(Chromosome=factor(7,levels = chrlevels),pos=55086725,fdr=2,gene="7p11.2/EGFR")
)






p_amp <- p_amp +geom_text_repel(
  data = focal_cnv_amp,
  aes(x=pos,y=fdr,label=gene),
  nudge_y      = 2,
  direction    = "x",
  angle        = 90,
  vjust        = 0.5,
  segment.size = 0,
  segment.color = 'black',
  force = 10,
  size=3,
  fontface="bold",
  family="Roboto Condensed"
)

#p_amp


# Deletion ----------------------------------------------------------------

deldata <- bind_rows(
  gcp$gis.scores %>% filter(Variant_Classification=="Del") %>% select(Chromosome,pos=Start_Position,G_Score,fdr),
  gcp$gis.scores %>% filter(Variant_Classification=="Del") %>% select(Chromosome,pos=End_Position ,G_Score,fdr),
  hg19centro %>% mutate(fdr=0,G_Score=0) %>% select(Chromosome,pos=Centro,fdr,G_Score),
  hg19backgroup %>% mutate(fdr=0,G_Score=0) %>%  select(Chromosome,pos=Start,fdr,G_Score),
  hg19backgroup %>% mutate(fdr=0,G_Score=0) %>%  select(Chromosome,pos=End,fdr,G_Score)
) %>% arrange(pos) %>% 
  filter(Chromosome %in% c(1:22)) 

fdrmax <- max(deldata$fdr)
deldata_limit <- deldata %>% group_by(Chromosome) %>% summarise(min=min(pos),max=max(pos)) %>% mutate(fdr_max=35,pos=NA_integer_) %>% mutate(Chromosome=factor(Chromosome,levels = chrlevels)) %>% arrange(Chromosome)

ampdata_limit$col <- rep(c('white','gray90'),11)

secbreacks <- NULL
secbreacks[1] <- gcp$gis.scores %>% filter(Variant_Classification=="Del") %>% arrange(G_Score) %>% filter(G_Score>=0.9) %>% slice(1) %>% pull(fdr)
secbreacks[2] <- gcp$gis.scores %>% filter(Variant_Classification=="Del") %>% arrange(G_Score) %>% filter(G_Score>=0.8) %>% slice(1) %>% pull(fdr)
secbreacks[3] <-gcp$gis.scores %>% filter(Variant_Classification=="Del") %>% arrange(G_Score) %>% filter(G_Score>=0.6) %>% slice(1) %>% pull(fdr)
secbreacks[4] <-gcp$gis.scores %>% filter(Variant_Classification=="Del") %>% arrange(G_Score) %>% filter(G_Score>=0.4) %>% slice(1) %>% pull(fdr)
secbreacks[5] <-gcp$gis.scores %>% filter(Variant_Classification=="Del") %>% arrange(G_Score) %>% filter(G_Score>=0.2) %>% slice(1) %>% pull(fdr)
secbreacks[6] <- -log10(0.25)
names(secbreacks) <- c(0.9,0.8,0.6,0.4,0.2,0)
names(secbreacks)[6] <- gcp$gis.scores %>% filter(Variant_Classification=="Del") %>% arrange(fdr) %>% filter(fdr>=-log10(0.25)) %>% slice(1) %>% pull(G_Score) %>% round(digits = 2)



p_del <- deldata %>% 
  ggplot()+
  geom_rect(data=deldata_limit,aes(xmin=min,xmax=max,ymin=0,ymax=64),fill=ampdata_limit$col)+
  geom_step(aes(pos,fdr),size=0.5,colour="blue",na.rm = TRUE,direction = 'vh')+
  facet_grid(~factor(Chromosome,levels = chrlevels),scales = "free_x",space = 'free')+
  geom_vline(data=hg19centro,aes(xintercept=Centro),size=0.2,colour="gray20",linetype=4)+
  geom_hline(data=hg19centro %>% mutate(fdrcutoff=2),aes(yintercept=fdrcutoff),size=0.2,colour="green4",linetype=1)+
  scale_x_continuous(breaks = pretty_breaks(),expand = c(0,0),sec.axis = dup_axis())+
  scale_y_continuous(breaks = c(1,2,4,8,16,32),expand = expand_scale(add = c(0, 0)), sec.axis = dup_axis(breaks = secbreacks,labels = names(secbreacks),name = 'G-Score'),labels = math_format(10^-.x),trans=reverselog_trans(2))+
  labs(x="",y="q-value")+
  theme(text = element_text(family = "Roboto Condensed"),panel.spacing = unit(0, "lines"),panel.grid.major.y = element_line(colour = 'gray50',linetype = 5),panel.background = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 14),strip.placement = "inside",strip.switch.pad.wrap = unit(0, "cm"),strip.text = element_blank(),strip.background = element_blank(),axis.line= element_line(colour = 'black'),axis.title.y.right = element_text(size=14),axis.text.y.right = element_text(size = 12),axis.title.x.top = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank() )+
  coord_cartesian(ylim = c(0.2,64))+panel_border(colour = 'gray80',size = 0.2)


focal_cnv_del <- focal_cnv %>% filter(Variant_Classification=="Del") %>% mutate(pos=(Start_Position+End_Position)/2,fdr=-log10(qvalues))%>% select(Chromosome,pos,fdr,Cytoband,Drive_Gene,TotalGene) %>% mutate(gene=Cytoband) %>% mutate(Chromosome=factor(Chromosome,levels = chrlevels))
#focal_cnv_del %>% View()

focal_cnv_del$gene[focal_cnv_del$Cytoband=="9p21.3"]="\n9p21.3/CDKN2A"
focal_cnv_del$gene[focal_cnv_del$Cytoband=="8p23.3"]="8p23.3/CSMD1"
focal_cnv_del$gene[focal_cnv_del$Cytoband=="19p13.3"]="19p13.3/STK11\n"

focal_cnv_del$gene[focal_cnv_del$Cytoband=="3q27.3"]=""

focal_cnv_del <- bind_rows(
  focal_cnv_del,
  tibble(Chromosome=factor(12,levels = chrlevels),pos=6055366,fdr=3,gene="12p13.1")
)

p_del <- p_del +geom_text_repel(
  data = focal_cnv_del,
  aes(x=pos,y=fdr,label=gene),
  nudge_y      = 2,
  direction    = "x",
  angle        = 90,
  vjust        = 3,
  segment.size = 0,
  segment.color = 'black',
  force = 10,
  size=3,
  fontface="bold",
  family="Roboto Condensed"
)

#p_del




# Cytoband plot -----------------------------------------------------------

p_cytoband <- gcp$chr.tbl %>%
  mutate(pos=(start+end)/2) %>% 
  filter(chr %in% c(1:22)) %>% 
  ggplot()+
  geom_rect(aes(xmin=start,xmax=end,ymin=0.2,ymax=0.8,fill=color),color="black")+
  geom_text(aes(x=pos,y = 0.5,label=chr),size=3,fontface = "bold",family="Roboto Condensed")+
  theme_void()+
  scale_fill_manual(values = c('white','gray50'))+
  theme(legend.position = 'none',text = element_text(family = "Roboto Condensed"))+
  scale_x_continuous(breaks = pretty_breaks(),expand = c(0,0))




# Combined Plots ----------------------------------------------------------
pdfhr()
p_combined_plot <- plot_grid(
  p_amp+theme(plot.margin=margin(b=-0.45,unit="cm")),
  p_cytoband+theme(plot.margin=margin(b=0,unit="cm")),
  p_del+theme(plot.margin=margin(t=0,unit="cm")),
  align = 'v',axis = "lr",rel_heights = c(2,0.1,2),ncol = 1)
ggsave(file="Gistic2_ztw.pdf",plot = p_combined_plot,width = 12,height = 6,device = cairo_pdf)




#### compare with lung ADC
library(readxl)
library(ggrepel)
library(hrbrthemes)
library(ggsci)
amp_adc <- read_excel('~/NIH-Work/EAGLE_NSLC/LUAD_paper/ng.3564-S2.xlsx',sheet = 17,skip = 2)
amp_adc <- left_join(
  amp_adc %>% slice(c(4)) %>% select(-cytoband) %>% gather(cytoband,gene),
  amp_adc %>% slice(c(1)) %>% select(-cytoband) %>% gather(cytoband,q2)
)

# compare with lung ADC ---------------------------------------------------
amp <- amp %>% filter(`wide peak boundaries` %in% focal_cnv$Wide_Peak)


ampall <- amp %>% select(cytoband,q1=`q value`) %>% full_join(amp_adc) %>%
  mutate(q2=as.numeric(q2)) %>% filter(q1<0.25|q2<0.25) %>% 
  mutate(gene=if_else(is.na(gene),cytoband,gene)) %>% 
  mutate(q1=if_else(is.na(q1),1,q1),q2=if_else(is.na(q2),1,q2))

ampall2 <- ampall %>% filter((q1<q2 & q1<10^-3)|q1<10^-10|q2<10^-10)

ampall %>% ggplot(aes(-log10(q1),-log10(q2)))+geom_point(color="black",pch=21,fill="gray20")+geom_label_repel(aes(label=gene,col=q1<q2),data = ampall2,label.size = NA,fill=NA)+geom_hline(yintercept = -log10(0.25),col="red",linetype=2)+geom_vline(xintercept = -log10(0.25),col="red",linetype=2)+theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,axis = "XYxy",axis_col = "black")+xlab("NCI Sherlock Lung NSLC q value (log2(-log10))")+ylab("Lung ADC q value log2((-log10))")+guides(color=guide_legend(title="q(NSLC)<q(ADC)"))+scale_color_aaas()+scale_color_aaas()+scale_x_continuous(trans='log2')+scale_y_continuous(trans='log2')+coord_cartesian(clip = 'off')


ggsave('gistic2_amp_comparsion.pdf',width = 8,height = 6,device = cairo_pdf)


del <- del %>% filter(`wide peak boundaries` %in% focal_cnv$Wide_Peak)

del_adc <- read_excel('~/NIH-Work/EAGLE_NSLC/LUAD_paper/ng.3564-S2.xlsx',sheet = 20,skip = 2)
del_adc <- left_join(
  del_adc %>% slice(c(4)) %>% select(-cytoband) %>% gather(cytoband,gene),
  del_adc %>% slice(c(1)) %>% select(-cytoband) %>% gather(cytoband,q2)
)

del_adc <- del_adc %>% mutate(cytoband=str_remove(cytoband,"\\.\\..*")) %>% group_by(cytoband) %>% arrange(gene) %>% slice(1) %>% ungroup()


delall <- del %>% select(cytoband,q1=`q value`) %>% full_join(del_adc) %>%
  mutate(q2=as.numeric(q2)) %>% filter(q1<0.25|q2<0.25) %>% 
  mutate(gene=if_else(is.na(gene),cytoband,gene)) %>% 
  mutate(q1=if_else(is.na(q1),1,q1),q2=if_else(is.na(q2),1,q2))

delall2 <- delall %>% filter((q1<q2 & q1<10^-3)|q1<10^-10|q2<10^-10)

delall %>% ggplot(aes(-log10(q1),-log10(q2)))+geom_point(color="black",pch=21,fill="gray20")+geom_label_repel(aes(label=gene,col=q1<q2),data = delall2,label.size = NA,fill=NA)+geom_hline(yintercept = -log10(0.25),col="red",linetype=2)+geom_vline(xintercept = -log10(0.25),col="red",linetype=2)+theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,axis = "xy",axis_col = "black")+xlab("NCI Sherlock Lung NSLC q value (log2(-log10))")+ylab("Lung ADC q value (log2(-log10))")+guides(color=guide_legend(title="q(NSLC)<q(ADC)"))+scale_color_aaas()+scale_x_continuous(trans='log2')+scale_y_continuous(trans='log2')+coord_cartesian(clip = 'off')

#delall %>% ggplot(aes(log2(-log10(q1)+1),log2(-log10(q2)+1)))+geom_point(color="black",pch=21,fill="gray20")+geom_label_repel(aes(label=gene,col=q1<q2),data = delall2,label.size = NA,fill=NA)+geom_hline(yintercept = log2(-log10(0.25)+1),col="red",linetype=2)+geom_vline(xintercept = log2(-log10(0.25)+1),col="red",linetype=2)+theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,axis = "xy",axis_col = "black")+xlab("NCI Sherlock Lung NSLC q value (log2(-log10))")+ylab("Lung ADC q value (log2(-log10))")+guides(color=guide_legend(title="q(NSLC)<q(ADC)"))+scale_color_aaas()

ggsave('gistic2_del_comparsion.pdf',width = 8,height = 6,device = cairo_pdf)



save.image('Gistic2.RData')


# gene level  -------------------------------------------------------------

genedata <- read_delim('result/all_thresholded.by_genes.txt',delim = '\t',col_names = T)

genedata <- genedata %>% 
  select(GeneSymbol=`Gene Symbol`,Cytoband,contains('NSLC')) %>% 
  pivot_longer(cols = -c(GeneSymbol,Cytoband),names_to = "Tumor_Barcode",values_to = 'CNV') 

save(genedata,file = 'all_thresholded.by_genes.RData')





# Function ----------------------------------------------------------------



reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}






GistChrom_ztw <-  function (gistic = NULL, fdrCutOff = 0.1, markBands = NULL, color = NULL, ref.build = "hg19", cytobandOffset = 0.01,  txtSize = 0.8, cytobandTxtSize = 0.6)  {
  g = getCytobandSummary(gistic)
  g = g[qvalues < fdrCutOff]
  g[, `:=`(Chromosome, sapply(strsplit(x = g$Wide_Peak_Limits, split = ":"), "[", 1))]
  g[, `:=`(loc, sapply(strsplit(x = g$Wide_Peak_Limits, split = ":"), "[", 2))]
  g[, `:=`(Start_Position, sapply(strsplit(x = g$loc, split = "-"),  "[", 1))]
  g[, `:=`(End_Position, sapply(strsplit(x = g$loc, split = "-"),  "[", 2))]
  g.lin = maftools:::transformSegments(segmentedData = g[, .(Chromosome,  Start_Position, End_Position, qvalues, Cytoband, Variant_Classification)])
  if (is.null(color)) {
    color = c(Amp = "red", Del = "blue")
  }
  gis.scores = maftools:::transformSegments(segmentedData = gistic@gis.scores,  build = ref.build)
  gis.scores$amp = ifelse(test = gis.scores$Variant_Classification ==  "Del", yes = -gis.scores$G_Score, no = gis.scores$G_Score)
  gis.scores$ystart = ifelse(test = gis.scores$Variant_Classification ==  "Del", yes = -cytobandOffset, no = cytobandOffset)
  fdrCutOff = -log10(fdrCutOff)
  gis.scores$Variant_Classification = ifelse(test = as.numeric(gis.scores$fdr) >  fdrCutOff, yes = gis.scores$Variant_Classification,  no = "neutral")
  gis.scores$Variant_Classification = factor(gis.scores$Variant_Classification,  levels = c("neutral", "Amp", "Del"))
  if (ref.build == "hg19") {
    chr.lens = c(249250621, 243199373, 198022430, 191154276, 
                 180915260, 171115067, 159138663, 146364022, 141213431, 
                 135534747, 135006516, 133851895, 115169878, 107349540, 
                 102531392, 90354753, 81195210, 78077248, 59128983, 
                 63025520, 48129895, 51304566, 155270560, 59373566)
  }
  else if (ref.build == "hg18") {
    chr.lens = c(247249719, 242951149, 199501827, 191273063, 
                 180857866, 170899992, 158821424, 146274826, 140273252, 
                 135374737, 134452384, 132349534, 114142980, 106368585, 
                 100338915, 88827254, 78774742, 76117153, 63811651, 
                 62435964, 46944323, 49691432, 154913754, 57772954)
  }
  else if (ref.build == "hg38") {
    chr.lens = c(248956422, 242193529, 198295559, 190214555, 
                 181538259, 170805979, 159345973, 145138636, 138394717, 
                 133797422, 135086622, 133275309, 114364328, 107043718, 
                 101991189, 90338345, 83257441, 80373285, 58617616, 
                 64444167, 46709983, 50818468, 156040895, 57227415)
  }
  else {
    stop("ref.build can only be hg18, hg19 or hg38")
  }
  chr.lens.cumsum = cumsum(chr.lens)
  nchrs = length(unique(gis.scores$Chromosome))
  chr.labels = c(1:22, "X", "Y")
  chr.tbl = data.table::data.table(chr = chr.labels, start = c(1,  chr.lens.cumsum[1:length(chr.lens.cumsum) - 1]), end = chr.lens.cumsum)
  chr.tbl$color = rep(c("black", "white"), length = nrow(chr.tbl))
  y_lims = pretty(gis.scores[, amp], na.rm = TRUE)
  gis.scores$Variant_Classification = factor(x = gis.scores$Variant_Classification,  levels = c("neutral", "Amp", "Del"))
  gis.scores = split(gis.scores, as.factor(as.character(gis.scores$Variant_Classification)))
  gis.scores = data.table::rbindlist(l = gis.scores, use.names = TRUE,fill = TRUE)
  return(list(gis.scores=gis.scores,g.lin=g.lin,chr.tbl=chr.tbl))
}








# CNV arm level data ------------------------------------------------------
# keep consistence with bb_heatmap

armdata <- read_delim('result/broad_values_by_arm.txt',delim = '\t',col_names = T)

armdata <- armdata %>% 
  pivot_longer(cols = -`Chromosome Arm`,names_to = 'Tumor_Barcode',values_to = 'tvalue') %>%
  mutate(calling=case_when(
    tvalue < -1.3 ~ "-2",
    tvalue < -0.25 & tvalue >= -1.3 ~ "-1",
    tvalue <= 0.25 & tvalue >= -0.25 ~ NA_character_,
    tvalue >0.25 & tvalue <= 0.9 ~ "1",
    tvalue > 0.9 ~ "2"
  )) %>% 
  #filter(value!=0) %>% 
  mutate(calling=factor(calling,levels = c("-2","-1","2","1")))%>%
  mutate(type=if_else(calling %in% c("1","2"),"Amp",if_else(calling %in% c("-1","-2"),"Del",NA_character_))) %>% 
  select(Tumor_Barcode,Arm=`Chromosome Arm`,everything())

armcut <- read_delim('result/broad_significance_results.txt',delim = '\t',col_names = T)

armcut %>% filter(`Amp q-value`<0.01)
armcut %>% filter(`Del q-value`<0.01)

save(armdata,armcut,file='broad_significance_results.RData')





