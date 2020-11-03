set_wd()
library(tidyverse)
#library(tidylog)
library(CNTools)
library(ggsci)
library(scales)
library(hrbrthemes)
library(cowplot)
library(valr)

load('sherlock_samples_unique.RData')
load('sherlock_wgd.RData')

load('cytoband_hg19.RData')
hg19centro <- cyto %>% filter(Centro!=0,Chromosome %in% c(1:22)) %>% select(chrom=Chromosome,start=Start,end=End,Centro) 



# read subclone data  -----------------------------------------------------

cnv <- read_delim('bb_subclone.txt',delim = '\t',col_names = T) %>% filter(Tumor_Barcode %in% sherlock_samples_unique$Tumor_Barcode)
cnv_lusc <- read_delim('../../FirstPaper/Bin/LUSC/bb_subclone.tsv',delim = '\t',col_names = T)
cnv_luad <- read_delim('../../FirstPaper/Bin/LUAD/bb_subclone.tsv',delim = '\t',col_names = T)

cnv <- bind_rows(
  cnv,
  cnv_luad
) 


cnv <- cnv %>% filter(chr %in% c(1:22,"X","Y"))
cnv <- cnv %>% 
  mutate(frac1_A=if_else(is.na(frac1_A),0,frac1_A)) %>% 
  mutate(frac2_A=if_else(is.na(frac2_A),0,frac2_A)) %>% 
  mutate(
    clone_frac=if_else(frac1_A>frac2_A,frac1_A,frac2_A),
    clone_nMaj=if_else(frac1_A>frac2_A,nMaj1_A,nMaj2_A),
    clone_nMin=if_else(frac1_A>frac2_A,nMin1_A,nMin2_A),
    subclone_frac=if_else(frac1_A<frac2_A,frac1_A,frac2_A),
    subclone_nMaj=if_else(frac1_A<frac2_A,nMaj1_A,nMaj2_A),
    subclone_nMin=if_else(frac1_A<frac2_A,nMin1_A,nMin2_A)
  ) %>% 
  select(Tumor_Barcode,chr:ntot,contains('clone'))

sherlock_bb_segmentation <- cnv 

# phenotype data 
pcawg_sample <- read_delim('~/NIH-Work/MutationSignature/mSigPortal/CBIIT/reference_paper_data/PCAWG/sample.all_projects.tsv',delim = '\t')
sbs <- read_csv('~/NIH-Work/MutationSignature/mSigPortal/CBIIT/reference_paper_data/The repertoire of mutational signatures in human cancer/Signature_in_Samples/PCAWG_sigProfiler_SBS_signatures_in_samples.csv',col_names = T)
sbs <- sbs %>% select(icgc_specimen_id=`Sample Names`,SBS4) %>% mutate(SBS4=if_else(SBS4>0,"Y","N")) %>% 
  left_join(pcawg_sample %>% select(icgc_specimen_id,icgc_sample_id))
pcawg <- read_delim('../../FirstPaper/Bin/icgc_sample_annotations_summary_table_LUNG.txt',delim = '\t')
pcawg <- pcawg %>% left_join(sbs)%>%
  #mutate(SBS4=if_else(is.na(SBS4),"N",SBS4)) %>% 
  select(Tumor_Barcode=tumour_aliquot_id,Study=histology_abbreviation,WGD_Status=wgd_status,Tumor_Purity=purity,SBS4) %>% 
  mutate(WGD_Status=if_else(WGD_Status=="wgd","Y","N"))

phodata <- sherlock_wgd %>% filter(Tumor_Barcode %in% sherlock_samples_unique$Tumor_Barcode) %>% 
  mutate(Study="Sherlock-Lung") %>% left_join(sherlock_samples_unique %>% select(Tumor_Barcode,Tumor_Purity)) %>% 
  select(Tumor_Barcode,Study,WGD_Status,Tumor_Purity) %>% 
  mutate(SBS4="N") %>% 
  bind_rows(pcawg)

# save result for SCNA signature ------------------------------------------
#cnv %>% 
#  select(Tumor_Barcode,chr,startpos,endpos,BAF,pval,LogR,clone_frac,nMajor=clone_nMaj,nMinor=clone_nMin) %>% 
#  write_delim(path = 'bb_signature_input_allsegments.txt',delim = '\t',na = '',col_names = T)

#cnv %>% 
#  filter(clone_frac==1) %>% 
#  select(Tumor_Barcode,chr,startpos,endpos,BAF,pval,LogR,clone_frac,nMajor=clone_nMaj,nMinor=clone_nMin) %>% 
#  write_delim(path = 'bb_signature_input_clonalsegments.txt',delim = '\t',na = '',col_names = T)


# library(biovizBase)    
# cytoband_color <- getBioColor("CYTOBAND") 
# #cytoband_color[1]<- 'grey80'
# p1 <- cytoGeno %>% filter(genome=='hg19') %>% 
#   mutate(chrom=str_remove(chrom,"chr")) %>% 
#   filter(chrom %in% c(1:22)) %>% 
#   ggplot(aes(fill=gieStain))+
#   geom_rect(aes(xmin=chromStart,xmax=chromEnd,ymin=0,ymax=1))+
#   facet_wrap(~chrom,nrow = 1,scales = 'free_x')+
#   scale_y_continuous(expand = c(0,0))+
#   scale_fill_manual(values = cytoband_color )+
#   theme_void()+
#   theme(panel.spacing = unit(0, "lines"),legend.position = "none",strip.text = element_blank())+
#   cowplot::panel_border(color="black")



# Define the value --------------------------------------------------------
totalsamples <- cnv %>% count(Tumor_Barcode) %>% pull(Tumor_Barcode) %>% length()
allsamples <- cnv %>% count(Tumor_Barcode) %>% pull(Tumor_Barcode) %>% unique()

# #version1 
# #filter(!is.na(clone_frac),clone_frac<1)
# cnvdata <- cnv %>%
#   mutate(clone_total=clone_nMin+clone_nMaj) %>% 
#   left_join(
#     sherlock_samples %>% select(Tumor_Barcode,Tumor_Purity,Tumor_Ploidy)
#   ) %>% 
#   left_join(
#     sherlock_wgd %>% select(Tumor_Barcode,WGD_Status)
#   ) %>% 
#   #mutate(clone_total=if_else(clone_total==0,0.4,clone_total)) %>% 
#   #mutate(clone_total=if_else(clone_total>5,5,clone_total)) %>% 
#   mutate(relative_copy2=clone_total-if_else(WGD_Status=="Y",4,2)) %>% 
#   mutate(relative_copy=log2(clone_total/Tumor_Ploidy)) %>% 
#   mutate(relative_copy=if_else(is.infinite(relative_copy)|relative_copy< -2,-2,relative_copy)) %>% 
#   mutate(relative_copy=if_else(relative_copy>2,2,relative_copy)) %>% 
#   mutate(relative_copy=if_else(relative_copy> -0.2 & relative_copy < 0.2,0,relative_copy))



# for sherlock version 2
# force copy neutral status for chr19HD region. 
load('sherlock_cnvinfo.RData')
chr19list <- sherlock_cnvinfo %>% filter(Chr19loss_Clust=="Y") %>% pull(Tumor_Barcode)
cnv <- cnv %>% mutate(clone_nMaj=if_else(Tumor_Barcode %in% chr19list,1,clone_nMaj),clone_nMin=if_else(Tumor_Barcode %in% chr19list,1,clone_nMin))

#version2
#filter(!is.na(clone_frac),clone_frac<1)
cnvdata <- cnv %>%
  filter(chr %in% c(1:22)) %>% 
  mutate(clone_total=clone_nMin+clone_nMaj) %>% 
  left_join( phodata) %>% 
  mutate(relative_copy=clone_total-if_else(WGD_Status=="Y",4,2)) %>% 
  mutate(relative_copy=if_else(relative_copy>4,4,relative_copy)) %>% 
  mutate(relative_copy=if_else(relative_copy< -4,-4,relative_copy)) %>% 
  mutate(relative_copy=if_else(clone_nMaj==0,-4,relative_copy))


# Dendrograms -------------------------------------------------------------
source('ggdendro.R')
segdata <- cnvdata %>% 
  mutate(num.mark=1000) %>% 
  select(ID=Tumor_Barcode,chrom=chr,loc.start=startpos,loc.end=endpos,num.mark,seg.mean=relative_copy) %>%
  mutate(seg.mean=if_else(is.na(seg.mean),0,seg.mean)) %>% 
  as.data.frame()
segdata %>% filter_all(any_vars(is.na(.)))
cnseg <- CNSeg(segdata)
rdseg <- getRS(cnseg, by = "region", imput = FALSE, XY = FALSE, what = "mean")
reducedseg <- rs(rdseg)

msegdata <- as.matrix(reducedseg[,-(1:3)])

msegdata <- apply(msegdata, 2, as.numeric)

### add samples without any change
if(dim(msegdata)[2]<totalsamples){
  numbr <- dim(msegdata)[1]
  numbc <- totalsamples-dim(msegdata)[2]
  cnvtmp <- matrix(rep(0,numbr*numbc),nrow = numbr,ncol =numbc )
  extrasamples <- allsamples[!(allsamples %in% colnames(msegdata))]
  colnames(cnvtmp) <- extrasamples
  msegdata <- cbind(msegdata,cnvtmp)
  
  cnvdata <- bind_rows(
    cnvdata,
    tibble(Tumor_Barcode=rep(extrasamples,each=22),chr=rep(1:22,numbc)) %>% mutate(chr=as.character(chr))
  )
  #change cnv data 
}


hc <- hclust(dist(t(msegdata),method = 'euclidean'),method = 'ward.D')
hcdata <- dendro_data_k(hc, 3)

cols <- c("#a9a9a9", "#1f77b4", "#ff7f0e", "#2ca02c")

p_dend <- plot_ggdendro(hcdata,
                        direction   = "lr",
                        scale.color = cols,
                        label.size  = 2.5,
                        branch.size = 0.5,
                        expand.y    = 0,nudge.label = 0)

p_dend <- p_dend+
  scale_x_continuous(breaks = NULL,limits = c(1, nrow(label(hcdata))),expand = c(0,0))+
  theme_void()

sample_order <- as_tibble(hcdata$labels) %>% select(Tumor_Barcode=label) %>% mutate(Seq=seq_along(Tumor_Barcode)) 

#chr19loss <- sample_order %>% tail(n=11) %>% mutate(chr19loss="Y") %>% select(-Seq)
#chr19loss <- cnvdata %>% filter(chr=="19",clone_nMaj==0) %>% mutate(size=endpos-startpos) %>% filter(size>30000000) %>% mutate(chr19loss="Y") %>% select(Tumor_Barcode,chr19loss)

cnvclust <- as_tibble(hcdata$labels) %>% select(Tumor_Barcode=label,CNV_Clust=clust) %>% mutate(CNV_Clust=paste0('NC',CNV_Clust)) %>% left_join(phodata)

cnvclust %>% write_csv('tmp_cluster.csv',col_names = T)

#save(cnvclust,file='cnvclust_LUAD.RData')

cnvclust <- cnvclust %>% filter(Study=="Sherlock-Lung") %>% select(Tumor_Barcode,CNV_Clust) %>% left_join(tibble(CNV_Clust2=c("C1","C2","C3"),CNV_Clust=c("NC3","NC2","NC1"))) %>% select(-CNV_Clust) %>% rename(CNV_Clust=CNV_Clust2)
#save(cnvclust,file='cnvclust.RData')

# Heatmap plot ------------------------------------------------------------

chrlevels <- seq(1:22)
p_heatmap <- cnvdata %>% 
  mutate(chr=factor(chr,levels = chrlevels)) %>% 
  left_join(sample_order) %>% 
  arrange(Seq,chr) %>% 
  #filter(Seq %in% c(1:20)) %>% 
  ggplot(aes(fill=relative_copy))+
  #geom_segment(aes(x=startpos,xend=endpos,y=Tumor_Barcode,yend=Tumor_Barcode),size=3)+
  geom_rect(aes(x=NULL,xmin=startpos,xmax=endpos,ymin=Seq-0.5,ymax=Seq+0.5))+
  scale_fill_gradient2(low = "blue",high = "red",mid = "white",midpoint = 0)+
  scale_x_continuous(breaks = pretty_breaks(),expand = c(0,0),sec.axis = dup_axis())+
  labs(fill="Relative Copy Number")+
  scale_y_continuous(breaks = sample_order$Seq,labels = sample_order$Tumor_Barcode,expand = expand_scale(add = c(0, 0)), sec.axis = dup_axis())+
  #facet_wrap(~chr,nrow = 1,scales = 'free_x',strip.position = "bottom")+
  facet_grid(~chr,scales = 'free_x',space = 'free',switch="both")+
  labs(x="",y="")+
  #theme_ipsum_rc()+
  theme(text = element_text(family = "Roboto Condensed"),panel.spacing = unit(0, "lines"),axis.text.x = element_blank(),axis.ticks = element_blank(),axis.text.y = element_text(size = 6),strip.placement = "outside",strip.switch.pad.wrap = unit(0, "cm"),strip.text = element_text(face = "bold",size=10),strip.background = element_rect(color = "gray50",fill="white"),panel.background = element_rect(fill = 'white'),axis.line= element_line(colour = 'black'),axis.title.y.right = element_blank(),axis.text.y.right = element_blank(),axis.ticks.y.right = element_blank(),axis.title.x.top = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank() )+
  coord_cartesian(clip="off")
#  coord_fixed(ratio = 1)
#panel_border(color = "black") 

#panel.background = element_rect(fill = 'green')
## the following code to change the color for strip
# g <- ggplot_gtable(ggplot_build(p_heatmap))
# stripr <- which(grepl('strip-b', g$layout$name))
# fills <- c(rep(c("gray50","gray80"),11))
# k <- 1
# for (i in stripr) {
#   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }
# p_heatmap <- ggdraw(g,)
# p_heatmap


p_heatmap_legend <- get_legend(
  # create some space to the left of the legend
  p_heatmap + theme(legend.box.margin = margin(8, 0, 5, 0))
)

p_heatmap <- p_heatmap+theme(legend.position = "none")



# Sample Annotation --------------------------------------------------------------
annodata <- sample_order %>% left_join(phodata)

# purity
annodata_1 <- annodata %>% select(Seq,Tumor_Purity) 
p_a1 <- annodata_1 %>% 
  ggplot(aes("anno1",Seq,fill=Tumor_Purity))+
  geom_tile()+
  scale_fill_viridis_c(option = "D")+
  labs(fill="Purity")+
  theme_void()+
  scale_y_discrete(breaks = NULL,limits = c(1, nrow(label(hcdata))),expand = c(0,0))

p_a1_legend <- get_legend(
  # create some space to the left of the legend
  p_a1 + theme(legend.box.margin = margin(8, 0, 5, 0))
)

p_a1 <- p_a1+theme(legend.position = "none")

# WGD
wgdcolor <- c('Y'='#f03b20','N'='#fee5d9')
annodata_2 <- annodata %>% select(Seq,WGD_Status) 
p_a2 <- annodata_2 %>% 
  ggplot(aes("anno2",Seq,fill=WGD_Status))+
  geom_tile()+
  scale_fill_manual(values = wgdcolor)+
  theme_void()+
  scale_y_continuous(breaks = NULL,limits = c(1, nrow(label(hcdata))),expand = c(0,0))


p_a2_legend <- get_legend(
  # create some space to the left of the legend
  p_a2 + theme(legend.box.margin = margin(8, 0, 5, 0))
)

p_a2 <- p_a2+theme(legend.position = "none")

#study 
annodata_3 <- annodata %>% select(Seq,Study) 
p_a3 <- annodata_3 %>% 
  ggplot(aes("anno3",Seq,fill=Study))+
  geom_tile()+
  scale_fill_manual(values = pal_jama()(3))+
  theme_void()+
  scale_y_continuous(breaks = NULL,limits = c(1, nrow(label(hcdata))),expand = c(0,0))


p_a3_legend <- get_legend(
  # create some space to the left of the legend
  p_a3 + theme(legend.box.margin = margin(8, 0, 5, 0))
)

p_a3 <- p_a3+theme(legend.position = "none")

#SBS4 
annodata_4 <- annodata %>% select(Seq,SBS4) 
p_a4 <- annodata_4 %>% 
  ggplot(aes("anno4",Seq,fill=SBS4))+
  geom_tile()+
  scale_fill_manual(values = pal_aaas()(2))+
  theme_void()+
  scale_y_continuous(breaks = NULL,limits = c(1, nrow(label(hcdata))),expand = c(0,0))


p_a4_legend <- get_legend(
  # create some space to the left of the legend
  p_a4 + theme(legend.box.margin = margin(8, 0, 5, 0))
)

p_a4 <- p_a4+theme(legend.position = "none")


#Histology
load('sherlock_histology.RData')
histcolor <- c('#a6611a','#f1b6da','#d01c8b')
names(histcolor) <- unique(sherlock_histology$Histology)

annodata_5 <- annodata %>% left_join(sherlock_histology) %>% mutate(Histology=if_else(Study=="Lung-AdenoCA","Adenocarcinomas", Histology)) %>% select(Seq,Histology) 
p_a5 <- annodata_5 %>% 
  ggplot(aes("anno4",Seq,fill=Histology))+
  geom_tile()+
  scale_fill_manual(values = histcolor)+
  theme_void()+
  scale_y_continuous(breaks = NULL,limits = c(1, nrow(label(hcdata))),expand = c(0,0))


p_a5_legend <- get_legend(
  # create some space to the left of the legend
  p_a5 + theme(legend.box.margin = margin(8, 0, 5, 0))
)

p_a5 <- p_a5+theme(legend.position = "none")

# Test chr19 loss  ---------------------------------------------------------
#annodata %>% left_join(chr19loss) %>% mutate(chr19loss=if_else(is.na(chr19loss),"N",chr19loss) ) %>% select(RTK_Altered_Status,chr19loss) %>% mutate(chr19loss=factor(chr19loss,levels = c('Y','N'))) %>% table() %>% fisher.test(alternative = 'greater')



# Combined figures --------------------------------------------------------
p_dend <- p_dend+theme(plot.margin=margin(r=0.6,unit="cm"))
p_heatmap <- p_heatmap+theme(plot.margin=margin(l=-1,unit="cm"))+panel_border(colour = 'gray80',size = 0.2)+theme(axis.text.y = element_blank())
p_combined_plot <- plot_grid(p_dend,p_heatmap,p_a1,p_a2,p_a3,p_a4,p_a5,align = 'h',axis = "tb",rel_widths = c(0.5,1.5,0.02,0.02,0.02,0.02,0.02),nrow = 1)
p_combined_legend <- plot_grid(NULL,NULL,NULL,p_heatmap_legend,p_a1_legend,p_a2_legend,p_a3_legend,p_a4_legend,p_a5_legend,NULL,NULL,nrow = 1,align = 'h')
p_combined <- plot_grid(p_combined_plot,p_combined_legend,nrow = 2,align = "v",rel_heights = c(12,1.2))

ggsave(file="Figures/heatmap_all_segements_LUAD_tmp.pdf",plot = p_combined,width = 22,height = 17,device = cairo_pdf)

