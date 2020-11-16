set_wd()
# load library ------------------------------------------------------------
library(tidyverse)
#library(tidylog)
library(scales)
library(ggsci)
library(cowplot)
library(hrbrthemes)
library(data.table)
library(broom)


# load data ---------------------------------------------------------------
load('sherlock_samples.RData')
load('sherlock_samples_unique.RData')
load('Sherlock_histology.RData')

sherlock_samples_unique
#samplesize <- dim(sherlock_samples_unique)[1]


# Color theme -------------------------------------------------------------
get_vcColors = function(alpha = 1){
  col = c(RColorBrewer::brewer.pal(12,name = "Paired"), RColorBrewer::brewer.pal(11,name = "Spectral")[1:3],'black', 'violet', 'royalblue', '#7b7060','#8dd3c7','#01665e','#2b8cbe','#c51b8a','#1c9099','#542788','#006d2c','#e41a1c','#253494', "#c994c7", "#dd1c77","#f03b20","#1F77B4FF", "#FF7F0EFF", "#2CA02CFF",'#993404','#ff7f00','#54278f','#253494','#35978f','#374E55FF')
  col = grDevices::adjustcolor(col = col, alpha.f = alpha)
  names(col) = names = c('Nonstop_Mutation','Frame_Shift_Del','IGR','Missense_Mutation','Silent','Nonsense_Mutation',
                         'RNA','Splice_Site','Intron','Frame_Shift_Ins','Nonstop_Mutation','In_Frame_Del','ITD','In_Frame_Ins',
                         'Translation_Start_Site',"Multi_Hit", 'Amp', 'Del', 'Complex_Event','Promoter','Fusion','INS','DEL','SNV','Chr19_Loss','HLA_LOH',"Amplification","Deletion",'MSI-L','MSI-H','WGD','C1','C2','C3','Kataegis','SV','synonymous_variant','p53_deficiency','RTK-RAS+','LOH')
  #col <- c(col,c('DEL'='red','INS'='blue','SNP'='gray30'))
  col
}

landscape_colors <- get_vcColors()

landscape_colors['In_Frame_Ins'] <- '#df65b0'
landscape_colors['Frame_Shift_Del'] <- '#1F78B4FF'
landscape_colors['Nonsense_Mutation'] <- '#a50f15'


histcolor <- c('#a6611a','#f1b6da','#d01c8b')
names(histcolor) <- unique(sherlock_histology$Histology)

landscape_colors <- c(landscape_colors,histcolor)


#nbreaks=c(0,50,100,150,200,250)
# Define the oncoplot funciton --------------------------------------------------------------
oncoplot <- function(data,data_clone = NULL, data_highlight = NULL, Gene_Sig=NULL,landscape_colors=NULL,gene_level=NULL,sample_level=NULL,sample_level0=sample_level0,GeneSortOnly=FALSE,namemax=NULL,tmar=0,bmar=0,nbreaks=c(0,20,40,60,80),p2_axis_hidden=FALSE,p2_hidden=FALSE,cell_height=0.95,legend_adjust=FALSE){
  
  # limited to sample level
  data <- data %>% filter(Tumor_Barcode %in% sample_level0)
  
  ## add blank gene 
  blankgene <- gene_level[!(gene_level %in% data$Gene)]
  if(length(blankgene)>0){
    n = length(blankgene)
    blankdata <- tibble(Subject=rep(data$Subject[1],n),Tumor_Barcode=rep(data$Tumor_Barcode[1],n),Gene=blankgene,Alteration=rep(NA_character_,n),Type=rep(data$Type[1],n))
    data <- bind_rows(data,blankdata)
  }
  
  samplesize <- length(sample_level0)
  altertype <- data$Type[1]
  #aggresate multiple hits
  # check Generating split-color rectangles from ggplot2 
  #https://stackoverflow.com/questions/22107666/generating-split-color-rectangles-from-ggplot2-geom-raster
  
  data0 <- data %>%
    group_by(Subject,Tumor_Barcode,Gene,Type) %>%
    summarise(Alteration_info=paste0(Alteration,collapse = ',')) %>%
    ungroup() %>%
    mutate(Alteration=if_else(Alteration_info=="NA",NA_character_,if_else(str_detect(Alteration_info,","),'Multi_Hit',Alteration_info))) %>%
    mutate(Alteration_info=if_else(str_detect(Alteration_info,","),Alteration_info,'')) %>%
    select(Subject,Tumor_Barcode, Gene, Type,Alteration, Alteration_info )
  
  data <- data %>% 
    group_by(Subject,Tumor_Barcode,Gene,Type) %>% 
    summarise(Alteration_info=paste0(Alteration,collapse = '/')) %>% 
    ungroup() %>% 
    mutate(Alteration=if_else(Alteration_info=="NA",NA_character_,Alteration_info)) %>% 
    select(Subject,Tumor_Barcode, Gene, Type,Alteration, Alteration_info )
  
  data <- as.data.table(data)
  data <- data[, strsplit(as.character(Alteration), "/"), by=list(Tumor_Barcode, Gene)]  # this expands "X/Y/Z" into three rows
  data[, shift:=(1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(Tumor_Barcode, Gene)]
  data[, height:=cell_height/.N, by=list(Tumor_Barcode, Gene)]
  
  data <- as_tibble(data) %>% rename(Alteration=V1) %>% mutate(shift=shift*cell_height)
  
  
  # define the order #
  Gene_Freq <- data0 %>% select(Gene,Tumor_Barcode) %>% unique() %>% count(Gene,sort=T) %>% mutate(n=if_else(Gene %in% blankgene,0L,n)) %>%  mutate(Freq=n/samplesize) %>% arrange(Freq) %>%  mutate(Freq=percent(Freq,accuracy = 0.1))
  print(Gene_Freq)
  if(is.null(gene_level)){
    gene_level <- Gene_Freq %>% pull(Gene)
  }else{
    Gene_Freq <- Gene_Freq %>% mutate(Gene=factor(Gene,levels = gene_level)) %>% arrange(Gene)
  }
  
  Gene_Freq_alte <- data0 %>% count(Gene,Alteration) %>%  mutate(Gene=factor(Gene,levels = gene_level)) %>% arrange(Gene)
  
  #sample_level <- databg$Tumor_Barcode
  input_sample_level=TRUE
  
  if(is.null(sample_level)){
    input_sample_level=FALSE
    sampleorder <- data0 %>% 
      mutate(Gene=factor(Gene,levels = rev(gene_level))) %>% 
      select(Tumor_Barcode,Gene,Alteration) %>% 
      pivot_wider(id_cols = Tumor_Barcode,names_from = 'Gene',values_from = 'Alteration') %>% 
      arrange_at(vars(one_of(rev(gene_level)))) %>% 
      pull(Tumor_Barcode) %>% 
      as.character()
    
    sample_level <- c(sampleorder,sample_level0[!c(sample_level0 %in% sampleorder)])
  }
  
  ## resort the data 
  databg <- crossing(Tumor_Barcode=sample_level,Gene=gene_level) %>% 
    mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level)) %>% 
    mutate(Gene=factor(Gene,levels = gene_level)) 
  
  data <- data %>% mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level)) %>% 
    mutate(Gene=factor(Gene,levels = gene_level))
  
  data0 <- data0 %>% mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level)) %>% 
    mutate(Gene=factor(Gene,levels = gene_level))
  
  ## same gene name length 
  if(!is.null(namemax)){
    gene_level2 <- str_pad(gene_level,namemax,pad = " ")
  }else{
    gene_level2 <- gene_level
  }
  
  # nameindex <- which.max(str_length(gene_level))
  # namelength <- max(str_length(gene_level))
  # namepad <- namemax-namelength
  # gene_level2 <- gene_level
  # if(namepad>0){
  #   gene_level2[nameindex] <- paste0(str_dup("X", times = namepad),gene_level[nameindex])
  # }
  # 
  
  # add subclone information
  # print(levels(databg$Tumor_Barcode))
  # print(levels(data$Tumor_Barcode))
  # print(levels(data0$Tumor_Barcode))
  
  p1 <- databg %>% 
    ggplot(aes(factor(Tumor_Barcode,levels = sample_level),as.integer(Gene))) +
    geom_tile(height=cell_height,fill="gray95")+
    theme_ipsum_rc(base_size = 10,axis_text_size = 10,grid = '',axis = '')+
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.title.y.right = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(face = "italic"),panel.background = element_blank(),legend.position = 'top',legend.justification = "top",legend.direction = 'vertical')+
    scale_y_continuous(expand = c(0,0),breaks = 1:length(gene_level2),labels = gene_level2,sec.axis = dup_axis(breaks = 1:length(gene_level2),labels = Gene_Freq$Freq))+
    ## data
    geom_tile(data = data,aes(factor(Tumor_Barcode,levels = sample_level),as.integer(Gene)+shift,fill=Alteration,height=height),size=0)+
    scale_fill_manual(values = landscape_colors)+
    geom_vline(xintercept = 1:samplesize-0.5,color="white",size=0.2)+
    panel_border(size = 0.3,color = 'gray70')+
    guides(fill = guide_legend(ncol = 1,title.position = "top",title=altertype))
  
  
  ## highlight multiple samples according to data_highlight
  ## highlight data can included both Tumor_Barcode and Gene, or just Tumor_Barcode
  
  if(!is.null(data_highlight)){
    if(dim(data_highlight)[2]==1){
      data_highlight <- data_highlight %>% 
        left_join(
          tibble(Tumor_Barcode=sample_level) %>% mutate(Seq=seq_along(Tumor_Barcode)) %>% mutate(Seq1=Seq-0.5,Seq2=Seq+0.5)
        )
      
      # merged near-by samples
      data_highlight <- data_highlight %>% arrange(Seq) %>% group_by(g = cumsum(cummax(lag(Seq2, default = first(Seq2))) < Seq1)) %>% 
        summarise(Seq1 = first(Seq1), Seq2 = max(Seq2),Tumor_Barcode=first(Tumor_Barcode)) %>% mutate(Gene=NA) %>% 
        mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level)) %>% 
        mutate(Gene=factor(Gene,levels = gene_level)) 
      
      p1 <- p1+geom_rect(data = data_highlight,aes(xmin=Seq1,xmax=Seq2,ymin=-Inf,ymax=Inf),fill=NA,col="black",size=0.3)
      
    } else {
      
      data_highlight <- data_highlight %>%
        mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level)) %>% 
        mutate(Gene=factor(Gene,levels = gene_level)) 
      
      p1 <- p1+geom_tile(data = data_highlight,aes(Tumor_Barcode,as.integer(Gene)),fill=NA,col="black",size=0.6)
      
    }
  }
  
  
  
  if(legend_adjust){
    p1 <- p1+guides(fill = guide_legend(nrow = 1,title.position = 'left',title.hjust = 0.5,title.vjust = 0.5))+theme(legend.position = "bottom")
  }
  
  
  if(!is.null(data_clone)){
    data_clone <- data_clone %>% filter(Tumor_Barcode %in% sample_level0)
    data_clone <- data_clone %>% mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level)) %>% 
      mutate(Gene=factor(Gene,levels = gene_level))
    p1 <- p1+geom_point(data=data_clone,pch=16,col="gray20",size=0.7)
  }
  
  oncoplot_legend <- get_legend(p1+theme(plot.margin=margin(b=-12,l=0,r=0,unit="cm"),legend.key.size = unit(0.18, "cm")))
  #legend.box.margin = margin(t = 0, b = 0, l = 0, r = 0,unit = 'cm')
  p1 <- p1+theme(legend.position="none")
  
  #  
  # gp1<- ggplot_gtable(ggplot_build(p1))
  # if(!is.null(maxWidth)){
  #   gp1$widths[1:3] <- maxWidth
  #   p1 <- as.ggplot(gp1)
  # } else{
  #   maxWidth <- grid:::unit.list(gp1$widths)
  # }
  # 
  #https://stackoverflow.com/questions/35822268/setting-width-in-gtable-object-collapses-plot-this-used-to-work-but-it-doesnt
  # extract the first three widths, 
  # corresponding to left margin, y lab, and y axis
  # gp1 <- ggplotGrob(p1) # convert to gtable
  # gp1$widths=grid:::unit.list(gp1$widths)
  # if(is.null(maxWidth)){ 
  #   maxWidth <- gp1$widths[1:3]
  # }else{
  #   gp1$widths[1:3] <- maxWidth
  #   p1 <- as.ggplot(gp1)
  # }
  
  if(p2_hidden){
    Gene_Freq_alte <-  Gene_Freq_alte %>% mutate(Alteration=factor(NA,levels =unique(Gene_Freq_alte$Alteration)))
  }
  
  p2 <- Gene_Freq_alte %>% 
    ggplot(aes(Gene,y=n,fill=fct_rev(Alteration)))+geom_bar(stat="identity",width = 0.5,size=0)+ 
    scale_x_discrete(expand = expand_scale(add=c(0.5,0.5)))+
    scale_y_continuous(breaks=nbreaks,expand = c(0,0),position = 'right',limits = c(0,max(nbreaks)))+
    coord_flip()+
    scale_fill_manual(values = landscape_colors)+
    theme(legend.position = 'none',panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.title.y.right = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.text.x=element_text(size = 6),axis.line.x = element_line(size = 0.2),axis.ticks.x=element_line(size=0.2))
  #p2 <- flush_ticks(p2,plot = FALSE)
  
  if(p2_axis_hidden){
    p2 <- p2+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line.x = element_blank())
  }
  
  
  
  if(!is.null(Gene_Sig)){
    p3 <- Gene_Sig %>% 
      #Gene_Freq %>% mutate(pvalue=n) %>% 
      ggplot(aes(Gene,y=pvalue))+geom_bar(stat="identity",width = 0.5,size=0,fill="gray60")+ 
      scale_x_discrete(expand = expand_scale(add=c(0.5,0.5)),position = 'top')+
      scale_y_reverse(breaks=pretty_breaks(n = 2),expand = c(0,0),position = 'right')+
      coord_flip()+
      scale_fill_manual(values = landscape_colors)+
      theme(legend.position = 'none',panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.title.y.right = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.text.x=element_text(size = 6),axis.line.x = element_line(size = 0.2),axis.ticks.x=element_line(size=0.2))
  }
  
  if(!is.null(Gene_Sig)){
    p_com <- align_plots(
      p3+theme(plot.margin=margin(l=0.2,r=0,t=tmar,b=bmar,unit="cm")),
      p1+theme(plot.margin=margin(r=0,t=tmar,b=bmar,unit="cm")),
      p2+theme(plot.margin=margin(l=-0.1,r=0.1,t=tmar,b=bmar,unit="cm")),
      align = 'h',
      axis = 'tb'
    )
    
  } else {
    p_com <- align_plots(p1+theme(plot.margin=margin(r=0,l=0.2,t=tmar,b=bmar,unit="cm")),
                         p2+theme(plot.margin=margin(l=-0.1,r=0.1,t=tmar,b=bmar,unit="cm")),
                         align = 'h',
                         axis = 'tb')
  }
  
  ## output order ##
  data <- data0
  if(!input_sample_level){
    
    if(GeneSortOnly){
      #use ony the Gene as order
      rankinfo <- data %>%
        select(Gene) %>%
        unique() %>%
        arrange(desc(Gene)) %>%
        mutate(Order=seq_along(Gene))
    }else{
      rankinfo <- data %>%
        select(Gene,Alteration) %>%
        unique() %>% 
        arrange(desc(Gene),Alteration) %>% 
        mutate(Order=seq_along(Gene))
    }
    
    maxrnak <- max(rankinfo$Order)+1
    
    rankinfo <- tibble(Tumor_Barcode=sample_level) %>% 
      left_join(
        data %>% 
          select(Tumor_Barcode,Gene,Alteration) %>% 
          left_join(rankinfo) %>% 
          group_by(Tumor_Barcode) %>% 
          arrange(Order) %>% 
          slice(1) %>% 
          ungroup()
      ) %>%
      mutate(Order=if_else(is.na(Order),as.integer(maxrnak),as.integer(Order)))
    
    rankinfo <- rankinfo %>% select(Tumor_Barcode,Order)
    colnames(rankinfo)[2] <- altertype
  } else{
    
    rankinfo <- tibble(Tumor_Barcode=sample_level) %>%
      mutate(Order=as.integer(seq_along(Tumor_Barcode))) 
    colnames(rankinfo)[2] <- altertype
  }
  
  return(list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level))
  
}




oncoplot2 <- function(data,sample_level0,Gene_Sig=NULL, scale_fill_ztw=scale_fill_viridis_c(),gene_level=NULL,sample_level=NULL,bmar=0,tmar=0,p2_axis_hidden=FALSE,p2_hidden=FALSE){
  
  # limited to sample level
  data <- data %>% filter(Tumor_Barcode %in% sample_level0)
  
  
  samplesize <- length(sample_level0)
  altertype <- data$Type[1]
  
  input_sample_level=TRUE
  
  if(is.null(sample_level)){
    input_sample_level=FALSE
    sampleorder <- data %>% 
      select(Tumor_Barcode,Gene,Alteration) %>% 
      arrange(desc(Alteration)) %>% 
      pull(Tumor_Barcode)
    
    sample_level <- c(sampleorder,sample_level0[!c(sample_level0 %in% sampleorder)])
  }
  gene_level <- data$Gene[1]
  ## resort the data 
  databg <- crossing(Tumor_Barcode=sample_level,Gene=gene_level) %>% 
    mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level)) %>% 
    mutate(Gene=factor(Gene,levels = gene_level))
  data <- data %>% mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level)) %>% 
    mutate(Gene=factor(Gene,levels = gene_level))
  
  alteration_mean <- number(mean(data$Alteration),accuracy = 0.01)
  
  p1 <- databg %>% 
    ggplot(aes(Tumor_Barcode,as.integer(Gene))) +
    geom_tile(height=0.95,fill="gray95")+
    theme_ipsum_rc(base_size = 10,axis_text_size = 10,grid = '',axis = '')+
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.title.y.right = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(face = "italic"),panel.background = element_blank(),legend.position = 'top',legend.justification = "top",legend.direction = 'vertical',legend.key.height =unit(1,"line") )+
    scale_y_continuous(expand = c(0,0),breaks = 1:length(gene_level),labels = gene_level,sec.axis = dup_axis(labels = alteration_mean))+
    ## data
    geom_tile(data = data,aes(Tumor_Barcode,as.integer(Gene),fill=Alteration),height=0.95)+scale_fill_ztw+
    geom_vline(xintercept = 1:samplesize-0.5,color="white",size=0.2)+
    panel_border(size = 0.3,color = 'gray70')+
    labs(fill=altertype)

  
  oncoplot_legend <- get_legend(p1+theme(plot.margin=margin(b=-12,l=0,r=0,unit="cm"),legend.key.size = unit(0.18, "cm")))
  #legend.box.margin = margin(t = 0, b = 0, l = 0, r = 0,unit = 'cm')
  p1 <- p1+theme(legend.position="none")
  
  datap2 <- data
  if(p2_hidden){
    datap2 <- data %>% mutate(Alteration=NA_real_)
  }
  
  p2 <- datap2 %>% ggplot(aes(y=Alteration))+
    geom_boxplot(outlier.color = NA)+
    scale_x_continuous(expand = expand_scale(add=c(0.5,0.5)))+
    scale_y_continuous(expand = c(0,0),position = 'right')+
    coord_flip(clip = 'off')+
    theme(legend.position = 'none',panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.title.y.right = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.text.x=element_text(size = 3),axis.line.x = element_line(size = 0.2),axis.ticks.x=element_line(size=0.2)) #axis.line.y = element_line(size=0.2)
  
  if(p2_axis_hidden){
    p2 <- p2+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line.x = element_blank())
  }
  
  
  if(!is.null(Gene_Sig)){
    p3 <- ggplot()+theme(panel.background = element_blank())+scale_y_discrete(expand = expand_scale(add=c(0.5,0.5)))+scale_x_continuous(expand = c(0,0))
  }
  
  if(!is.null(Gene_Sig)){
    p_com <- align_plots(
      p3+theme(plot.margin=margin(l=0.2,r=0,t=tmar,b=bmar,unit="cm")),
      p1+theme(plot.margin=margin(r=0,b=bmar,t=tmar,unit="cm")),
      p2+theme(plot.margin=margin(l=-0.1,r=0.1,t=tmar,b=bmar,unit="cm")),
      align = 'h',
      axis = 'tb'
    )
    
  } else {
    p_com <- align_plots(p1+theme(plot.margin=margin(r=0,l=0.2,b=bmar,t=tmar,unit="cm")),
                         p2+theme(plot.margin=margin(l=-0.1,r=0.1,b=bmar,t=tmar,unit="cm")),
                         align = 'h',
                         axis = 'tb')
  }
  
  ## output order ##
  
  if(!input_sample_level){
    rankinfo <- tibble(Tumor_Barcode=sample_level) %>% 
      left_join(
        data %>% 
          select(Tumor_Barcode,Alteration)
      ) %>% 
      arrange(desc(Alteration)) %>% 
      mutate(Order=as.integer(seq_along(Tumor_Barcode)))
    
    rankinfo <- rankinfo %>% select(Tumor_Barcode,Order)
    colnames(rankinfo)[2] <- altertype
  } else{
    
    rankinfo <- tibble(Tumor_Barcode=sample_level) %>%
      mutate(Order=as.integer(seq_along(Tumor_Barcode))) 
    colnames(rankinfo)[2] <- altertype
  }
  
  return(list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level))
  
}


oncoplot3 <- function(data,sample_level0,Gene_Sig=NULL, landscape_colors=NULL,gene_level=NULL,sample_level=NULL,bmar=0,tmar=0,height=10){
  # limited to sample level
  data <- data %>% filter(Tumor_Barcode %in% sample_level0)
  
  samplesize <- length(sample_level0)
  altertype <- data$Type[1]
  
  input_sample_level=TRUE
  
  if(is.null(sample_level)){
    input_sample_level=FALSE
    sampleorder <- data %>% 
      select(Tumor_Barcode,Gene,Alteration) %>% 
      arrange(desc(Alteration)) %>% 
      pull(Tumor_Barcode)
    
    sample_level <- c(sampleorder,sample_level0[!c(sample_level0 %in% sampleorder)])
  }
  
  
  if(is.null(gene_level)){
    gene_level <- data %>% group_by(Gene) %>% summarise(mean=mean(Alteration)) %>% arrange(mean) %>% pull(Gene)
  }
  
  ## resort the data 
  
  data <- data %>% mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level)) %>% 
    mutate(Gene=factor(Gene,levels = gene_level))
  
  
  p1 <- data %>% ggplot(aes(Tumor_Barcode,Alteration,fill=Gene))+
    geom_bar(position="stack", stat="identity",size=0)+
    #scale_x_discrete(expand =c(0,0))+
    scale_fill_manual(values = landscape_colors)+
    theme_ipsum_rc(base_size = 10,axis_text_size = 6,grid = 'Y',axis = 'XY',axis_col = 'black')+
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.title.y.right = element_blank(),axis.text.x = element_blank(),panel.background = element_blank(),legend.position = 'top',legend.justification = "top",legend.direction = 'vertical',legend.key.height =unit(1,"line"),axis.ticks.y = element_line(color = 'black') )+
    scale_y_continuous(expand = c(0,0),breaks = pretty_breaks(n = 3),labels = comma,sec.axis = dup_axis())+
    #panel_border(size = 0.3,color = 'black')+
    labs(fill=altertype)
  
  oncoplot_legend <- get_legend(p1+theme(plot.margin=margin(b=-12,l=0,r=0,unit="cm"),legend.key.size = unit(0.18, "cm")))
  #legend.box.margin = margin(t = 0, b = 0, l = 0, r = 0,unit = 'cm')
  p1 <- p1+theme(legend.position="none")
  
  p2 <- ggplot()+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line.x = element_blank(),panel.background = element_blank())
  
  if(!is.null(Gene_Sig)){
    p3 <- ggplot()+theme(panel.background = element_blank())+scale_y_discrete(expand = expand_scale(add=c(0.5,0.5)))+scale_x_continuous(expand = c(0,0))
  }
  
  if(!is.null(Gene_Sig)){
    p_com <- align_plots(
      p3+theme(plot.margin=margin(l=0.2,r=0,t=tmar,b=bmar,unit="cm")),
      p1+theme(plot.margin=margin(r=0,b=bmar,t=tmar,unit="cm")),
      p2+theme(plot.margin=margin(l=-0.1,r=0.1,t=tmar,b=bmar,unit="cm")),
      align = 'h',
      axis = 'tb'
    )
  } else {
    p_com <- align_plots(p1+theme(plot.margin=margin(r=0,l=0.2,b=bmar,t=tmar,unit="cm")),
                         p2+theme(plot.margin=margin(l=-0.1,r=0.1,b=bmar,t=tmar,unit="cm")),
                         align = 'h',
                         axis = 'tb')
  }
  
  ## output order ##
  
  if(!input_sample_level){
    rankinfo <- tibble(Tumor_Barcode=sample_level) %>% 
      left_join(
        data %>% 
          select(Tumor_Barcode,Alteration)
      ) %>% 
      arrange(desc(Alteration)) %>% 
      mutate(Order=as.integer(seq_along(Tumor_Barcode)))
    
    rankinfo <- rankinfo %>% select(Tumor_Barcode,Order)
    colnames(rankinfo)[2] <- altertype
  } else{
    
    rankinfo <- tibble(Tumor_Barcode=sample_level) %>%
      mutate(Order=as.integer(seq_along(Tumor_Barcode))) 
    colnames(rankinfo)[2] <- altertype
  }
  
  return(list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level,height=height))
  
}


oncoplot4 <- function(sample_level,bmar=0,tmar=0,height=10,labelsize=4,Gene_Sig=NULL){
  
  data <- tibble(Tumor_Barcode=sample_level,value=1) %>% mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level))
  
  p1 <- data %>% 
    ggplot(aes(Tumor_Barcode,value,label=Tumor_Barcode))+
    geom_text(angle=90,size=labelsize,hjust=1,vjust=0.5)+
    scale_y_continuous(limits = c(0,1),expand = c(0,0))+
    theme_void(base_family = "Roboto Condensed")
  #    theme(plot.margin=margin(t=0,b=-5,unit="cm"))
  
  
  p2 <- ggplot()+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line.x = element_blank(),panel.background = element_blank())
  
  if(!is.null(Gene_Sig)){
    p3 <- ggplot()+theme(panel.background = element_blank())+scale_y_discrete(expand = expand_scale(add=c(0.5,0.5)))+scale_x_continuous(expand = c(0,0))
  }
  
  if(!is.null(Gene_Sig)){
    p_com <- align_plots(
      p3+theme(plot.margin=margin(l=0.2,r=0,t=tmar,b=bmar,unit="cm")),
      p1+theme(plot.margin=margin(r=0,b=bmar,t=tmar,unit="cm")),
      p2+theme(plot.margin=margin(l=-0.1,r=0.1,t=tmar,b=bmar,unit="cm")),
      align = 'h',
      axis = 'tb'
    )
  } else {
    p_com <- align_plots(p1+theme(plot.margin=margin(r=0,l=0.2,b=bmar,t=tmar,unit="cm")),
                         p2+theme(plot.margin=margin(l=-0.1,r=0.1,b=bmar,t=tmar,unit="cm")),
                         align = 'h',
                         axis = 'tb')
  }
  
  ## output order ##
  return(list(oncoplot=p_com,oncoplot_legend=NULL,sample_level=sample_level,gene_level=NULL,height=height))
}



# oncoplot_combined <- function(...,oncolist=NULL){
#   oncolist <- c(list(...), oncolist)
#   num_plots <- length(oncolist)
#   plotlist <- list()
#   leglist <- list()
#   sizelist <- vector(length = num_plots)
#   for(i in 1:num_plots){
#     plotlist[[i]] <- oncolist[[i]]$oncoplot
#     leglist[[i]] <- oncolist[[i]]$oncoplot_legend
#     sizelist[i] <- length(oncolist[[i]]$gene_level)
#   }
#   #print(length(plotlist))
#   #print(sizelist)
#  
#   
#   
#   
#   
#   oncoplot_com <- plot_grid(plotlist = plotlist,
#                             align = 'v',
#                             axis = 'lr',
#                             rel_widths = rep(1,num_plots),
#                             rel_heights = sizelist,
#                             ncol= 1)
#   oncoplot_legend <- plot_grid(plotlist = c(leglist,list(NULL)),
#                                align = 'h',
#                                rel_widths = c(rep(1,num_plots),5),
#                                nrow = 1)
#   oncoplot_final <- plot_grid(
#     oncoplot_com,
#     oncoplot_legend,
#     align = 'v',axis = 'l',ncol = 1)
#   return(oncoplot_final)
# }



oncoplot_combined <- function(...,oncolist=NULL){
  oncolist <- c(list(...), oncolist)
  num_plots <- length(oncolist)
  plotlist1 <- list()
  plotlist2 <- list()
  leglist <- list()
  sizelist <- vector(length = num_plots)
  for(i in 1:num_plots){
    plotlist1[[i]] <- (oncolist[[i]]$oncoplot)[[1]]
    plotlist2[[i]] <- (oncolist[[i]]$oncoplot)[[2]]
    leglist[[i]] <- oncolist[[i]]$oncoplot_legend
    if(is.null(oncolist[[i]]$height)){
      sizelist[i] <- length(oncolist[[i]]$gene_level)
    }else {
      sizelist[i] <- oncolist[[i]]$height
    }
    
  }
  #print(length(plotlist))
  #print(sizelist)
  
  oncoplot_com1 <- plot_grid(plotlist = plotlist1,
                             align = 'v',
                             axis = 'lr',
                             rel_widths = rep(1,num_plots),
                             rel_heights = sizelist,
                             ncol= 1)
  oncoplot_com2 <- plot_grid(plotlist = plotlist2,
                             align = 'v',
                             axis = 'lr',
                             rel_widths = rep(1,num_plots),
                             rel_heights = sizelist,
                             ncol= 1)
  
  
  
  oncoplot_com <- plot_grid(oncoplot_com1,
                            oncoplot_com2,
                            align = 'h',
                            axis = 'tb',
                            rel_widths = c(10,1),
                            rel_heights = c(1,1),
                            nrow = 1)
  oncoplot_legend <- plot_grid(plotlist = c(leglist,list(NULL)),
                               align = 'h',
                               rel_widths = c(rep(1,num_plots),1),
                               nrow = 1)
  oncoplot_final <- plot_grid(
    oncoplot_com,
    oncoplot_legend,rel_heights = c(1,0.2),
    align = 'v',axis = 'l',ncol = 1)
  return(oncoplot_final)
}






## most frequent mutated drive genes ##
load('sherlock_altered_maf.RData')
load('sherlock_drivegene_list.RData')
load('sherlock_drive.RData')

genetmp <- sherlock_altered_maf %>% 
  filter(Hugo_Symbol %in% sherlock_drivegene_list$Gene,Variant_Type %in% c("SNP","INS","DEL","Promoter")) %>% 
  count(Subject,Hugo_Symbol) %>% 
  count(Hugo_Symbol,sort = T) %>% 
  mutate(freq=n/232) %>% 
  filter(freq>0.03) %>% 
  pull(Hugo_Symbol)

genetmp <- unique(c(genetmp, "TP53","RBM10","EGFR","KRAS","ARID1A","CDKN2A","NKX2-1","SMAD4","UBA1","TERT",sherlock_drive$SYMBOL))

data_top <- sherlock_altered_maf %>% 
  filter(Hugo_Symbol %in% genetmp,Variant_Type %in% c("SNP","INS","DEL","Promoter")) %>% 
  select(Subject,Gene=Hugo_Symbol,Alteration=Variant_Classification) %>%
  mutate(Type="Mutation_Driver") %>% 
  left_join(sherlock_samples_unique %>% select(Subject,Tumor_Barcode)) %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type)


# Clone and subclone information #
load('sherlock_clonality.RData')
clone_info <- sherlock_clonality %>% mutate(Variant_Classification=if_else(Variant_Classification=="5'Flank","Promoter",Variant_Classification)) %>% 
  select(Tumor_Barcode,Gene=Gene_Name,Alteration=Variant_Classification,Clone) %>% 
  filter(Gene %in% data_top$Gene,Alteration %in% data_top$Alteration,Clone=="Y") %>% 
  unique()

data_clone <- data_top %>% 
  left_join(clone_info) %>% 
  filter(Clone=="Y")


# SFTB nocoding  --------------------------------------------------------------------
load('sherlock_sftp.RData')

sherlock_sftp_info <- sherlock_sftp_info %>% 
  mutate(Chromosome=as.character(Chromosome)) %>% 
  select(Subject,CHROM=Chromosome,POS=Start_Position,REF=Reference_Allele,ALT=Tumor_Seq_Allele2,Gene=Hugo_Symbol,Alteration=Variant_Type) %>% 
  left_join(sherlock_clonality %>% mutate(Subject=str_sub(Tumor_Barcode,1,9)) %>% select(Subject,CHROM,POS,REF,ALT,Clone ) )


data_sftb <- sherlock_sftp_info %>% 
  select(Subject,Gene,Alteration,Clone) %>%
  mutate(Type="Mutation_SFTP") %>% 
  left_join(sherlock_samples_unique %>% select(Subject,Tumor_Barcode)) %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type,Clone) %>% 
  mutate(Alteration=if_else(Alteration=="SNP","SNV",Alteration))

data_sftb_clone <- data_sftb %>% filter(Clone=="Y")
data_sftb <- data_sftb %>% select(-Clone)


# CNV bar -----------------------------------------------------------------
#loh hla
# load('sherlock_lohhla.RData')
# data_hla <- sherlock_lohhla %>%
#   mutate(Gene="HLA",Alteration="HLA_LOH",Type="SCNA") %>% 
#   left_join(sherlock_samples_unique %>% select(Subject,Tumor_Barcode)) %>% 
#   select(Subject,Tumor_Barcode,Gene,Alteration,Type) 
# 
# load('sherlock_cnvinfo.RData')
# data_chr19loss <-  sherlock_cnvinfo %>% 
#   filter(Chr19loss_Clust=="Y") %>% 
#   select(Tumor_Barcode,Chr19loss_Clust) %>%
#   mutate(Subject=str_sub(Tumor_Barcode,1,9),Gene="Chr19",Alteration="Chr19_Loss",Type="SCNA") %>% 
#   select(Subject,Tumor_Barcode,Gene,Alteration,Type) 

load('sherlock_gistic_gene.RData')
data_gistic <- sherlock_gistic_gene %>% filter(GeneSymbol %in% c('CDKN2A','MDM2','MCL1','TERT','EGFR','BRAF','MYC','MDM2','NKX2-1','GRIN2A','CSMD1','STK11'))
data_gistic <- data_gistic %>% 
  filter(CNV %in% c(-2,2)) %>% 
  mutate(CNV=if_else(CNV==2,"Amplification","Deletion")) %>% 
  mutate(Subject=str_sub(Tumor_Barcode,1,9),Gene=GeneSymbol,Alteration=CNV,Type="SCNA_Focal") %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type) 

#data_scna <- bind_rows(data_chr19loss,data_hla,data_gistic)
data_scna <- data_gistic

# CNV ARM bar -------------------------------------------------------------
load('sherlock_cnv_arm.RData')
amparm <- sherlock_cnv_arm_info %>% filter(`Amp q-value`<0.01) %>% pull(Arm)
delarm <- sherlock_cnv_arm_info %>% filter(`Del q-value`<0.01) %>% pull(Arm)

## use heatmap arm instead of gistic arm
data_arm <- sherlock_cnv_arm_heatmap %>% 
  filter((calling=="2" & Arm %in% amparm)|(calling=="-2" & Arm %in% delarm)) %>% 
  mutate(Subject=str_sub(Tumor_Barcode,1,9),Gene=Arm,Alteration=type,Type="SCNA_Arm") %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type) %>% 
  mutate(Alteration=if_else(Alteration=="Amp","Amplification","Deletion")) %>% 
  mutate(Gene=paste0("Chr",Gene))

#data_arm <- bind_rows(data_arm,data_chr19loss,data_hla)

## remove chr19 loss subject
load('sherlock_cnvinfo.RData')
tmp <- sherlock_cnvinfo %>% filter(Chr19loss_Clust=="Y") %>% pull(Tumor_Barcode)
data_arm <- data_arm %>% filter(!(str_detect(Gene,'Chr19') & (Tumor_Barcode %in% tmp)))


# Oncogenic Fusion --------------------------------------------------------
load('sherlock_drive_fusion.RData')
#fgeneset1 <- c('RET','ALK','ROS1','FGFR2','MET','NRG1','AXL')
#fgeneset2 <- sherlock_drive_fusion %>% count(Gene,Subject) %>% count(Gene,sort = T) %>% filter(n>1) %>% pull(Hugo_Symbol)
#fgeneset <- unique(c(fgeneset1,fgeneset2))
data_fusion <- sherlock_drive_fusion %>%
  #  filter(Gene %in% fgeneset) %>% 
  mutate(Alteration="Fusion",Type="Fusion") %>% 
  left_join(sherlock_samples_unique %>% select(Subject,Tumor_Barcode)) %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type) 


# Germline mutations ------------------------------------------------------

load('sherlock_germline.RData')
load('sherlock_drivegene_list.RData')

genelist1 <- sherlock_germline_rare %>% select(Subject,HUGO_Symbol) %>% unique() %>% count(HUGO_Symbol,sort = T) %>% filter(n>2) %>% pull(HUGO_Symbol) 

genelist2 <- c('EGFR','ALK','ATM','BARD1','BRCA2','BRCA1','BRIP1','CHEK2','MRE11','NRN','PALB2','RAD51C','RAD51D','MSH6','TP53','ATM','BRCA1','XPC','OGG1','MUTYH','RAD51D','CHEK2','LIG4','POLG','GJB2','MET')

genelist <- c(genelist1,genelist2)

genelist <- genelist[! genelist %in% c("MET",'CYLD')]

data_germline <- 
  sherlock_germline_rare %>% 
  filter(HUGO_Symbol %in% genelist | HUGO_Symbol %in% sherlock_drivegene_list$Gene ) %>% 
  mutate(Type="Germline_Mutation") %>% 
  left_join(sherlock_samples_unique %>% select(Subject,Tumor_Barcode)) %>% 
  select(Subject,Tumor_Barcode,Gene=HUGO_Symbol,Alteration=VEP_Most_Severe_Consequence,Type)

data_germline <- data_germline %>% mutate(
  Alteration=recode(
    Alteration,
    "splice_donor_variant" = "Splice_Site",
    "stop_gained"='Nonsense_Mutation',
    "missense_variant"='Missense_Mutation',
    "frameshift_variant"='Frame_Shift_Del',
    "splice_acceptor_variant"="Splice_Site",
    "start_lost"='Nonsense_Mutation'
  )
)



# Catagoloy value ---------------------------------------------------------

load('sherlock_wgd.RData')
load('sherlock_msi.RData')
load('sherlock_cnvinfo.RData')
load('cnvclust.RData')
load('sherlock_kataegis.RData')

data_feature0 <- sherlock_histology %>%
  mutate(Subject,Tumor_Barcode,Gene="Histology",Alteration=Histology,Type="Other_feature") %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type) 

data_feature1 <- sherlock_wgd %>% 
  filter(WGD_Status=="Y") %>% 
  mutate(Gene="WGD",Alteration="WGD",Type="Other_feature") %>% 
  filter(Tumor_Barcode %in% sherlock_samples_unique$Tumor_Barcode) %>% 
  mutate(Subject=str_sub(Tumor_Barcode,1,9)) %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type) %>% 
  bind_rows(data_feature0)

# data_feature1 <- sherlock_msi %>% 
#   filter(MSI_Status !="MSS") %>% 
#   mutate(Gene="MSI",Alteration=MSI_Status,Type="Other_feature") %>% 
#   filter(Tumor_Barcode %in% sherlock_samples_unique$Tumor_Barcode) %>% 
#   mutate(Subject=str_sub(Tumor_Barcode,1,9)) %>% 
#   select(Subject,Tumor_Barcode,Gene,Alteration,Type) %>% 
#   bind_rows(data_feature1)

#sherlock_cnvinfo
data_feature1 <- cnvclust %>% 
  mutate(Gene="SCNA_Cluster",Alteration=CNV_Clust,Type="Other_feature") %>% 
  mutate(Subject=str_sub(Tumor_Barcode,1,9)) %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type) %>% 
  bind_rows(data_feature1)

data_feature1 <- sherlock_kataegis %>% 
  mutate(Gene="Kataegis",Alteration='Kataegis',Type="Other_feature") %>% 
  filter(Tumor_Barcode %in% sherlock_samples_unique$Tumor_Barcode) %>% 
  mutate(Subject=str_sub(Tumor_Barcode,1,9)) %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type) %>% 
  bind_rows(data_feature1)


load('sherlock_lohhla.RData')
data_hla <- sherlock_lohhla %>%
  mutate(Gene="HLA_LOH",Alteration="HLA_LOH",Type="Other_feature") %>% 
  left_join(sherlock_samples_unique %>% select(Subject,Tumor_Barcode)) %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type) 

# load('sherlock_cnvinfo.RData')
# data_chr19loss <-  sherlock_cnvinfo %>% 
#   filter(Chr19loss_Clust=="Y") %>% 
#   select(Tumor_Barcode,Chr19loss_Clust) %>%
#   mutate(Subject=str_sub(Tumor_Barcode,1,9),Gene="Chr19_HD",Alteration="Chr19_Loss",Type="Other_feature") %>% 
#   select(Subject,Tumor_Barcode,Gene,Alteration,Type) 

load('sherlock_p53.RData')
data_p53 <- sherlock_p53 %>% 
  filter(MDM2_TP53=="Y") %>% 
  mutate(Gene="p53_deficiency",Alteration='p53_deficiency',Type="Other_feature") %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type) 

load('sherlock_rtk.RData')
data_rtk <- sherlock_rtk %>% 
  filter(RTK_Altered_Status=='RTK-RAS+') %>% 
  mutate(Gene="RTK-RAS+",Alteration='RTK-RAS+',Type="Other_feature",Subject=str_sub(Tumor_Barcode,1,9)) %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type) 


data_feature1 <- bind_rows(data_feature1,data_hla) %>% 
  # bind_rows(data_chr19loss) %>% 
  bind_rows(data_p53) %>% 
  bind_rows(data_rtk)

featurelevel <- rev(c("Histology","SCNA_Cluster","RTK-RAS+","Kataegis" ,"WGD","p53_deficiency","HLA_LOH"))                                  
#data_feature1 <- data_feature1 %>% mutate(Gene=factor(Gene,levels = c()))


# Continues value ---------------------------------------------------------
load('sherlock_samples_unique.RData')

data_purity <- sherlock_samples_unique %>%
  select(Subject,Tumor_Barcode,Tumor_Purity) %>% 
  mutate(Gene="Tumor_Purity",Alteration=Tumor_Purity,Type="Tumor_Purity") %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type)

data_ratio <- sherlock_samples_unique %>%
  select(Subject,Tumor_Barcode,Subclonal_Mutation_Ratio) %>% 
  mutate(Gene="Subclonal_Ratio",Alteration=Subclonal_Mutation_Ratio,Type="Subclonal_Ratio") %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type)

data_nrpcc <- sherlock_samples_unique %>%
  select(Subject,Tumor_Barcode,NRPCC) %>% 
  mutate(Gene="NRPCC",Alteration=NRPCC,Type="NRPCC") %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type)

data_cnvratio <- sherlock_cnvinfo %>% 
  mutate(Gene="SCNA_Ratio",Alteration=CNV_Coverage,Type="SCNA_Ratio") %>% 
  mutate(Subject=str_sub(Tumor_Barcode,1,9)) %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type)






# Backgroud data 
load('sherlock_sis_numbers.RData')

data_sis <- sherlock_sis_numbers %>% 
  pivot_longer(cols = -Tumor_Barcode) %>% 
  mutate(Subject=str_sub(Tumor_Barcode,1,9),Gene=name,Alteration=value,Type='Mutations') %>% 
  mutate(Gene=str_remove(Gene,"_Num")) %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type) %>% 
  filter(Tumor_Barcode %in% sherlock_samples_unique$Tumor_Barcode)


# make oncoplot  ----------------------------------------------------------
sample_level0 <- sherlock_samples_unique$Tumor_Barcode
result_top <- oncoplot(data = data_top,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
tmplevel <- result_top$gene_level
tmplevel[26:27] <-tmplevel[c(27:26)]
result_top <- oncoplot(data = data_top,landscape_colors = landscape_colors,sample_level0 = sample_level0,gene_level = tmplevel,GeneSortOnly = TRUE)


result_sftb <- oncoplot(data = data_sftb,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_scna <- oncoplot(data = data_scna,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_arm <- oncoplot(data = data_arm,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_fusion <- oncoplot(data = data_fusion,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_germline <- oncoplot(data=data_germline,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_feature1 <- oncoplot(data = data_feature1,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)


# tmp <- sherlock_histology %>% select(Tumor_Barcode,Histology) %>% 
#   left_join(result_top$sample_level) %>%
#   left_join(
#     result_sftb$sample_level
#   ) %>% left_join(
#     result_fusion$sample_level
#   ) %>% left_join(
#     result_scna$sample_level
#   ) %>% left_join(
#     result_arm$sample_level
#   )    %>% left_join(
#     result_germline$sample_level
#   ) %>% left_join(
#     result_feature1$sample_level
#   ) %>% arrange(Histology,Mutation_Driver,Fusion,Mutation_SFTP,SCNA_Focal,SCNA_Arm,Germline_Mutation,Other_feature) 
# tmp %>% write_csv(tmp,path = 'sample_new_level.csv',col_names = T,na = '')
tmp <- read_csv('sample_new_level.csv',col_names = T)
sample_new_level <- tmp %>%   pull(Tumor_Barcode)


result_top <- oncoplot(data = data_top,data_clone = data_clone,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,tmar = -0.2,bmar = -0.05,gene_level = tmplevel)
result_sftb <- oncoplot(data = data_sftb,data_clone = data_sftb_clone,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
result_fusion <- oncoplot(data = data_fusion,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
result_scna <- oncoplot(data = data_scna,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
result_arm <- oncoplot(data = data_arm,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
resul_germline <- oncoplot(data=data_germline,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
result_feature1 <- oncoplot(data = data_feature1,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,gene_level = featurelevel,bmar = 0,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_purity <- oncoplot2(data = data_purity,sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_viridis_c(),tmar=0,bmar = -0.1,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_ratio <- oncoplot2(data = data_ratio,sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_viridis_c(option = "C"),tmar=-0.2,bmar = -0.1,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_nrpcc <- oncoplot2(data = data_nrpcc,sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_material(palette = 'green'),tmar=-0.2,bmar =0,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_cnvratio <- oncoplot2(data = data_cnvratio,sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_material(palette = 'blue'),tmar=-0.2,bmar =0,p2_axis_hidden = TRUE,p2_hidden = TRUE)


result_sis <- oncoplot3(data = data_sis,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,tmar=0.2,bmar=-0.2,height = 6)

oncoplot_final <- oncoplot_combined(result_sis,result_top,result_fusion,result_sftb,result_scna,result_arm,resul_germline,result_feature1,result_purity,result_ratio,result_nrpcc,result_cnvratio)
save_plot(filename = 'Figures/Genome_landscape2.pdf',plot = oncoplot_final,base_height = 18,base_width = 16,device=cairo_pdf)

data_all1 <- bind_rows(
  data_top,
  data_sftb,
  data_scna,
  data_fusion,
  data_arm,
  data_germline,
  data_feature1
)

data_all2 <- bind_rows(
  data_sis,
  data_purity,
  data_ratio,
  data_nrpcc,
  data_cnvratio
)

save.image('sherlock_landscape_v2.RData')


# add LOH data # 
load('HRD_LOH.RData')

data_hrd <- HRD_LOH %>% mutate(Type="HRD_LOH")
result_hrd <- oncoplot(data = data_hrd,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
oncoplot_final2 <- oncoplot_combined(result_sis,result_top,result_fusion,result_sftb,result_scna,result_arm,resul_germline,result_hrd,result_feature1,result_purity,result_ratio,result_nrpcc,result_cnvratio)
save_plot(filename = 'Figures/Genome_landscape2.1.pdf',plot = oncoplot_final2,base_height = 20,base_width = 15,device=cairo_pdf)





# TP53 + MDM2 -------------------------------------------------------------
data_tmp <- bind_rows(
  data_top %>% filter(Gene %in% c('TP53','MDM2')),
  data_gistic %>% filter(Gene %in% c('TP53','MDM2'))
)

result_top <- oncoplot(data = data_tmp,landscape_colors = landscape_colors,sample_level0 = sample_level0,tmar = 0,bmar = -0.05)
samplenewlevel <- sherlock_histology %>% left_join(result_top$sample_level) %>% arrange(Mutation_Driver,Histology) %>% pull(Tumor_Barcode)

result_top <- oncoplot(data = data_tmp,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = samplenewlevel, tmar = 0,bmar = -0.05)
result_feature0 <- oncoplot(data = data_feature0,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = samplenewlevel,bmar = 0,p2_axis_hidden = TRUE,p2_hidden = TRUE)
oncoplot_tmp<- oncoplot_combined(result_top,result_feature0)
save_plot(filename = 'Figures/MDM2_TP53_tmp.pdf',plot = oncoplot_tmp,base_height = 1,base_width = 15,device=cairo_pdf)




# TOP genes  --------------------------------------------------------------
tmp <- sherlock_histology %>% count(Histology) %>% rename(total=n)
genetmp <- sherlock_altered_maf %>% 
  left_join(sherlock_histology) %>% 
  filter(Variant_Type %in% c("SNP","INS","DEL","Promoter")) %>% 
  group_by(Histology) %>% 
  count(Subject,Hugo_Symbol) %>% 
  select(-n) %>% 
  count(Hugo_Symbol,sort = T) %>% 
  left_join(tmp) %>% 
  mutate(freq=n/total) %>% 
  filter(freq>0.03,n>=3) %>% 
  arrange(desc(freq)) %>% 
  slice(1:30) %>% 
  pull(Hugo_Symbol)

#genetmp <- c(genetmp, "TP53","RBM10","EGFR","KRAS","ARID1A","CDKN2A","NKX2-1","SMAD4","UBA1","TERT")


data_high <- sherlock_altered_maf %>% 
  filter(Hugo_Symbol %in% genetmp,Variant_Type %in% c("SNP","INS","DEL","Promoter")) %>% 
  select(Subject,Gene=Hugo_Symbol,Alteration=Variant_Classification) %>%
  mutate(Type="Mutation_Driver") %>% 
  left_join(sherlock_samples_unique %>% select(Subject,Tumor_Barcode)) %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type)

load('sherlock_clonality.RData')
clone_info <- sherlock_clonality %>% mutate(Variant_Classification=if_else(Variant_Classification=="5'Flank","Promoter",Variant_Classification)) %>% 
  select(Tumor_Barcode,Gene=Gene_Name,Alteration=Variant_Classification,Clone) %>% 
  filter(Gene %in% data_top$Gene,Alteration %in% data_top$Alteration,Clone=="Y") %>% 
  unique()

data_high_clone <- data_high %>% 
  left_join(clone_info) %>% 
  filter(Clone=="Y")

result_high <- oncoplot(data = data_high,data_clone = data_high_clone,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_feature0 <- oncoplot(data = data_feature0,landscape_colors = landscape_colors,sample_level0 = sample_level0,bmar = 0,p2_axis_hidden = TRUE,p2_hidden = TRUE)

tmp <- sherlock_histology %>% 
  select(Tumor_Barcode,Histology) %>%
  left_join(result_high$sample_level) %>% 
  arrange(Histology,Mutation_Driver) %>% 
  pull(Tumor_Barcode)

result_high <- oncoplot(data = data_high,data_clone = data_high_clone,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = tmp,GeneSortOnly = TRUE)
result_feature0 <- oncoplot(data = data_feature0,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = tmp,bmar = 0,p2_axis_hidden = TRUE,p2_hidden = TRUE)
tmp<- oncoplot_combined(result_high,result_feature0)
save_plot(filename = 'Figures/top_mutated_genes.pdf',plot = tmp,base_height = 5,base_width = 12,device=cairo_pdf)



# Hotspot mutation frequency ----------------------------------------------

tmp <- sherlock_altered_maf %>% 
  filter(Subject %in% sherlock_histology$Subject) %>% 
  left_join(sherlock_histology) %>% 
  filter(Variant_Classification!="Fusion") %>% 
  mutate(AAChange=if_else(is.na(AAChange)|AAChange=="",paste(Chromosome,Start_Position,Tumor_Seq_Allele1,Tumor_Seq_Allele2,sep = ":"),AAChange)) %>% 
  group_by(Histology) %>% 
  count(Subject,Hugo_Symbol,AAChange) %>% 
  select(-n) %>% 
  count(Hugo_Symbol,AAChange) %>% 
  arrange(desc(n)) %>% 
  ungroup()

tmpp <- tmp %>% ungroup() %>% count(Hugo_Symbol,AAChange) %>% filter(n>1) %>% arrange(desc(n)) %>% select(-n)

tmpp %>% left_join(tmp)%>%
  mutate(AAChange=paste0(Hugo_Symbol,", ",AAChange)) %>% 
  mutate(AAChange=fct_reorder(AAChange,n,.desc = T)) %>% 
  mutate(Histology=factor(Histology,levels = sort(names(histcolor)))) %>% 
  ggplot(aes(AAChange,n,fill=Histology))+
  geom_bar(position = "stack",stat="identity")+
  labs(y="Number of mutated samples",x="")+
  theme_ipsum_rc(axis_title_size = 14,axis_title_just = "m",axis = "xy",ticks = T,grid = "Y")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),legend.position = c(0.92,0.75))+
  scale_fill_manual(values = histcolor,drop=FALSE)

ggsave('Figures/Mutaiton_hotspot.pdf',width = 10,height = 5,device = cairo_pdf)



# Chr19_HD_alteration -----------------------------------------------------

sample_level0 <- sherlock_cnvinfo %>% filter(Chr19loss_Clust=="Y") %>% pull(Tumor_Barcode)
sample_new_level <- data_sis %>% filter(Tumor_Barcode %in% sample_level0,Gene=="SNV") %>% arrange(desc(Alteration)) %>% pull(Tumor_Barcode)

datatmp1 <- bind_rows(
  data_top %>% filter(Tumor_Barcode %in% sample_level0) %>% mutate(Type="Mutation")
  #data_sftb %>% filter(Tumor_Barcode %in% sample_level0) %>% mutate(Type="Mutation")
)

datatmp2 <- bind_rows(
  data_scna %>% filter(Tumor_Barcode %in% sample_level0) %>% mutate(Type="SCNA"),
  data_arm %>% filter(Tumor_Barcode %in% sample_level0) %>% filter(!str_detect(Gene,'Chr19')) %>% mutate(Type="SCNA")
)

datatmp3 <- data_feature1 %>% 
  filter(Tumor_Barcode %in% sample_level0) %>%
  filter(!str_detect(Gene,'Chr19'),!str_detect(Gene,'p53'),!str_detect(Gene,'RTK'),!str_detect(Gene,'SCNA'))


result_tmp1 <- oncoplot(data = datatmp1,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,tmar = -0.2,bmar = -0.05,p2_hidden = TRUE,p2_axis_hidden = TRUE)
result_tmp2 <- oncoplot(data = datatmp2,landscape_colors = landscape_colors,gene_level = rev(c("TERT","MDM2","MYC","Chr7p","Chr9q","Chr15q","Chr17p","Chr18q")),sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_hidden = TRUE,p2_axis_hidden = TRUE)
result_tmp3 <- oncoplot(data = datatmp3,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = 0,p2_hidden = TRUE,p2_axis_hidden = TRUE)


chr19snv <- read_csv('chr19_mutations.csv',col_names = T)

data_tmp0 <- chr19snv %>% 
  filter(Tumor_Barcode %in% sample_level0) %>% 
  bind_rows(tibble(Tumor_Barcode="NSLC-0019-T01",chr19N=0,Ntotal=12)) %>% 
  select(Tumor_Barcode,chr19N,Ntotal) %>% 
  mutate(Ntotal=Ntotal-chr19N) %>% 
  pivot_longer(cols = -Tumor_Barcode,names_to = 'Gene',values_to = 'Alteration') %>% 
  mutate(Subject=str_sub(Tumor_Barcode,1,9),Type="Count") %>% 
  mutate(Gene=if_else(Gene=="chr19N","Mutations in chr19","Mutations in other chromosomes")) %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type)


data_tmp00 <- data_tmp0 %>% pivot_wider(names_from = Gene,values_from = Alteration) %>% mutate(Alteration=`Mutations in chr19`/(`Mutations in chr19`+`Mutations in other chromosomes`)) %>% mutate(Gene="Ratio",Type=Gene)%>% select(Subject,Tumor_Barcode,Gene,Alteration,Type)

data_tmp00 <- data_tmp00 %>% mutate(Alteration=Alteration/(59128983/3074258928))



# number of mutations in chr19 HD not due to the purity

load('sherlock_samples_unique.RData')
load('sherlock_cnvinfo.RData')
chr19snv <- read_csv('chr19_mutations.csv',col_names = T)


tmp <- sherlock_samples_unique %>% left_join(chr19snv) %>% left_join(sherlock_cnvinfo)

tmp <- tmp %>% filter(Tumor_Purity<0.3|Chr19loss_Clust=="Y") %>%
  mutate(Chr19loss_Clust=if_else(is.na(Chr19loss_Clust),"N",Chr19loss_Clust))


tmp %>% ggplot(aes(Chr19loss_Clust,chr19N/Ntotal))+geom_boxplot()

ggsave('tmp.pdf',height = 6,width = 2)

wilcox.test(tmp$chr19N/tmp$Ntotal~tmp$Chr19loss_Clust)




## add telomere length ## 
load('sherlock_tl.RData')
data_tmp000 <- data_tmp00 %>% select(Subject,Tumor_Barcode) %>% 
  left_join(
    sherlock_samples_unique %>% select(Tumor_Barcode) %>% left_join(sherlock_tl) %>% mutate(TLratio=Tumor_TL/Normal_TL) %>% arrange(desc(TLratio)) %>% mutate(rank=seq_along(TLratio))
  ) %>% 
  mutate(Gene="TL",Alteration=TLratio,Type="TL") %>% 
  select(Subject,Tumor_Barcode,Tumor_Barcode,Gene,Alteration,Type) 




tmpcol <- c('#ff7f00','#6a3d9a')
names(tmpcol) <- unique(data_tmp0$Gene)
landscape_colors <- c(landscape_colors,tmpcol)

result_tmp0 <- oncoplot3(data = data_tmp0,landscape_colors = tmpcol,sample_level0 = sample_level0,sample_level = sample_new_level,tmar=0.1,bmar=0.15,height = 6)
result_tmp00 <- oncoplot2(data = data_tmp00,sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_viridis_c(),tmar=-0.2,bmar = 0.1,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_tmp000 <- oncoplot2(data = data_tmp000,sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_material("purple"),tmar=-0.2,bmar = 0.1,p2_axis_hidden = TRUE,p2_hidden = TRUE)


result_sample <- oncoplot4(sample_level = sample_new_level,height = 4,labelsize = 2)


oncoplot_final <- oncoplot_combined(result_tmp0,result_tmp00,result_tmp000,result_tmp1,result_tmp2,result_tmp3,result_sample)


save_plot(filename = 'tmp.pdf',plot = oncoplot_final,base_height = 6,base_width = 2.5,device=cairo_pdf)











# Telomere length  --------------------------------------------------------
library(broom)
library(ggrepel)
load('sherlock_samples_unique.RData')
load('sherlock_clinical.RData')
load('sherlock_tl.RData')
data_tl <- sherlock_tl %>% mutate(TLratio=Tumor_TL/Normal_TL) %>%  select(Tumor_Barcode,Tumor_TL,Normal_TL,TLratio)

load('sherlock_landscape_v2.RData')
data_all <- bind_rows(
  data_top,
  data_sftb,
  data_scna,
  data_arm,
  data_fusion,
  data_feature1 %>% mutate(Type="@") %>% filter(!(Gene %in% c( "SCNA_Cluster" ,"MSI", "p53_deficiency", "RTK-RAS+")))
)
data_all <- data_all %>% 
  select(Tumor_Barcode,Gene,Type) %>%
  unique() %>% 
  mutate(Key=paste0(Gene,"_",Type)) %>% mutate(Key=str_remove(Key,"_@")) %>% 
  select(Tumor_Barcode,Key) 

tmpkey <- data_all %>% count(Key) %>% filter(n>5,n<232-5) %>% pull(Key)

data_all <- data_all %>% 
  filter(Key %in% tmpkey) %>% 
  mutate(Status="Y") %>% 
  pivot_wider(names_from = 'Key',values_from = 'Status') %>% 
  right_join(sherlock_samples_unique %>% select(Subject,Tumor_Barcode)) %>% 
  mutate_all(~replace_na(., 'N')) %>% 
  pivot_longer(cols = -c(Subject,Tumor_Barcode),names_to = 'Alternation',values_to = 'Status') %>% 
  mutate(Status=factor(Status,levels = c('N','Y'))) %>% 
  left_join(sherlock_clinical %>% select(Subject,age_at_diagnosis))

data_all <- data_all %>% left_join(data_tl)

meantl <- data_all %>% group_by(Alternation,Status) %>% summarise(mean=mean(TLratio)) %>% ungroup() %>% pivot_wider(names_from = "Status",values_from = 'mean') %>% mutate(TLdiff=Y-N)

data_all %>% 
  group_by(Alternation) %>% 
  do(tidy(wilcox.test(TLratio~Status,data=.))) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  mutate(fdr=p.adjust(p.value)) 

result_all <- data_all %>% 
  group_by(Alternation) %>% 
  do(tidy(t.test(TLratio~Status,data=.))) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  mutate(fdr=p.adjust(p.value,method = 'BH')) %>% 
  left_join(meantl)



# data_all %>% 
#   group_by(Alternation) %>% 
#   do(tidy(lm(TLratio~Status,data=.))) %>% 
#   ungroup() %>% 
#   arrange(p.value) %>% 
#   filter(term=="StatusY") %>% 
#   mutate(fdr=p.adjust(p.value)) %>% 
#   left_join(meantl)

result_all2 <- result_all %>% filter(fdr<0.1|Y>1.1|Y<0.9)

result_all2$Alternation[1] <- 'Chr9q Loss'
result_all2$Alternation[2] <- 'Chr22q Loss'
result_all2$Alternation[3] <- 'HLA LOH'
result_all2$Alternation[4] <- 'Chr9p Loss'
result_all2$Alternation[5] <- 'Chr18q Loss'
result_all2$Alternation[6] <- 'TERT Amplificaiton'
result_all2$Alternation[7] <- 'Chr7q Gain'
result_all2$Alternation[8] <- 'Chr19 HD'
result_all2$Alternation[9] <- 'SMAD4 Mutation'
result_all2$Alternation[10] <- 'STK11 Loss'
result_all2$Alternation[11] <- 'CDKN2A Mutation'
result_all2$Alternation[12] <- 'MET Mutation'




result_all %>% 
  ggplot(aes(Y,-log10(p.value)))+
  geom_point(pch=21,fill="gray50",size=3)+
  geom_text_repel(data=result_all2,aes(label=Alternation),size=3,force = 10)+
  labs(x="T/N TL Ratio",y='-log10(P-value)')+
  scale_x_continuous(breaks = pretty_breaks(n = 6),limits = c(0.7,1.3),expand = c(0,0))+
  scale_y_continuous(breaks = pretty_breaks(n=6),expand = c(0.02,0))+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 12,base_size = 10,axis = "XY",grid = "XYy",ticks = TRUE)+
  panel_border(color = 'black')+
  #geom_vline(xintercept = 0,linetype=2,color="grey50")+
  geom_hline(yintercept = -log10(0.0055),linetype=2,color="grey50")+
  geom_vline(xintercept = 1,linetype=2,color="grey50")

ggsave('Figures/TL_ratio_alteration.pdf',width = 7,height = 6,device = cairo_pdf)



## top candidate vs porportion of C1/C2/C3

bind_rows(
  data_arm %>% filter(Gene=="Chr9p"|Gene=="Chr9q") %>% mutate(Gene="Chr9p/q Loss"),
  data_arm %>% filter(Gene=="Chr22q") %>% mutate(Gene="Chr22q Loss"),
  data_hla %>% mutate(Gene="HLA LOH")
) %>% 
  left_join(sherlock_cnvinfo) %>% 
  mutate(CNV_Clust=factor(CNV_Clust,levels = c('C2','C1','C3'),labels = c('C1','C2','C3'))) %>% 
  ggplot(aes(Gene,fill=CNV_Clust))+geom_bar(position = "fill",width = 0.7)+
  scale_fill_manual("SCNA_Cluster",values = c("#FF7F0EFF","#2CA02CFF","#1F77B4FF"))+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  labs(x='',y='Freqeuncy (%)')+
  scale_y_percent()+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))

ggsave('Figures/TL_ratio_alteration2.pdf',width = 5,height = 7,device = cairo_pdf)




# mutual exclusivity ------------------------------------------------------

load('sherlock_landscape_v2.RData')
load('sherlock_clinical.RData')

data_all <- bind_rows(
  data_top,
  data_sftb,
  data_scna,
  data_arm,
  data_hla,
  #data_chr19loss,
  data_fusion,
  data_fusion %>% mutate(Gene="All_Fusion"),
  data_germline,
  data_feature1 %>% filter(Gene!="SCNA_Cluster"),
  data_feature1 %>% filter(Gene=="SCNA_Cluster") %>% mutate(Gene=Alteration),
)


tmp <- sherlock_samples_unique %>% 
  left_join(sherlock_clinical) %>% 
  mutate(passive_smoking=factor(passive_smoking,levels = c('Y','N'))) %>% 
  select(Tumor_Barcode,passive_smoking) %>% 
  left_join(
    data_all %>% mutate(Gene=paste0(Gene,"|",Type)) %>% select(Gene,Tumor_Barcode) %>% unique()
  ) %>%
  mutate(Mutation=factor("Y",levels = c("Y","N"))) %>% 
  pivot_wider(names_from = "Gene",values_from = "Mutation",values_fill = list(Mutation="N"))%>%
  select(-Tumor_Barcode)
# select(-`NA`) 

result <- tmp %>% select(c(1,2)) %>% table() %>% fisher.test() %>% tidy() %>% mutate(gene1='geneA',gene2='geneB') %>% slice(0)
size <- dim(tmp)[2]
for(i in 1:(size-1)){
  for(j in (i+1):size){
    geneA <- colnames(tmp)[i]
    geneB <- colnames(tmp)[j]
    result <- tmp %>% select(c(i,j)) %>% table() %>% fisher.test() %>% tidy() %>% mutate(gene1=geneA,gene2=geneB) %>% bind_rows(result)
  }
}

result <- result %>% select(gene1,gene2,estimate:conf.high)

View(result)

save(result,file = 'Mutual_exlusive.RData')



data_tmp <- data_all %>% filter((Gene=="MDM2" & Type=="SCNA_Focal") | (Gene=="TP53" & Type== "Mutation_Driver")) %>% mutate(Type="MDM2_TP53")
tmp<- oncoplot(data = data_tmp,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
plot(tmp$oncoplot[[1]])
tmp<- oncoplot_combined(tmp)
save_plot(filename = 'Figures/MDM2_TP53_KRAS.pdf',plot = tmp,base_height = 1,base_width = 12,device=cairo_pdf)

sherlock_p53 <- sherlock_samples_unique %>% select(Subject,Tumor_Barcode) %>% left_join(data_tmp) %>% mutate(MDM2_TP53=if_else(is.na(Type),"N","Y")) %>% unique() %>% select(Subject,Tumor_Barcode,MDM2_TP53,Gene,Alteration)
save(sherlock_p53,file='sherlock_p53.RData')



data_all2 <- bind_rows(
  data_top,
  data_sftb,
  data_scna,
  data_arm,
  data_fusion,
  data_germline,
  data_feature1 %>% filter(Gene!="SCNA_Cluster"),
  data_tmp %>% mutate(Gene="TP53_MDM2")
)



# chr1q amplificaiton vs WGD 
tmp1 <- data_arm %>% filter(Gene=='Chr1q') %>% mutate(Type="tmp")
tmp2 <- data_feature1 %>% filter(Gene=="WGD") %>% mutate(Type="tmp")
tmp <- bind_rows(tmp1,tmp2)
tmprs <- oncoplot(data = tmp,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE,p2_hidden = TRUE)
tmp<- oncoplot_combined(tmprs)
save_plot(filename = 'Figures/WGD_chr1q.pdf',plot = tmp,base_height = 1,base_width = 12,device=cairo_pdf)



# Mutational signature ----------------------------------------------------
load('sherlock_landscape_v2.RData')

data_all <- bind_rows(
  data_top,
  data_sftb,
  data_scna,
  data_arm,
  data_hla,
  data_chr19loss,
  data_fusion,
  data_fusion %>% mutate(Gene="All_Fusion"),
  data_germline,
  data_feature1 %>% filter(Gene!="SCNA_Cluster")
)


load('sherlock_signature.RData')
load('sherlock_clinical.RData')

tmp <- sherlock_samples_unique %>% 
  select(Tumor_Barcode) %>% 
  left_join(
    data_all %>% mutate(Gene=paste0(Gene,"|",Type)) %>% select(Gene,Tumor_Barcode) %>% unique()
  ) %>%
  mutate(Mutation=factor("Y",levels = c("N","Y"))) %>% 
  pivot_wider(names_from = "Gene",values_from = "Mutation",values_fill = list(Mutation="N")) %>% 
  select(-`NA`) %>% 
  pivot_longer(col=-Tumor_Barcode,names_to="Gene",values_to="Mutation") %>% 
  left_join(sherlock_SBS) %>% 
  left_join(sherlock_ID) %>% 
  left_join(sherlock_DBS) %>% 
  pivot_longer(col=-c(Tumor_Barcode,Gene,Mutation))


tmp2 <- tmp %>% 
  group_by(Gene,name) %>% 
  do(tidy(lm(value~Mutation,data=.))) %>% 
  filter(term=="MutationY") %>% 
  ungroup() %>% 
  arrange(p.value)


ttmp <- tmp %>% filter(value>0) %>% count(Gene,Mutation,name) %>% filter(n>5) %>% count(Gene,name) %>% filter(n>1)


tmp2 <- tmp %>% 
  filter(value>0) %>% 
  left_join(ttmp) %>% 
  filter(n==2) %>% 
  group_by(Gene,name) %>% 
  do(tidy(lm(value~Mutation,data=.))) %>% 
  filter(term=="MutationY") %>% 
  ungroup() %>% 
  arrange(p.value)




# Quiety tumor ------------------------------------------------------------
# define quiety tumor on the right 33 patients
load('sherlock_clinical.RData')
tmp <- tibble(Tumor_Barcode=sample_new_level) %>% left_join(sherlock_cnvinfo) %>% tail(33) %>%mutate(Subject=str_sub(Tumor_Barcode,1,9)) %>% pull(Subject)
tmp2 <- sherlock_clinical %>% mutate(group=if_else(Subject %in% tmp,"Quiet tumors (N=33)","Other tumors (N=199)"))
tmp2 %>% mutate(group2=if_else(age_at_diagnosis<45,"Y","N")) %>% select(group,group2) %>% table()%>% fisher.test()
# p-value= 0.01014  OR= 3.67
tibble(Subject=tmp) %>% left_join(sherlock_samples_unique) %>% write_csv(path = 'quiet_tumor.csv',col_names = T)

px <- tmp2 %>% mutate(group2=if_else(age_at_diagnosis<45,"Y","N")) %>% select(group,group2) %>% table() %>% as.data.frame() %>% group_by(group) %>% mutate(perc = round(Freq/sum(Freq),2)) %>% ggplot(aes(group,perc,fill=group2))+geom_col()+geom_text(aes(label=paste0(perc*100,"%")),position = position_stack(vjust = 0.5))+ theme_ipsum_rc(base_size = 12,axis_title_size = 14,axis = 'y',axis_title_just = "m",ticks = TRUE,axis_col = "black")+theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),panel.grid.major.y = element_line())+labs(y="Percentage of patients")+scale_fill_d3()+scale_y_continuous(breaks = pretty_breaks(),label=percent,expand = c(0.01,0.01))+labs(fill="Patients with age bellow 45")+theme(legend.position = "top")


tmp2 %>% mutate(stage=if_else(stage %in% c('IA','IB'),'I','II'))%>% group_by(group) %>% count(stage)
py <- tmp2 %>% mutate(stage=if_else(stage %in% c('IA','IB'),'I','I+')) %>% mutate(stage=fct_rev(fct_inorder(stage)))%>% select(stage,group) %>% table() %>% as.data.frame() %>% group_by(group) %>% mutate(perc = round(Freq/sum(Freq),2)) %>% ggplot(aes(group,perc,fill=stage))+geom_col()+geom_text(aes(label=paste0(perc*100,"%")),position = position_stack(vjust = 0.5))+ theme_ipsum_rc(base_size = 12,axis_title_size = 14,axis = 'y',axis_title_just = "m",ticks = TRUE,axis_col = "black")+theme(axis.title.x = element_blank(),axis.ticks.y = element_line(colour = "black"),legend.text = element_text(size = 10),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),panel.grid.major.y = element_line(),axis.title.y = element_text(margin = margin(r = 5)))+labs(y="Percentage of patients")+scale_fill_d3()+scale_y_continuous(breaks = pretty_breaks(),label=percent,expand = c(0.01,0.01))+labs(fill="Tumor stage")+theme(legend.position = "top")

plot_grid(px,py,align = "h",nrow = 1)
ggsave('Figures/Quiet_Tumor_features2.pdf',width = 6,height = 8,device = cairo_pdf)




# mutation signature
quiet_tumor <- tmp2 %>% mutate(group2=if_else(age_at_diagnosis<50,"Y","N")) %>% select(Subject,group)

quiet_tumor %>% left_join(sherlock_clinical) %>% select(group,gender) %>% table() %>% fisher.test()
quiet_tumor %>% left_join(sherlock_clinical) %>% select(group,passive_smoking) %>% table() %>% fisher.test()
tmp <- quiet_tumor %>% left_join(sherlock_sis_numbers %>% mutate(Subject=str_sub(Tumor_Barcode,1,9)))
wilcox.test(SNV_Num~group,data=tmp)
tmp %>% mutate(Mutation=SNV_Num+INS_Num) %>% ggplot(aes(group,log2(Mutation)))+geom_boxplot()



load('sherlock_signature.RData')
source('~/NIH-Work/MutationSignature/mSigPortal/Visualization/Sigvisualfunc.R')
#source('Sigvisualfunc.R')

quiet_tumor %>% 
  left_join(sherlock_SBS %>% filter(Tumor_Barcode %in% sherlock_samples_unique$Tumor_Barcode ) %>% mutate(Tumor_Barcode=str_sub(Tumor_Barcode,1,9)) %>% rename(Subject=Tumor_Barcode)) %>% 
  group_by(group) %>% 
  summarise_at(vars(contains("SBS")),sum) %>% 
  gather(Signature,Weight,contains("SBS")) %>% 
  ggplot(aes(group,Weight,fill=factor(Signature,levels = names(SBScolor))))+geom_bar(stat="identity",position="fill",col="gray95",size=0)+theme_ipsum_rc(base_size = 12,axis_title_just = "m",axis_title_size = 14,axis = "y",axis_col = "black")+theme(axis.ticks.x = element_blank(), axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),panel.grid.major=element_line(),legend.position = "right",legend.box.background = element_blank(),legend.box.spacing = unit(0.6,"cm"),legend.key = element_rect(size =0),legend.key.height = unit(0.8,"cm"),axis.ticks.y = element_line(colour = "black"),legend.key.size = unit(0.5, "cm"),legend.text = element_text(size = 10),legend.box.margin=margin(l=-8))+scale_y_continuous(expand = c(0.01,0.01),breaks = pretty_breaks(),label=percent)+xlab("")+scale_fill_manual("Signature",values = SBScolor)+ylab("Mutation contribution (%)")


# individual
SBS96_Clustering(sigdata = sherlock_SBS,clustern = 2 )


quiet_tumor %>% 
  left_join(sherlock_ID %>% filter(Tumor_Barcode %in% sherlock_samples_unique$Tumor_Barcode ) %>% mutate(Tumor_Barcode=str_sub(Tumor_Barcode,1,9)) %>% rename(Subject=Tumor_Barcode)) %>% 
  group_by(group) %>% 
  summarise_at(vars(contains("ID")),sum) %>% 
  gather(Signature,Weight,contains("ID")) %>% 
  ggplot(aes(group,Weight,fill=Signature))+geom_bar(stat="identity",position="fill",col="gray95",size=0)+theme_ipsum_rc(base_size = 12,axis_title_just = "m",axis_title_size = 14,axis = "y",axis_col = "black")+theme(axis.ticks.x = element_blank(), axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),panel.grid.major=element_line(),legend.position = "right",legend.box.background = element_blank(),legend.box.spacing = unit(0.6,"cm"),legend.key = element_rect(size =0),legend.key.height = unit(0.8,"cm"),axis.ticks.y = element_line(colour = "black"),legend.key.size = unit(0.5, "cm"),legend.text = element_text(size = 10),legend.box.margin=margin(l=-8))+scale_y_continuous(expand = c(0.01,0.01),breaks = pretty_breaks(),label=percent)+xlab("")+scale_fill_d3("category10")+ylab("Mutation contribution (%)")


#MUTATIONS
data_quiet <- sherlock_altered_maf %>% filter(Subject %in% quiet_tumor$Subject[quiet_tumor$group=="Quiet tumors (N=33)"]) %>% 
  filter(Variant_Type %in% c("SNP","INS","DEL","Promoter")) %>% 
  select(Subject,Gene=Hugo_Symbol,Alteration=Variant_Classification) %>%
  mutate(Type="Mutation_Driver") %>% 
  left_join(sherlock_samples_unique %>% select(Subject,Tumor_Barcode)) %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration,Type)

geneset <- data_quiet %>% count(Subject,Gene) %>% count(Gene,sort=T) %>% filter(n>1) %>% pull(Gene)
data_quiet <- data_quiet %>% filter(Gene %in% geneset) %>% mutate(Tumor_Barcode=Subject)

result_quiet <- oncoplot(data = data_quiet,landscape_colors = landscape_colors,sample_level0 = quiet_tumor$Subject[quiet_tumor$group=="Quiet tumors (N=33)"],GeneSortOnly = TRUE)
plot(result_quiet$oncoplot[[1]])


# telomere length
load('sherlock_tl.RData')
quiet_tumor %>% left_join(sherlock_tl) %>% mutate(TLratio=Tumor_TL/Normal_TL) %>% 
  ggplot(aes(group,TLratio))+geom_boxplot()

tmp <- quiet_tumor %>% left_join(sherlock_tl) %>% mutate(TLratio=Tumor_TL/Normal_TL)
wilcox.test(TLratio~group,data=tmp)

save(quiet_tumor,file = 'quiet_tumor.RData')




# landscape sorted by age -------------------------------------------------
set_wd()
load('sherlock_landscape_v2.RData')
load('sherlock_clinical.RData')
libztw()
library(data.table)
sample_new_level <- sherlock_samples_unique %>% left_join(sherlock_clinical) %>% arrange(age_at_diagnosis) %>% pull(Tumor_Barcode)

result_top <- oncoplot(data = data_top,data_clone = data_clone,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,tmar = -0.2,bmar = -0.05,gene_level = tmplevel)
result_sftb <- oncoplot(data = data_sftb,data_clone = data_sftb_clone,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
result_fusion <- oncoplot(data = data_fusion,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
result_scna <- oncoplot(data = data_scna,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
result_arm <- oncoplot(data = data_arm,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
resul_germline <- oncoplot(data=data_germline,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
result_feature1 <- oncoplot(data = data_feature1,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = 0,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_purity <- oncoplot2(data = data_purity,sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_viridis_c(),tmar=0,bmar = -0.1,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_ratio <- oncoplot2(data = data_ratio,sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_viridis_c(option = "C"),tmar=-0.1,bmar = 0,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_nrpcc <- oncoplot2(data = data_nrpcc,sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_material(palette = 'green'),tmar=-0.2,bmar =0,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_cnvratio <- oncoplot2(data = data_cnvratio,sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_material(palette = 'blue'),tmar=-0.2,bmar =0.2,p2_axis_hidden = TRUE,p2_hidden = TRUE)

result_sis <- oncoplot3(data = data_sis,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,tmar=0.2,bmar=-0.2,height = 6)

oncoplot_final <- oncoplot_combined(result_sis,result_top,result_fusion,result_sftb,result_scna,result_arm,resul_germline,result_feature1,result_purity,result_ratio,result_nrpcc,result_cnvratio)
save_plot(filename = 'Figures/Genome_landscape_sortedBy_age.pdf',plot = oncoplot_final,base_height = 20,base_width = 15,device=cairo_pdf)





## STK11
tmp1 <- data_gistic %>% filter(Gene=="STK11")
tmp2 <- sherlock_altered_maf %>% filter(Hugo_Symbol=="STK11") %>% left_join(sherlock_samples_unique) %>% select(Subject,Tumor_Barcode,Gene=Hugo_Symbol,  Alteration=Variant_Classification, Type=Variant_Type )

sherlock_STK11 <- bind_rows(tmp1,tmp2) 

save(sherlock_STK11,file = 'sherlock_STK11.RData')








# TMB/TL doesn't correlated with tumor purity -----------------------------
load('sherlock_tl.RData')
load('sherlock_samples_unique.RData')
load('sherlock_tmb.RData')




sherlock_tl %>% 
  mutate(TLratio=Tumor_TL/Normal_TL) %>% 
  right_join(sherlock_samples_unique) %>% 
  ggplot(aes(Tumor_Purity,TLratio))+
  geom_point(pch=21,fill=gray(0.5),size=3)+
  geom_smooth(method = 'lm')+
  labs(x="Tumor purity",y="T/N TL ratio")+
  theme_ipsum_rc(base_size = 12,axis_title_just = "m",axis_title_size = 14,ticks = T)+
  scale_x_continuous(breaks = pretty_breaks(),expand = c(0,0),limits = c(0.1,1))+
  scale_y_continuous(breaks = pretty_breaks(),expand = c(0,0),limits = c(0,2))+
  panel_border()

ggsave(filename = 'Purity_TL.pdf',width = 6,height = 5.5,device = cairo_pdf)


sherlock_samples_unique %>% left_join(sherlock_tmb) %>% 
  ggplot(aes(Tumor_Purity,log10(TMB)))+
  geom_point(pch=21,fill=gray(0.5),size=3)+
  geom_smooth(method = 'lm')+
  labs(x="Tumor purity",y="TMB (log10)")+
  theme_ipsum_rc(base_size = 12,axis_title_just = "m",axis_title_size = 14,ticks = T)+
  scale_x_continuous(breaks = pretty_breaks(),expand = c(0,0),limits = c(0.1,1))+
  scale_y_continuous(breaks = pretty_breaks(),expand = c(0,0),limits = c(-2,2))+
  panel_border()
ggsave(filename = 'Purity_TMB.pdf',width = 6,height = 5.5,device = cairo_pdf)





# APOBEC by Domitry and EGFR ----------------------------------------------
load('sherlock_samples_unique.RData')
load('sherlock_rtk.RData')
apobec <- read_delim('~/NIH-Work/MutationSignature/Dmitry_Gordenin/software/NSLC/output/NSLC_APOBEC_MAF_sorted_sum_all_fisher_Pcorr.txt',delim="\t",col_names = TRUE) %>% filter(Sample %in% sherlock_samples_unique$Tumor_Barcode)

apobec <- apobec %>%
  select(Sample,APOBEC_enrich,qvalue=`BH_Fisher_p-value_tCw`) %>%
  mutate(APOBEC=if_else(qvalue<0.05,"APOBEC_high","APOBEC_low")) %>% 
  mutate(Subject=str_sub(Sample,1,9)) %>% 
  left_join(
    sherlock_rtk_info %>% filter(Hugo_Symbol=="EGFR") %>% select(Subject) %>% mutate(EGFR="MT") %>% unique()
  ) %>% 
  mutate(EGFR=if_else(is.na(EGFR),"WT",EGFR))




# SCNA_C3 group -----------------------------------------------------------
load('sherlock_cnvinfo.RData')
load('sherlock_landscape_v2.RData')
library(data.table)
sample_level0 <- sherlock_cnvinfo %>% filter(CNV_Clust=="C3") %>% pull(Tumor_Barcode)

result_top <- oncoplot(data = data_top,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
#tmplevel <- result_top$gene_level
#tmplevel[13:14] <-tmplevel[c(14:13)]
#result_top <- oncoplot(data = data_top,landscape_colors = landscape_colors,sample_level0 = sample_level0,gene_level = tmplevel,GeneSortOnly = TRUE)


result_sftb <- oncoplot(data = data_sftb,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_scna <- oncoplot(data = data_scna,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_arm <- oncoplot(data = data_arm,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_fusion <- oncoplot(data = data_fusion,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_germline <- oncoplot(data=data_germline,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_feature1 <- oncoplot(data = data_feature1,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)


sample_new_level <- result_top$sample_level %>%
  left_join(
    result_sftb$sample_level
  ) %>% left_join(
    result_fusion$sample_level
  ) %>% left_join(
    result_scna$sample_level
  ) %>% left_join(
    result_arm$sample_level
  )    %>% left_join(
    result_germline$sample_level
  ) %>% left_join(
    result_feature1$sample_level
  ) %>% arrange(Mutation_Driver,Fusion,Mutation_SFTP,SCNA_Focal,SCNA_Arm,Germline_Mutation,Other_feature) %>% 
  pull(Tumor_Barcode)


result_top <- oncoplot(data = data_top,data_clone = data_clone,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,tmar = -0.2,bmar = -0.05)
result_sftb <- oncoplot(data = data_sftb,data_clone = data_sftb_clone,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
result_fusion <- oncoplot(data = data_fusion,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
result_scna <- oncoplot(data = data_scna,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
result_arm <- oncoplot(data = data_arm,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
resul_germline <- oncoplot(data=data_germline,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
result_feature1 <- oncoplot(data = data_feature1,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = 0,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_purity <- oncoplot2(data = data_purity,sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_viridis_c(),tmar=0,bmar = -0.1,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_ratio <- oncoplot2(data = data_ratio,sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_viridis_c(option = "C"),tmar=-0.1,bmar = 0,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_nrpcc <- oncoplot2(data = data_nrpcc,sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_material(palette = 'green'),tmar=-0.2,bmar =0,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_cnvratio <- oncoplot2(data = data_cnvratio,sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_material(palette = 'blue'),tmar=-0.2,bmar =0.2,p2_axis_hidden = TRUE,p2_hidden = TRUE)


result_sis <- oncoplot3(data = data_sis,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,tmar=0.2,bmar=-0.2,height = 6)

oncoplot_final <- oncoplot_combined(result_sis,result_top,result_fusion,result_sftb,result_scna,result_arm,resul_germline,result_feature1,result_purity,result_ratio,result_nrpcc,result_cnvratio)
save_plot(filename = 'Figures/Genome_landscape2_C3.pdf',plot = oncoplot_final,base_height = 16,base_width = 15,device=cairo_pdf)

