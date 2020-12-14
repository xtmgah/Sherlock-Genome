oncoplot_all_fcn <- function(file1=NULL,file2=NULL,data_clone=NULL,data_highlight=NULL,Gene_Sig=NULL,gene_level_input=list(),landscape_colors=NULL,sample_level=NULL,sample_level0=NULL,GeneSortOnly=FALSE,namemax=NULL,cell_height=0.95,tmar=0,bmar=0,nbreaks=c(0,20,40,60,80),p2_axis_hidden=TRUE,onco1_p2_hidden=FALSE,onco2_p2_hidden=TRUE,legend_adjust=FALSE,scale_fill_ztw=NULL,height=6){
  # Read in file1- contains the categorical and continuous variables in the main part of the plot
  # input_data1 <- read.csv("oncoplot_combined_practice.csv", header = TRUE, sep = ",")
  # input_data2 <- read.csv("mutations_practice_file.csv", header = TRUE, sep = ",")
  input_data1 <- read.csv(file1, header = TRUE, sep =",")
  input_data1 <- tibble(input_data1)
  
  # Read in file2- contains the continuous mutation data for the bar chart at the top of the plot
  input_data2 <- read.csv(file2, header = TRUE, sep = ",")
  input_data2 <- tibble(input_data2)
  input_data2$Alteration <- as.double(input_data2$Alteration)
  
  # Split by Type
  unique_data_types <- unique(input_data1$Type)
  # put all categorical variables into one tibble
  categorical_data <- tibble(Subject=character(), Tumor_Barcode= character(), Gene = character(), Alteration = factor(), Type = character())
  # put all continuous variables into one tibble
  continuous_data <- tibble(Subject=character(), Tumor_Barcode= character(), Gene = character(), Alteration = double(), Type = character())
  #subset by each type
  i <- 1
  for(each in unique_data_types){
    type <- unique_data_types[[i]]
    # look at the first alteration value for a subset
    x <- str_sub(assign(paste0("subset_",type), subset(input_data1, input_data1$Type== unique_data_types[[i]]))[1, "Alteration"],1,1)
    # "M" (Missense mutation)
    range <- seq(0,9,1)
    # [1] 0 1 2 3 4 5 6 7 8 9
    x %in% range
    if (x%in% range){
      continuous_data <- rbind(continuous_data, assign(paste0("subset_",type), subset(input_data1, input_data1$Type== unique_data_types[[i]])))
      continuous_data$Alteration <- as.double(continuous_data$Alteration)
      # print(continuous_data)
      # typeof(continuous_data$Alteration)    
      continuous_data <- as_tibble(continuous_data)
    }else{
      # put all categorical variables into one tibble
      categorical_data <- rbind(categorical_data, assign(paste0("subset_",type), subset(input_data1, input_data1$Type== unique_data_types[[i]])))
      categorical_data$Alteration <- as.factor(categorical_data$Alteration)
      categorical_data <- as_tibble(categorical_data)
      # print(categorical_data)
      # typeof(categorical_data$Alteration)
    }
    i <- i + 1
  }
  
  print(categorical_data)
  print(continuous_data)
  
  oncolist <- list()

  #############################################################################
  # oncoplot 2 equivalent 
  continuous_data_type <- c(unique(continuous_data$Type))
  
  i <- 1
  # for loop to go through each item in the continuous data type list 
  for(each in 1:length(continuous_data_type)){
    
    data <- subset(continuous_data, Type == continuous_data_type[[i]])
    scale_fill_ztw <- list(scale_fill_viridis_c(), scale_fill_viridis_c(option = "C"),scale_fill_material(palette = 'green'), scale_fill_material(palette = 'blue'),scale_fill_viridis_b(option = "D"), scale_fill_material(palette= 'pink'), scale_fill_viridis_b(option = "A"), scale_fill_viridis_c(option = "A"), scale_fill_material(palette = "orange"), scale_fill_material(palette = "indigo"))
    #colors_cont<- c("pink", "cyan", "deep-purple", "orange")
    # scale_fill_ztw <- scale_fill_material(palette =  cont_colors[[i]])
    # data <- subset(datasets_input, Gene == data_list[[i]])
    # limited to sample level
    data <- data %>% filter(Tumor_Barcode %in% sample_level0)
    
    samplesize <- length(sample_level0)
    #altertype <- data$Type[1]
    altertype <- data$Type[1]
    print(altertype)
    
    input_sample_level=TRUE
    
    if(is.null(sample_level)){
      input_sample_level=FALSE
      sampleorder <- data %>% 
        select(Tumor_Barcode,Gene,Alteration) %>% 
        arrange(Gene,desc(Alteration)) %>% 
        pull(Tumor_Barcode)
      
      sample_level <- c(sampleorder,sample_level0[!c(sample_level0 %in% sampleorder)])
    }
    
    gene_level <- data$Gene[1]
    
    ## resort the data 
    databg <- crossing(Tumor_Barcode=sample_level,Gene=gene_level) %>% 
      mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = unique(sample_level))) %>% 
      mutate(Gene=factor(Gene,levels = gene_level))
    data <- data %>% mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = unique(sample_level))) %>% 
      mutate(Gene=factor(Gene,levels = gene_level))
    
    # number(mean(data$Alteration[which(data$Type == "Tumor_Purity")]),accuracy = 0.01)
    alteration_mean <- number(mean(data$Alteration),accuracy = 0.01)
    print(alteration_mean)
    
    p1 <- databg %>% 
      ggplot(aes(Tumor_Barcode,as.integer(Gene))) +
      geom_tile(height=cell_height,fill="gray95")+
      theme_ipsum_rc(base_size = 10,axis_text_size = 10,grid = '',axis = '')+
      theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.title.y.right = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(face = "italic"),panel.background = element_blank(),legend.position = 'top',legend.justification = "top",legend.direction = 'vertical',legend.key.height =unit(1,"line") )+
      scale_y_continuous(expand = c(0,0),breaks = 1:length(gene_level),labels = gene_level,sec.axis = dup_axis(labels = alteration_mean))+
      ## data
      geom_tile(data = data,aes(Tumor_Barcode,as.integer(Gene),fill=Alteration),height=cell_height)+scale_fill_ztw[[i]]+
      geom_vline(xintercept = 1:samplesize-0.5,color="white",size=0.2)+
      panel_border(size = 0.3,color = 'gray70')+
      labs(fill=altertype)
    
    oncoplot_legend <- get_legend(p1+theme(plot.margin=margin(b=-12,l=0,r=0,unit="cm"),legend.key.size = unit(0.1, "cm")))
    #legend.box.margin = margin(t = 0, b = 0, l = 0, r = 0,unit = 'cm')
    p1 <- p1+theme(legend.position="none")
    
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
          mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = unique(sample_level))) #%>% 
        #mutate(Gene=factor(Gene,levels = gene_level)) 
        
        p1 <-p1+geom_rect(data = data_highlight,aes(xmin=Seq1,xmax=Seq2,ymin=-Inf,ymax=Inf),fill=NA,col="black",size=0.3)
        
      } else {
        
        data_highlight <- data_highlight %>%
          mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = unique(sample_level))) %>% 
          mutate(Gene=factor(Gene,levels = gene_level)) 
        
        p1 <-p1+geom_tile(data = data_highlight,aes(Tumor_Barcode,as.integer(Gene)),fill=NA,col="black",size=0.6)
        
      }
    }
    
    datap2 <- data
    if(onco2_p2_hidden){
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
    # have to add in the gene level/Gene somehow
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
      # colnames(rankinfo)[2] <- altertype #Warning message:In colnames(rankinfo)[2] <- altertype:number of items to replace is not a multiple of replacement length
    }
    
    print(rankinfo)
    
    result <- list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level)
    #return(list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level))
    oncolist <- list.append(oncolist, result)
    
    if(i < length(continuous_data_type)){
      i <- i + 1
    }else{
      break
    }
  }
  #####################################################################
  
  #oncoplot 1 equivalent
  # categorical data sorting,formatting, and plotting
  categorical_data_type <- c(unique(categorical_data$Type))
  # Mutation Driver, SCNA_Arm
  # loop through each Type
  i <- 1
  # for loop to go through each item in the input list 
  for(each in 1:length(categorical_data_type)){
    data <- subset(categorical_data, Type == categorical_data_type[[i]])
    # limited to sample level
    data <- data %>% filter(Tumor_Barcode %in% sample_level0)
    
    gene_level <- gene_level_input[[i]]
    
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
    
    print(altertype)
    
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
    
    print("data 0")
    print(data0)
    print("data")
    print(data)
    
    data <- as.data.table(data)
    data <- data[, strsplit(as.character(Alteration), "/"), by=list(Tumor_Barcode, Gene)]  # this expands "X/Y/Z" into three rows
    data[, shift:=(1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(Tumor_Barcode, Gene)]
    data[, height:=0.95/.N, by=list(Tumor_Barcode, Gene)]
    
    data <- as_tibble(data) %>% rename(Alteration=V1) %>% mutate(shift=shift*0.95)
    print("data")
    print(data)
    
    # define the order #
    Gene_Freq <- data0 %>% select(Gene,Tumor_Barcode) %>% unique() %>% count(Gene,sort=T) %>% mutate(n=if_else(Gene %in% blankgene,0L,n)) %>%  mutate(Freq=n/samplesize) %>% arrange(Freq) %>%  mutate(Freq=percent(Freq,accuracy = 0.1))
    print(Gene_Freq)
    
    # how can gene_level ever ve null if it is needed for creating blankgene?
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
    
    data <- data %>% mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = unique(sample_level))) %>% 
      mutate(Gene=factor(Gene,levels = gene_level))
    
    data0 <- data0 %>% mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = unique(sample_level))) %>% 
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
      geom_tile(data = data,aes(factor(Tumor_Barcode,levels = sample_level),as.integer(Gene) + shift,fill=Alteration,height=height),size=0)+ 
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
        
        p1 <-p1+geom_rect(data = data_highlight,aes(xmin=Seq1,xmax=Seq2,ymin=-Inf,ymax=Inf),fill=NA,col="black",size=0.3)
        
      } else {
        
        data_highlight <- data_highlight %>%
          mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level)) %>% 
          mutate(Gene=factor(Gene,levels = gene_level)) 
        
        p1 <-p1+geom_tile(data = data_highlight,aes(Tumor_Barcode,as.integer(Gene)),fill=NA,col="black",size=0.6, height = cell_height)
        
      }
    }
    
    if(legend_adjust){
      p1 <- p1+guides(fill = guide_legend(nrow = 1,title.position = 'left',title.hjust = 0.5,title.vjust = 0.5))+theme(legend.position = "bottom")
    }
    
    
    if(!is.null(data_clone)){
      data_clone <- data_clone %>% filter(Tumor_Barcode %in% sample_level0)
      data_clone <- data_clone %>% mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = unique(sample_level))) %>% 
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
    
    if(onco1_p2_hidden){
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
    
    if(p2_axis_hidden & i!= 1){
      p2 <- p2+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line.x = element_blank())
    }#else if (i == 1){
    #  p2_axis_hidden <- FALSE
    #}#else{
    #  p2_axis_hidden <- p2_axis_hidden
    #}
    
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
        #use only the Gene as order
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
    
    result <- list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level)
    #return(list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level))
    oncolist <- list.append(oncolist, result)
    
    print(gene_level2)
    print(oncolist)
    
    if(i < length(continuous_data_type)){
      i <- i + 1
    }else{
      break
    }
    
    #print(i)
  }  

##################################################################
  
  #oncoplot3 <- function(data,sample_level0,Gene_Sig=NULL, landscape_colors=NULL,gene_level=NULL,sample_level=NULL,bmar=0,tmar=0,height=10){
    # limited to sample level
    data <- input_data2
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
    
    result <- list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level)
    #return(list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level))
    oncolist <- list.append(oncolist, result)
    
#############################################################################
  # oncoplot combined equivalent  
  # intialized oncolist at the top so that it could be appended to as each of the plot are generated
  #oncoplot_combined <- function(...,oncolist=NULL){
  #  oncolist <- c(list(...), oncolist)
  #oncolist <- rev(oncolist)
  #print(oncolist)
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
  #return(oncoplot_final)
  
### Might need to edit this part in the future
  tmp <- sherlock_histology %>% select(Tumor_Barcode,Histology) %>% 
    left_join(result_top$sample_level) %>%
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
    ) %>% arrange(Histology,Mutation_Driver,Fusion,Mutation_SFTP,SCNA_Focal,SCNA_Arm,Germline_Mutation,Other_feature) 
  tmp %>% write_csv(tmp,path = 'sample_new_level_test.csv',col_names = T,na = '')
  #tmp <- read_csv('sample_new_level_test.csv',col_names = T)
  #sample_new_level <- tmp %>%   pull(Tumor_Barcode)
  
  
  return(oncoplot_final)
}

#Error in oncoplot_all_fcn(file1 = file1, file2 = file2, data_clone = NULL,  : 
#                            formal argument "data_clone" matched by multiple actual arguments
oncoplot_combined_plot <-oncoplot_all_fcn(file1="oncoplot_combined_practice.csv", file2="mutations_practice_file.csv",gene_level_input=gene_level_input,sample_level0=sherlock_samples_unique$Tumor_Barcode,sample_level=NULL,namemax=NULL,cell_height=0.95,data_highlight=NULL,legend_adjust=FALSE,landscape_colors=landscape_colors,Gene_Sig=NULL,GeneSortOnly = TRUE)  
tmp <- read_csv('sample_new_level_test.csv', col_names = T)
sample_new_level <- tmp %>% pull(Tumor_Barcode)
oncoplot_combined_plot <-oncoplot_all_fcn(file1="oncoplot_combined_practice.csv", file2="mutations_practice_file.csv",data_clone=data_clone,gene_level_input=gene_level_input,sample_level0=sherlock_samples_unique$Tumor_Barcode,sample_level=sample_new_level,namemax=NULL,cell_height=0.95,data_highlight=NULL,legend_adjust=FALSE,landscape_colors=landscape_colors,Gene_Sig=NULL,GeneSortOnly = TRUE)  

save_plot(filename = 'Genome_landscapeTest2b.pdf',plot = oncoplot_combined_plot,base_height = 18,base_width = 16) #device=cairo_pdf)

