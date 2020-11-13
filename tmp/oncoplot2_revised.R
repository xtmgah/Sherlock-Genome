oncoplot2 <- function(data_purity= c(),data_ratio= c(),data_nrpcc= c(),data_cnvratio= c(),data_highlight = NULL,sample_level0,Gene_Sig=NULL, scale_fill_ztw=scale_fill_viridis_c(option = "C"),gene_level=NULL,sample_level=NULL,bmar=0,tmar=0,p2_axis_hidden=FALSE,p2_hidden=FALSE){
  
  #initialize list for data sets that were input
  datasets_input <-tibble(Subject=character(), Tumor_Barcode= character(), Gene = character(), Alteration = double(), Type = character(), stringsAsFactors = FALSE)
  #initialize list for loop to take in each dataset
  data_list <- list()
  
  if(!is.null(data_purity)){
    data_list <- append(data_list,  'Tumor_Purity')
    #append the data to the tibble
    datasets_input <- rbind(datasets_input, data_purity)
    
  }
  
  if(!is.null(data_ratio)){
    data_list <- append(data_list,  'Subclonal_Ratio')
    #append the data to the tibble
    datasets_input <- rbind(datasets_input, data_ratio)
    
  }
  
  if(!is.null(data_nrpcc)){
    data_list <- append(data_list,  'NRPCC')
    #append the data to the tibble
    datasets_input <- rbind(datasets_input, data_nrpcc)
    
  }
  
  if(!is.null(data_cnvratio)){
    data_list <- append(data_list,  'SCNA_Ratio')
    #append the data to the tibble
    datasets_input <- rbind(datasets_input, data_cnvratio)
    
  }
  
  i <- 1
  # for loop to go through each item in the input list 
  for(each in 1:length(data_list)){
    data <- subset(datasets_input, Gene == data_list[[i]])
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
          mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level)) %>% 
          mutate(Gene=factor(Gene,levels = gene_level)) 
        
        p1 <-p1+geom_rect(data = data_highlight,aes(xmin=Seq1,xmax=Seq2,ymin=-Inf,ymax=Inf),fill=NA,col="black",size=0.3)
        
      } else {
        
        data_highlight <- data_highlight %>%
          mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level)) %>% 
          mutate(Gene=factor(Gene,levels = gene_level)) 
        
        p1 <-p1+geom_tile(data = data_highlight,aes(Tumor_Barcode,as.integer(Gene)),fill=NA,col="black",size=0.6)
        
      }
    }
    
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
    
    # store each output to then put in result to return
    if(i == 1){
      result_one <- list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level)
      #print(result_one)
    }else if(i == 2){
      result_two <- list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level)
      #print(result_two)
    }else if(i == 3){
      result_three <- list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level)
      #print(result_three)
    }else{
      result_four <- list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level)
      #print(result_four)
    }
    #return(list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level))
    i <- i + 1
    #print(i)# store each output to then put in result to return
  }
  
  result_onco2_all <- list(result_one, result_two, result_three, result_four)
  #result_one
  #result_two
  #result_three
  #result_four
  
  
}

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


#tmp <- sherlock_histology %>% select(Tumor_Barcode,Histology) %>% 
#  left_join(result_top$sample_level) %>%
#  left_join(
#    result_sftb$sample_level
#  ) %>% left_join(
#    result_fusion$sample_level
#  ) %>% left_join(
#    result_scna$sample_level
#  ) %>% left_join(
#    result_arm$sample_level
#  )    %>% left_join(
#    result_germline$sample_level
#  ) %>% left_join(
#    result_feature1$sample_level
#  ) %>% arrange(Histology,Mutation_Driver,Fusion,Mutation_SFTP,SCNA_Focal,SCNA_Arm,Germline_Mutation,Other_feature) 
tmp %>% write_csv(tmp,file = 'sample_new_level.csv',col_names = T,na = '')
tmp <- read_csv('sample_new_level.csv',col_names = T)
sample_new_level <- tmp %>%   pull(Tumor_Barcode)


result_top <- oncoplot(data = data_top,data_clone = data_clone, data_highlight = data_highlight,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,tmar = -0.2,bmar = -0.05,gene_level = tmplevel)
result_sftb <- oncoplot(data = data_sftb,data_clone = data_sftb_clone,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
result_fusion <- oncoplot(data = data_fusion,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
result_scna <- oncoplot(data = data_scna,data_highlight = data_highlight,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
result_arm <- oncoplot(data = data_arm,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
resul_germline <- oncoplot(data=data_germline,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,p2_axis_hidden = TRUE)
result_feature1 <- oncoplot(data = data_feature1,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,gene_level = featurelevel,bmar = 0,p2_axis_hidden = TRUE,p2_hidden = TRUE)
#result_purity <- oncoplot2(data = data_purity,data_highlight = data_highlight, sample_level0 = sample_level0,sample_level = sample_new_level,tmar=0,bmar = -0.1,p2_axis_hidden = TRUE,p2_hidden = TRUE)
#result_ratio <- oncoplot2(data = data_ratio,sample_level0 = sample_level0,sample_level = sample_new_level,tmar=-0.2,bmar = -0.1,p2_axis_hidden = TRUE,p2_hidden = TRUE)
#result_nrpcc <- oncoplot2(data = data_nrpcc,sample_level0 = sample_level0,sample_level = sample_new_level,tmar=-0.2,bmar =0,p2_axis_hidden = TRUE,p2_hidden = TRUE)
#result_cnvratio <- oncoplot2(data = data_cnvratio,sample_level0 = sample_level0,sample_level = sample_new_level,tmar=-0.2,bmar =0,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_oncoplot2_revised <- oncoplot2(data_purity = data_purity, data_ratio= data_ratio, data_nrpcc = data_nrpcc, data_cnvratio = data_cnvratio, sample_level0 = sample_level0,sample_level = sample_new_level,tmar=0,bmar = -0.1,p2_axis_hidden = TRUE,p2_hidden = TRUE)

result_sis <- oncoplot3(data = data_sis, data_highlight = data_highlight, landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,tmar=0.2,bmar=-0.2,height = 6)

oncoplot_final <- oncoplot_combined(result_sis,result_top,result_fusion,result_sftb,result_scna,result_arm,resul_germline,result_feature1,result_oncoplot2_revised[[1]], result_oncoplot2_revised[[2]], result_oncoplot2_revised[[3]], result_oncoplot2_revised[[4]])#result_purity,result_ratio,result_nrpcc,result_cnvratio)
save_plot(filename = 'Genome_landscape_test.pdf',plot = oncoplot_final,base_height = 18,base_width = 16) #device = cairo_pdf
