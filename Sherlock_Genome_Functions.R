
load('./www/Genomic Data/Sherlock/Integrative_Analysis/sherlock_data_all.RData')
load('./www/Genomic Data/Sherlock/Survival_Analysis/suvdata.RData')
#load('./test.RData')
load('covdata0.RDS')
load('./www/Genomic Data/Sherlock/Mutations/sherlock_maf.RData')

mdata0 <- sherlock_data_full %>% 
  mutate(Gene=paste0(Gene,"|",Type)) %>% 
  select(Tumor_Barcode,Gene,Alteration) %>% 
  pivot_wider(names_from = "Gene",values_from = "Alteration")

oncoplot_colors <- read_csv('oncoplot_colors.csv')

pdfhra <- function(){
  d <- read.csv(extrafont:::fonttable_file(), stringsAsFactors = FALSE)
  d[grepl("Light", d$FontName), ]$FamilyName <- font_rc_light
  write.csv(d, extrafont:::fonttable_file(), row.names = FALSE)
}

roboto_font_import <- function(){
  hrbrthemes::import_roboto_condensed()
  library(hrbrthemes)
  pdfhra()
  extrafont::loadfonts()
}


# pdfhr2 <- function (...) 
# {
#   require(showtext)
#   font_add_google(name = "Roboto Condensed", family = "Roboto Condensed", 
#                   regular.wt = 400, bold.wt = 700)
#   showtext_auto()
#   showtext_opts(dpi = 112)
# }


# os_detect function ----------------------------------------------------------------
# detect os to apply correct file paths
os_detect <- function(){
  os_name <- Sys.info()[['sysname']]
  return(os_name)
}

# read_colnames function ----------------------------------------------------------------
# read in column names
read_colnames <- function(filename, sep="\t",header=TRUE){
  file <- read.delim(file=filename)
  return(colnames(file))
}

# read_in_file function ----------------------------------------------------------------
read_in_file <- function(filename, sep="\t",header=TRUE){
  file <- read.delim(file=filename,sep="\t", header=TRUE)
  return(file)
}

# files_list function ----------------------------------------------------------------
files_list <- function(...){
  dirs <- list.dirs(full.names=FALSE,recursive = FALSE)
  i <- 1
  files_each_dir <- list()
  files_ind_dir <- list()
  for(each in dirs){
    files_ind_dir[[i]] <- list()
    sub_dir <- list.files(list.dirs(full.names=FALSE,recursive = FALSE)[[i]])
    sub_dir_files <- list.files(paste0("./",dirs[[i]],"/",sub_dir))
    # sub_dir <- dirs[[i]] %>% list.files(recursive = TRUE, full.names = FALSE)
    if(length(sub_dir_files) >0){
      sub_dir <- c(sub_dir, sub_dir_files)
    }
    
    files_ind_dir[[i]] <- append(files_ind_dir[[i]], sub_dir)
    files_each_dir <- append(files_each_dir, files_ind_dir[[i]])
    i <- i + 1
  }
  
  # show the user only the files that are actually loaded (for example, no png files, pdf, etc); add to list as files are added
  module_files <- c("all_ngspurity_output.txt","MutationTime_Proportion.txt","QC_sample_level.txt","QC_subject_level.txt","Study_Overview.Rmd","genomePlot_list.txt")
  i <- 1
  for(each in files_ind_dir){
    if(length(files_ind_dir[[i]]) != 0){
      files_ind_dir[[i]] <- files_ind_dir[[i]][which(files_ind_dir[[i]] %in% module_files)]
      files_ind_dir[[i]]
    }else{
      # print("No files available for that module at this time.")
    }
    
    i <- i + 1
    # print(i)
  }
  
  return(files_ind_dir)
}

# sherlock_genome_filter function ----------------------------------------------------------------
# filter a dataframe given user input conditions
sherlock_genome_filter <- function(data, conditions, column_selection){ #sherlock_genome_filter(data=data, conditions) conditions <- "MCN_WGD == 'nWGD'"
  library(stringr)
  
  #conditions <- "MCN_WGD == 'nWGD'"
  #conditions <- c("MCN_WGD == 'nWGD'" , "PGA >= 0.52")
  # conditions <- c("MCN_WGD == 'nWGD'" "PGA >= 0.52",head(100)) # might want to somehow add this?
  
  # print(conditions)
  if(str_detect(conditions, ",")){
    condition_split <- str_split(conditions, ",")
    # print(condition_split)
    i <- 1
    filter_list <- c()
    #filter_condition <- data %>% filter( PGA >= 0.52) %>% filter(MCN_WGD == 'nWGD')
    for(each in condition_split[[1]]){
      # print(condition_split[[1]][i])
      #test[i] <- paste0("filter(", condition_split[[1]][i], ")")
      filter_list <- append(filter_list,paste0("filter(", condition_split[[1]][i], ")"))
      #print(test[i])
      #filter_list <- append(filter_list, test[i])
      i <- i + 1
      # print(i)
    }
    # print(filter_list)
    filter_list <- paste0(filter_list, collapse= " %>% ") # puts %>% between each filter
    filter_df <- paste0("filter_condition <- data %>% ",filter_list)
    
  }else{
    filter_df <- paste0("filter_condition <- data %>% filter(",conditions, ")") # %>% select(",column_selection, ")")
    print(filter_df)
  }
  
  # print(column_selection)
  # if(str_detect(column_selection, ",")){
  #   column_select_split <- str_split(column_selection, ",'")
  #   print(column_select_split)
  if(length(column_selection > 1)){
    i <- 1
    column_list <- c()
    for(each in column_selection){
      print(column_selection[i])
      column_list <- append(column_list,column_selection[i])
      i <- i + 1
    }
    
    column_list <- paste0(column_list, collapse= ",")
    column_list <- (paste0("select(",column_list,")"))
    print(column_list)
  }else{
    column_list <- (paste0("select(",column_selection,")"))
  }
  
  
  
  filter_df <- paste0("filter_condition <- data %>% filter( ",conditions, ") %>%", column_list)
  #print(filter_df)
  eval(parse(text=filter_df))
  
  return(filter_condition)
  
}

# change_data_type function ----------------------------------------------------------------
change_data_type <- function(data){
  if(is.numeric(data)){
    data <- as.character(data)
  }else if(is.character(data)){
    data <- as.numeric(data)
  }
  return(data)
}

# validate_vardf function ----------------------------------------------------------------
validate_vardf <- function(data, forces=NULL, nachars=c("","na","Na","NA","nan","NAN","Nan"), Nmin=5, Nmin_drop=FALSE,excludes=NULL, lump=FALSE){
  #names_oringal <- colnames(data)
  # force data type
  if(!is.null(excludes)){
    names_all <- colnames(data)
    data_excludes <- data  %>% dplyr::select(one_of(excludes))
    data <- data %>% dplyr::select(-one_of(excludes))
  }
  
  if(!is.null(forces)){ 
    data <- data %>% mutate(across(one_of(forces),change_data_type))  
  }
  
  # process na
  data <- data %>% mutate(across(where(is.character), ~  if_else(.x %in% nachars, NA_character_, .x)))
  
  # process the integer
  data <- data %>% mutate(across(where(is.integer), as.numeric))
  
  # convert numbers to characters if less than xx unique value
  #data %>% mutate(across(where(is.numeric), ~ if_else(n_distinct(.) < Nmin, .x, .x)))
  Nmin_names <- data %>% summarise(across(where(is.numeric),n_distinct)) %>% dplyr::select_if(function(x) x<Nmin) %>% colnames()
  if(length(Nmin_names)>0){
    if(Nmin_drop){
      data <- data %>% dplyr::select(-one_of(Nmin_names))
      #names_oringal <- names_oringal[!(names_oringal %in% Nmin_names)]
    }
    else{
      data <- data %>% mutate(across(one_of(Nmin_names),change_data_type)) 
    }
  }
  
  ## make the factors and reorder the levels
  ## exclude for the all unique value columns
  #data <- data %>% mutate(across(where(is.character),~ fct_lump(fct_infreq(as.factor(.x)),prop = 0.2)))
  # if(lump){
  #   data <- data %>% mutate(across(where(~ is.character(.) & (n_distinct(.) < 0.8*n())),~ fct_lump(fct_infreq(as.factor(.x)),prop = 0.2)))
  # }
  if(lump){
    nrows <- dim(data)[1]
    data <- data %>% mutate(across(where(~ is.character(.) & (n_distinct(.)>1) & (n_distinct(.) < 0.8*nrows)),~ fct_lump(fct_infreq(as.factor(.x)),prop = 0.2)))
  }
  # if change the order
  #data <- data %>% select(names_oringal)
  
  if(!is.null(excludes)){
    names_keep <- c(colnames(data),colnames(data_excludes))
    names_all <- names_all[names_all %in% names_keep]
    data <- bind_cols(data_excludes,data) %>% select(one_of(names_all))
  }
  
  return(as_tibble(data))
}

# df_columns_subset_function(dataframe, col_names))
df_columns_subset_function <- function(dataframe, col_names= NULL){
  if("All columns" %in% col_names){
    data <- as.data.frame(dataframe)
  }
  if(!("All columns" %in% col_names)){
    data <- as.data.frame(dataframe[,col_names])
  }
  data <- validate_vardf(data)
  
  return(data)
  
}


# inspect_data_function function ----------------------------------------------------------------
# use inspectdf package to explore data
inspect_data_function <- function(dataframe, type_of_inspection=NULL,column_name=NULL ){
  if(column_name == "All columns"){
    data <- as.data.frame(dataframe)
  }
  
  if(column_name != "All columns"){
    data <- as.data.frame(dataframe[,column_name])
  }
  
  data <- validate_vardf(data)
  
  # inspect_na()
  if(type_of_inspection =="na"){
    x <- inspect_na(data)
    x <- show_plot(x)
  }
  # inspect_cat()
  if(type_of_inspection == "cat"){
    x <- inspect_cat(data)
    x <- show_plot(x)
  } 
  # inspect_imb()
  if(type_of_inspection == "cat_levels"){
    x <- inspect_imb(data)
    x <- show_plot(x)
  }
  # inspect_num()
  if(type_of_inspection == "num"){
    #dataframe <- as.data.frame(as.numeric(unlist(dataframe)))
    x <- inspect_num(data)
    x <- show_plot(x)
  }
  # inspect_types()
  if(type_of_inspection == "types"){
    x <- inspect_types(data)
    x <- show_plot(x)
  }
  return(x)
}

# figure_display_ngspurity function ----------------------------------------------------------------
# output figures for ngspurity data -- needs to be for others in the future?
figure_display_ngspurity <- function(tumor_barcode=NULL,battenberg=NULL,type=NULL,project_code=NULL){
  ngspurity_qc <- read_in_file("NGSpurity/ngspurity_qc_file.txt")
  # pull the filename from the ngspurity_qc_file table
  file_path <- ngspurity_qc$File[which(ngspurity_qc$Tumor_Barcode==tumor_barcode & ngspurity_qc$Battenberg==battenberg & ngspurity_qc$Type==type)]

  if(os_detect() %in% c("Linux","Darwin")){
    filename <- sub(".","NGSpurity", file_path)
    print(getwd())
    print(filename)
    list(src = filename, alt = paste0(tumor_barcode, "_", battenberg, "_", type))
  }else{
    filename <- sub(".",paste0("Genomic Data/", project_code,"/NGSpurity"), file_path)
    print(getwd())
    print(filename)
    return(filename)
  }
}


figure_display_mutationTime <- function(tumor_barcode=NULL){
  mutationTime_file <- read_in_file("Clonal_Evolution/MutationTime/MutationTime_Proportion.txt")
  # pull the filename from the ngspurity_qc_file table
  mut_barcode <- mutationTime_file$Tumor_Barcode[which(mutationTime_file$Tumor_Barcode==tumor_barcode)]
  print(mut_barcode)
  file_path <- paste0("Clonal_Evolution/MutationTime/", mut_barcode, "_MTime.pdf")
  print(file_path)
  
  if(os_detect() %in% c("Linux","Darwin")){
    filename <- file_path
    print(getwd())
    # print(filename)
    list(src = filename, alt = paste0(tumor_barcode, "_MTime.pdf"))
  }else{
    filename <- sub(".",paste0("Genomic Data/", project_code,"/Clonal_Evolution/"), file_path)
    print(getwd())
    print(filename)
    return(filename)
  }
}

#fisher_result_test_run <- fishergroup(vartmp=vartmp, sp_group=NULL, samplelist=NULL, var2name="WGD_Status|Overall_Feature", excludes= NULL, excludes_cat= NULL, keeps= NULL, keeps_cat= NULL, minfreq= 0.03, freq_column= 'Freq', method= "fisher.test", glm_formula= "Var1 ~ Var2 + Gender", covdata= covdata0, subfold="TMP", anndata=NULL, fdrcutoff=0.1)
#fisher_result_test_run_glm <- fishergroup(vartmp=vartmp, sp_group=NULL, samplelist=NULL, var2name="WGD_Status|Overall_Feature", excludes= NULL, excludes_cat= NULL, keeps= NULL, keeps_cat= NULL, minfreq= 0.03, freq_column= 'Freq', method= "glm", glm_formula= "Var1 ~ Var2 + Gender", covdata= covdata0, subfold="TMP", anndata=NULL, fdrcutoff=0.1)

# fishergroup function ----------------------------------------------------------------
fishergroup <- function(vartmp, sp_group=NULL, samplelist=NULL, var2name = NULL, excludes = NULL,excludes_cat = NULL, keeps = NULL,keeps_cat = NULL, minfreq=0.03, freq_column='Freq',  method = "fisher.test", glm_formula = "Var1 ~ Var2 + Gender", covdata = covdata0, subfold="TMP",anndata=NULL, fdrcutoff=0.1){
  
  # vartmp <- c('Tumor_Barcode', 'HR_Status|Overall_Feature')
  # sp_group = NULL
  # samplelist = NULL
  # var2name = "WGD_Status|Overall_Feature"
  # excludes = NULL
  # excludes_cat =NULL
  # keeps =NULL
  # keeps_cat = NULL
  # minfreq-0.03
  # minfreq=0.03
  # freq_column= 'Freq'
  # method = "fisher.test" method = "glm"
  # glm_formula = "Var1 ~ Var2 + Gender"
  # fdrcutoff=0.1
  # covdata = covdata0
  # subfold = 'TMP'
  # anndata= NULL

  #pdfhr2()
  
  vartmp <- c('Tumor_Barcode',vartmp)
  
  mdata <- mdata0
  
  if(!is.null(sp_group)){
    tmpxx <- sherlock_overall %>% filter(SP_Group %in% sp_group) %>% pull(Tumor_Barcode)
    mdata <- mdata0 %>% filter(Tumor_Barcode %in% tmpxx)
  }
  
  if(!is.null(samplelist)){
    mdata <- mdata0 %>% filter(Tumor_Barcode %in% samplelist)
  }
  
  mdata <- mdata %>% 
    pivot_longer(cols = -one_of(vartmp)) %>% 
    drop_na() 
  
  colnames(mdata)[2] <- 'Var1'
  colnames(mdata)[4] <- 'Var2'
  
  if(!is.null(excludes)){
    mdata <- mdata %>% filter(!name %in% excludes)
  }
  
  if(!is.null(excludes_cat)){
    mdata <- mdata %>% mutate(TMP=name) %>% separate(col = TMP,into = c('TMP1','TMP2'),sep = '\\|') %>% filter(!TMP2 %in% excludes_cat) %>% select(-TMP1,-TMP2)
  }
  
  if(!is.null(keeps)){
    mdata <- mdata %>% filter(name %in% keeps)
  }
  
  if(!is.null(keeps_cat)){
    mdata <- mdata %>% mutate(TMP=name) %>% separate(col = TMP,into = c('TMP1','TMP2'),sep = '\\|') %>% filter(TMP2 %in% keeps_cat) %>% select(-TMP1,-TMP2)
  }
  
  tmp <- mdata %>% select(Var1,name,Var2) %>% unique() %>% count(name) %>% filter(n<3) %>% pull(name)
  
  if(method == "fisher.test") {
    ## for fisher exact test 
    result <- mdata %>% 
      filter(!(name %in% tmp)) %>% 
      group_by(name) %>% 
      do(tidy(fisher.test(.$Var1,.$Var2))) %>% 
      ungroup() %>% 
      select(name:conf.high) %>% 
      filter(!is.na(p.value)) %>% 
      arrange(p.value) %>% 
      mutate(estimate = log2(estimate))
    
    xlabtmp <- "Fisher's Exact Test, log2(OR)"
  }
  
  if(method == "glm") {
    ## for logistic linear regression, adjusted for age, gender, histology, 
    glm_formula = as.formula(glm_formula)
    result <-  mdata %>%
      filter(!(name %in% tmp)) %>% 
      left_join(covdata) %>%
      group_by(name) %>%  
      mutate(Var1= as.integer(as.factor(Var1))-1) %>% 
      do(tresult = safely(glm)(glm_formula, data=. )) %>% 
      mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
      filter(!tresult_null) %>% 
      mutate(fit = list(tidy(tresult[['result']]))) %>% 
      select(name,fit) %>% 
      unnest(cols = c(fit)) %>% 
      ungroup() %>% 
      filter(str_detect(term,'^Var2')) %>% 
      arrange(p.value)
    
    xlabtmp <- 'Logistic Regression Coefficient'
    
  }
  
  if(is.null(sp_group)|length(sp_group)>1){
    tmp_freq <- sherlock_freq %>% select(name,one_of(freq_column))
  }else{
    tmp_freq <- sherlock_freq %>% select(name,contains(sp_group))
  }
  colnames(tmp_freq)[2] <- 'Freq'
  
  result <- result %>% mutate(fdr = p.adjust(p.value,method = 'BH')) %>% left_join(tmp_freq)

  ## fdr cutoff
  tmp <- result %>% filter(fdr<fdrcutoff) %>% dim() %>% .[[1]]
  if(tmp>0){
    fdrline <- mean(result$p.value[c(tmp,tmp+1)])
  }else{
    fdrline <- 0
  }
  
  #tmpsize <- mdata %>% count(Tumor_Barcode) %>% dim() %>% .[[1]]
  # tmpsize <- mdata %>% count(Tumor_Barcode,name)  %>% count(name) %>% rename(size=n)
  # tmpfreq <- mdata %>% count(name,Var2) %>% filter(!is.na(Var2)) %>% left_join(tmpsize)%>% mutate(Freq=n/size) %>% group_by(name) %>% arrange(Freq) %>% slice(1) %>% ungroup() %>% select(name,Freq) %>% mutate(Freq=if_else(Freq==1,0,Freq))
  
  tmp <- result %>% mutate(xx=name) %>% separate(col = xx,into = c('Gene','Type'),sep = '\\|') %>% filter(Freq>minfreq)
  
  anncoloms <- c('name','Gene','Type','estimate','p.value','Freq')
  
  if(!is.null(anndata)){
    tmp2 <- anndata
  }else {
    tmp2 <- tmp %>% filter(is.finite(estimate),Freq>0.3 & p.value<0.05) %>% arrange(p.value) %>% slice(1:20) %>% ungroup()%>% select(one_of(anncoloms))
    tmp2 <- tmp %>% filter(is.finite(estimate)) %>% group_by(estimate>0) %>% arrange(p.value) %>% slice(1:20) %>% ungroup() %>% select(one_of(anncoloms)) %>% bind_rows(tmp2) %>% unique()
    tmp2 <- tmp %>% filter(is.finite(estimate)) %>% group_by(estimate>0) %>% arrange(desc(abs(estimate))) %>% slice(1:5) %>% ungroup() %>% select(one_of(anncoloms))%>% bind_rows(tmp2) %>% unique()
    tmp2 <- tmp %>% filter(is.infinite(estimate)) %>% group_by(estimate>0|is.infinite(estimate)) %>% arrange(p.value) %>% slice(1:5) %>% ungroup() %>% select(one_of(anncoloms)) %>% bind_rows(tmp2) %>% unique() %>% arrange(p.value)
    
  }
  
  
  # output prefix
  prefix <- str_replace_all(vartmp[2],"[/ ]","_")
  prefix <- paste0(subfold,'/',str_remove(prefix,'\\|.*'))

  if(!is.null(sp_group)){
    prefix <- paste0(prefix,"_",paste0(sp_group,collapse = '_'))
  }
  
  tmp_plot1 <- tmp %>% 
    ggplot(aes((estimate),-log10(p.value),fill=Type))+
    geom_hline(yintercept = -log10(0.05),linetype=2,col='red',size=0.5)+
    geom_hline(yintercept = -log10(fdrline),linetype=2,col='#006d2c',size=0.5)+
    geom_point(aes(size=Freq),pch=21,stroke=0.2)+
    scale_fill_manual(values = sherlock_type_colors[sort(unique(tmp$Type))])+
    scale_size_binned()+
    ggrepel::geom_text_repel(data=tmp2,aes(label=Gene),max.overlaps = 30)+
    labs(x=xlabtmp,y='-log10(P-value)')+
    guides(fill = guide_legend(override.aes = list(size=3.5)))+
    theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid="XY",ticks = T)+
    panel_border(color = 'black',size = 0.5)+
    coord_cartesian(clip = 'off')
  
  # all_results <- list.append(all_results, tmp_plot1)
  #ggsave(paste0(prefix,'_enriched0.pdf'),width = 14,height = 10,device = cairo_pdf)
  tmp_plot2 <- tmp %>% 
    filter(estimate!=0,is.finite(estimate)) %>% 
    ggplot(aes((estimate),-log10(p.value),fill=Type))+
    geom_hline(yintercept = -log10(0.05),linetype=2,col='red',size=0.5)+
    geom_hline(yintercept = -log10(fdrline),linetype=2,col='#006d2c',size=0.5)+
    geom_point(aes(size=Freq),pch=21,stroke=0.2)+
    scale_fill_manual(values = sherlock_type_colors[sort(unique(tmp$Type))])+
    scale_size_binned()+
    ggrepel::geom_text_repel(data=tmp2 %>% filter(estimate!=0,is.finite(estimate)) ,aes(label=Gene),max.overlaps = 30)+
    labs(x=xlabtmp,y='-log10(P-value)')+
    guides(fill = guide_legend(override.aes = list(size=3.5)))+
    theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid="XY",ticks = T)+
    panel_border(color = 'black',size = 0.5)+
    coord_cartesian(clip = 'off')
  
  # all_results <- list.append(all_results, tmp_plot2)
  #ggsave(paste0(prefix,'_enriched.pdf'),width = 14,height = 10,device = cairo_pdf)

 ### proportion plot ### 
 
  if(!is.null(var2name)){
    
    var1name = vartmp[2]
    # var2name = 'Kataegis|Overall_Feature'
    # var2name = '12q15|SCNA_Focal_Cytoband'
    # var2name = 'SBS288L|Signature_Denovo'
    # var2name = 'EGFR|Mutation_Driver'
    # var2name = 'SBS288I|Signature_Denovo'
    mdata <- mdata %>% filter(name == var2name)
    vartmp1 <- mdata %>% group_by(Var1) %>% tally() %>% dplyr::mutate(percent=percent_format(accuracy = 0.1)(n/sum(n))) %>% mutate(Var1_Lab=paste0(Var1,' (',percent,')')) %>% select(Var1,Var1_Lab)
    vartmp2 <- mdata %>% group_by(Var2) %>% tally() %>% dplyr::mutate(percent=percent_format(accuracy = 0.1)(n/sum(n))) %>% mutate(Var2_Lab=paste0(Var2,' (',percent,')')) %>% select(Var2,Var2_Lab)
    var_prop_plot <- mdata %>% 
      group_by(Var1,Var2) %>% 
      dplyr::tally()%>%
      dplyr::mutate(percent=n/sum(n)) %>% 
      ungroup() %>% 
      left_join(vartmp1) %>% 
      left_join(vartmp2) %>% 
      ggplot(aes(x=Var1_Lab, y=n, fill=Var2_Lab))+
      geom_bar(stat="identity", position ="fill",width = 0.8)+
      geom_text(aes(label=paste0(n,'\n',sprintf("%1.1f", percent*100),"%")), position=position_fill(vjust=0.5), colour="white")+
      scale_y_continuous(breaks = pretty_breaks(),labels = percent_format(),expand = c(0,0))+
      scale_fill_jama()+
      theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,axis_text_size = 10,grid = 'Yy',ticks = FALSE)+
      #panel_border(color = 'black')+
      labs(fill=var2name,x=paste0('\n',var1name), y = 'Percentage')
    
    var2name <- str_replace_all(var2name,"[/ ]","_")
    #ggsave(filename = paste0(prefix,'_enriched_',str_remove(var2name,'\\|.*'),'.pdf'),width = 5,height = 8,device = cairo_pdf )
  }
  
   # return(all_results)
   if(!is.null(var2name)){
     return(list(result,tmp_plot1,var_prop_plot))
   }else{
     return(list(result,tmp_plot1))
   }
}

# fisherbarplot function ----------------------------------------------------------------
fisherbarplot <- function(vartmp, sp_group=NULL, samplelist=NULL, var2name = NULL,subfold="TMP"){
  
  #pdfhr2()
  mdata <- mdata0
  
  if(!is.null(sp_group)){
    tmpxx <- sherlock_overall %>% filter(SP_Group %in% sp_group) %>% pull(Tumor_Barcode)
    mdata <- mdata0 %>% filter(Tumor_Barcode %in% tmpxx)
  }
  
  if(!is.null(samplelist)){
    mdata <- mdata0 %>% filter(Tumor_Barcode %in% samplelist)
  }
  
  
  mdata <- mdata %>% select(Tumor_Barcode,one_of(c(vartmp,var2name))) %>% drop_na() 
  
  colnames(mdata)[2] <- 'Var1'
  colnames(mdata)[3] <- 'Var2'
  
  
  # output prefix
  prefix <- str_replace_all(vartmp,"[/ ]","_")
  prefix <- paste0(subfold,'/',str_remove(prefix,'\\|.*'))
  
  print(prefix)
  
  if(!is.null(sp_group)){
    prefix <- paste0(prefix,"_",paste0(sp_group,collapse = '_'))
  }
  
  
  ### proportion plot ### 
  var1name = vartmp
  # var2name = 'Kataegis|Overall_Feature'
  # var2name = '12q15|SCNA_Focal_Cytoband'
  # var2name = 'SBS288L|Signature_Denovo'
  # var2name = 'EGFR|Mutation_Driver'
  # var2name = 'SBS288I|Signature_Denovo'
  #mdata <- mdata %>% filter(name == var2name)
  vartmp1 <- mdata %>% group_by(Var1) %>% tally() %>% dplyr::mutate(percent=percent_format(accuracy = 0.1)(n/sum(n))) %>% mutate(Var1_Lab=paste0(Var1,' (',percent,')')) %>% select(Var1,Var1_Lab)
  vartmp2 <- mdata %>% group_by(Var2) %>% tally() %>% dplyr::mutate(percent=percent_format(accuracy = 0.1)(n/sum(n))) %>% mutate(Var2_Lab=paste0(Var2,' (',percent,')')) %>% select(Var2,Var2_Lab)
  fisher_barplot <- mdata %>% 
    group_by(Var1,Var2) %>% 
    dplyr::tally()%>%
    dplyr::mutate(percent=n/sum(n)) %>% 
    ungroup() %>% 
    left_join(vartmp1) %>% 
    left_join(vartmp2) %>% 
    ggplot(aes(x=Var1_Lab, y=n, fill=Var2_Lab))+
    geom_bar(stat="identity", position ="fill",width = 0.8)+
    geom_text(aes(label=paste0(n,'\n',sprintf("%1.1f", percent*100),"%")), position=position_fill(vjust=0.5), colour="white")+
    scale_y_continuous(breaks = pretty_breaks(),labels = percent_format(),expand = c(0,0))+
    scale_fill_jama()+
    theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,axis_text_size = 10,grid = 'Yy',ticks = FALSE)+
    #panel_border(color = 'black')+
    labs(fill=var2name,x=paste0('\n',var1name), y = 'Percentage')
  
  var2name <- str_replace_all(var2name,"[/ ]","_")
  #ggsave(filename = paste0(prefix,'_enriched_',str_remove(var2name,'\\|.*'),'.pdf'),width = 5,height = 8,device = cairo_pdf )
  
  return(fisher_barplot)
}

join_assoc_data <- function(df_list= list()){ #list()
  # only works for data sets that contain 'Tumor_Barcode' at the moment
  # data_qc_sample and data --> join by "Tumor_Barcode"
  # data_qc_subject --> Barcode (includes Tumor_Barcodes and Normal_Barcode in the same column)
  
  #joined_df <- list(data_qc_sample,data) %>% reduce(left_join, by= "Tumor_Barcode")
  # i <- 1
  # df_list_tum_barcode <- list()
  # df_list_no_tum_barcode <- list()
  # for(each in df_list){
  #   if("Tumor_Barcode" %in% colnames(data.frame(df_list[i]))){
  #     # print(paste0("Yes Tumor_Barcode",i))
  #     # df_list_tum_barcode <- list.append(df_list[i], df_list_tum_barcode)
  #     df_list_tum_barcode <- append(df_list_tum_barcode, df_list[i], after = length(df_list_tum_barcode))
  #     # print(df_list_tum_barcode)
  #   }else{
  #     # print(paste0("No Tumor_Barcode",i))
  #     # df_list_no_tum_barcode <- list.append(df_list[i], df_list_no_tum_barcode)
  #     df_list_no_tum_barcode <- append(df_list_no_tum_barcode, df_list[i], after = length(df_list_no_tum_barcode))
  #     # print(df_list_no_tum_barcode)
  #   }
  #   i <- i + 1
  #   # print(i)
  # }
  joined_df <- df_list %>% reduce(left_join, by= "Tumor_Barcode")
  # print(dim(joined_df))
  return(joined_df)
  
}

sherlock_genome_association <- function(data, Var1, Var2, regression, formula, filter_zero1=FALSE, filter_zero2=FALSE, log_var1=FALSE, log_var2=FALSE, type="parametric", collapse_var1=NULL, collapse_var2=NULL, xlab=xlab, ylab=ylab, output_plot, file_ext = "png",plot_height,plot_width) {
  
  data <- validate_vardf(data)
  
  if(regression){
    
    data_types <- data %>% inspectdf::inspect_types()
    factor_vars <- as.vector(data_types$col_name[["factor"]])
    character_vars <- as.vector(data_types$col_name[["character"]])
    numeric_vars <- as.vector(data_types$col_name[["numeric"]])
    
    ## for regression module
    supported_types <- c("lm", "glm")
    
    # use different family depending on if the response variable (categorical family= binomial, continuous family=gaussian)
    # or tell user the response needs to be categorical since it is glm (logistic regression)
    
    # if(!str_detect(formula,"~")){
    #   stop("Please check your formula for regression, for example, lm( mpg ~ vs + gear")
    # }
    # 
    if(startsWith(formula, "lm")){
      # if(str_detect(formula, fixed("lm"))){
      print(paste0("str_detect(formula, fixed(lm)",str_detect(formula, fixed("lm"))))
      str_spl_formula <- str_split(formula, "\\(")
      print(str_spl_formula)
      formula <- paste0(str_spl_formula[[1]][1],"(formula=",str_spl_formula[[1]][2])
      print(formula)
    }
    
    if(startsWith(formula, "glm")){
      # if(str_detect(formula, fixed("glm"))){
      print(paste0("str_detect(formula, fixed(glm)",str_detect(formula, fixed("glm"))))
      str_spl_formula <- str_split(formula, "\\(")
      print(str_spl_formula)
      print(str_spl_formula[[1]][2])
      # get the first variable that is in the formula model (i.e. before the ~)
      response_var <- str_split(str_spl_formula[[1]][2],"~")
      print(response_var)
      # trim whitespace
      response_var2 <- trimws(response_var[[1]][1])
      print(response_var2)
      if(response_var2 %in% factor_vars | response_var2 %in% character_vars){
        formula <- paste0(str_spl_formula[[1]][1],"(family=binomial,formula=", str_spl_formula[[1]][2]) # glm(MCN_WGD ~ PGA) 
        print(formula)
      }
      if(response_var2 %in% numeric_vars){
        formula <- paste0(str_spl_formula[[1]][1],"(family=gaussian,formula=", str_spl_formula[[1]][2]) # glm(PGA ~ MCN_WGD)
        print(formula)
      }
      
    }
    
    # input_formula <- paste0("mod <- data %>% ",type, "(", formula,", data=.)")
    input_formula <- paste0("mod <- data %>% ", formula)
    print(input_formula)
    eval(parse(text=input_formula))
    
    p <- ggstatsplot::ggcoefstats(
      x= mod,
      point.args = list(color = "red", size = 3, shape = 15),
      exclude.intercept = TRUE,
      #title = formula,
      ggtheme = hrbrthemes::theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)
    ) + # note the order in which the labels are entered
      ggplot2::labs(x = "Regression Coefficient", y = NULL)
    
    filename <- paste0("multivariable_", formula, ".", file_ext)
    print(filename)
  
  }else{
    
    # print(Var1)
    # print(Var2)
    ## subset data
    data <- data %>% select(one_of(c(Var1,Var2)))
    colnames(data) <- c("Var1","Var2")
    var1_type <- if_else(is.factor(data[[1]]),"categorical", if_else(is.numeric(data[[1]]),"continuous",NA_character_))
    var2_type <- if_else(is.factor(data[[2]]),"categorical", if_else(is.numeric(data[[2]]),"continuous",NA_character_))
    
    # print(var1_type)
    # print(var2_type)
    # 
    # print(filter_zero1)
    # print(filter_zero2)
    # 
    # print(log_var1)
    # print(log_var2)
    # 
    # print(collapse_var1)
    # print(collapse_var2)
    # 
    if(is.na(var1_type)|is.na(var2_type)){
      stop("Please check your data type of these two selected variables")
    }
    
    # # process data or filtering data
    
    if(filter_zero1 & var1_type == 'continuous') {
      data <- data %>% filter(Var1 != 0)
    }else{
      print("none")
    }
    
    if(filter_zero2 & var2_type == 'continuous') {
      data <- data %>% filter(Var2 != 0)
    }else{
      print("none")
    }
    
    
    if(log_var1 & var1_type == 'continuous') {
      data <- data %>% filter(Var1>0) %>% mutate(Var1 = log2(Var1))
    }else{
      print("none")
    }
    
    if(log_var2 & var2_type == 'continuous') {
      data <- data %>% filter(Var2>0) %>% mutate(Var2 = log2(Var2))
    }else{
      print("none")
    }
    
    if(var1_type =="categorical" && !is.null(collapse_var1)){
      if(! (collapse_var1 %in% data$Var1)){
        print("Warning: categorical value does not exist in data, please input the correct level of the categorical variable for variable1.")
      }else{
        data$Var1 <- fct_other(data$Var1,keep = collapse_var1)
      }
    }
    
    if(var2_type =="categorical" && !is.null(collapse_var2)){
      if(! (collapse_var2 %in% data$Var2)){
        print("Warning: categorical value does not exist in data, please input the correct level of the categorical variable for variable2.")
      }else{
        data$Var2 <- fct_other(data$Var2,keep = collapse_var2)
      }
    }
    
    ## association test based on the types
    
    # for continues vs continues
    if(var1_type == 'continuous' & var2_type == 'continuous'){
      # supported types: "parametric", "nonparametric", "robust", "bayes", "skit
      supported_types <- c("parametric", "nonparametric", "robust", "bayes", "skit")
      if(!(type %in% supported_types)){
        print("Warning: selected statistic parameter type does not supported for selected data types, use the default parameteric type")
        type = "parametric"
      }
      
      if(type == "skit"){
        skit_res <- SKIT::skit(data$Var1,data$Var2,nboot = 1000)
        skit_res <- skit_res$pvalues[1]
        skit_lab <- paste0("P-value by SKIT test: ",if_else(skit_res <= 1/1000, "<1e-03",as.character(scientific(skit_res,digits = 3))))
        xlab=paste0(xlab,"\n",skit_lab)
        type = "parametric"
      }
      
      p <-  ggstatsplot::ggscatterstats(
        data = data,
        x = Var1,
        y = Var2, 
        xlab= xlab,
        ylab = ylab,
        marginal.type = "boxplot",
        xfill = "#009E73",
        yfill = "#D55E00",
        ggtheme = hrbrthemes::theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14),
        type=type
      )
      
    }
    
    # for categorical vs categorical
    if(var1_type == 'categorical' & var2_type == 'categorical'){
      # supported types: "parametric", "nonparametric", "robust", "bayes"
      supported_types <- c("parametric", "nonparametric", "robust", "bayes","fisher")
      if(!(type %in% supported_types)){
        print("Warning: selected statistic parameter type does not supported for selected data types, use the default parameteric type")
        type = "parametric"
      }
      
      if(type == "fisher"){
        fisher_res <- tidy(fisher.test(data$Var1,data$Var2))
        if("estimate" %in% colnames(fisher_res)){
          fisher_lab <- paste0("Fisher Exact Test: P-value =", scientific(fisher_res$p.value,digits = 3),", OR = ",number_format(accuracy = 0.01)(fisher_res$estimate))
        }else{
          fisher_lab <- paste0("Fisher Exact Test: P-value =", scientific(fisher_res$p.value,digits = 3))
        }
        xlab=paste0(xlab,"\n",fisher_lab)
        type = "parametric"
      }
      
      p <-  ggstatsplot::ggbarstats(
        data = data,
        x = Var1,
        y = Var2, 
        xlab= ylab,
        legend.title = xlab,
        ggtheme = hrbrthemes::theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14),
        type=type
      )
    }
    
    # for categorical vs continues 
    if(var1_type != var2_type ){
      # supported types: "parametric", "nonparametric", "robust", "bayes"
      # switch the name if Var2 is categorical
      if(var2_type == 'categorical'){
        p <-  ggstatsplot::ggbetweenstats(
          data = data,
          x = Var2,
          y = Var1, 
          xlab= ylab,
          ylab = xlab,
          ggtheme = hrbrthemes::theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14),
          type=type
        )
        
      }else {
        p <-  ggstatsplot::ggbetweenstats(
          data = data,
          x = Var1,
          y = Var2, 
          xlab= xlab,
          ylab = ylab,
          ggtheme = hrbrthemes::theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14),
          type=type
        )
        
      } 
      
    }
      filename <- paste0("bivariable_", Var1 ,"_", Var2, ".", file_ext)
      # print(filename)
    
    }

      if(output_plot == FALSE){
        return(p)
      }else{
        # filename <- paste0("output_plot.", file_ext)
        ggsave(filename = filename ,plot = p,width = plot_width,height = plot_height)
        return(p)
      }
}
    

sherlock_genome_association_group <- function(data, Var1, Var2, Group_Var, regression=FALSE, formula=NULL, filter_zero1=NULL, filter_zero2=NULL,log_var1=FALSE,log_var2=FALSE,type="parametric", collapse_var1=NULL, collapse_var2=NULL) {

  data <- validate_vardf(data,excludes = Group_Var)
  if(regression){
    supported_types <- c("lm", "glm")
    
    # if(is.null(formula)|!str_detect(formula,"~")){
    #   stop("Please check your formula for regression, for example, lm( mpg ~ vs + gear")
    # }
    
    colnames(data)[colnames(data) == Group_Var] <- 'Group'
    
    if(str_detect(formula, "lm")==TRUE){
      str_spl_formula <- str_split(formula, "\\(") 
      formula <- paste0("(", str_spl_formula[[1]][2])
      formula <- gsub("\\(", "", formula)
      formula <- gsub("\\)", "", formula)
      type <- "lm"
      
    }
    
    if(str_detect(formula, "glm")==TRUE){
      str_spl_formula <- str_split(formula, "\\(") 
      formula <- paste0("(", str_spl_formula[[1]][2]) # need to add in family somehow
      formula <- gsub("\\(", "", formula)
      formula <- gsub("\\)", "", formula)
      print(formula)
      type <- "glm"
      
    }
    
    if(type =="lm"){
      result <- data %>% group_by(Group) %>% do(tidy(lm(formula,data=.))) %>% ungroup() %>% filter(term!="(Intercept)") %>% filter(!is.na(p.value))
    }
    
    if(type == "glm"){
      result <- data %>% group_by(Group) %>% do(tidy(glm(formula,data=.))) %>% ungroup() %>% filter(term!="(Intercept)") %>% filter(!is.na(p.value))
    }
    
    colnames(result)[1] <- tolower(Group_Var)
    result <- result %>% ungroup() %>% group_by(term) %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH") %>% mutate(formula=formula)
    
    # transpose the result output to make it fit better and easier to read in shiny app
    result <- t(result)
    colnames(result) <- (1:ncol(result))
    
  }else{

    ## subset data
    data <- data %>% select(one_of(c(Group_Var,Var1,Var2)))
    colnames(data) <- c("Group","Var1","Var2")
    var1_type <- if_else(is.factor(data[["Var1"]]),"categorical", if_else(is.numeric(data[["Var1"]]),"continuous",NA_character_))
    var2_type <- if_else(is.factor(data[["Var2"]]),"categorical", if_else(is.numeric(data[["Var2"]]),"continuous",NA_character_))
    
    if(is.na(var1_type)|is.na(var2_type)){
      stop("Please check your data type of these two selected variables")
    }

    # process data or filtering data
    if(!is.null(filter_zero1) & var1_type == 'continuous') {
      filter_zero1 <-  as.numeric(filter_zero1)
      if(!is.na(filter_zero1)){
        data <- data %>% filter(Var1 > filter_zero1)
      }

    }

    if(!is.null(filter_zero2) & var2_type == 'continuous') {
      filter2 <-  as.numeric(filter_zero2)
      if(!is.na(filter_zero2)){
        data <- data %>% filter(Var2 > filter_zero2)
      }
    }

    if(log_var1 & var1_type == 'continuous') {
      data <- data %>% filter(Var1>0) %>% mutate(Var1 = log2(Var1))
    }

    if(log_var2 & var2_type == 'continuous') {
      data <- data %>% filter(Var2>0) %>% mutate(Var2 = log2(Var2))
    }

    if(var1_type =="categorical" && !is.null(collapse_var1)){
      if(! (collapse_var1 %in% data$Var1)){
        print("Warning: categorical value does not exist in data, please input the correct level of the categorical variables for variable1.")
      }else{
        data$Var1 <- fct_other(data$Var1,keep = collapse_var1)
      }
    }

    if(var2_type =="categorical" && !is.null(collapse_var2)){
      if(! (collapse_var2 %in% data$Var2)){
        print("Warning: categorical value does not exist in data, please input the correct level of the categorical variables for variable1.")
      }else{
        data$Var2 <- fct_other(data$Var2,keep = collapse_var2)
      }
    }

    ## association test based on the types

    # for continues vs continues
    if(var1_type == 'continuous' & var2_type == 'continuous'){
      # supported types: "parametric", "nonparametric", "robust", "bayes", "skit
      supported_types <- c("parametric", "nonparametric", "robust", "bayes", "skit")
      if(!(type %in% supported_types)){
        print("Warning: selected statistic parameter type does not supported for selected data types, use the default parameteric type")
        type = "parametric"
      }

      if(type == "skit"){

        result <- tibble(Group=character(),parameter1=character(),parameter2=character(),p.value=numeric(), method=character(), n.obs=integer())
        for(sig in unique(data$Group)){
          tmp <- data %>% filter(Group==sig)
          skit_res <- SKIT::skit(tmp$Var1,tmp$Var2,nboot = 1000)
          skit_res <- as.numeric(skit_res$pvalues[1])
          nobs <- tmp %>% filter(!is.na(Var1),!is.na(Var2)) %>% dim() %>% .[[1]]
          result <- tibble(Group=sig,parameter1="Var1",parameter2="Var2",p.value=skit_res, method="SKIT test", n.obs=nobs) %>%
            bind_rows(result)
        }

      }else{
        result <- data %>% group_by(Group) %>% group_modify(~statsExpressions::corr_test(data = .,x=Var1,y=Var2,type=type) %>% select(-expression)) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
      }

      result$parameter1 <- Var1
      result$parameter2 <- Var2
      colnames(result)[1:3] <- c(tolower(Group_Var),"variable_name1","varible_name2")

    }

    # for categorical vs categorical
    if(var1_type == 'categorical' & var2_type == 'categorical'){
      # supported types: "parametric", "nonparametric", "robust", "bayes"
      supported_types <- c("parametric", "nonparametric", "robust", "bayes","fisher")
      if(!(type %in% supported_types)){
        print("Warning: selected statistic parameter type does not supported for selected data types, use the default parameteric type")
        type = "parametric"
      }

      tmp <- data %>% count(Group,Var1,Var2) %>% count(Group) %>% filter(n>2) %>% pull(Group)

      if(type == "fisher"){
        result <-  data %>% filter(Group %in% tmp) %>%  nest_by(Group) %>% mutate(test=list(fisher.test(data$Var1,data$Var2))) %>% summarise(tidy(test)) %>% arrange(p.value) %>% ungroup() %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% ungroup()
      }else{
        result <- data %>% filter(Group %in% tmp) %>%  group_by(Group) %>% group_modify(~tryCatch(expr = statsExpressions::contingency_table(data = .,x=Var1,y=Var2,type=type), error = function(e) NULL)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
      }
      result <- result %>% mutate(variable_name1=Var1, variable_name2 = Var2) %>% select(Group,variable_name1,variable_name2,everything())
      colnames(result)[1] <- c(tolower(Group_Var))
    }

    # for categorical vs continues
    if(var1_type != var2_type ){
      # supported types: "parametric", "nonparametric", "robust", "bayes"
      # switch the name if Var2 is categorical
      ## remove unique value

      ## decide two sample test or oneway_annovar

      if(var1_type=="categorical"){
        if(length(levels(data$Var1))==2){

          tmp <- data %>% group_by(Group) %>% filter(!is.na(Var1),!is.na(Var2))%>%  summarise(n1=n_distinct(Var1),n2=n_distinct(Var2)) %>% filter(n1!=2|n2==1) %>% pull(Group)
          result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::two_sample_test(data = .,x=Var1,y=Var2,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
        }

        if(length(levels(data$Var1))>2){
          tmp <- data %>% filter(!is.na(Var1),!is.na(Var2)) %>% group_by(Group,Var1) %>% summarise(SD=sd(Var2)) %>% filter(SD==0) %>% select(Group) %>% unique() %>% pull(Group) %>% as.character()
          result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::oneway_anova(data = .,x=Var1,y=Var2,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
        }

        result$parameter1 <- Var1
        result$parameter2 <- Var2
        colnames(result)[1:3] <- c(tolower(Group_Var),"variable_name1","varible_name2")
      }else{

        #
        # result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::two_sample_test(data = .,x=Var2,y=Var1,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
        #
        if(length(levels(data$Var2))==2){

          tmp <- data %>% group_by(Group) %>% filter(!is.na(Var1),!is.na(Var2))%>%  summarise(n1=n_distinct(Var2),n2=n_distinct(Var1)) %>% filter(n1!=2|n2==1) %>% pull(Group)
          result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::two_sample_test(data = .,x=Var2,y=Var1,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
        }

        if(length(levels(data$Var2))>2){
          tmp <- data %>% filter(!is.na(Var1),!is.na(Var2)) %>% group_by(Group,Var2) %>% summarise(SD=sd(Var1)) %>% filter(SD==0) %>% select(Group) %>% unique() %>% pull(Group) %>% as.character()
          result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::oneway_anova(data = .,x=Var2,y=Var1,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
        }

        result$parameter1 <- Var2
        result$parameter2 <- Var1
        colnames(result)[1:3] <- c(tolower(Group_Var),"variable_name1","varible_name2")
      }
    }

  }

  return(result)
}


# Function ----------------------------------------------------------------


SurvZTW <- function(suvdata,plot=TRUE,legend.labs=NULL,keyname="Strata",filename=NULL,width = 7,height = 5.5,keyterm="KeyY"){
  fit <- coxph(Surv(Survival_Month, Death) ~ Key+Age+Gender+Stage+Smoking+Histology, data = suvdata) ### overall 
  suvpvalue <- broom::tidy(fit,conf.int=T) %>% filter(term==keyterm) %>% mutate(lab=paste0("P = ",format(p.value,digits = 2,scientific = FALSE),"\nHR = ",format(exp(estimate),digits = 2,scientific = FALSE)," (95% CI: ",format(exp(conf.low),digits = 2,scientific = FALSE),"-",format(exp(conf.high),digits = 2,scientific = FALSE),")")) %>% select(term,p.value,lab)
  if(plot){
    fit2 <- coxph(Surv(Survival_Month, Death) ~ strata(Key)+Age+Gender+Stage+Smoking+Histology, data = suvdata) ## strata curve ##
    ##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
    suvfit <- survfit(fit2)
    #plot(suvfit)
    suvfit$call$formula <- as.formula("Surv(Survival_Month, Death) ~ strata(Key)+Age+Gender+Stage+Smoking+Histology")
    
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

# enviroment function -----------------------------------------------------
Survgroup <- function(vartmp, sp_group, reference, keyname, filename){
  #vartmp <- "CR_Overall|Signature_CR"
  if(is.null(sp_group)){
    slists <- sherlock_overall %>% pull(Tumor_Barcode)
  }else{
    slists <- sherlock_overall %>% filter(SP_Group==sp_group) %>% pull(Tumor_Barcode)
  }
  
  # print(paste0("slists:", slists))
  suvdata_tmp <- 
    sherlock_data_full %>% 
    mutate(Key=paste0(Gene,"|",Type)) %>% 
    filter(Key==vartmp) %>% 
    filter(!is.na(Alteration)) %>% 
    select(Tumor_Barcode,Key=Alteration) %>% 
    left_join(suvdata) %>% 
    mutate(Key=factor(Key)) %>% 
    filter(Tumor_Barcode %in% slists)
  
  print(head(suvdata_tmp))
  
  if(!is.null(reference)){
    olevel <- levels(suvdata_tmp$Key)
    nlevel <- c(olevel[olevel == reference],olevel[olevel != reference])
    levels(suvdata_tmp$Key) <- nlevel
  }
  
  if(is.null(keyname)){
    keyname <-  vartmp
    print(keyname)
  }
  
  if(is.null(filename)){
    print(is.null(filename))
    filename <- paste0(str_replace_all(vartmp,'[/ \\|]','_'),'.pdf')
    print(filename)
  }else{
    print(filename)
  }
  
  SurvZTWm(suvdata = suvdata_tmp,plot = TRUE,keyname=str_remove(vartmp,'.*\\|'),pvalsize = 3,filename = filename)
}

SurvZTWm <- function(suvdata,plot=TRUE,legend.labs=NULL,keyname="Strata",filename=NULL,width = 7,height = 5.5,pvalsize=4,pvalab=NULL){
  
  fit <- coxph(Surv(Survival_Month, Death) ~ Key+Age+Gender+Stage+Smoking+Histology, data = suvdata) ### overall 
  suvpvalue <- tidy(fit,conf.int = T) %>% filter(str_detect(term,"Key")) %>% mutate(lab=paste0("P = ",format(p.value,digits = 2,scientific = FALSE),"\nHR = ",format(exp(estimate),digits = 2,scientific = FALSE)," (95% CI: ",format(exp(conf.low),digits = 2,scientific = FALSE),"-",format(exp(conf.high),digits = 2,scientific = FALSE),")")) %>% select(term,p.value,lab) %>% 
    mutate(term=str_remove(term,"Key")) %>% mutate(lab=paste0(term,":\n",lab))
  # print(suvpvalue)
  if(is.null(pvalab)){
    pvalab <- paste0(suvpvalue$lab,collapse = "\n")
  }
  
  # print("Fit complete")
  # print(fit)
  fit2 <- coxph(Surv(Survival_Month, Death) ~ strata(Key)+Age+Gender+Stage+Smoking+Histology, data = suvdata) ## strata curve ##
  ##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
  
  # print("Fit2 complete")
  # print(fit2)
  
  suvfit <- survfit(fit2)
  # print("Suvfit complete")
  # print(suvfit)
  #plot(suvfit)
  suvfit$call$formula <- as.formula("Surv(Survival_Month, Death) ~ strata(Key)+Age+Gender+Stage+Smoking+Histology")
  
  if(is.null(legend.labs)) 
  {
    legend.labs <- paste0(levels(suvdata$Key)," (", suvfit$n,")")
    print(legend.labs)
    print(length(legend.labs))
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
  # ggsurv$plot
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
  
  # ggsurv
  
  if(is.null(filename)){
    filename=paste0(keyname,"-Survival.pdf")
  }
  
  print(filename)

  # ggsave(filename,width = width,height = height,device = cairo_pdf)
  
  # print(suvpvalue)
  # return(suvpvalue)
  return(list(suvpvalue, ggsurv))
}

SurvZTWms <- function(suvdata,plot=TRUE,legend.labs=NULL,keyname="Strata",filename=NULL,width = 7,height = 5.5,pvalsize=4,pvalab=NULL){
  
  fit <- coxph(Surv(Survival_Month, Death) ~ Key+Age+Gender+Stage+Histology, data = suvdata) ### overall 
  suvpvalue <- tidy(fit,conf.int = T) %>% filter(str_detect(term,"Key")) %>% mutate(lab=paste0("P = ",format(p.value,digits = 2,scientific = FALSE),"\nHR = ",format(exp(estimate),digits = 2,scientific = FALSE)," (95% CI: ",format(exp(conf.low),digits = 2,scientific = FALSE),"-",format(exp(conf.high),digits = 2,scientific = FALSE),")")) %>% select(term,p.value,lab) %>% 
    mutate(term=str_remove(term,"Key")) %>% mutate(lab=paste0(term,":\n",lab))
  if(is.null(pvalab)){
    pvalab <- paste0(suvpvalue$lab,collapse = "\n")
  }
  
  fit2 <- coxph(Surv(Survival_Month, Death) ~ strata(Key)+Age+Gender+Stage+Histology, data = suvdata) ## strata curve ##
  ##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
  suvfit <- survfit(fit2)
  #plot(suvfit)
  suvfit$call$formula <- as.formula("Surv(Survival_Month, Death) ~ strata(Key)+Age+Gender+Stage+Histology")
  
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

# oncoplot_data_prep <- function(genomic_alt= list()){
oncoplot_data_prep <- function(genomic_alts = c(), opt_four_freq = c(), freq_table = NULL){ 
  # genomic_alt <- list("Gender|Overall_Feature", "RYR3|Mutation_Driver")
  # genomic_alt <- list("Gender|Overall_Feature","Smoking|Overall_Feature","WGD_Status|Overall_Feature","MMR_Status|Overall_Feature","HR_Status|Overall_Feature",
       # "HRDetect_Status|Overall_Feature","HLA_LOH|Overall_Feature","Kataegis|Overall_Feature","EBV|Overall_Feature","EUR|Overall_Feature","EAS|Overall_Feature",
       # "nonSmoking_Comparision|Overall_Feature","EAS_Comparision|Overall_Feature","Piano_Forte|Overall_Feature","Piano|Overall_Feature","Forte|Overall_Feature",
       # "Mezzo_forte|Overall_Feature","Smoking_SBS4_others|Overall_Feature","Passive_Smoking|Overall_Feature","Passive_Smoking_non|Overall_Feature")

  negative_values <- c("No", "N", "NA", "NULL", "Wild-type","WT", NA)
  # negative_values <- c("No", "N", "NA", "NULL", "Wild-type","WT", NA, "Non-Smoker", "MMR_deficient", "HR_deficient", "nWGD")
  # filtered_df <- data.frame()
  
  # option 2- make input into individual entries
  if(str_detect(genomic_alts, "\\n")){
    genomic_alts <- str_split(genomic_alts, "\\n")
    genomic_alts <- unlist(genomic_alts)
    if(str_detect(genomic_alts , ",")){
      
    }
       
    i <- 1
    for(each in genomic_alts){
      genomic_alts[i] <- trimws(genomic_alts[i])
      print(genomic_alts[i])
      i <- i + 1
    }
    print(genomic_alts)
  }

  i <- 1
  # option 1 and 2 (individual genomic alterations)
  filtered_df_ind_genalt <- list()
  filtered_df_source_cat <- list()
  for(each in genomic_alts){
    print(i)
    if(str_detect(genomic_alts[1], "\\|")){
      print(str_detect(genomic_alts[1], "\\|"))
      ind_gen_alt <- str_split(genomic_alts[i], "\\|")
      print(ind_gen_alt)
      filtered_gen_alt <- sherlock_data_full %>% filter(Gene==ind_gen_alt[[1]][1] & Type==ind_gen_alt[[1]][2])
      filtered_gen_alt <- filtered_gen_alt %>% filter(!(Alteration %in% negative_values))
      filtered_df_ind_genalt <- append(filtered_df_ind_genalt,list(filtered_gen_alt))
    
      # i <- i + 1
    }else{ # option 3 and 4 (source category selection)
      filtered_gen_alt <- sherlock_data_full %>% filter(Type == genomic_alts[i])
      
      if(!is.null(opt_four_freq)){
        print(opt_four_freq)
        feat_type <- filtered_gen_alt %>% pull(Gene) %>% unique() %>% paste0("|", genomic_alts[i]) #paste0("|Overall_Feature")
        cat_freq_filter <- freq_table %>% filter(name %in% feat_type) %>% select(name, Freq) %>% arrange(desc(Freq)) %>% filter(Freq > opt_four_freq) %>%
          separate(col = name, into=c('Gene','Type'), sep= "\\|")
        filtered_gen_alt <- filtered_gen_alt %>% filter(Type %in% cat_freq_filter$Type) %>% filter(Gene %in% cat_freq_filter$Gene)
      }
      filtered_gen_alt <- filtered_gen_alt %>% filter(!(Alteration %in% negative_values))
      filtered_df_source_cat <- append(filtered_df_source_cat,list(filtered_gen_alt))
      # i <- i + 1
    }
     i <- i + 1
     print(i)
  }
  
  if(length(filtered_df_ind_genalt) != 0){
    option <- "ind_genalt"
    ind_gen_alt_bind1 <- rep(paste0("filtered_df_ind_genalt[[",1:length(filtered_df_ind_genalt),"]]"))
    ind_gen_alt_bind2 <- paste0(ind_gen_alt_bind1, collapse= " , ")
    ind_gen_alt_bind_list <- paste0("filtered_df_ind_genalt <- rbind(", ind_gen_alt_bind2, ")")
    eval(parse(text=ind_gen_alt_bind_list))
    return(list(option, filtered_df_ind_genalt))
  }
 
  if(length(filtered_df_source_cat) !=0){
    option <- "source_cat"
    return(list(option, filtered_df_source_cat))
  }


}


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

# landscape_colors <- get_vcColors()
# 
# landscape_colors['In_Frame_Ins'] <- '#df65b0'
# landscape_colors['Frame_Shift_Del'] <- '#1F78B4FF'
# landscape_colors['Nonsense_Mutation'] <- '#a50f15'
# 
# 
# histcolor <- c('#a6611a','#f1b6da','#d01c8b')
# names(histcolor) <- unique(sherlock_histology$Histology)
# 
# landscape_colors <- c(landscape_colors,histcolor)

# after separate legend function is created, get rid of legend adjust
oncoplot <- function(data,data_clone = NULL, Gene_Sig=NULL,landscape_colors=NULL,gene_level=NULL,sample_level=NULL,sample_level0=sample_level0,GeneSortOnly=FALSE,namemax=NULL,tmar=0,bmar=0,nbreaks=c(0,20,40,60,80),p2_axis_hidden = TRUE,p2_hidden =TRUE,cell_height=0.95,legend_adjust=FALSE, frequency= 0.1, gen_alt_second= FALSE, order_by_input= FALSE){

  # data_alts <- data %>% pull(Alteration) %>% unique()
  # color_list <- read.delim("../landscape_colors.csv", header = TRUE, sep = ",")
  # i <- 1
  # for(each in data_alts){
  #   if(!(data_alts[i] %in% names(landscape_colors))){
  #     # assign color to each alteration without a color already set 
  #     landscape_colors[data_alts[i]] <- color_list$Value[i]
  #   }
  #   
  #   i <- i + 1
  # }

  # limited to sample level
  data <- data %>% filter(Tumor_Barcode %in% sample_level0)
  # sample_level0 <- sherlock_samples_unique$Tumor_Barcode
  

  ## add blank gene
  blankgene <- gene_level[!(gene_level %in% data$Gene)]
  # print(blankgene)
  if(length(blankgene)>0){
    n = length(blankgene)
    blankdata <- tibble(Subject=rep(data$Subject[1],n),Tumor_Barcode=rep(data$Tumor_Barcode[1],n),Gene=blankgene,Alteration=rep(NA_character_,n),Type=rep(data$Type[1],n))
    data <- bind_rows(data,blankdata)
  }

  samplesize <- length(sample_level0)
  altertype <- data$Type[1]
  altergene <- data$Gene[1]
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
  # Gene_Freq <- data0 %>% select(Gene,Tumor_Barcode) %>% unique() %>% count(Gene,sort=T) %>% mutate(n=if_else(Gene %in% blankgene,0L,n)) %>%  mutate(Freq=n/samplesize) %>% arrange(Freq) %>% mutate(Freq=percent(Freq,accuracy = 0.1)) #filter(Freq > .1) %>% mutate(Freq=percent(Freq,accuracy = 0.1))
  # Gene_Freq <- data0 %>% select(Gene,Tumor_Barcode) %>% unique() %>% count(Gene,sort=T) %>% mutate(n=if_else(Gene %in% blankgene,0L,n)) %>%  mutate(Freq=n/samplesize) %>% arrange(Freq) %>% filter(Freq >= frequency) %>% mutate(Freq=percent(Freq,accuracy = 0.1))
  Gene_Freq <- data0 %>% select(Gene,Tumor_Barcode) %>% unique() %>% count(Gene,sort=T) %>% mutate(n=if_else(Gene %in% blankgene,0L,n)) %>%  mutate(Freq=n/samplesize) %>% arrange(Freq) %>% filter(Freq >= frequency) %>% mutate(Freq=percent(Freq,accuracy = 0.1))
  
  #Gene_Freq <- data0 %>% select(Gene,Tumor_Barcode) %>% unique() %>% count(Gene,sort=T) %>% mutate(n=if_else(Gene %in% blankgene,0L,n)) %>%  mutate(Freq=n/samplesize) %>% arrange(Freq) %>% mutate(Freq=percent(Freq,accuracy = 0.1))
  Gene_Freq2 <- data0 %>% select(Gene,Tumor_Barcode) %>% unique() %>% count(Gene,sort=T) %>% mutate(n=if_else(Gene %in% blankgene,0L,n)) %>%  mutate(Freq=n/samplesize) %>% arrange(Freq)

  if(is.null(gene_level)){
    gene_level <- Gene_Freq %>% pull(Gene)
  }else{
    Gene_Freq <- Gene_Freq %>% mutate(Gene=factor(Gene,levels = gene_level)) %>% arrange(Gene)
  }

  Gene_Freq_alte <- data0 %>% count(Gene,Alteration) %>% mutate(Gene=factor(Gene,levels = gene_level)) %>% arrange(Gene) %>% filter(Gene %in% Gene_Freq$Gene)
  print(Gene_Freq)
  print(Gene_Freq_alte)
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

  print(gene_level)
  print(gene_level2)
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

  print(length(gene_level2))
  print(length(Gene_Freq$Freq))
  
  legend_size_cat <- theme(legend.box.background = element_blank(),legend.key.size = unit(0.5, "cm"), legend.key.height = unit(0.3, "cm"), legend.key.width =unit(0.4, "cm"), legend.position=c(0.5,1),legend.justification = c(0.5,1))
  # guides(fill=guide_legend(ncol=legendnrow,byrow=FALSE))+legend_size_cat

  # main plot
  p1 <- databg %>%
    ggplot(aes(factor(Tumor_Barcode,levels = sample_level),as.integer(Gene))) +
    geom_tile(height=cell_height,fill="gray95")+
    theme_ipsum_rc(base_size = 10,axis_text_size = 10,grid = '',axis = '')+
    # theme(plot.background = element_rect(fill = "#ECF0F5"))+ # panel.border = element_rect(fill = "#ECF0F5"))+
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.title.y.right = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(face = "italic"),panel.background = element_blank(),legend.position = 'bottom', legend.justification = 'bottom',legend.direction = 'vertical')+ #legend.justification = "top",
    scale_y_continuous(expand = c(0,0),breaks = 1:length(gene_level2),labels = gene_level2,sec.axis = dup_axis(breaks = 1:length(gene_level2),labels = Gene_Freq$Freq))+
    ## data
    geom_tile(data = data,aes(factor(Tumor_Barcode,levels = sample_level),as.integer(Gene)+shift,fill=Alteration,height=height),size=0)+
    scale_fill_manual(values = landscape_colors)+
    #scale_fill_manual(values= landscape_colors[which(names(landscape_colors) %in% data$Alteration)]) + 
    # scale_fill_manual(values = c("#d01c8b","#f1b6da"))+
    geom_vline(xintercept = 1:samplesize-0.5,color="white",size=0.2)+
    panel_border(size = 0.3,color = 'gray70')+
    guides(fill = guide_legend(ncol = 1,title.position = "top",title=altertype))+legend_size_cat
  
  if(gen_alt_second){
    onco_legend <- list()
    ind_gen_alts_plot_list <- list()
    #i <-1 
    i <- length(Gene_Freq$Gene)
    for(each in 1:length(Gene_Freq$Gene)){
      databg1 <- databg %>% filter(Gene== Gene_Freq$Gene[i]) 
      data1 <- data %>% filter(Gene== Gene_Freq$Gene[i])
      
    #for(each in 1:length(unique(databg$Gene))){
      #if(order_by_input){
      #   databg1 <- databg %>% filter(Gene==gene_level[i])
      #   data1 <- data %>% filter(Gene == gene_level[i])
      # }else{
      #  databg1 <- databg %>% filter(Gene==unique(databg$Gene)[i]) 
      #  data1 <- data %>% filter(Gene == unique(databg$Gene)[i])
      #}
      
      # make legend for individual genomic alterations
      ind_gen_alts_plot <- databg1 %>%
        ggplot(aes(factor(Tumor_Barcode,levels = sample_level),as.integer(Gene))) +
        geom_tile(height=cell_height,fill="gray95")+
        theme_ipsum_rc(base_size = 10,axis_text_size = 10,grid = '',axis = '')+
        # theme(plot.background = element_rect(fill = "#ECF0F5"))+ # panel.border = element_rect(fill = "#ECF0F5"))+
        theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.title.y.right = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(face = "italic"),panel.background = element_blank(),legend.position = 'bottom', legend.justification = 'bottom',legend.direction = 'vertical')+ #legend.justification = "top",
        scale_y_continuous(expand = c(0,0),breaks = 1:length(gene_level2),labels = gene_level2,sec.axis = dup_axis(breaks = 1:length(gene_level2),labels = Gene_Freq$Freq))+
        ## data
        geom_tile(data = data1 ,aes(factor(Tumor_Barcode,levels = sample_level),as.integer(Gene)+shift,fill=Alteration,height=height),size=0)+
        scale_fill_manual(values = landscape_colors[which(names(landscape_colors) %in% data1$Alteration)], name = Gene_Freq$Gene[i])+
        # scale_fill_manual(values = c("#d01c8b","#f1b6da"))+
        geom_vline(xintercept = 1:samplesize-0.5,color="white",size=0.2)+
        panel_border(size = 0.3,color = 'gray70')#+
      #guides(fill = guide_legend(ncol = 1,title.position = "top",title=altertype))+legend_size_cat
        pleg <- get_legend(ind_gen_alts_plot + theme(legend.key.size = unit(0.18, "cm")))
        ind_gen_alts_plot <- ind_gen_alts_plot +theme(legend.position='none')
      #pleg <- get_legend(pleg+theme(plot.margin=margin(l=0,r=-0.5,unit="cm"),legend.key.size = unit(0.18, "cm")))
        onco_legend <- list.append(onco_legend, pleg)
        ind_gen_alts_plot_list <- list.append(ind_gen_alts_plot_list, ind_gen_alts_plot)
        i <- i - 1
    } 
  }
  
  # if(order_by_input){
  #   reps <- rep(paste0("ind_gen_alts_plot_list[[",1:length(ind_gen_alts_plot_list),"]]"))
  #   reps <- paste0(reps, collapse= " , ")
  #   p1_init <- paste0("p1_by_input <- align_plots(",reps, ")")
  #   eval(parse(text=p1_init))
  # }
  # 
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
  #oncoplot_legend <- get_legend(p1)
  legend.box.margin = margin(t = 0, b = 0, l = 0, r = 0,unit = 'cm')
  p1 <- p1+theme(legend.position='none')

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
    # scale_fill_manual(values = c("#d01c8b","#f1b6da"))
    scale_fill_manual(values = landscape_colors)+
    theme(plot.background = element_rect(fill = "#ECF0F5"))+
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

  # if(gen_alt_second){
  #   reps <- rep(paste0("onco_legend[[",1:length(onco_legend),"]]"))
  #   #reps <- rep(paste0("onco_legend[[",1:length(onco_legend),"]]+theme(plot.margin=margin(r=0,l=0.2,t=tmar,b=bmar,unit='cm'))"))
  #   reps <- paste0(reps, collapse= " , ")
  #   leg_com_1 <- paste0("leg_com <- align_plots(",reps, ",align= 'h', axis= 'tb')")
  #   eval(parse(text=leg_com_1))
  # }
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
    # colnames(rankinfo)[2] <- altergene
  } else{

    rankinfo <- tibble(Tumor_Barcode=sample_level) %>%
      mutate(Order=as.integer(seq_along(Tumor_Barcode)))
    colnames(rankinfo)[2] <- altertype
    # colnames(rankinfo)[2] <- altergene
  }

  if(gen_alt_second){
    return(list(oncoplot=p_com,oncoplot_legend=onco_legend,sample_level=rankinfo,gene_level=gene_level,gene_freq=Gene_Freq2))
  }else{
    return(list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level,gene_freq=Gene_Freq2))
  }
  #return(list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level,gene_freq=Gene_Freq2))

}

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
  # oncoplot_com2 <- plot_grid(plotlist = plotlist2,
  #                            align = 'v',
  #                            axis = 'lr',
  #                            rel_widths = rep(1,num_plots),
  #                            rel_heights = sizelist,
  #                            ncol= 1)
  # 
  
  
  oncoplot_com <- plot_grid(oncoplot_com1,
                            # oncoplot_com2,
                            align = 'h',
                            axis = 'tb',
                            rel_widths = c(10,1),
                            rel_heights = c(1,1),
                            nrow = 1)
  oncoplot_legend <- plot_grid(plotlist = c(leglist,list(NULL)),
                            align = 'h',
                            rel_widths = c(rep(1,num_plots),1),
                            nrow = 1)
  # oncoplot_legend <- plegends
  oncoplot_final <- plot_grid(
    oncoplot_com,
    oncoplot_legend,
    rel_heights = c(1,0.2),
    align = 'v',axis = 'l',ncol = 1)
  return(oncoplot_final)
}

lolliplot_setup <- function(gene, group){
  # load and prepare data---------------------------------------------------------------
  #load('../RDS/sherlock_data_all.RData')
  #load('sherlock_maf.RData')
  
  #tmp <- read_csv('../RDS/oncoplot_colors.csv')
  #tmp <- read_csv('oncoplot_colors.csv')
  # tmp <- oncoplot_colors
  # landscape_colors <- tmp$Color
  # names(landscape_colors) <- tmp$Name
  # 
  gff0 = system.file('extdata', 'protein_domains.RDs', package = 'maftools')
  gff0 = readRDS(file = gff0)
  gff0 <- as_tibble(gff0)
  
  
  ## input paramters
  # gene <- 'TP53'
  # group <- NULL # N_A, NU,S_U
  # minN <- 5
  
  
  if(is.null(group)){
    samplelist <- sherlock_overall$Tumor_Barcode
  }else{
    samplelist <- sherlock_overall %>% filter(SP_Group %in% group)
  }
  
  tdata1 <- sherlock_maf %>%
    filter(Hugo_Symbol %in% gene, Tumor_Barcode %in% samplelist,Variant_Classification != "Splice_Site") %>% 
    mutate(AAChange.refGene = if_else(is.na(AAChange.refGene),'Splice_Site',AAChange.refGene)) %>% 
    select(Tumor_Barcode,Variant_Classification,AAChange.refGene) %>% 
    separate_rows(AAChange.refGene,sep = ',') %>% 
    separate(AAChange.refGene,into = c('Gene','TS','Exon','cDNA','AAChange'),sep = ':') %>% 
    mutate(AAPosition = parse_number(str_remove(AAChange,'p.'))) 
  
  
  tdata2 <- sherlock_maf %>%
    filter(Hugo_Symbol %in% gene, Tumor_Barcode %in% samplelist,Variant_Classification == "Splice_Site") %>% 
    select(Tumor_Barcode,Variant_Classification,GeneDetail.refGene) %>% 
    separate_rows(GeneDetail.refGene,sep = ';') %>% 
    separate(GeneDetail.refGene,into = c('TS','Exon','cDNA'),sep = ':') %>% 
    mutate(Gene = gene) %>% 
    mutate(AAPosition = ceiling(parse_number(str_remove(cDNA,'^c.'))/3)) %>% 
    mutate(AAChange = paste0('X',AAPosition,'_splice'))
  
  tdata0 <- bind_rows(tdata1,tdata2)
  
  tslist <- tdata0 %>% filter(!is.na(TS),TS %in% gff0$refseq.ID) %>% pull(TS) %>% unique()
  
  #tslist_select <- tslist[1]
  #tslist_select <- "NM_201283"
  
  return(list(tdata0, tslist, samplelist))
  
}

lolliplot_plot <- function(gene, tdata0, tslist, tslist_input, samplelist, minN){
  
  tslist_select <- tslist[which(tslist==tslist_input)] %>% unlist()
  print(tslist_select)
  
  tmp <- oncoplot_colors
  landscape_colors <- tmp$Color
  names(landscape_colors) <- tmp$Name
  
  gff0 = system.file('extdata', 'protein_domains.RDs', package = 'maftools')
  gff0 = readRDS(file = gff0)
  gff0 <- as_tibble(gff0)
  
  #tslist_select <- tslist[1]
  #tslist_select <- "NM_201283"
  
  print(tdata0)
  tdata <- tdata0 %>% filter(TS==tslist_select | is.na(TS)) %>% unique() 
  #tdata <- tdata0 %>% filter(TS=="NM_001126115" | is.na(TS)) %>% unique() 
  #print(head(tdata))
  
  ## label 
  
  genefreq <- percent_format(accuracy = 0.1)(length(unique(tdata$Tumor_Barcode))/length(samplelist))
  
  genetitle <- paste0(gene,' mutated in ',genefreq, ' tumors (N=',length(samplelist),')')
  
  
  
  # Check the domain --------------------------------------------------------
  gff <- gff0 %>% filter(HGNC==gene,refseq.ID==tslist_select) %>% mutate(AAPosition=(End+Start)/2) %>% mutate(n=-0.4)
  aasize <- gff$aa.length[1]
  dcolors <- c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f')
  gff$color <- dcolors[1:dim(gff)[1]]
  
  tdata_count <- tdata %>% count(AAPosition,AAChange,Variant_Classification)
  ymax <- max(tdata_count$n)*1.05
  
  # double check on this part that it is supposed to be minN
  #tdata_count2 <- tdata_count %>% filter(n>5)
  tdata_count2 <- tdata_count %>% filter(n>minN)
  
  lolli_plot <- tdata_count %>% 
    ggplot(aes(x=AAPosition,y=n))+
    geom_segment(aes(x=AAPosition,xend=AAPosition,y=-0.03*ymax,yend=n),size=0.2)+
    geom_point(aes(fill=Variant_Classification),pch=21,size=4)+
    ggrepel::geom_text_repel(
      data = tdata_count2,
      aes(label=AAChange,color=Variant_Classification),
      force_pull   = 0, # do not pull toward data points
      nudge_y      = 0.2,
      direction    = "x",
      angle        = 90,
      vjust        = 0.5,
      hjust = -1,
      segment.size = 0.1,
      max.iter = 1e4, max.time = 1
    ) +
    geom_rect(aes(xmin=0,xmax=aasize,ymin=-0.06*ymax,ymax=-0.03*ymax),fill="#cccccc",show.legend=FALSE)+
    geom_rect(data=gff,aes(xmin=Start,xmax=End,ymin=-0.07*ymax,ymax=-0.02*ymax),fill=gff$color,show.legend=FALSE)+
    geom_text(data=gff,aes(x=AAPosition,y=-0.045*ymax,label=Label),family = 'Roboto Condensed',fontface = 'plain',size=3.5)+
    scale_fill_manual(values =landscape_colors,breaks = tdata_count$Variant_Classification)+
    scale_color_manual(values =landscape_colors,breaks = tdata_count$Variant_Classification,guide='none')+
    scale_x_continuous(breaks = pretty_breaks(n = 8),limits = c(0,aasize),expand = expansion(mult = c(0.01,0.01)))+
    scale_y_continuous(breaks = pretty_breaks(),limits = c(-0.08*ymax,ymax),expand = expansion(mult = c(0,0)))+
    guides(color='none',fill=guide_legend(nrow=1))+
    theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid=FALSE,axis = T,axis_col = 'black',ticks = T)+
    theme(legend.position = 'bottom',legend.text=element_text(size=12))+
    labs(x='',y=paste0('# ',gene,' Mutations'),fill="Mutation Type",subtitle = genetitle)
  
  #ggsave(filename = 'tmp.pdf',width = 18,height = 6,device = cairo_pdf)
  
  #tdata
  return(list(lolli_plot, tdata))
  
}
  

# # oncoplot <- function(data,data_clone = NULL, Gene_Sig=NULL,landscape_colors=NULL,gene_level=NULL,sample_level=NULL,sample_level0=sample_level0,GeneSortOnly=FALSE,namemax=NULL,tmar=0,bmar=0,nbreaks=c(0,20,40,60,80),p2_axis_hidden=FALSE,p2_hidden=FALSE,cell_height=0.95,legend_adjust=FALSE)
# # # oncoplot test using filtered_gen_alt ("Gender|Overall_Feature")
# # # data = filtered_gen_alt
# # data = sherlock_data_full %>% filter(Type == "Mutation_Driver")
# # data = sherlock_data_full %>% filter(Type == "Mutation_Driver") %>% filter(Alteration == "Yes")
# # data = sherlock_data_full %>% filter(Type == "Mutation_Driver") %>% filter(Alteration == "Yes")
# # data = sherlock_data_full %>% filter(Gene== "Gender",Type == "Overall_Feature")
# # # data = test[[1]]
# # landscape_colors = landscape_colors
# # sample_level0 = sample_level0
# # gene_level = NULL
# # GeneSortOnly = TRUE
# # sample_level = sample_new_level
# # tmar = 0
# # bmar = 0
# # 
# # oncoplot(data,data_clone=NULL,landscape_colors,sample_level0,gene_level,GeneSortOnly,sample_level,tmar,bmar)

# i <- 1
# for(each in unique(test$Gene)){
#   # print(test$Gene[i])
#   alt_vals <- test %>% filter(Gene == test$Gene[i]) %>% pull(Alteration) %>% unique()
#   # print(alt_vals)
#   print(paste0(test$Gene[i], "|", "Overall_Feature"))
#   # print(paste0(test$Gene[i], ":", alt_vals))
#   i <- i + 1
# }

# [1] "Gender"
# [1] "Male"   "Female"
# [1] "Smoking"
# [1] "Smoker"     "Non-Smoker"
# [1] "WGD_Status"
# [1] "WGD"  "nWGD"
# [1] "MMR_Status"
# [1] "MMR_proficient" "MMR_deficient" 
# [1] "HR_Status"
# [1] "HR_proficient" "HR_deficient" 
# [1] "HRDetect_Status"
# [1] "No"  "Yes"
# [1] "HLA_LOH"
# [1] "N" "Y"
# [1] "Kataegis"
# [1] "Yes" "No" 
# [1] "EBV"
# [1] "No"  "Yes"
# [1] "EUR"
# [1] "Yes" "No" 
# [1] "EAS"
# [1] "No"  "Yes"
# [1] "nonSmoking_Comparision"
# [1] "No"  NA    "Yes"
# [1] "EAS_Comparision"
# [1] NA    "No"  "Yes"
# [1] "Piano_Forte"
# [1] "No"  "Yes" NA   
# [1] "Piano"
# [1] "No"  "Yes"
# [1] "Forte"
# [1] "Yes" "No" 
# [1] "Mezzo_forte"
# [1] "No"  "Yes"
# [1] "Smoking_SBS4_others"
# [1] "No"  "Yes"
# [1] "Passive_Smoking"
# [1] "Yes" NA    "No" 
# [1] "Passive_Smoking_non"
# [1] NA    "Yes" "No" 

