
parent <- getwd()

# onStop(function() {
#   setwd(parent)
#   print(getwd())
#   unlink("www/Genomic Data/User_project", recursive = TRUE)
# })

load('./www/Genomic Data/Sherlock_TCGA/Integrative_Analysis/sherlock_data_all.RData')
load('./www/Genomic Data/Sherlock_TCGA/Integrative_Analysis/sherlock_type_colors.RData')

#load('./www/Genomic Data/Sherlock/SCNA/bb_heatmap_inputs.RData')
#load('./www/Genomic Data/Sherlock/Integrative_Analysis/sherlock_variable.RData')
#load('./www/Genomic Data/Sherlock/Survival_Analysis/suvdata.RData')
#load('./test.RData')
#load('./www/Genomic Data/Sherlock/Integrative_Analysis/covdata0.RDS') # fisher enrichment glm only?
#load('./www/Genomic Data/Sherlock/Mutations/sherlock_maf.RData')
#color_list <- read.delim("/Users/kleinam/Sherlock-Genome/www/Genomic Data/landscape_colors.csv", header = TRUE, sep = ",")
#color_list <- read.delim("./www/Genomic Data/landscape_colors.csv", header = TRUE, sep = ",")
#load('~/Documents/Sherlock_Genome/www/Genomic Data/Sherlock/SCNA/bb_heatmap_inputs.RData') # scna


# mdata0 <- sherlock_data_full %>% 
#   mutate(Gene=paste0(Gene,"|",Type)) %>% 
#   select(Tumor_Barcode,Gene,Alteration) %>% 
#   pivot_wider(names_from = "Gene",values_from = "Alteration")



#oncoplot_colors <- read_csv(paste0(parent, '/oncoplot_colors.csv'))
#color_list <- read_delim("landscape_colors_all.csv", delim = ",") # same list as oncoplot, but without values assigned to a color (more general color list)

# pdfhra <- function(){
#   d <- read.csv(extrafont:::fonttable_file(), stringsAsFactors = FALSE)
# }
# 

# display debugging messages in R (if local) 
# and in the console log (if running in shiny)
# debug_msg <- function(...) {
#   is_local <- Sys.getenv('SHINY_PORT') == ""
#   in_shiny <- !is.null(shiny::getDefaultReactiveDomain())
#   txt <- toString(list(...))
#   if (is_local) message(txt)
#   if (in_shiny) shinyjs::runjs(sprintf("console.debug(\"%s\")", txt))
# }

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

# group_ui_list_build <- function(data){
#   # format params for selectizeGroupUI funcitonality for filtering and subsetting
#   data_cols <- colnames(data)
#   param_list <- list()
#   for(each in data_cols){
#     curr_col <- paste0(each, "=list(inputId ='", each, "',title ='", each, ":')")
#     #eval(parse(text=curr_col))
#     param_list <- append(param_list,curr_col)
#     param_list <- noquote(param_list)
#     #print(param_list)
#     param_list2 <- paste0(param_list[1:20], collapse = ",")
#     param_list2 <- paste0("list(", param_list2, ")")
#     #print(param_list2)
#   }
#   
#   print(noquote(param_list2))
#   return(param_list2)
# }

# os_detect function ----------------------------------------------------------------
# detect os to apply correct file paths
os_detect <- function(){
  os_name <- Sys.info()[['sysname']]
  return(os_name)
}

# while_test <- function(){
#   while (!is.null(dev.list()))  dev.off()
# }
# read_colnames function ----------------------------------------------------------------
# read in column names
# read_colnames <- function(filename, sep="\t",header=TRUE){
#   file <- read.delim(file=filename)
#   return(colnames(file))
# }

# # read_delim function ----------------------------------------------------------------
# read_delim <- function(filename, sep="\t",header=TRUE){
#   file <- read.delim(file=filename,sep="\t", header=TRUE)
#   return(file)
# }

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

# check filter input or lm/glm input (add later)
check_inputs <- function(data, inputs=NULL, formula=NULL){
  # example with one condition
  #MCN_WGD == 'nWGD'

  data <- validate_vardf(data)
  
  print(formula)
  print(inputs)
  if(!is.null(formula)){
    if(startsWith(formula, "lm (") | startsWith(formula, "lm(") | startsWith(formula, "glm(") | startsWith(formula, "glm (")){
      str_spl_formula <- str_split(formula, "\\(") 
      formula <- paste0("(", str_spl_formula[[1]][2])
      formula <- gsub("\\(", "", formula)
      formula <- gsub("\\)", "", formula)
      
      inputs_split <- str_split(formula, '[~+]')
      #string_list <- c()
      
      tmp <- str_trim(inputs_split[[1]])
      
      if(all(tmp %in% colnames(data))){
        var_info <- 'All variables are found in the data.'
      }else{
        print(paste0('One or more variables are not found in the data.'))
        var_info <- paste0('One or more variables are not found in the data.')
      }
      # for(each in inputs_split[[1]]){
      #   print(each)
      #   #tmp <- str_split(each, '[==>=<=]')[[1]][1]
      #   tmp <- each %>% str_trim()
      #   print(paste0('tmp:',tmp))
      #   #if(!(tolower(tmp) %in% tolower(colnames(data)))){
      #   if(!(tmp %in% colnames(data))){
      #     print(paste0('One or more variables are not found in the data.'))
      #     var_info <- paste0('One or more variables are not found in the data.')
      #   }else{
      #     var_info <- 'All variables are found in the data.'
      #   }
      # }
      
    }else{
      var_info <- 'The model type entered is not one of the available options (lm or glm).'
    }
  }else{
    if(all(grepl("[\\|&,]", inputs))){
      inputs_split <- str_split(inputs, '[\\|&,]')
      string_list <- c()
    }else{
      inputs_split <- inputs
      print(paste0('inputs_split',inputs_split))
    }
    
    for(each in inputs_split[[1]]){
      print(each)
      tmp <- str_split(each, '[==>=<=]')[[1]][1]
      tmp <- tmp %>% str_trim()
      print(paste0('tmp:',tmp))
      #if(!(tolower(tmp) %in% tolower(colnames(data)))){
      if(!(tmp %in% colnames(data))){
        print(paste0('One or more variables are not found in the data.'))
        var_info <- paste0('One or more variables are not found in the data.')
      }else{
        var_info <- 'All variables are found in the data.'
      }
    }
    #filter_condition <- data %>% filter( PGA >= 0.52) %>% filter(MCN_WGD == 'nWGD')
   
  }
  
  return(var_info)
}

check_gene <- function(data, gene){
  if(gene %in% data$Hugo_Symbol){
    gene_info <- 'Gene entered exists in the data.'
  }else{
    gene_info <- 'Warning: gene entered does not exist in data. Please change gene input and try again.'
  }
  
  return(gene_info)
}

# for survival analysis reference level input check
check_factor_levels <- function(data, variable, reference_level){
  
  print(paste0('statement:',reference_level %in% variable))
  #if (! (collapse_level %in% data[,variable])){
  if( ! (reference_level %in% variable)){
    level_info <- "Warning: Reference level does not exist in data, please input one of the levels output above."
  }else{
    level_info <- "Reference level exists in data."
  }
  
  #print(level_info)
  return(level_info)
}

check_levels <- function(data, variable, collapse_level){ # varlidate_vardf being run when evaluating variable 1 and variable 2
  
  print(paste0('statement:',collapse_level %in% data[,variable]))
  if (! (collapse_level %in% data[,variable])){
    level_info <- "Warning: categorical value does not exist in data, please input the correct level of the categorical variables for the variable."
  }else{
    level_info <- "Categorical value exists in data."
  }
  
  #print(level_info)

  # if(collapse_var1 != ''){
  #   if(!(collapse_var1 %in% data$variable_one)){
  #     level_info_var1 <- "Warning: categorical value does not exist in data, please input the correct level of the categorical variables for variable one."
  #   }else{
  #     level_info_var1 <- "Categorical value exists in data."
  #   }
  # }else{
  #   level_info_var1 <- NULL
  # }
  # 
  # if(collapse_var2 != ''){
  #   if(!(collapse_var2 %in% data$variable_two)){
  #     level_info_var2 <- "Warning: categorical value does not exist in data, please input the correct level of the categorical variables for variable two."
  #   }else{
  #     level_info_var2 <- "Categorical value exists in data."
  #   }
  # }else{
  #   level_info_var2 <- NULL
  # }
  # 
  return(level_info)
}


# sherlock_genome_filter function ----------------------------------------------------------------
# filter a dataframe given user input conditions
sherlock_genome_filter <- function(data, conditions = NULL, col_names){ #sherlock_genome_filter(data=data, conditions) conditions <- "MCN_WGD == 'nWGD'"
  #library(stringr)
  
  data <- validate_vardf(data)
  #column_selection <- col_names
  
  #conditions <- "MCN_WGD == 'nWGD'"
  #conditions <- c("MCN_WGD == 'nWGD'" , "PGA >= 0.52") MCN_WGD == 'nWGD',PGA >= 0.52 Wave =='W3', Gender== 'Female'
  # conditions <- c("MCN_WGD == 'nWGD'" "PGA >= 0.52",head(100)) # might want to somehow add this?
  if(!is.null(conditions)){
    #var_check <- check_inputs(data, inputs = conditions)
    #if(var_check == 'All variables are found in the data.'){
    #str_split(conditions, [==><%in%])
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
    #  }else{
    #   stop('One or more variables entered are not in the current data.')
    #   #message_to_console <- 'One or more variables entered are not in the current data.'
    # }
  }

  # print(column_selection)
  # if(str_detect(column_selection, ",")){
  #   column_select_split <- str_split(column_selection, ",'")
  #   print(column_select_split)
  # }
  # # if("All columns" %in% col_names){
  # #   data <- as.data.frame(data)
  # # }
  # # if(!("All columns" %in% col_names)){
  # #   data <- as.data.frame(data[,col_names])
  # # }
  # # 
  # # data <- validate_vardf(data)
  # 

  if(length(col_names > 1)){
    i <- 1
    # if("All columns" %in% col_names){
    #   column_list <- col_names[col_names != "All columns"]
    #   column_list <- paste0(column_list, collapse= ",")
    # }else{
      column_list <- paste0(col_names, collapse= ",")
    # }
    column_list <- (paste0("select(",column_list,")"))
  }else{
    column_list <- (paste0("select(",col_names,")"))
  }
  
  #filter_df <- paste0("filter_condition <- data %>% filter( ",conditions, ") %>%", column_list)
  if(!is.null(conditions)){
    filter_df <- paste0("filter_condition <- data %>% filter( ",conditions, ") %>%", column_list)
  }else{
    filter_df <- paste0("filter_condition <- data %>%", column_list)
  }
  
  #print(filter_df)
  #print(paste("098",eval(parse(text=filter_df))))
  eval(parse(text=filter_df))
  
  # test <- isTRUE(filter_condition)
  # print(test)
  # 
  #return(filter_condition)
  # if(!is.null(conditions)){
  #   if(var_check == 'One or more variables are not found in the data.'){
  #     print("HERE HERE")
  #     filter_df <- paste0("filter_condition <- data %>%", column_list)
  #     eval(parse(text=filter_df))
  #     message_to_console <- 'One or more variables entered are not in the current data.'
  #   }else{
  #     message_to_console <- 'All variables are found in the data.'
  #   }
  #   return(list(filter_condition, message_to_console))
  # }else{
  #   message_to_console <- NULL
  #   #print(message_to_console)
  #   return(list(filter_condition, message_to_console))
  # 
  # }
  
  return(filter_condition)
  
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

# inspect_data_function function ----------------------------------------------------------------
# use inspectdf package to explore data
inspect_data_function <- function(dataframe, type_of_inspection=NULL,column_name=NULL ){
  if(column_name == "All columns"){
    data <- as.data.frame(dataframe)
  }
  
  if(column_name != "All columns"){
    data <- as.data.frame(dataframe[,column_name])
    colnames(data) <- column_name
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
# figure_display_ngspurity <- function(tumor_barcode=NULL,battenberg=NULL,type=NULL,project_code=NULL){
#   ngspurity_qc <- read_delim("NGSpurity/ngspurity_qc_file.txt")
#   # pull the filename from the ngspurity_qc_file table
#   file_path <- ngspurity_qc$File[which(ngspurity_qc$Tumor_Barcode==tumor_barcode & ngspurity_qc$Battenberg==battenberg & ngspurity_qc$Type==type)]
# 
#   if(os_detect() %in% c("Linux","Darwin")){
#     filename <- sub(".","NGSpurity", file_path)
#     print(getwd())
#     print(filename)
#     list(src = filename, alt = paste0(tumor_barcode, "_", battenberg, "_", type))
#   }else{
#     filename <- sub(".",paste0("Genomic Data/", project_code,"/NGSpurity"), file_path)
#     print(getwd())
#     print(filename)
#     return(filename)
#   }
# }


#figure_display_mutationTime <- function(tumor_barcode=NULL){
# figure_display_mutationTime <- function(tumor_barcode=NULL, project_code=NULL){
#   #mutationTime_file <- read_delim("Clonal_Evolution/MutationTime/MutationTime_Proportion.txt")
#   # pull the filename from the ngspurity_qc_file table
#   #mut_barcode <- mutationTime_file$Tumor_Barcode[which(mutationTime_file$Tumor_Barcode==tumor_barcode)]
#   mut_barcode <- tumor_barcode
#   print(mut_barcode)
#   file_path <- paste0("Clonal_Evolution/MutationTime/", mut_barcode, "_MTime.pdf")
#   print(file_path)
#   
#   if(os_detect() %in% c("Linux","Darwin")){
#     filename <- file_path
#     print(getwd())
#     # print(filename)
#     list(src = filename, alt = paste0(tumor_barcode, "_MTime.pdf"))
#   }else{
#     filename <- sub(".",paste0("Genomic Data/", project_code,"/Clonal_Evolution/"), file_path)
#     print(getwd())
#     print(filename)
#     return(filename)
#   }
# }

# figure_display_clustered_mut <- function(tumor_barcode=NULL, project_code=NULL){
#   #mutationTime_file <- read_delim("Clonal_Evolution/MutationTime/MutationTime_Proportion.txt")
#   # pull the filename from the ngspurity_qc_file table
#   #mut_barcode <- mutationTime_file$Tumor_Barcode[which(mutationTime_file$Tumor_Barcode==tumor_barcode)]
#   mut_barcode <- tumor_barcode
#   print(mut_barcode)
#   file_path <- paste0("Mutational_Signatures/Clustered_Mutations/", mut_barcode, "_Mutation_Clustering.pdf")
#   print(file_path)
#   
#   if(os_detect() %in% c("Linux","Darwin")){
#     filename <- file_path
#     print(getwd())
#     # print(filename)
#     list(src = filename, alt = paste0(tumor_barcode, "_Mutation_Clustering.pdf"))
#   }else{
#     filename <- sub(".",paste0("Genomic Data/", project_code,"/"), file_path)
#     print(getwd())
#     print(filename)
#     return(filename)
#   }
# }

figure_display_tmb <- function(filename){
  file_path <- paste0("./", filename)
  print(file_path)
  
  if(os_detect() %in% c("Linux","Darwin")){
    filename <- file_path
    print(getwd())
    # print(filename)
    list(src = filename, alt = filename)
  }else{
    filename <- file_path
    print(getwd())
    print(filename)
    return(filename)
  }
}


# figure_display_scna_gistic <- function(amp_or_del){
#   #mutationTime_file <- read_delim("Clonal_Evolution/MutationTime/MutationTime_Proportion.txt")
#   # pull the filename from the ngspurity_qc_file table
#   #mut_barcode <- mutationTime_file$Tumor_Barcode[which(mutationTime_file$Tumor_Barcode==tumor_barcode)]
#   
#   print(amp_or_del)
#   file_path <- paste0("SCNA/", amp_or_del, "_qplot.pdf")
#   print(file_path)
#   
#   if(os_detect() %in% c("Linux","Darwin")){
#     filename <- file_path
#     print(getwd())
#     # print(filename)
#     list(src = filename, alt = paste0(amp_or_del, "_qplot.pdf"))
#   }else{
#     filename <- sub(".",paste0("Genomic Data/", project_code,"/"), file_path)
#     print(getwd())
#     print(filename)
#     return(filename)
#   }
# }

#fisher_result_test_run <- fishergroup(vartmp=vartmp, sp_group=NULL, samplelist=NULL, var2name="WGD_Status|Overall_Feature", excludes= NULL, excludes_cat= NULL, keeps= NULL, keeps_cat= NULL, minfreq= 0.03, freq_column= 'Freq', method= "fisher.test", glm_formula= "Var1 ~ Var2 + Gender", covdata= covdata0, subfold="TMP", anndata=NULL, fdrcutoff=0.1)
#fisher_result_test_run_glm <- fishergroup(vartmp=vartmp, sp_group=NULL, samplelist=NULL, var2name="WGD_Status|Overall_Feature", excludes= NULL, excludes_cat= NULL, keeps= NULL, keeps_cat= NULL, minfreq= 0.03, freq_column= 'Freq', method= "glm", glm_formula= "Var1 ~ Var2 + Gender", covdata= covdata0, subfold="TMP", anndata=NULL, fdrcutoff=0.1)

# fishergroup function ----------------------------------------------------------------
fishergroup <- function(mdata0,freq_data, group_data, vartmp, sp_group=NULL, samplelist=NULL, var2name = NULL, excludes = NULL,excludes_cat = NULL, keeps = NULL,keeps_cat = NULL, minfreq=0.03, freq_column='Freq',  method = "fisher.test", glm_formula = "Var1 ~ Var2 + Gender", covdata = covdata0, subfold="TMP",anndata=NULL, fdrcutoff=0.1){
  
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
  
  print('in fisher group function')
  
  vartmp <- c('Tumor_Barcode',vartmp)
  
  mdata <- mdata0
  
  if(!is.null(sp_group)){
    #tmpxx <- sherlock_overall %>% filter(SP_Group %in% sp_group) %>% pull(Tumor_Barcode)
    tmpxx <- group_data %>% filter(SP_Group %in% sp_group) %>% pull(Tumor_Barcode)
    
    mdata <- mdata0 %>% filter(Tumor_Barcode %in% tmpxx)
  }
  
  if(!is.null(samplelist)){
    mdata <- mdata0 %>% filter(Tumor_Barcode %in% samplelist)
  }
  
  mdata <- mdata %>% 
    pivot_longer(cols = -one_of(vartmp)) %>% 
    drop_na() 
  
  print(head(mdata))
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
  
  print('passed fisher.test if statement')
  
  if(method == "glm") {
    ## for logistic linear regression can be adjusted for assigned population, age, gender, smoking, histology
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
    #tmp_freq <- sherlock_freq %>% select(name,one_of(freq_column))
    tmp_freq <- freq_data %>% select(name,one_of(freq_column))
  }else{
    #tmp_freq <- sherlock_freq %>% select(name,contains(sp_group))
    tmp_freq <- freq_data %>% select(name,contains(sp_group))
    
  }
  colnames(tmp_freq)[2] <- 'Freq'
  
  result <- result %>% mutate(fdr = p.adjust(p.value,method = 'BH')) %>% left_join(tmp_freq)
  print(paste0('result$fdr',head(result$fdr)))

  ## fdr cutoff
  tmp <- result %>% filter(fdr<fdrcutoff) %>% dim() %>% .[[1]]
  if(tmp>0){
    fdrline <- mean(result$p.value[c(tmp,tmp+1)])
  }else{
    fdrline <- 0
  }
  
  print('passed fdr cutoff point')
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
  
  print('generating plots')
  
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
     #return(list(tmp,tmp_plot1,var_prop_plot))
   }else{
     return(list(result,tmp_plot1))
     #return(list(tmp,tmp_plot1))
   }
}

# fisherbarplot function ----------------------------------------------------------------
fisherbarplot <- function(mdata0, group_data, vartmp, sp_group=NULL, samplelist=NULL, var2name = NULL,subfold="TMP"){
  
  #pdfhr2()
  mdata <- mdata0
  
  if(!is.null(sp_group)){
    #tmpxx <- sherlock_overall %>% filter(SP_Group %in% sp_group) %>% pull(Tumor_Barcode)
    tmpxx <- group_data %>% filter(SP_Group %in% sp_group) %>% pull(Tumor_Barcode)
    
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
    #labs(fill=var2name,x='EAS_Comparison|Overall_Feature', y = 'Percentage')
  
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
  #joined_df <- df_list %>% reduce(left_join, by= "Tumor_Barcode")
  
  # look at column names in each dataframe in the input list (assumption that all have Tumor Barcode to start with)
  colnames_data_input <- c()
  print(paste0('df list length:', length(df_list)))
  for(each in df_list){
    colnames_data_input <- append(colnames_data_input, colnames(each), after=length(colnames(each)))
  }
  
  all_join_columns <- unique(colnames_data_input[which(duplicated(colnames_data_input))])
  print(all_join_columns)
  # need to check column names and that they are the same data type- or change to the same data type- ex. age from character to double or whatever the type is in 'x' dataframe
  # join by other column names that are the same depending on data input selected
  # ex. sample and survival data can also be joined by Subject, Gender, and Age in addition to Tumor Barcode
  #joined_df <- df_list %>% reduce(left_join, by = c("Tumor_Barcode", "Subject","Gender","Age"))
  joined_df <- df_list %>% reduce(left_join, by = all_join_columns)
  # print(dim(joined_df))
  return(list(joined_df, all_join_columns))
  
}

# sherlock_genome_association <- function(data, Var1, Var2, regression, formula, filter_zero1=FALSE, filter_zero2=FALSE, log_var1=FALSE, log_var2=FALSE, type="parametric", collapse_var1=NULL, collapse_var2=NULL, xlab=xlab, ylab=ylab, output_plot, file_ext = "png",plot_height,plot_width) {
#   
#   data <- validate_vardf(data)
#   
#   if(regression){
#     
#     data_types <- data %>% inspectdf::inspect_types()
#     factor_vars <- as.vector(data_types$col_name[["factor"]])
#     character_vars <- as.vector(data_types$col_name[["character"]])
#     numeric_vars <- as.vector(data_types$col_name[["numeric"]])
#     
#     ## for regression module
#     supported_types <- c("lm", "glm")
#     
#     # use different family depending on if the response variable (categorical family= binomial, continuous family=gaussian)
#     # or tell user the response needs to be categorical since it is glm (logistic regression)
#     
#     # if(!str_detect(formula,"~")){
#     #   stop("Please check your formula for regression, for example, lm( mpg ~ vs + gear")
#     # }
#     # 
#     if(startsWith(formula, "lm")){
#       # if(str_detect(formula, fixed("lm"))){
#       print(paste0("str_detect(formula, fixed(lm)",str_detect(formula, fixed("lm"))))
#       str_spl_formula <- str_split(formula, "\\(")
#       print(str_spl_formula)
#       formula <- paste0(str_spl_formula[[1]][1],"(formula=",str_spl_formula[[1]][2])
#       print(formula)
#     }
#     
#     if(startsWith(formula, "glm")){
#       # if(str_detect(formula, fixed("glm"))){
#       print(paste0("str_detect(formula, fixed(glm)",str_detect(formula, fixed("glm"))))
#       str_spl_formula <- str_split(formula, "\\(")
#       print(str_spl_formula)
#       print(str_spl_formula[[1]][2])
#       # get the first variable that is in the formula model (i.e. before the ~)
#       response_var <- str_split(str_spl_formula[[1]][2],"~")
#       print(response_var)
#       # trim whitespace
#       response_var2 <- trimws(response_var[[1]][1])
#       print(response_var2)
#       if(response_var2 %in% factor_vars | response_var2 %in% character_vars){
#         formula <- paste0(str_spl_formula[[1]][1],"(family=binomial,formula=", str_spl_formula[[1]][2]) # glm(MCN_WGD ~ PGA) 
#         print(formula)
#       }
#       if(response_var2 %in% numeric_vars){
#         formula <- paste0(str_spl_formula[[1]][1],"(family=gaussian,formula=", str_spl_formula[[1]][2]) # glm(PGA ~ MCN_WGD)
#         print(formula)
#       }
#       
#     }
#     
#     # input_formula <- paste0("mod <- data %>% ",type, "(", formula,", data=.)")
#     input_formula <- paste0("mod <- data %>% ", formula)
#     print(input_formula)
#     eval(parse(text=input_formula))
#     
#     p <- ggstatsplot::ggcoefstats(
#       x= mod,
#       point.args = list(color = "red", size = 3, shape = 15),
#       exclude.intercept = TRUE,
#       #title = formula,
#       ggtheme = hrbrthemes::theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)
#     ) + # note the order in which the labels are entered
#       ggplot2::labs(x = "Regression Coefficient", y = NULL)
#     
#     filename <- paste0("multivariable_", formula, ".", file_ext)
#     print(filename)
#   
#   }else{
#     
#     # print(Var1)
#     # print(Var2)
#     ## subset data
#     data <- data %>% select(one_of(c(Var1,Var2)))
#     print(head(data))
#     colnames(data) <- c("Var1","Var2")
#     var1_type <- if_else(is.factor(data[[1]]) | is.character(data[[1]]),"categorical", if_else(is.numeric(data[[1]]),"continuous",NA_character_))
#     var2_type <- if_else(is.factor(data[[2]]) | is.character(data[[2]]),"categorical", if_else(is.numeric(data[[2]]),"continuous",NA_character_))
#     
#     print(typeof(data[[1]]))
#     print(typeof(data[[2]]))
#     print(is.factor(data[[1]]))
#     print(is.factor(data[[2]]))
#     print(var1_type)
#     print(var2_type)
#     # 
#     # print(filter_zero1)
#     # print(filter_zero2)
#     # 
#     # print(log_var1)
#     # print(log_var2)
#     # 
#     # print(collapse_var1)
#     # print(collapse_var2)
#     # 
#     if(is.na(var1_type)|is.na(var2_type)){
#       stop("Please check your data type of these two selected variables")
#     }
#     
#     # # process data or filtering data
#     
#     if(filter_zero1 & var1_type == 'continuous') {
#       data <- data %>% filter(Var1 != 0)
#     }else{
#       print("none")
#     }
#     
#     if(filter_zero2 & var2_type == 'continuous') {
#       data <- data %>% filter(Var2 != 0)
#     }else{
#       print("none")
#     }
#     
#     
#     if(log_var1 & var1_type == 'continuous') {
#       data <- data %>% filter(Var1>0) %>% mutate(Var1 = log2(Var1))
#     }else{
#       print("none")
#     }
#     
#     if(log_var2 & var2_type == 'continuous') {
#       data <- data %>% filter(Var2>0) %>% mutate(Var2 = log2(Var2))
#     }else{
#       print("none")
#     }
#     
#     if(var1_type =="categorical" && !is.null(collapse_var1)){
#       if(! (collapse_var1 %in% data$Var1)){
#         print("Warning: categorical value does not exist in data, please input the correct level of the categorical variable for variable1.")
#       }else{
#         data$Var1 <- fct_other(data$Var1,keep = collapse_var1)
#       }
#     }
#     
#     if(var2_type =="categorical" && !is.null(collapse_var2)){
#       if(! (collapse_var2 %in% data$Var2)){
#         print("Warning: categorical value does not exist in data, please input the correct level of the categorical variable for variable2.")
#       }else{
#         data$Var2 <- fct_other(data$Var2,keep = collapse_var2)
#       }
#     }
#     
#     ## association test based on the types
#     
#     # for continues vs continues
#     if(var1_type == 'continuous' & var2_type == 'continuous'){
#       # supported types: "parametric", "nonparametric", "robust", "bayes", "skit
#       supported_types <- c("parametric", "nonparametric", "robust", "bayes", "skit")
#       if(!(type %in% supported_types)){
#         print("Warning: selected statistic parameter type does not supported for selected data types, use the default parameteric type")
#         type = "parametric"
#       }
#       
#       if(type == "skit"){
#         skit_res <- SKIT::skit(data$Var1,data$Var2,nboot = 1000)
#         skit_res <- skit_res$pvalues[1]
#         skit_lab <- paste0("P-value by SKIT test: ",if_else(skit_res <= 1/1000, "<1e-03",as.character(scientific(skit_res,digits = 3))))
#         xlab=paste0(xlab,"\n",skit_lab)
#         type = "parametric"
#       }
#       
#       p <-  ggstatsplot::ggscatterstats(
#         data = data,
#         x = Var1,
#         y = Var2, 
#         xlab= xlab,
#         ylab = ylab,
#         marginal.type = "boxplot",
#         xfill = "#009E73",
#         yfill = "#D55E00",
#         ggtheme = hrbrthemes::theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14),
#         type=type
#       )
#       
#     }
#     
#     # for categorical vs categorical
#     if(var1_type == 'categorical' & var2_type == 'categorical'){
#       # supported types: "parametric", "nonparametric", "robust", "bayes"
#       supported_types <- c("parametric", "nonparametric", "robust", "bayes","fisher")
#       if(!(type %in% supported_types)){
#         print("Warning: selected statistic parameter type does not supported for selected data types, use the default parameteric type")
#         type = "parametric"
#       }
#       
#       if(type == "fisher"){
#         fisher_res <- tidy(fisher.test(data$Var1,data$Var2))
#         if("estimate" %in% colnames(fisher_res)){
#           fisher_lab <- paste0("Fisher Exact Test: P-value =", scientific(fisher_res$p.value,digits = 3),", OR = ",number_format(accuracy = 0.01)(fisher_res$estimate))
#         }else{
#           fisher_lab <- paste0("Fisher Exact Test: P-value =", scientific(fisher_res$p.value,digits = 3))
#         }
#         xlab=paste0(xlab,"\n",fisher_lab)
#         type = "parametric"
#       }
#       
#       p <-  ggstatsplot::ggbarstats(
#         data = data,
#         x = Var1,
#         y = Var2, 
#         xlab= ylab,
#         legend.title = xlab,
#         ggtheme = hrbrthemes::theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14),
#         type=type
#       )
#     }
#     
#     # for categorical vs continues 
#     if(var1_type != var2_type ){
#       # supported types: "parametric", "nonparametric", "robust", "bayes"
#       # switch the name if Var2 is categorical
#       if(var2_type == 'categorical'){
#         p <-  ggstatsplot::ggbetweenstats(
#           data = data,
#           x = Var2,
#           y = Var1, 
#           xlab= ylab,
#           ylab = xlab,
#           ggtheme = hrbrthemes::theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14),
#           type=type
#         )
#         
#       }else {
#         p <-  ggstatsplot::ggbetweenstats(
#           data = data,
#           x = Var1,
#           y = Var2, 
#           xlab= xlab,
#           ylab = ylab,
#           ggtheme = hrbrthemes::theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14),
#           type=type
#         )
#         
#       } 
#       
#     }
#       filename <- paste0("bivariable_", Var1 ,"_", Var2, ".", file_ext)
#       # print(filename)
#     
#     }
# 
#       if(output_plot == FALSE){
#         return(p)
#       }else{
#         # filename <- paste0("output_plot.", file_ext)
#         ggsave(filename = filename ,plot = p,width = plot_width,height = plot_height)
#         return(p)
#       }
# }
#     
# 
# sherlock_genome_association_group <- function(data, Var1, Var2, Group_Var, regression=FALSE, formula=NULL, filter_zero1=NULL, filter_zero2=NULL,log_var1=FALSE,log_var2=FALSE,type="parametric", collapse_var1=NULL, collapse_var2=NULL) {
# 
#   data <- validate_vardf(data,excludes = Group_Var)
#   if(regression){
#     supported_types <- c("lm", "glm")
#     
#     # if(is.null(formula)|!str_detect(formula,"~")){
#     #   stop("Please check your formula for regression, for example, lm( mpg ~ vs + gear")
#     # }
#     
#     colnames(data)[colnames(data) == Group_Var] <- 'Group'
#     
#     if(str_detect(formula, "lm")==TRUE){
#       str_spl_formula <- str_split(formula, "\\(") 
#       formula <- paste0("(", str_spl_formula[[1]][2])
#       formula <- gsub("\\(", "", formula)
#       formula <- gsub("\\)", "", formula)
#       type <- "lm"
#       
#     }
#     
#     if(str_detect(formula, "glm")==TRUE){
#       str_spl_formula <- str_split(formula, "\\(") 
#       formula <- paste0("(", str_spl_formula[[1]][2]) # need to add in family somehow
#       formula <- gsub("\\(", "", formula)
#       formula <- gsub("\\)", "", formula)
#       print(formula)
#       type <- "glm"
#       
#     }
#     
#     if(type =="lm"){
#       result <- data %>% group_by(Group) %>% do(tidy(lm(formula,data=.))) %>% ungroup() %>% filter(term!="(Intercept)") %>% filter(!is.na(p.value))
#     }
#     
#     if(type == "glm"){
#       result <- data %>% group_by(Group) %>% do(tidy(glm(formula,data=.))) %>% ungroup() %>% filter(term!="(Intercept)") %>% filter(!is.na(p.value))
#     }
#     
#     colnames(result)[1] <- tolower(Group_Var)
#     result <- result %>% ungroup() %>% group_by(term) %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH") %>% mutate(formula=formula)
#     
#     # transpose the result output to make it fit better and easier to read in shiny app
#     result <- t(result)
#     colnames(result) <- (1:ncol(result))
#     
#   }else{
# 
#     ## subset data
#     data <- data %>% select(one_of(c(Group_Var,Var1,Var2)))
#     colnames(data) <- c("Group","Var1","Var2")
#     # var1_type <- if_else(is.factor(data[["Var1"]]),"categorical", if_else(is.numeric(data[["Var1"]]),"continuous",NA_character_))
#     # var2_type <- if_else(is.factor(data[["Var2"]]),"categorical", if_else(is.numeric(data[["Var2"]]),"continuous",NA_character_))
#     var1_type <- if_else(is.factor(data[["Var1"]]) | is.character(data[["Var1"]]),"categorical", if_else(is.numeric(data[["Var1"]]),"continuous",NA_character_))
#     var2_type <- if_else(is.factor(data[["Var2"]]) | is.character(data[["Var2"]]),"categorical", if_else(is.numeric(data[["Var2"]]),"continuous",NA_character_))
#     
#     if(is.na(var1_type)|is.na(var2_type)){
#       stop("Please check your data type of these two selected variables")
#     }
# 
#     # process data or filtering data
#     if(!is.null(filter_zero1) & var1_type == 'continuous') {
#       filter_zero1 <-  as.numeric(filter_zero1)
#       if(!is.na(filter_zero1)){
#         data <- data %>% filter(Var1 > filter_zero1)
#       }
# 
#     }
# 
#     if(!is.null(filter_zero2) & var2_type == 'continuous') {
#       filter2 <-  as.numeric(filter_zero2)
#       if(!is.na(filter_zero2)){
#         data <- data %>% filter(Var2 > filter_zero2)
#       }
#     }
# 
#     if(log_var1 & var1_type == 'continuous') {
#       data <- data %>% filter(Var1>0) %>% mutate(Var1 = log2(Var1))
#     }
# 
#     if(log_var2 & var2_type == 'continuous') {
#       data <- data %>% filter(Var2>0) %>% mutate(Var2 = log2(Var2))
#     }
# 
#     if(var1_type =="categorical" && !is.null(collapse_var1)){
#       if(! (collapse_var1 %in% data$Var1)){
#         print("Warning: categorical value does not exist in data, please input the correct level of the categorical variables for variable1.")
#       }else{
#         data$Var1 <- fct_other(data$Var1,keep = collapse_var1)
#       }
#     }
# 
#     if(var2_type =="categorical" && !is.null(collapse_var2)){
#       if(! (collapse_var2 %in% data$Var2)){
#         print("Warning: categorical value does not exist in data, please input the correct level of the categorical variables for variable1.")
#       }else{
#         data$Var2 <- fct_other(data$Var2,keep = collapse_var2)
#       }
#     }
# 
#     ## association test based on the types
# 
#     # for continues vs continues
#     if(var1_type == 'continuous' & var2_type == 'continuous'){
#       # supported types: "parametric", "nonparametric", "robust", "bayes", "skit
#       supported_types <- c("parametric", "nonparametric", "robust", "bayes", "skit")
#       if(!(type %in% supported_types)){
#         print("Warning: selected statistic parameter type does not supported for selected data types, use the default parameteric type")
#         type = "parametric"
#       }
# 
#       if(type == "skit"){
# 
#         result <- tibble(Group=character(),parameter1=character(),parameter2=character(),p.value=numeric(), method=character(), n.obs=integer())
#         for(sig in unique(data$Group)){
#           tmp <- data %>% filter(Group==sig)
#           skit_res <- SKIT::skit(tmp$Var1,tmp$Var2,nboot = 1000)
#           skit_res <- as.numeric(skit_res$pvalues[1])
#           nobs <- tmp %>% filter(!is.na(Var1),!is.na(Var2)) %>% dim() %>% .[[1]]
#           result <- tibble(Group=sig,parameter1="Var1",parameter2="Var2",p.value=skit_res, method="SKIT test", n.obs=nobs) %>%
#             bind_rows(result)
#         }
# 
#       }else{
#         result <- data %>% group_by(Group) %>% group_modify(~statsExpressions::corr_test(data = .,x=Var1,y=Var2,type=type) %>% select(-expression)) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
#       }
# 
#       result$parameter1 <- Var1
#       result$parameter2 <- Var2
#       colnames(result)[1:3] <- c(tolower(Group_Var),"variable_name1","variable_name2")
# 
#     }
# 
#     # for categorical vs categorical
#     if(var1_type == 'categorical' & var2_type == 'categorical'){
#       # supported types: "parametric", "nonparametric", "robust", "bayes"
#       supported_types <- c("parametric", "nonparametric", "robust", "bayes","fisher")
#       if(!(type %in% supported_types)){
#         print("Warning: selected statistic parameter type does not supported for selected data types, use the default parameteric type")
#         type = "parametric"
#       }
# 
#       tmp <- data %>% count(Group,Var1,Var2) %>% count(Group) %>% filter(n>2) %>% pull(Group)
# 
#       if(type == "fisher"){
#         result <-  data %>% filter(Group %in% tmp) %>%  nest_by(Group) %>% mutate(test=list(fisher.test(data$Var1,data$Var2))) %>% summarise(tidy(test)) %>% arrange(p.value) %>% ungroup() %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% ungroup()
#       }else{
#         result <- data %>% filter(Group %in% tmp) %>%  group_by(Group) %>% group_modify(~tryCatch(expr = statsExpressions::contingency_table(data = .,x=Var1,y=Var2,type=type), error = function(e) NULL)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
#       }
#       result <- result %>% mutate(variable_name1=Var1, variable_name2 = Var2) %>% select(Group,variable_name1,variable_name2,everything())
#       colnames(result)[1] <- c(tolower(Group_Var))
#     }
# 
#     # for categorical vs continues
#     if(var1_type != var2_type ){
#       # supported types: "parametric", "nonparametric", "robust", "bayes"
#       # switch the name if Var2 is categorical
#       ## remove unique value
# 
#       ## decide two sample test or oneway_annovar
# 
#       if(var1_type=="categorical"){
#         if(length(levels(data$Var1))==2){
# 
#           tmp <- data %>% group_by(Group) %>% filter(!is.na(Var1),!is.na(Var2))%>%  summarise(n1=n_distinct(Var1),n2=n_distinct(Var2)) %>% filter(n1!=2|n2==1) %>% pull(Group)
#           result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::two_sample_test(data = .,x=Var1,y=Var2,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
#         }
# 
#         if(length(levels(data$Var1))>2){
#           tmp <- data %>% filter(!is.na(Var1),!is.na(Var2)) %>% group_by(Group,Var1) %>% summarise(SD=sd(Var2)) %>% filter(SD==0) %>% select(Group) %>% unique() %>% pull(Group) %>% as.character()
#           result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::oneway_anova(data = .,x=Var1,y=Var2,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
#         }
# 
#         result$parameter1 <- Var1
#         result$parameter2 <- Var2
#         colnames(result)[1:3] <- c(tolower(Group_Var),"variable_name1","variable_name2")
#       }else{
# 
#         #
#         # result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::two_sample_test(data = .,x=Var2,y=Var1,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
#         #
#         if(length(levels(data$Var2))==2){
# 
#           tmp <- data %>% group_by(Group) %>% filter(!is.na(Var1),!is.na(Var2))%>%  summarise(n1=n_distinct(Var2),n2=n_distinct(Var1)) %>% filter(n1!=2|n2==1) %>% pull(Group)
#           result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::two_sample_test(data = .,x=Var2,y=Var1,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
#         }
# 
#         if(length(levels(data$Var2))>2){
#           tmp <- data %>% filter(!is.na(Var1),!is.na(Var2)) %>% group_by(Group,Var2) %>% summarise(SD=sd(Var1)) %>% filter(SD==0) %>% select(Group) %>% unique() %>% pull(Group) %>% as.character()
#           result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::oneway_anova(data = .,x=Var2,y=Var1,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
#         }
# 
#         result$parameter1 <- Var2
#         result$parameter2 <- Var1
#         colnames(result)[1:3] <- c(tolower(Group_Var),"variable_name1","variable_name2")
#       }
#     }
# 
#   }
# 
#   return(result)
# }

sherlock_genome_association_test <- function(data, Var1, Var2, Group_Var=NULL, regression=FALSE, formula=NULL, filter_zero1=FALSE, filter_zero2=FALSE,log_var1=FALSE,log_var2=FALSE,type="parametric", collapse_var1=NULL, collapse_var2=NULL,xlab=NULL,ylab=NULL) {
  
  # glm warning: Warning message: glm.fit: fitted probabilities numerically 0 or 1 occurred 3/22/23
  # Warning: Error in contrasts<-: contrasts can be applied only to factors with 2 or more levels 3/22/23
  
  if(is.null(Group_Var)){
    data <- validate_vardf(data)
  }else{
    data <- validate_vardf(data,excludes = Group_Var)
  }
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
      # if(startsWith(formula, "lm")){
      #   # if(str_detect(formula, fixed("lm"))){
      #   print(paste0("str_detect(formula, fixed(lm)",str_detect(formula, fixed("lm"))))
      #   str_spl_formula <- str_split(formula, "\\(")
      #   print(str_spl_formula)
      #   formula <- paste0(str_spl_formula[[1]][1],"(formula=",str_spl_formula[[1]][2])
      #   print(formula)
      # }
     
      if(!is.null(Group_Var)){ 
        colnames(data)[colnames(data) == Group_Var] <- 'Group'
      }
      
      if(startsWith(formula, "lm")==TRUE){
        str_spl_formula <- str_split(formula, "\\(") 
        formula <- paste0("(", str_spl_formula[[1]][2])
        formula <- gsub("\\(", "", formula)
        formula <- gsub("\\)", "", formula)
        type <- "lm"
        
      }
      
      if(startsWith(formula, "glm")==TRUE){
        print('entering glm section')
        str_spl_formula <- str_split(formula, "\\(") 
        formula <- paste0("(", str_spl_formula[[1]][2]) # need to add in family somehow
        formula <- gsub("\\(", "", formula)
        formula <- gsub("\\)", "", formula)
        print(formula)
        
        response_var <- str_split(str_spl_formula[[1]][2],"~")
        print(response_var)
        # trim whitespace
        response_var2 <- trimws(response_var[[1]][1])
        print(response_var2)
        if(response_var2 %in% factor_vars | response_var2 %in% character_vars){
          #formula <- paste0(str_spl_formula[[1]][1],"(family=binomial,formula=", str_spl_formula[[1]][2]) # glm(MCN_WGD ~ PGA)
          formula <- str_split(formula, "~")
          formula <- paste0("as.factor(",formula[[1]][1],") ~", formula[[1]][2])
          print(formula)
          family<- "binomial"
        }
        if(response_var2 %in% numeric_vars){
          #formula <- paste0(str_spl_formula[[1]][1],"(family=gaussian,formula=", str_spl_formula[[1]][2]) # glm(PGA ~ MCN_WGD)
          family<- "gaussian"
        }
        
        type <- "glm"
        print(type)
      }

      if(is.null(Group_Var)){
        # input_formula <- paste0("mod <- data %>% ",type, "(", formula,", data=.)")
        if(type == "lm"){
          input_formula <- paste0("mod <-", type, "(", formula, ",data = data)") 
        }
        if(type == "glm"){
          input_formula <- paste0("mod <-", type, "(family=",family, ", formula=", formula, ",data = data)") 
        }
        
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
      
        result <- p
        
      }else{
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
      }
    
    }else{ # if not regression (bivariable instead)
    
    ## subset data
      print(Group_Var)
      if(!is.null(Group_Var)){
        #print(head(data))
        data <- data %>% select(one_of(c(Group_Var,Var1,Var2)))
        colnames(data) <- c("Group","Var1","Var2")
        var1_type <- if_else(is.factor(data[["Var1"]]) | is.character(data[["Var1"]]),"categorical", if_else(is.numeric(data[["Var1"]]),"continuous",NA_character_))
        var2_type <- if_else(is.factor(data[["Var2"]]) | is.character(data[["Var2"]]),"categorical", if_else(is.numeric(data[["Var2"]]),"continuous",NA_character_))
        # print(paste0('var1_type:',var1_type))
        # print(paste0('var2_type:',var2_type))
        # print(head(data))
      }else{
        data <- data %>% select(one_of(c(Var1,Var2)))
        colnames(data) <- c("Var1","Var2")
        # var1_type <- if_else(is.factor(data[["Var1"]]),"categorical", if_else(is.numeric(data[["Var1"]]),"continuous",NA_character_))
        # var2_type <- if_else(is.factor(data[["Var2"]]),"categorical", if_else(is.numeric(data[["Var2"]]),"continuous",NA_character_))
        var1_type <- if_else(is.factor(data[[1]]) | is.character(data[[1]]),"categorical", if_else(is.numeric(data[[1]]),"continuous",NA_character_))
        var2_type <- if_else(is.factor(data[[2]]) | is.character(data[[2]]),"categorical", if_else(is.numeric(data[[2]]),"continuous",NA_character_))
        print(paste0('var1_type:',var1_type))
        print(paste0('var2_type:',var2_type))
        print(head(data))
      }
      
          if(is.na(var1_type)|is.na(var2_type)){
            stop("Please check your data type of these two selected variables")
          }

          # process data or filtering data
          if(!is.null(filter_zero1) & var1_type == 'continuous') {
            #filter_zero1 <-  as.numeric(filter_zero1)
            data <- data %>% filter(Var1 > 0)
          }
            # if(!is.na(filter_zero1)){
            #   data <- data %>% filter(Var1 > filter_zero1)
            # }

          if(!is.null(filter_zero2) & var2_type == 'continuous') {
            #filter2 <-  as.numeric(filter_zero2)
            data <- data %>% filter(Var2 > 0)
          }
            # if(!is.na(filter_zero2)){
            #   data <- data %>% filter(Var2 > filter_zero2)
            # }

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
              print("Warning: categorical value does not exist in data, please input the correct level of the categorical variables for variable2.")
              #collapse_var2_warning <- "Warning: categorical value does not exist in data, please input the correct level of the categorical variables for variable2."
            }else{
              data$Var2 <- fct_other(data$Var2,keep = collapse_var2)
            }
          }

          ## association test based on the types

          # for continues vs continues
        if("Group" %in% colnames(data)){
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
            colnames(result)[1:3] <- c(tolower(Group_Var),"variable_name1","variable_name2")
            
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
            print('var1_type is not the same as var2_type')
            print(var1_type)
            print(length(levels(data$Var1)))
            print(length(levels(as.factor(data$Var1))))
            print(length(unique(data$Var1)))
            # supported types: "parametric", "nonparametric", "robust", "bayes"
            # switch the name if Var2 is categorical
            ## remove unique value
            
            ## decide two sample test or oneway_annovar
            
            if(var1_type=="categorical"){
              #if(length(levels(data$Var1))==2){
              if(length(levels(as.factor(data$Var1)))==2){
                
                tmp <- data %>% group_by(Group) %>% filter(!is.na(Var1),!is.na(Var2))%>%  summarise(n1=n_distinct(Var1),n2=n_distinct(Var2)) %>% filter(n1!=2|n2==1) %>% pull(Group)
                result <- data  %>% filter(!Group %in% tmp) %>% group_by(Group) %>% group_modify(~statsExpressions::two_sample_test(data = .,x=Var1,y=Var2,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
              }
              
              #if(length(levels(data$Var1))>2){
              if(length(levels(as.factor(data$Var1)))>2){
                print(paste0('levels >2'))
                tmp <- data %>% filter(!is.na(Var1),!is.na(Var2)) %>% group_by(Group,Var1) %>% summarise(SD=sd(Var2)) %>% filter(SD==0) %>% select(Group) %>% unique() %>% pull(Group) %>% as.character()
                result <- data  %>% filter(!Group %in% tmp) %>% group_by(Group) %>% group_modify(~statsExpressions::oneway_anova(data = .,x=Var1,y=Var2,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
                #print(result)
              }
              
              result$parameter1 <- Var1
              result$parameter2 <- Var2
              colnames(result)[1:3] <- c(tolower(Group_Var),"variable_name1","variable_name2")
              print(head(result))
            }else{
              
              #
              # result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::two_sample_test(data = .,x=Var2,y=Var1,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
              #
              if(length(levels(as.factor(data$Var2)))==2){
                
                print("here")
                tmp <- data %>% group_by(Group) %>% filter(!is.na(Var1),!is.na(Var2))%>%  summarise(n1=n_distinct(Var2),n2=n_distinct(Var1)) %>% filter(n1!=2|n2==1) %>% pull(Group)
                result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::two_sample_test(data = .,x=Var2,y=Var1,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
              }
              
              if(length(levels(as.factor(data$Var2)))>2){
                print("here")
                tmp <- data %>% filter(!is.na(Var1),!is.na(Var2)) %>% group_by(Group,Var2) %>% summarise(SD=sd(Var1)) %>% filter(SD==0) %>% select(Group) %>% unique() %>% pull(Group) %>% as.character()
                result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::oneway_anova(data = .,x=Var2,y=Var1,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
              }
              
              result$parameter1 <- Var2
              result$parameter2 <- Var1
              colnames(result)[1:3] <- c(tolower(Group_Var),"variable_name1","variable_name2")
            }
          }
          
        }else{
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
                print('now I am here')
                # supported types: "parametric", "nonparametric", "robust", "bayes"
                supported_types <- c("parametric", "nonparametric", "robust", "bayes","fisher")
                if(!(type %in% supported_types)){
                  print("Warning: selected statistic parameter type does not supported for selected data types, use the default parameteric type")
                  type = "parametric"
                }

                if(type == "fisher"){
                  #fisher_res <- tidy(fisher.test(data$Var1,data$Var2))
                  #Warning: Error in fisher.test: FEXACT error 7(location). LDSTP=18300 is too small for this problem, (pastp=18.7841, ipn_0:=ipoin[itp=6]=4679, stp[ipn_0]=21.6643). Increase workspace or consider using 'simulate.p.value=TRUE'
                  fisher_res <- tidy(fisher.test(data$Var1,data$Var2, simulate.p.value = TRUE))
                  #fisher_res <- tidy(fisher.test(data$Var1,data$Var2, simulate.p.value = FALSE))
                  if("estimate" %in% colnames(fisher_res)){
                    fisher_lab <- paste0("Fisher Exact Test: P-value =", scientific(fisher_res$p.value,digits = 3),", OR = ",number_format(accuracy = 0.01)(fisher_res$estimate))
                  }else{
                    fisher_lab <- paste0("Fisher Exact Test: P-value =", scientific(fisher_res$p.value,digits = 3))
                  }
                  xlab=paste0(xlab,"\n",fisher_lab)
                  type = "parametric"
                }
                
                # generate palette when levels of variable are greater than default palette (Dark2, 8 colors)
                n <- length(unique(data$Var1))
                print(paste0('n:',n))
                palette <- distinctColorPalette(n)
                
                if(n > 8){
                  print('now I am trying to generate the plot for more than 8 colors')
                  p <-  ggstatsplot::ggbarstats(
                    data = data,
                    x = Var1,
                    y = Var2,
                    xlab= ylab,
                    legend.title = xlab,
                    ggtheme = hrbrthemes::theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14),
                    type=type,
                    ggplot.component = theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
                  ) + ggplot2::scale_fill_manual(values = palette)
                }else{
                  print('now I am trying to generate the plot for 8 colors or less')
                  p <-  ggstatsplot::ggbarstats(
                    data = data,
                    x = Var1,
                    y = Var2,
                    xlab= ylab,
                    legend.title = xlab,
                    ggtheme = hrbrthemes::theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14),
                    type=type,
                    ggplot.component = theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
                  )
                }
               
              
                p
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
          #       filename <- paste0("bivariable_", Var1 ,"_", Var2, ".", file_ext)
          #       # print(filename)
          #     
          #     }
          # 
          #       if(output_plot == FALSE){
          #         return(p)
          #       }else{
          #         # filename <- paste0("output_plot.", file_ext)
          #         ggsave(filename = filename ,plot = p,width = plot_width,height = plot_height)
          #         return(p)
          #       }
          # }
          result <- p
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

survprep <- function(data, suvdata, vartmp, sp_group, reference, slists){
  
  print('in prep function for survival analysis')
  
  suvdata_tmp <- 
    #sherlock_data_full %>% 
    data %>% 
    mutate(Key=paste0(Gene,"|",Type)) %>% 
    filter(Key==vartmp) %>% 
    filter(!is.na(Alteration)) %>% 
    select(Tumor_Barcode,Key=Alteration) %>% 
    left_join(suvdata) %>% 
    mutate(Key=factor(Key)) %>% 
    filter(Tumor_Barcode %in% slists)
  
  print(head(suvdata_tmp))
  
  # not originally in function. will remove
  # mdata0 <- sherlock_data_full %>% 
  #   mutate(Gene=paste0(Gene,"|",Type)) %>% 
  #   select(Tumor_Barcode,Gene,Alteration) %>% 
  #   pivot_wider(names_from = "Gene",values_from = "Alteration")
  
  # if(!is.null(reference)){

  olevel <- levels(suvdata_tmp$Key)
    # nlevel <- c(olevel[olevel == reference],olevel[olevel != reference])
    # levels(suvdata_tmp$Key) <- nlevel
  # }
  
  return(list(suvdata_tmp, olevel))
  #return(suvdata_tmp)
}

# environment function -----------------------------------------------------
#Survgroup <- function(data, suvdata, vartmp, sp_group, reference, keyname, filename, slists){
Survgroup <- function(data, vartmp, reference, keyname, filename){
  #vartmp <- "CR_Overall|Signature_CR"
  # if(is.null(sp_group)){
  #   slists <- sherlock_overall %>% pull(Tumor_Barcode)
  # }else{
  #   slists <- sherlock_overall %>% filter(SP_Group==sp_group) %>% pull(Tumor_Barcode)
  # }
  
  # print(paste0("slists:", slists))
  
  print('in first function for survival analysis')
  
  # suvdata_tmp <-
  #   #sherlock_data_full %>%
  #   data %>%
  #   mutate(Key=paste0(Gene,"|",Type)) %>%
  #   filter(Key==vartmp) %>%
  #   filter(!is.na(Alteration)) %>%
  #   select(Tumor_Barcode,Key=Alteration) %>%
  #   left_join(suvdata) %>%
  #   mutate(Key=factor(Key)) #%>%
  #   filter(Tumor_Barcode %in% slists)

  # print(head(suvdata_tmp))
  # 
  # # not originally in function. will remove
  # # mdata0 <- sherlock_data_full %>% 
  # #   mutate(Gene=paste0(Gene,"|",Type)) %>% 
  # #   select(Tumor_Barcode,Gene,Alteration) %>% 
  # #   pivot_wider(names_from = "Gene",values_from = "Alteration")
  # 
  # if(!is.null(reference)){
  #   olevel <- levels(suvdata_tmp$Key)
  #   nlevel <- c(olevel[olevel == reference],olevel[olevel != reference])
  #   levels(suvdata_tmp$Key) <- nlevel
  # }
  
  suvdata_tmp <- data
  
  olevel <- levels(suvdata_tmp$Key)
  nlevel <- c(olevel[olevel == reference],olevel[olevel != reference])
  #levels(suvdata_tmp$Key) <- nlevel
  
  suvdata_tmp$Key <- factor(suvdata_tmp$Key, levels = nlevel)
  
  print(olevel)
  print(nlevel)
  print(levels(suvdata_tmp$Key))
  
  suvdata <- suvdata_tmp
  
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
  
  SurvZTWm(suvdata = suvdata,plot = TRUE,keyname=str_remove(vartmp,'.*\\|'),pvalsize = 3,filename = filename)
}

SurvZTWm <- function(suvdata,plot=TRUE,legend.labs=NULL,keyname="Strata",filename=NULL,width = 7,height = 5.5,pvalsize=4,pvalab=NULL){
  
  print('in SurvZTWm function now')
  fit <- coxph(Surv(Survival_Month, Death) ~ Key+Age+Gender+Stage+Smoking+Histology, data = suvdata) ### overall 
  suvpvalue <- tidy(fit,conf.int = T) %>% filter(str_detect(term,"Key")) %>% mutate(lab=paste0("P = ",format(p.value,digits = 2,scientific = FALSE),"\nHR = ",format(exp(estimate),digits = 2,scientific = FALSE)," (95% CI: ",format(exp(conf.low),digits = 2,scientific = FALSE),"-",format(exp(conf.high),digits = 2,scientific = FALSE),")")) %>% select(term,p.value,lab) %>% 
    mutate(term=str_remove(term,"Key")) %>% mutate(lab=paste0(term,":\n",lab))
  # print(suvpvalue)
  if(is.null(pvalab)){
    pvalab <- paste0(suvpvalue$lab,collapse = "\n")
  }
  
  print("Fit complete")
  # print(fit)
  fit2 <- coxph(Surv(Survival_Month, Death) ~ strata(Key)+Age+Gender+Stage+Smoking+Histology, data = suvdata) ## strata curve ##
  ##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata)  ## all strata, wrong way
  
  print("Fit2 complete")
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
oncoplot_data_prep <- function(genomic_alts = c(), opt_four_freq = c(), freq_table = NULL, input_opt){ 
  # genomic_alt <- list("Gender|Overall_Feature", "RYR3|Mutation_Driver")
  # genomic_alt <- list("Gender|Overall_Feature","Smoking|Overall_Feature","WGD_Status|Overall_Feature","MMR_Status|Overall_Feature","HR_Status|Overall_Feature",
       # "HRDetect_Status|Overall_Feature","HLA_LOH|Overall_Feature","Kataegis|Overall_Feature","EBV|Overall_Feature","EUR|Overall_Feature","EAS|Overall_Feature",
       # "nonSmoking_Comparision|Overall_Feature","EAS_Comparision|Overall_Feature","Piano_Forte|Overall_Feature","Piano|Overall_Feature","Forte|Overall_Feature",
       # "Mezzo_forte|Overall_Feature","Smoking_SBS4_others|Overall_Feature","Passive_Smoking|Overall_Feature","Passive_Smoking_non|Overall_Feature")

  print(paste0('genomic alts being used in function call:', genomic_alts))
  negative_values <- c("No", "N", "NA", "NULL", "Wild-type","WT", NA)
  # negative_values <- c("No", "N", "NA", "NULL", "Wild-type","WT", NA, "Non-Smoker", "MMR_deficient", "HR_deficient", "nWGD")
  # filtered_df <- data.frame()
  
  # option 2- make input into individual entries
  if(input_opt == 2){
    if(str_detect(genomic_alts, "[\\n,]")){
      if(str_detect(genomic_alts, ",")){
        genomic_alts <- str_split(genomic_alts, ",")
      }
      if(str_detect(genomic_alts, "\\n")){
        genomic_alts <- str_split(genomic_alts, "\\n")
      }

      genomic_alts <- unlist(genomic_alts)

      i <- 1
      for(each in genomic_alts){
        genomic_alts[i] <- trimws(genomic_alts[i])
        print(genomic_alts[i])
        i <- i + 1
      }

      print(genomic_alts)

    }
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
        print(paste0('opt_four_freq', opt_four_freq))
        feat_type <- filtered_gen_alt %>% pull(Gene) %>% unique() %>% paste0("|", genomic_alts[i]) #paste0("|Overall_Feature")
        print(feat_type)
        cat_freq_filter <- freq_table %>% filter(name %in% feat_type) %>% select(name, Freq) %>% arrange(desc(Freq)) %>% filter(Freq > opt_four_freq) %>%
          separate(col = name, into=c('Gene','Type'), sep= "\\|")
        print(cat_freq_filter)
        filtered_gen_alt <- filtered_gen_alt %>% filter(Type %in% cat_freq_filter$Type) %>% filter(Gene %in% cat_freq_filter$Gene)
        print(filtered_gen_alt)
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
# get_vcColors = function(alpha = 1){
#   col = c(RColorBrewer::brewer.pal(12,name = "Paired"), RColorBrewer::brewer.pal(11,name = "Spectral")[1:3],'black', 'violet', 'royalblue', '#7b7060','#8dd3c7','#01665e','#2b8cbe','#c51b8a','#1c9099','#542788','#006d2c','#e41a1c','#253494', "#c994c7", "#dd1c77","#f03b20","#1F77B4FF", "#FF7F0EFF", "#2CA02CFF",'#993404','#ff7f00','#54278f','#253494','#35978f','#374E55FF')
#   col = grDevices::adjustcolor(col = col, alpha.f = alpha)
#   names(col) = names = c('Nonstop_Mutation','Frame_Shift_Del','IGR','Missense_Mutation','Silent','Nonsense_Mutation',
#                          'RNA','Splice_Site','Intron','Frame_Shift_Ins','Nonstop_Mutation','In_Frame_Del','ITD','In_Frame_Ins',
#                          'Translation_Start_Site',"Multi_Hit", 'Amp', 'Del', 'Complex_Event','Promoter','Fusion','INS','DEL','SNV','Chr19_Loss','HLA_LOH',"Amplification","Deletion",'MSI-L','MSI-H','WGD','C1','C2','C3','Kataegis','SV','synonymous_variant','p53_deficiency','RTK-RAS+','LOH')
#   #col <- c(col,c('DEL'='red','INS'='blue','SNP'='gray30'))
#   col
# }

get_vcColors <- function(){
  
  oncoplot_colors <- read_delim(paste0(parent, '/oncoplot_colors.csv'), delim = ',')
  #oncoplot_colors <- read_delim("~/Documents/Sherlock_Genome/oncoplot_colors.csv", delim = ",")
  onco_colors <- oncoplot_colors$Color
  names(onco_colors) = oncoplot_colors$Name
  onco_colors
}

landscape_colors_fcn <- function(){
  landscape_colors <- get_vcColors()
  landscape_colors['In_Frame_Ins'] <- '#df65b0'
  landscape_colors['Frame_Shift_Del'] <- '#1F78B4FF'
  landscape_colors['Nonsense_Mutation'] <- '#a50f15'
  landscape_colors['C>A'] <- '#2EBAED'
  landscape_colors['C>G'] <- '#000000'
  landscape_colors['C>T'] <- '#DE1C14'
  landscape_colors['T>A'] <- '#D4D2D2'
  landscape_colors['T>C'] <- '#ADCC54'
  landscape_colors['T>G'] <- '#F0D0CE'
  landscape_colors['SNP'] <- '#F0D0CE'
  landscape_colors['Multiple'] <- 'black'
  
  landscape_colors
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
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.title.y.right = element_blank(),axis.text.x = element_blank(),panel.background = element_blank(),legend.position = 'bottom', legend.justification = 'bottom',legend.direction = 'vertical')+ #legend.justification = "top",
    #theme_ipsum_rc(base_size = 10,axis_text_size = 10,grid = '',axis = '')+
    scale_y_continuous(expand = c(0,0),breaks = 1:length(gene_level2),labels = gene_level2,sec.axis = dup_axis(breaks = 1:length(gene_level2),labels = Gene_Freq$Freq))+
    ## data
    geom_tile(data = data,aes(factor(Tumor_Barcode,levels = sample_level),as.integer(Gene)+shift,fill=Alteration,height=height),size=0)+
    scale_fill_manual(values = landscape_colors)+
    #scale_fill_manual(values= landscape_colors[which(names(landscape_colors) %in% data$Alteration)]) + 
    # scale_fill_manual(values = c("#d01c8b","#f1b6da"))+
    geom_vline(xintercept = 1:samplesize-0.5,color="white",size=0.2)+
    panel_border(size = 0.3,color = 'gray70')+
    guides(fill = guide_legend(ncol = 1,title.position = "top",title=altertype))+legend_size_cat+ theme(legend.position = 'top')
  
  # p1 <- p1 + theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.title.y.right = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(face = "italic"),panel.background = element_blank(),legend.position = 'bottom', legend.justification = 'bottom',legend.direction = 'vertical') #legend.justification = "top",
  # p1 <- p1 + theme_ipsum_rc(base_size = 10,axis_text_size = 10,grid = '',axis = '')
  # 
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
        panel_border(size = 0.3,color = 'gray70') + guides(fill = guide_legend(ncol = 1,title.position = "top",title=Gene_Freq$Gene[i]))+legend_size_cat+ theme(legend.position = 'top')
      #pleg <- get_legend(ind_gen_alts_plot + theme(legend.key.size = unit(0.18, "cm")))
      pleg <- get_plot_component(ind_gen_alts_plot, 'guide-box-top', return_all = TRUE) #+ theme(legend.key.size = unit(0.18, "cm"))
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


  #oncoplot_legend <- get_legend(p1+theme(plot.margin=margin(b=-12,l=0,r=0,unit="cm"),legend.key.size = unit(0.18, "cm")))
  oncoplot_legend <- get_plot_component(p1, 'guide-box-top')
  #legend.box.margin = margin(t = 0, b = 0, l = 0, r = 0,unit = 'cm')
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
    #rel_heights = c(1,0.2),
    rel_heights = c(3,0.2),
    align = 'v',axis = 'l',ncol = 1)
  return(oncoplot_final)
}

lolliplot_setup <- function(gene, group, data_overall, data_maf){
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
  gff0 = readRDS(file = gff0) %>% select(-12)
  gff0 <- as_tibble(gff0)
  print(head(gff0))
  
  
  ## input paramters
  # gene <- 'TP53'
  # group <- NULL # N_A, NU,S_U
  # minN <- 5
  
  if(is.null(group)){
    #samplelist <- sherlock_overall$Tumor_Barcode
    samplelist <- data_overall$Tumor_Barcode
  }else{
    #samplelist <- sherlock_overall %>% filter(SP_Group %in% group) %>% pull(Tumor_Barcode)
    samplelist <- data_overall %>% filter(SP_Group %in% group) %>% pull(Tumor_Barcode)
  }
  
  #tdata1 <- sherlock_maf %>%
  tdata1 <- data_maf %>%
    filter(Hugo_Symbol %in% gene, Tumor_Barcode %in% samplelist,Variant_Classification != "Splice_Site") %>% 
    mutate(AAChange.refGene = if_else(is.na(AAChange.refGene),'Splice_Site',AAChange.refGene)) %>% 
    select(Tumor_Barcode,Variant_Classification,AAChange.refGene) %>% 
    separate_rows(AAChange.refGene,sep = ',') %>% 
    separate(AAChange.refGene,into = c('Gene','TS','Exon','cDNA','AAChange'),sep = ':') %>% 
    mutate(AAPosition = parse_number(str_remove(AAChange,'p.'))) 
  
  
  #tdata2 <- sherlock_maf %>%
  tdata2 <- data_maf %>%
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

lolliplot_plot <- function(gene, tdata0, tslist, tslist_input, samplelist, minN, domain_annotation){
  
  #grouplist <- unique(sherlock_overall$SP_Group)
  
  ## input paramters
  # gene <- 'TP53'
  # group <- NULL # N_A, NU,S_U
  # minN <- 5
  # domain_annotation <- FALSE
  
  tslist_select <- tslist[which(tslist==tslist_input)] %>% unlist()
  print(tslist_select)
  
  # tmp <- oncoplot_colors
  # landscape_colors <- tmp$Color
  # names(landscape_colors) <- tmp$Name
  
  landscape_colors <- get_vcColors()
  # landscape_colors <- tmp$Color
  # names(landscape_colors) <- tmp$Name
  
  gff0 = system.file('extdata', 'protein_domains.RDs', package = 'maftools')
  gff0 = readRDS(file = gff0) %>% select(-12)
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
  tdata_count2 <- tdata_count %>% filter(n>=minN)
  
  lolli_plot <- tdata_count %>% 
    ggplot(aes(x=AAPosition,y=n))+
    geom_segment(aes(x=AAPosition,xend=AAPosition,y=-0.03*ymax,yend=n),size=0.2)+
    geom_point(aes(fill=Variant_Classification),pch=21,size=4)+
    ggrepel::geom_text_repel(
      data = tdata_count2,
      aes(label=AAChange,color=Variant_Classification),
      force_pull   = 0, # do not pull toward data points
      nudge_y      = 0.1,
      direction    = "x",
      angle        = 90,
      vjust        = 0.1,
      hjust = -1,
      segment.size = 0.1,
      max.iter = 1e4, max.time = 1
    ) +
    geom_rect(aes(xmin=0,xmax=aasize,ymin=-0.06*ymax,ymax=-0.03*ymax),fill="#cccccc",show.legend=FALSE)+
    geom_rect(data=gff,aes(xmin=Start,xmax=End,ymin=-0.07*ymax,ymax=-0.02*ymax),fill=gff$color,show.legend=FALSE)+
    #geom_text(data=gff,aes(x=AAPosition,y=-0.045*ymax,label=Label),family = 'Roboto Condensed',fontface = 'plain',size=3.5)+
    scale_fill_manual(values =landscape_colors,breaks = tdata_count$Variant_Classification)+
    scale_color_manual(values =landscape_colors,breaks = tdata_count$Variant_Classification,guide='none')+
    scale_x_continuous(breaks = pretty_breaks(n = 8),limits = c(0,aasize),expand = expansion(mult = c(0.01,0.01)))+
    scale_y_continuous(breaks = pretty_breaks(),limits = c(-0.08*ymax,ymax),expand = expansion(mult = c(0,0)))+
    guides(color='none',fill=guide_legend(nrow=1))+
    theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid=FALSE,axis = T,axis_col = 'black',ticks = T)+
    theme(legend.position = 'bottom',legend.text=element_text(size=12))+ 
    labs(x='',y=paste0('# ',gene,' Mutations'),fill="Mutation Type",subtitle = genetitle)
  
  
  
  if(domain_annotation){ lolli_plot <- lolli_plot + ggrepel::geom_text_repel(data=gff,aes(x=AAPosition,y=-0.045*ymax,label=Label),family = 'Roboto Condensed',fontface = 'plain',size=3.5,direction = 'x')}
  
  #ggsave(filename = 'tmp.pdf',width = 18,height = 6,device = cairo_pdf)
  #tdata
 
  return(list(lolli_plot, tdata))
  
}
  
genomic_features_wordcloud <- function(data, word_alts, tumor_barcode){
  
  #color_list <- read_delim("../../../landscape_colors_all.csv", delim = ",")
  landscape_colors <- landscape_colors_fcn()
  #print("pre_filtering")
  #print(head(data))
  data <- data %>% filter(Tumor_Barcode== tumor_barcode)
  #print("filtering by tumor barcode selected")
  
  data <- data %>% filter(!is.na(name))
  
  i <- 1
  new_colors <- c()
  for(each in word_alts){
    if(!(word_alts[i] %in% names(landscape_colors))){
      data_alt <- word_alts[i]
      # print(data_alt)
      # assign color to each alteration without a color already set
      # landscape_colors[data_alts[i]] <- sample(color_list$Value, 1)
      # new_colors[data_alt] <- sample(color_list$Color, 1)
      new_colors[data_alt] <- rand_color(n = 1)
      while(new_colors[data_alt] %in% landscape_colors){
        new_colors[data_alt] <- rand_color(n = 1)
      }
      # print(new_colors)
      # print(paste0(data_alts[i], landscape_colors[data_alts[i]]))
      # names(new_colors)[i] <- data_alt
    }
    i <- i + 1
    # print(landscape_colors)
  }

  landscape_colors <- c(landscape_colors, new_colors)
  #print(head(landscape_colors))
  landscape_colors <- landscape_colors[which(names(landscape_colors) %in% data$Alteration)] %>% c()

  #print(landscape_colors)
  
  # make into tibble to create quick plot to pull legend from to put under wordclouds generated
  color_plot <- tibble(Alteration=names(landscape_colors), color_code=landscape_colors, tmp=2)
  color_plot <- color_plot %>% filter(Alteration %in% unique(data$Alteration))
  p <- ggplot(color_plot, aes(x=color_code, y=tmp)) + geom_point(aes(fill=Alteration), pch=21, size=4) + scale_fill_manual(values=landscape_colors) + 
    guides(color='none',fill=guide_legend(nrow=2)) + theme(legend.title=element_text(family = 'Roboto Condensed'),legend.text=element_text(family = 'Roboto Condensed')) #+ theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid=FALSE,axis = T,axis_col = 'black',ticks = T) + theme(legend.position="bottom")
  #p_legend <- p + guides(color='none',fill=guide_legend(nrow=2))
  #p_legend <- get_plot_component(p, 'guide-box', return_all = TRUE)
  p_legend <- get_legend(p) 
  print(paste0('p_legend: ', p_legend))
  #_legend <- guide_legend(nrow=1)
  
  plotlist <- list()
  for(each in sort(unique(data$Wordcloud_cat))){ # Driver_Events, Others
    print(each)
    data_each <- data %>% filter(Wordcloud_cat==each)
    print(dim(data_each))
    set.seed(50)
    wordcloud_plot <- ggplot(data_each, aes(label= Gene, size= Freq, color = Alteration)) + geom_text_wordcloud_area() + scale_size_area(max_size= 18) + 
      scale_color_manual(values=landscape_colors) + theme_ipsum_rc()+ theme(panel.background = element_rect(fill="gray98")) +
      theme(plot.margin = unit(c(0,0,0,0), "cm"))

    plotlist <- list.append(plotlist, wordcloud_plot)
  }
  
  print(plotlist)
  
  # put plots together in a plot_grid
  num_plots <- length(plotlist) # should be 1 or 2 (driver events and others are the only possibilities)
  reps <- rep(paste0("plotlist[[",1:num_plots,"]]"))
  reps <- paste0(reps, collapse= ",")
  #plot_grid
  wordcloud_plotgrid <- paste0("wordcloud_plotgrid_final <- plot_grid(nrow=1,ncol=length(plotlist), rel_widths=0.5, labels=c('Driver Events','Others'), label_size= 20, label_fontface='bold',", reps, ")")
  eval(parse(text=wordcloud_plotgrid))
  
  wordcloud_and_legend <- plot_grid(wordcloud_plotgrid_final, p_legend, nrow=2, rel_heights=c(10,3))
  
  # set.seed(50)
  # wordcloud_plot <- ggplot(data, aes(label= Gene, size= Freq, color= Freq)) + 
  #   geom_text_wordcloud(show.legend = TRUE) + scale_color_gradient(low = "darkred", high = "red") + #scale_color_manual(values=landscape_colors) +
  #   #geom_text_wordcloud_area(area_corr_power= 1, eccentricity= 1, show.legend = TRUE) + scale_color_manual(values=landscape_colors) + 
  #   guides(color='none',fill=guide_legend(nrow=1))+
  #   theme_ipsum_rc(strip_text_size= 14, base_family = "Roboto Condensed") +
  #   #theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid=FALSE,axis = T, axis_col= 'black') +
  #   #facet_wrap(vars(Wordcloud_cat), nrow=2, ncol= 1) +
  #   theme(legend.position = 'right',legend.text=element_text(size=12)) + labs(title= "Driver_Events")#scale_size_area(max_size=8) +
    
  return(wordcloud_and_legend)

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

# mutation summary

mutation_summary <- function(data,samples){
  
  landscape_colors <- landscape_colors_fcn()
  #sample_input <- c('IGC-02-1001-T03','IGC-02-1016-T01') #c('IGC-02-1001-T03') c('IGC-02-1001-T03','IGC-02-1016-T01')
  # maf_sample <- data %>% filter(Tumor_Barcode %in% c(samples))
  # maf_sample <- data %>% filter(SP_Group %in% samples)
  maf_sample <- data
  #print(paste0('maf_sample:', maf_sample))
  
  tmp <- maf_sample %>% 
    count(Variant_Classification) %>% 
    #mutate(Variant_Classification = gsub('_',' ',Variant_Classification)) %>% 
    mutate(Variant_Classification = fct_reorder(Variant_Classification,desc(n)))
  
  var_class_order <- levels(tmp$Variant_Classification)
  
  # variant classification plot
  var_class_plot <- maf_sample %>% 
    count(Variant_Classification) %>% 
    #mutate(Variant_Classification = factor(Variant_Classification, levels = rev(var_class_order))) %>%
    #mutate(Variant_Classification = gsub('_',' ',Variant_Classification)) %>%
    mutate(Variant_Classification = fct_reorder(Variant_Classification,n)) %>%
    ggplot(aes(y=Variant_Classification,x=n,fill=Variant_Classification)) +
    geom_col(col = 'black', linewidth = 0.1, width = 0.8) +
    scale_x_continuous(expand = c(0,0)) +
    #geom_col(width = 0.7, col = 'black', linewidth = 0.4) + 
    ggtitle('Variant Classification') + 
    theme_ipsum_rc(plot_title_face = 'plain') + theme(legend.position = 'none',axis.title.x=element_blank(), axis.text.y = element_text(size = 8), axis.title.y=element_blank(), axis.ticks.y=element_blank(),panel.border = element_blank()) +
    scale_fill_manual(values = landscape_colors)
    #scale_y_discrete(labels = gsub('_',' ',unique(Variant_Classification)))
    #theme_bw(base_size = 15) +
    #theme(plot.title = element_text(size=14), axis.title.x=element_blank(), axis.text.y = element_text(size = 8), axis.title.y=element_blank(), axis.ticks.y=element_blank(), legend.position = 'none',panel.border = element_blank())
  
  var_type_plot <- maf_sample %>%
    count(Variant_Type) %>%
    mutate(Variant_Type = fct_reorder(Variant_Type,n)) %>%
    ggplot(aes(y=Variant_Type,x=n,fill=Variant_Type)) +
    geom_col(col = 'black', linewidth = 0.1, width = 0.8) +
    scale_x_continuous(expand = c(0,0)) +
    #geom_col(width = 0.7, col = 'black', linewidth = 0.4) + 
    ggtitle('Variant Type') + 
    theme_ipsum_rc(plot_title_face = 'plain') + theme(legend.position = 'none',axis.title.x=element_blank(), axis.text.y = element_text(size = 8), axis.title.y=element_blank(), axis.ticks.y=element_blank(),panel.border = element_blank()) +
    #scale_fill_manual(values = color_list$Color)
    scale_fill_manual(values = landscape_colors)
    #theme_bw(base_size = 15) +
    #theme(plot.title = element_text(size=14), axis.title.x=element_blank(), axis.text.y = element_text(size = 8), axis.title.y=element_blank(), axis.ticks.y=element_blank(), legend.position = 'none',panel.border = element_blank())

  #snv_colors <- oncoplot_colors$Color[which(oncoplot_colors$Name %in% unique(tmp$snv_class))] %>% c()
  
  snv_class_plot <- maf_sample %>%
    mutate(snv_class=paste0(Reference_Allele,">",Tumor_Seq_Allele2)) %>%
    count(snv_class) %>%
    filter(!nchar(snv_class) > 3) %>%
    mutate(snv_class=if_else(snv_class=='A>C', 'T>G',
                             if_else(snv_class=='A>G', 'T>C',
                                     if_else(snv_class=='A>T','T>A',
                                             if_else(snv_class=='G>A','C>T',
                                                     if_else(snv_class=='G>C','C>G',
                                                             if_else(snv_class=='G>T','C>A',snv_class))))))) %>%
    bind_rows() %>% group_by(snv_class) %>% summarise(n= sum(n))%>%
    mutate(n_total=sum(n)) %>% mutate(snv_class = factor(snv_class, levels = c('C>A','C>G','C>T','T>C','T>A','T>G'))) %>% 
    ggplot(aes(y=snv_class,x=n/n_total,fill=snv_class)) +
    geom_col(col = 'black', linewidth = 0.1, width = 0.8) + geom_text(aes(x= (n/n_total + 0.05),label=n), family = 'Roboto Condensed', fontface = 'bold') +
    scale_x_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0.0, 0.25, 0.50, 0.75, 1.00)) +
    #geom_col(width = 0.7, col = 'black', linewidth = 0.4) +
    ggtitle('SNV Class') + 
    theme_ipsum_rc(plot_title_face = 'plain') +
    theme(legend.position = 'none',axis.title.x=element_blank(), axis.text.y = element_text(size = 8), axis.title.y=element_blank(), axis.ticks.y=element_blank(),panel.border = element_blank()) +
    scale_fill_manual(values = landscape_colors)
    # theme_bw(base_size = 15) +
    # theme(plot.title = element_text(size=14), axis.title.x=element_blank(), axis.text.y = element_text(size = 8), axis.title.y=element_blank(), axis.ticks.y=element_blank(), legend.position = 'none',panel.border = element_blank())

  sample_tab <- maf_sample %>% group_by(Tumor_Barcode)%>% count() %>% arrange(desc(n))
  sample_order <- maf_sample %>% group_by(Tumor_Barcode)%>% count() %>% arrange(desc(n)) %>% pull(Tumor_Barcode)
  median_var_per_sample <- median(sample_tab$n)
  
  variants_per_sample_plot <- maf_sample %>% 
    group_by(Tumor_Barcode) %>% mutate(Tumor_Barcode = factor(Tumor_Barcode, levels = sample_order)) %>%
    mutate(Variant_Classification = factor(Variant_Classification, levels = var_class_order)) %>% 
    count(Variant_Classification) %>% 
    ggplot(aes(x=Tumor_Barcode,y=n,fill=Variant_Classification)) +
    geom_col(col = 'black', linewidth = 0.1, width = 0.8) +
    geom_hline(aes(yintercept = median_var_per_sample,color = 'red'), linetype = 2) +
    ggtitle(label = 'Variants per sample', subtitle = paste0('Median: ',median_var_per_sample)) +
    theme_ipsum_rc(plot_title_face = 'plain') +
    theme(legend.position = 'none',axis.text.x = element_blank(),axis.title.x=element_blank(), axis.text.y = element_text(size = 8), axis.title.y=element_blank(), axis.ticks.y=element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_manual(values = landscape_colors)
    # theme_bw(base_size = 15) +
    # theme(plot.title = element_text(size=14), 
    #     axis.text.x = element_blank(), axis.title.x=element_blank(), axis.text.y = element_text(size = 8), axis.title.y=element_blank(), axis.ticks.x=element_blank(), legend.position = 'none',panel.border = element_blank())
  
  variant_classification_summary_boxplot <- maf_sample %>%
    count(Tumor_Barcode,Variant_Classification) %>% 
    #mutate(Variant_Classification = fct_reorder(Variant_Classification,desc(n))) %>%
    mutate(Variant_Classification = factor(Variant_Classification, levels = var_class_order)) %>%
    ggplot(aes(x=Variant_Classification,y=n,color=Variant_Classification)) +
    geom_boxplot() +
    ggtitle('Variant Classification Summary') + 
    theme_ipsum_rc(plot_title_face = 'plain') +
    theme(legend.position = 'none',axis.text.x = element_blank(),axis.title.x=element_blank(), axis.text.y = element_text(size = 8), axis.title.y=element_blank(), axis.ticks.y=element_blank(),panel.border = element_blank()) +
    scale_color_manual(values = landscape_colors)
    # theme_bw(base_size = 15) +
    # theme(plot.title = element_text(size=14), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(size = 8), axis.title.y=element_blank(), axis.ticks.y=element_blank(), legend.position = 'none',panel.border = element_blank())
  
  #sherlock_maf_freq for top 10 mutated genes (symbol, n, frequency) or do calculation
  sizedata <- n_distinct(maf_sample$Tumor_Barcode)
  frequency_data <- maf_sample %>% 
    count(Tumor_Barcode, Hugo_Symbol) %>%
    select(-n) %>%
    count(Hugo_Symbol) %>% 
    arrange(desc(n)) %>%
    mutate(Freq=100*n/sizedata) %>% 
    mutate(nMutSample=n) %>%
    select(-n) %>% arrange(desc(Freq)) #%>% mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = Hugo_Symbol)) #%>% slice_min(n=10, order_by = desc(Freq))
  
  frequency_top10 <- frequency_data[1:10,] %>% arrange(nMutSample) %>% mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = Hugo_Symbol))
  
  # variant_classification_data_per_gene <- maf_sample %>%
  #   count(Hugo_Symbol,Variant_Classification) %>% group_by(Hugo_Symbol) %>% mutate(total=sum(n)) %>% arrange(desc(total)) %>% distinct(Hugo_Symbol)
  # 
  # variant_classification_data_per_gene <- variant_classification_data_per_gene$Hugo_Symbol[1:10]
  # 
  # x_lim_tmp <- maf_sample %>%
  #   count(Hugo_Symbol,Variant_Classification) %>% group_by(Hugo_Symbol) %>% mutate(total=sum(n)) %>% arrange(desc(total)) 
  # 
  # x_limit <- x_lim_tmp$total[1]
  
  # n_vals <- maf_sample %>% 
  #   filter(Hugo_Symbol %in% variant_classification_data_per_gene) %>%
  #   #group_by(Hugo_Symbol) %>%
  #   count(Hugo_Symbol, Variant_Classification) %>% 
  #   left_join(frequency_data, by = 'Hugo_Symbol') %>% group_by(Hugo_Symbol) 
  
  tmp <- maf_sample %>% filter(Hugo_Symbol %in% frequency_top10$Hugo_Symbol) %>% group_by(Tumor_Barcode) %>% count(Variant_Classification, Hugo_Symbol) 
  
  #test <- test[1:12,]
  
  #test <- test %>% filter(Tumor_Barcode =='17MH0022_4R')
  
  #test <- maf_sample %>% filter(Hugo_Symbol %in% frequency_top10$Hugo_Symbol) %>% group_by(Tumor_Barcode) %>% count(Variant_Classification, Hugo_Symbol)  %>% filter(Tumor_Barcode == '5V8LTAZB')
  
  # count tumor barcodes to make sure all are being processed
  new_df <- tibble('Tumor_Barcode' = character(), 'Variant_Classification' = factor(), 'Hugo_Symbol' = character(), 'n' = integer())
  for(each in unique(tmp$Tumor_Barcode)){
    print(each)
    tmp2 <- tmp %>% filter(Tumor_Barcode == each)
    #print(tmp)
    if(nrow(tmp2) == 1){
      #next
      print('all good here, going to the next one')
      new_df <- add_row(new_df, tmp2)
    }else{
      length_var_class <- length(unique(tmp2$Variant_Classification))
      if(length_var_class > 1){
        length_genes <- length(unique(tmp2$Hugo_Symbol))
        print(length_genes)
        for(each in unique(tmp2$Hugo_Symbol)){
          tmp3 <- tmp2 %>% filter(Hugo_Symbol == each) #%>% count(Variant_Classification)
          print(tmp3)
          if(nrow(tmp3) > 1){
            color <- 'black'
            tmp3 <- tmp3 %>% mutate(Variant_Classification = 'Multiple') #%>% mutate(n = 0) %>% distinct()
            new_df <- add_row(new_df, tmp3)
          }else{
            color <- 'original'
            new_df <- add_row(new_df, tmp3)
          }
          
          print(color)
          
        }
      }else{
        new_df <- add_row(new_df, tmp2)
      }
    }
  }
  
  new_df <- new_df %>% distinct()
  
  # new_df <- new_df %>% select(-n) %>%
  #   count(Tumor_Barcode, Hugo_Symbol) %>%
  #   select(-n) %>%
  #   count(Hugo_Symbol) %>%
  #   arrange(desc(n)) %>%
  #   mutate(Freq=100*n/sizedata)
  
  top10_mutated_genes <- new_df %>% group_by(Hugo_Symbol, Tumor_Barcode) %>% mutate(nMutSample = length(unique(Tumor_Barcode))) %>% mutate(Freq = nMutSample/sizedata*100) %>% mutate(Freq = sum(Freq))
  #top10_mutated_genes_plot$nMutSample[duplicated(top10_mutated_genes_plot$Hugo_Symbol)] <- 0
  top10_mutated_genes <- top10_mutated_genes %>% select(-n) %>% distinct() %>% group_by(Hugo_Symbol) %>% mutate(Freq = nMutSample/sizedata*100) %>% mutate(Freq = paste0(round(sum(Freq),1),"%"))
  
  max_x_limit <- top10_mutated_genes %>% group_by(Hugo_Symbol) %>% reframe(max_nMutSample = sum(nMutSample)) %>% reframe(max_x = max(max_nMutSample)) %>% pull(max_x)
  max_x_limit <- max_x_limit + 100
  
  top10_mutated_genes_plot <- top10_mutated_genes %>% group_by(Hugo_Symbol) %>%
    # mutate(Variant_Classification = factor(Variant_Classification, levels = var_class_order)) %>%
    # mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = rev(variant_classification_data_per_gene))) %>%
    ggplot(aes(x=nMutSample, y=Hugo_Symbol, fill = Variant_Classification)) + geom_col(col = 'black', linewidth = 0.1, width = 0.8)+
    #geom_text(aes(y = Hugo_Symbol, label = Freq, group = Hugo_Symbol))+
    scale_x_continuous(expand = c(0,0), limits = c(0,max_x_limit)) +
    scale_fill_manual(values = landscape_colors)+ theme_ipsum_rc(plot_title_face = 'plain') + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = 'none')
    
  # need labels for frequency % at the end of each bar
  mutSamples <- top10_mutated_genes %>% group_by(Hugo_Symbol) %>% mutate(Tumor_Barcode = 'Sample',Variant_Classification = 'Var_Class', nMutSample = sum(nMutSample)) %>% distinct()
  
  #top10_mutated_genes_plot <- top10_mutated_genes_plot + geom_text_repel(data = mutSamples, aes(x = nMutSample, label = Freq),family = 'Roboto Condensed', color = "black", size = 5) #, nudge_y = 0.0, nudge_x = 1)

  top10_mutated_genes_plot <- top10_mutated_genes_plot + geom_text(data = mutSamples, aes(x= nMutSample,label=Freq), nudge_x = 25, family = 'Roboto Condensed', fontface = 'bold') + ggtitle('Top 10 Mutated Genes') 
    # based on raw count instead of frequency
  # top10_mutated_genes_plot <- maf_sample %>% 
  #   filter(Hugo_Symbol %in% variant_classification_data_per_gene) %>%
  #   #group_by(Hugo_Symbol) %>%
  #   count(Hugo_Symbol, Variant_Classification) %>% 
  #   left_join(frequency_data, by = 'Hugo_Symbol') %>% group_by(Hugo_Symbol) %>%
  #   mutate(Variant_Classification = factor(Variant_Classification, levels = var_class_order)) %>%
  #   mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = rev(variant_classification_data_per_gene))) %>% 
  #   #filter(Hugo_Symbol %in% variant_classification_data_per_gene) %>%
  #   ggplot(aes(x=n, y=Hugo_Symbol, fill = Variant_Classification)) +
  #   geom_col(col = 'black', linewidth = 0.1, width = 0.8) + 
  #   geom_text(aes(x=max(n) + 125,label=paste0(round(Freq, digits = 1),"%")), family = 'Roboto Condensed', fontface = 'bold') +
  #   scale_x_continuous(expand = c(0,0), limits = c(0,x_limit + 100)) +
  #   ggtitle('Top 10 Mutated Genes') +
  #   theme_ipsum_rc(plot_title_face = 'plain') +
  #   theme(legend.position = 'none',axis.title.x=element_blank(), axis.text.y = element_text(size = 8), axis.title.y=element_blank(), axis.ticks.y=element_blank(),panel.border = element_blank()) +
  #   scale_fill_manual(values = oncoplot_colors$Color[9:1])

  # top10_mutated_genes_plot <- frequency_data %>% 
  #   ggplot(aes(x=Freq, y=Hugo_Symbol)) +
  #   geom_col(col = 'black', linewidth = 0.1, width = 0.8)
    
  mut_summary_all_plot <- plot_grid(var_class_plot, var_type_plot, snv_class_plot,variants_per_sample_plot,variant_classification_summary_boxplot,top10_mutated_genes_plot)
  return(mut_summary_all_plot)
}

tcgaCompareWGS2 <- function (project_code, maf.mutload, capture_size = NULL, tcga_capture_size = 1, cohortName = NULL, ytext="Somatic mutation prevalence\n(number mutations per megabase log10)",
                             tcga_cohorts = NULL, primarySite = FALSE, col = c("gray70", 
                                                                               "black"), bg_col = c("#EDF8B1", "#2C7FB8"), medianCol = "red", 
                             logscale = TRUE, rm_hyper = FALSE,tcga_cohorts_file=NULL,cohortsize=NULL) 
{
  par(family = "Roboto Condensed")
  if(is.null(tcga_cohorts_file)){
    top('Please provide TCGA cohort file!!!')
  }
  load(tcga_cohorts_file)
  
  if (primarySite) {
    tcga.cohort = tcga.cohort[, .(Tumor_Sample_Barcode, total, 
                                  site)]
    colnames(tcga.cohort)[3] = "cohort"
  }
  else {
    tcga.cohort = tcga.cohort[, .(Tumor_Sample_Barcode, total, 
                                  cohort)]
  }
  if (!is.null(tcga_cohorts)) {
    tcga.cohort = tcga.cohort[cohort %in% tcga_cohorts]
    if (nrow(tcga.cohort) == 0) {
      stop("Something went wrong. Provide correct names for 'tcga_cohorts' arguments")
    }
  }
  
  
  if (is.null(cohortName)) {
    cohortName = paste0("Input", seq_len(length(maf.mutload)))
  }
  else if (length(cohortName) != length(maf.mutload)) {
    stop("Please provide names for all input cohorts")
  }
  names(maf.mutload) = cohortName
  maf.mutload = data.table::rbindlist(l = maf.mutload, idcol = "cohort")
  tcga.cohort$total = as.numeric(as.character(tcga.cohort$total))
  maf.mutload$total = as.numeric(as.character(maf.mutload$total))
  if (rm_hyper) {
    tcga.cohort = split(tcga.cohort, as.factor(tcga.cohort$cohort))
    tcga.cohort = lapply(tcga.cohort, function(x) {
      xout = boxplot.stats(x = x$total)$out
      if (length(xout) > 0) {
        message(paste0("Removed ", length(xout), " outliers from ", 
                       x[1, cohort]))
        x = x[!total %in% xout]
      }
      x
    })
    tcga.cohort = data.table::rbindlist(l = tcga.cohort)
    xout = boxplot.stats(x = maf.mutload$total)$out
    if (length(xout) > 0) {
      message(paste0("Removed ", length(xout), " outliers from Input MAF"))
      maf.mutload = maf.mutload[!total %in% xout]
    }
  }
  if (!is.null(capture_size)) {
    maf.mutload[, `:=`(total, total/capture_size)]
    tcga.cohort[, `:=`(total, total/tcga_capture_size)]
  }
  tcga.cohort = rbind(tcga.cohort, maf.mutload)
  tcga.cohort.med = tcga.cohort[, .(.N, median(total)), cohort][order(V2, 
                                                                      decreasing = TRUE)]
  tcga.cohort$cohort = factor(x = tcga.cohort$cohort, levels = tcga.cohort.med$cohort)
  colnames(tcga.cohort.med) = c("Cohort", "Cohort_Size", "Median_Mutations")
  tcga.cohort$TCGA = ifelse(test = tcga.cohort$cohort %in% 
                              cohortName, yes = "Input", no = "TCGA")
  tcga.cohort = split(tcga.cohort, as.factor(tcga.cohort$cohort))
  plot.dat = lapply(seq_len(length(tcga.cohort)), function(i) {
    x = tcga.cohort[[i]]
    x = data.table::data.table(rev(seq(i - 1, i, length.out = nrow(x))), 
                               x[order(total, decreasing = T), total], x[, TCGA])
    x
  })
  names(plot.dat) = names(tcga.cohort)
  if (logscale) {
    y_lims = range(log10(data.table::rbindlist(l = plot.dat)[, 
                                                             V2]))
  }
  else {
    y_lims = range(data.table::rbindlist(l = plot.dat)[, 
                                                       V2])
  }
  
  y_lims[y_lims== -Inf] <- -3
  y_max = ceiling(max(y_lims))
  y_min = floor(min(y_lims))
  y_lims = c(y_min, y_max)
  y_at = pretty(y_lims)
  par(mar = c(6, 4, 3, 1))
  #dev.off()
  #print(dev.cur())
  pdf(file = paste0(project_code, '_', cohortName,'_tmb.pdf'),width = 14,height = 5,useDingbats = FALSE) #######
  plot(NA, NA, xlim = c(-0.15, length(plot.dat)+0.15), ylim = y_lims, xaxs='i',yaxs='i',
       axes = FALSE, xlab = NA, ylab = NA)
  rect(xleft = seq(0, length(plot.dat) - 1, 1), ybottom = min(y_lims), 
       xright = seq(1, length(plot.dat), 1), ytop = y_max, col = grDevices::adjustcolor(col = bg_col, 
                                                                                        alpha.f = 0.2), border = NA)
  abline(h = pretty(y_lims), lty = 2, col = "gray70")
  ## define corhort col
  cohortCol <- col[-1]
  names(cohortCol) <- cohortName
  plotCol <- rep("black",length(names(plot.dat)))
  names(plotCol) <- names(plot.dat)
  plotCol[names(cohortCol)] <- cohortCol
  
  lapply(seq_len(length(plot.dat)), function(i) {
    x = plot.dat[[i]]
    if (x[1, V3] == "Input") {
      if (logscale) {
        points(x$V1, log10(x$V2), pch = 16, cex = 0.4, 
               col = as.character(cohortCol[names(plot.dat)[i]]))
        #plotCol[i] <- as.character(cohortCol[names(plot.dat)[i]])
        # print(c(i,ccc,col[ccc+1]))
        # print(names(plot.dat))
      }
      else {
        points(x$V1, x$V2, pch = 16, cex = 0.4, col = as.character(cohortCol[names(plot.dat)[i]]))
        #plotCol[i] <- as.character(cohortCol[names(plot.dat)[i]])
        
        # print(c(i,col[ccc+1]))
        ccc <- ccc + 1 
        
      }
    }
    else {
      if (logscale) {
        #print(x$V1)
        points(x$V1, log10(x$V2), pch = 16, cex = 0.4, 
               col = col[1])
      }
      else {
        #print(x$V1)
        points(x$V1, x$V2, pch = 16, cex = 0.4, col = col[1])
      }
    }
  })
  axis(side = 2, at = y_at, las = 2, line = -1, tick = FALSE)
  samp_sizes = lapply(plot.dat, nrow)
  if(!is.null(cohortsize)){
    samp_sizes[cohortName] <- cohortsize 
  }
  
  for(i in seq_len(length(names(plot.dat)))){
    axis(side = 1, at = seq(0.5, length(plot.dat) - 0.5, 1)[i], labels = paste0(names(plot.dat)[i],"  "), las = 2, tick = FALSE, line = -1,col.axis = plotCol[names(plot.dat)[i]])
    axis(side = 3, at = seq(0.5, length(plot.dat) - 0.5, 1)[i], labels = paste0("  ",unlist(samp_sizes[i])), las = 2, tick = FALSE, line = -1.2, col.axis = plotCol[names(plot.dat)[i]], font = 3)
  }
  
  #print(plotCol)
  
  if (logscale) {
    if (is.null(capture_size)) {
      mtext(text = "TMB", side = 2, line = 1.2)
    }
    else {
      mtext(text = ytext, side = 2, line = 1.2)
    }
  }
  else {
    if (is.null(capture_size)) {
      mtext(text = "TMB", side = 2, line = 1.2)
    }
    else {
      mtext(text = "TMB log10", side = 2, line = 1.6)
    }
  }
  tcga.cohort.med[, `:=`(Median_Mutations_log10, log10(Median_Mutations))]
  if (logscale) {
    lapply(seq_len(nrow(tcga.cohort.med)), function(i) {
      segments(x0 = i - 1, x1 = i, y0 = tcga.cohort.med[i, 
                                                        Median_Mutations_log10], y1 = tcga.cohort.med[i, 
                                                                                                      Median_Mutations_log10], col = medianCol)
    })
  }
  else {
    lapply(seq_len(nrow(tcga.cohort.med)), function(i) {
      segments(x0 = i - 1, x1 = i, y0 = tcga.cohort.med[i, 
                                                        Median_Mutations], y1 = tcga.cohort.med[i, Median_Mutations], 
               col = medianCol)
    })
  }
  box()
  dev.off()
  #tcga.cohort.med
}


# generate necessary data and clustering output to use in scna_clustering_part_two, heatmap generation
scna_clustering_part_one <- function(data1, data2){
  
  data1 <- validate_vardf(data1)
  # specific to sherlock
  #sample_order <- BBsolution4 %>% select(Tumor_Barcode) %>% mutate(Seq=seq_along(Tumor_Barcode)) 
  
  # general across all projects
  sample_order <- data1 %>% select(Tumor_Barcode) %>% mutate(Seq=seq_along(Tumor_Barcode)) 

  # read subclone data  -----------------------------------------------------
  
  #cnv <- BBprofile %>% mutate(frac1_A=as.numeric(frac1_A),frac2_A=as.numeric(frac2_A),nMaj1_A=as.integer(nMaj1_A),nMaj2_A=as.integer(nMaj2_A), nMin1_A=as.integer(nMin1_A),nMin2_A=as.integer(nMin2_A))
  cnv <- data2 %>% mutate(frac1_A=as.numeric(frac1_A),frac2_A=as.numeric(frac2_A),nMaj1_A=as.integer(nMaj1_A),nMaj2_A=as.integer(nMaj2_A), nMin1_A=as.integer(nMin1_A),nMin2_A=as.integer(nMin2_A))
  
  #cnv %>% filter(chr %in% c(1:22,"X","Y"),Tumor_Barcode %in% sherlock_samples_unique$Tumor_Barcode) %>% select(Tumor_Barcode,everything()) %>% write_delim('bb_subclone_232.txt',delim = '\t',col_names = T)
  cnv <- cnv %>% filter(chr %in% c(1:22,"X","Y"),Tumor_Barcode %in% sample_order$Tumor_Barcode)
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
  
  
  # Define the value --------------------------------------------------------
  totalsamples <- cnv %>% count(Tumor_Barcode) %>% pull(Tumor_Barcode) %>% length()
  allsamples <- cnv %>% count(Tumor_Barcode) %>% pull(Tumor_Barcode) %>% unique()
  
  
  cnvdata <- cnv %>%
    filter(chr %in% c(1:22)) %>% 
    mutate(clone_total=clone_nMin+clone_nMaj) %>% 
    left_join(
      #BBsolution4 %>% select(Tumor_Barcode,WGD_Status,Tumor_Purity,Tumor_Ploidy=BB_Ploidy)
      data1 %>% select(Tumor_Barcode,WGD_Status,Tumor_Purity,Tumor_Ploidy=BB_Ploidy)
    ) %>% 
    mutate(WGD_Status=if_else(is.na(WGD_Status),'nWGD',WGD_Status)) %>% 
    mutate(relative_copy=clone_total-if_else(WGD_Status=="WGD",4,2)) %>% 
    mutate(relative_copy=if_else(relative_copy>4,4,relative_copy)) %>% 
    mutate(relative_copy=if_else(relative_copy< -4,-4,relative_copy)) %>% 
    mutate(relative_copy=if_else(clone_nMaj==0,-4,relative_copy)) %>% 
    mutate(startpos=as.integer(startpos),endpos=as.integer(endpos)) 
  
  
  # Dendrograms -------------------------------------------------------------
  #source('~/NIH-Work/EAGLE_NSLC/FirstPaper/Biowulf/ggdendro.R')
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
  
  cols <- c("#a9a9a9", "#2ca02c", "#ff7f0e","#1f77b4")
  
  p_dend <- plot_ggdendro(hcdata,
                          direction   = "lr",
                          scale.color = cols,
                          label.size  = 2.5,
                          branch.size = 0.5,
                          expand.y    = 0,nudge.label = 0)
  
  p_dend <- p_dend+
    scale_x_continuous(breaks = NULL,limits = c(1, nrow(label(hcdata))),expand = c(0,0))+
    theme_void()
  
  sample_order <- as_tibble(hcdata$labels) %>% select(Tumor_Barcode=label) %>% unique() %>% mutate(Seq=seq_along(Tumor_Barcode)) 
  
  #chr19loss <- sample_order %>% tail(n=11) %>% mutate(chr19loss="Y") %>% select(-Seq)
  #chr19loss <- cnvdata %>% filter(chr=="19",clone_nMaj==0) %>% mutate(size=endpos-startpos) %>% filter(size>30000000) %>% mutate(chr19loss="Y") %>% select(Tumor_Barcode,chr19loss)
  #save(chr19loss,file='chr19loss.RData')
  
  cnvclust <- as_tibble(hcdata$labels) %>% select(Tumor_Barcode=label,CNV_Clust=clust) %>% mutate(CNV_Clust=paste0('C',CNV_Clust))
  #cnvclust <- cnvclust %>% mutate(SCNA_Group=if_else(CNV_Clust == 'C1', 'Forte',if_else(CNV_Clust == 'C2', 'Piano', 'Mezzo-forte'))) 
  #save(cnvclust,file='cnvclust_wgs_all.RData')
  #save(cnvclust,file='cnvclust_wgs_all.RData')
  
  return(list(cnvdata, p_dend, sample_order, cnvclust,totalsamples,cnseg))
  
}
  
  # Heatmap plot (clustering output) ------------------------------------------------------------
scna_clustering_part_two <- function(cnvdata, p_dend, sample_order,totalsamples,cnseg, hg38centro, cyto){
  chrlevels <- seq(1:22)
  p_heatmap <- cnvdata %>% 
    mutate(chr=factor(chr,levels = chrlevels)) %>%
    mutate(startpos=as.integer(startpos),endpos=as.integer(endpos)) %>% 
    left_join(sample_order) %>% 
    arrange(Seq,chr) %>% 
    #filter(Seq %in% c(1:20)) %>% 
    ggplot(aes(fill=relative_copy))+
    #geom_segment(aes(x=startpos,xend=endpos,y=Tumor_Barcode,yend=Tumor_Barcode),size=3)+
    geom_rect(aes(xmin=startpos,xmax=endpos,ymin=Seq-0.5,ymax=Seq+0.5))+
    facet_grid(~chr,scales = 'free_x',space = 'free',switch="both")+
    scale_fill_gradient2(low = "blue",high = "red",mid = "white",midpoint = 0)+
    scale_x_continuous(breaks = pretty_breaks(),expand = c(0,0),sec.axis = dup_axis())+
    labs(fill="Relative Copy Number")+
    scale_y_continuous(breaks = sample_order$Seq,labels = sample_order$Tumor_Barcode,expand = expand_scale(add = c(0, 0)), sec.axis = dup_axis())+
    #facet_wrap(~chr,nrow = 1,scales = 'free_x',strip.position = "bottom")+
    labs(x="",y="")+
    #theme_ipsum_rc()+
    theme(text = element_text(family = "Roboto Condensed"),panel.spacing = unit(0, "lines"),axis.text.x = element_blank(),axis.ticks = element_blank(),axis.text.y = element_text(size = 6),strip.placement = "outside",strip.switch.pad.wrap = unit(0, "cm"),strip.text = element_text(face = "bold",size=10),strip.background = element_rect(color = "gray50",fill="white"),panel.background = element_rect(fill = 'white'),axis.line= element_line(colour = 'black'),axis.title.y.right = element_blank(),axis.text.y.right = element_blank(),axis.ticks.y.right = element_blank(),axis.title.x.top = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank() )+
    coord_cartesian(clip="off")
  
  #ggsave('heatmap_all_segements_raw.pdf',width = 14,height = 8,plot = p_heatmap,device = cairo_pdf)
  
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
  # annodata <- sample_order %>% 
  #   left_join(
  #     BBsolution4 %>% select(Tumor_Barcode,WGD_Status,Tumor_Purity,Tumor_Ploidy=BB_Ploidy) 
  #   ) %>% 
  #   left_join(wgs_groups_info %>% select(Tumor_Barcode,Study,Assigned_Population,Smoker=Smoking) %>% unique()) %>% 
  #   left_join(clinical_data)
  # 
  # # purity
  # annodata_1 <- annodata %>% select(Seq,Tumor_Purity) 
  # p_a1 <- annodata_1 %>% 
  #   ggplot(aes("anno1",Seq,fill=Tumor_Purity))+
  #   geom_tile()+
  #   scale_fill_viridis_c(option = "D")+
  #   labs(fill="Purity")+
  #   theme_void()+
  #   scale_y_discrete(breaks = NULL,limits = c(1, nrow(label(hcdata))),expand = c(0,0))
  # 
  # p_a1_legend <- get_legend(
  #   # create some space to the left of the legend
  #   p_a1 + theme(legend.box.margin = margin(8, 0, 5, 0))
  # )
  # 
  # p_a1 <- p_a1+theme(legend.position = "none")
  # 
  # # WGD
  # annodata_2 <- annodata %>% select(Seq,WGD_Status) 
  # p_a2 <- annodata_2 %>% 
  #   ggplot(aes("anno2",Seq,fill=WGD_Status))+
  #   geom_tile()+
  #   scale_fill_aaas()+
  #   theme_void()+
  #   scale_y_continuous(breaks = NULL,expand = c(0,0))
  # #limits = c(1, nrow(label(hcdata)))
  # 
  # 
  # p_a2_legend <- get_legend(
  #   # create some space to the left of the legend
  #   p_a2 + theme(legend.box.margin = margin(8, 0, 5, 0))
  # )
  # 
  # p_a2 <- p_a2+theme(legend.position = "none")
  # 
  # 
  # # Histology
  # annodata_3 <- annodata %>% select(Seq,Histology) 
  # p_a3 <- annodata_3 %>% 
  #   ggplot(aes("anno3",Seq,fill=Histology))+
  #   geom_tile()+
  #   scale_fill_d3(palette = "category20b")+
  #   theme_void()+
  #   scale_y_continuous(breaks = NULL,expand = c(0,0))
  # #limits = c(1, nrow(label(hcdata)))
  # 
  # 
  # p_a3_legend <- get_legend(
  #   # create some space to the left of the legend
  #   p_a3 + theme(legend.box.margin = margin(8, 0, 5, 0))
  # )
  # 
  # p_a3 <- p_a3+theme(legend.position = "none")
  # 
  # 
  # # Smoker
  # annodata_4 <- annodata %>% select(Seq,Smoker) 
  # p_a4 <- annodata_4 %>% 
  #   ggplot(aes("anno4",Seq,fill=Smoker))+
  #   geom_tile()+
  #   scale_fill_npg()+
  #   theme_void()+
  #   scale_y_continuous(breaks = NULL,expand = c(0,0))
  # #limits = c(1, nrow(label(hcdata)))
  # 
  # 
  # p_a4_legend <- get_legend(
  #   # create some space to the left of the legend
  #   p_a4 + theme(legend.box.margin = margin(8, 0, 5, 0))
  # )
  # 
  # p_a4 <- p_a4+theme(legend.position = "none")
  # 
  # # Population  
  # annodata_5 <- annodata %>% select(Seq,Assigned_Population) 
  # p_a5 <- annodata_5 %>% 
  #   ggplot(aes("anno4",Seq,fill=Assigned_Population))+
  #   geom_tile()+
  #   scale_fill_npg()+
  #   theme_void()+
  #   scale_y_continuous(breaks = NULL,expand = c(0,0))
  # #limits = c(1, nrow(label(hcdata)))
  # 
  # 
  # p_a5_legend <- get_legend(
  #   # create some space to the left of the legend
  #   p_a5 + theme(legend.box.margin = margin(8, 0, 5, 0))
  # )
  # 
  # p_a5 <- p_a5+theme(legend.position = "none")
  # Test chr19 loss  ---------------------------------------------------------
  #annodata %>% left_join(chr19loss) %>% mutate(chr19loss=if_else(is.na(chr19loss),"N",chr19loss) ) %>% select(RTK_Altered_Status,chr19loss) %>% mutate(chr19loss=factor(chr19loss,levels = c('Y','N'))) %>% table() %>% fisher.test(alternative = 'greater')
  
  
  # Frequency_Plot ----------------------------------------------------------
  # use bin window as gene id 
  
  hg38 <- left_join(
    cnvdata %>% group_by(chr) %>% arrange(startpos) %>% slice(1) %>% select(chr,startpos) %>% ungroup(),
    cnvdata %>% group_by(chr) %>% arrange(desc(endpos)) %>% slice(1) %>% select(chr,endpos) %>% ungroup()
  ) %>% mutate(len=endpos-startpos+1)
  
  bin=1000000
  hg38bins <- NULL
  for(i in 1:22){
    size=hg38$len[i]
    chr=hg38$chr[i]
    startx=hg38$startpos[i]
    endx=hg38$endpos[i]
    tmp <- tibble(chr=chr,start=seq(from = startx,to=endx,by = bin),end=c(seq(from = startx-1,to=endx,by = bin)[-1],size))
    lastval <- as.integer(tail(tmp,1)[3])
    if( lastval < endx){
      tmp <- bind_rows(tmp,tibble(chr=chr,start=lastval+1,end=endx))
    }
    #tmp <- tmp %>% mutate(start2=start+startx-1,end2=end+startx-1)
    hg38bins <- bind_rows(hg38bins,tmp)
  }
  
  hg38info <- 
    hg38bins %>% 
    mutate(geneid=paste(chr,start,end,sep="_"),genename=geneid) %>%
    select(chrom=chr,start,end,geneid,genename) %>% 
    mutate(start=as.integer(start),end=as.integer(end)) %>% 
    as.data.frame()
  
  
  
  # Freq_LOH ----------------------------------------------------------------
  segdata3 <- cnvdata %>% 
    mutate(LOH=if_else(clone_nMin==0 & clone_nMaj==2,1,0)) %>% 
    mutate(num.mark=1000) %>% 
    select(ID=Tumor_Barcode,chrom=chr,loc.start=startpos,loc.end=endpos,num.mark,seg.mean=LOH) %>%
    mutate(seg.mean=if_else(is.na(seg.mean),0,seg.mean)) %>% 
    as.data.frame()
  segdata3 %>% filter_all(any_vars(is.na(.)))
  cnseg3 <- CNSeg(segdata3)
  
  rdseg3 <- getRS(cnseg3, by = "gene", imput = FALSE, XY = FALSE, what = "mean",geneMap = hg38info)
  reducedseg3 <- rs(rdseg3)
  msegdata3 <- as.matrix(reducedseg3[,-(1:5)])
  msegdata3 <- apply(msegdata3, 2, as.numeric)
  
  
  freqdata_loh <- bind_cols(
    reducedseg3[,1:3],
    as_data_frame(msegdata3)
  ) %>% 
    pivot_longer(cols = -c(chrom,start,end)) %>% 
    as_tibble() %>% 
    mutate(start=as.integer(start),end=as.integer(end))
  
  freqdata_loh <- freqdata_loh %>% 
    filter(value!=0) %>% 
    mutate(calling="LOH")
  
  
  freqdata_loh <- freqdata_loh  %>% 
    count(chrom,start,end,calling) %>% 
    mutate(freq=-n/totalsamples)
  
  # remove cytoband centro
  hg38centro2 <- bed_intersect(hg38centro,freqdata_loh %>% mutate(calling=as.character(calling)),suffix = c('','.y')) %>% select(chrom:Centro)
  
  freqdata_loh <- bed_intersect(freqdata_loh %>% mutate(calling=as.character(calling)),hg38centro,invert = T)
  freqdata_loh <- bind_rows(
    freqdata_loh %>% mutate(Centro=0),
    hg38centro2 %>% mutate(calling="LOH",n=0,freq=-1e-36)
  ) %>% 
    arrange(chrom,start,end) %>% select(-Centro)
  
  
  freqdata_loh <- bind_rows(
    freqdata_loh %>% rename(pos=start) %>% select(-end),
    freqdata_loh %>% rename(pos=end) %>% select(-start)
  ) %>%
    filter(!is.na(calling)) %>% 
    unique()
  
  freqdata_loh <- freqdata_loh %>% arrange(chrom,pos) %>% mutate(type="Del",calling=NA_character_) %>% mutate(chrom=factor(chrom,levels = chrlevels))
  
  
  # Freq_CNV ----------------------------------------------------------------
  # segdata <- cnvdata %>% 
  #   mutate(num.mark=1000) %>% 
  #   select(ID=Tumor_Barcode,chrom=chr,loc.start=startpos,loc.end=endpos,num.mark,seg.mean=relative_copy) %>%
  #   mutate(seg.mean=if_else(is.na(seg.mean),0,seg.mean)) %>% 
  #   as.data.frame()
  # segdata %>% filter_all(any_vars(is.na(.)))
  # cnseg <- CNSeg(segdata)
  rdseg2 <- getRS(cnseg, by = "gene", imput = FALSE, XY = FALSE, what = "mean",geneMap = hg38info)
  reducedseg2 <- rs(rdseg2)
  msegdata2 <- as.matrix(reducedseg2[,-(1:5)])
  msegdata2 <- apply(msegdata2, 2, as.numeric)
  
  
  freqdata <- bind_cols(
    reducedseg2[,1:3],
    as_data_frame(msegdata2)
  ) %>% 
    pivot_longer(cols = -c(chrom,start,end)) %>% 
    as_tibble() %>% 
    mutate(start=as.integer(start),end=as.integer(end))
  
  #version1
  # freqdata <- freqdata %>% 
  #   mutate(calling=case_when(
  #     value < log2(1/2) ~ "-2",
  #     value < 0 & value >=log2(1/2) ~ "-1",
  #     value ==0 ~ NA_character_,
  #     value>0 & value < log2(3/2)  ~ "1",
  #     value > log2(3/2) ~ "2"
  #   )) %>% 
  #   filter(value!=0) %>% 
  #   mutate(calling=factor(calling,levels = c("-2","-1","2","1")))%>%
  #   mutate(type=if_else(calling %in% c("1","2"),"Amp","Del"))
  
  cplevels <- c("Amplification (Copy gain >=4)","Gain","Loss","Homozygous deletion")
  
  freqdata <- freqdata %>% 
    mutate(calling=case_when(
      value <= -4 ~ "Homozygous deletion",
      value < 0 & value > -4 ~ "Loss",
      value ==0 ~ NA_character_,
      value>0 & value < 4  ~ "Gain",
      value >= 4 ~ "Amplification (Copy gain >=4)"
    )) %>% 
    filter(value!=0) %>% 
    mutate(calling=factor(calling,levels = cplevels)) %>% 
    mutate(type=if_else(calling %in% cplevels[1:2],"Amp","Del"))
  
  
  freqdata <- freqdata  %>% 
    count(chrom,start,end,type,calling) %>% 
    mutate(freq=if_else(calling %in% cplevels[1:2],n/totalsamples,-n/totalsamples)) 
  
  # add boundary to line up heatmap 
  #heatscale <- cnvdata %>% select(chrom=chr,startpos,endpos) %>% group_by(chrom) %>% summarise(start=min(startpos),end=max(endpos)) %>% ungroup() %>% mutate(Centro=0)
  # remove cytoband centro
  hg38centro2 <- bed_intersect(hg38centro,freqdata %>% mutate(calling=as.character(calling)),suffix = c('','.y')) %>% select(chrom:Centro)
  
  freqdata <- bed_intersect(freqdata %>% mutate(calling=as.character(calling)),hg38centro,invert = T)
  freqdata <- bind_rows(
    freqdata %>% mutate(Centro=0),
    hg38centro2 %>% mutate(type="Amp",calling="Amplification (Copy gain >=4)",n=0,freq=1e-36),
    hg38centro2 %>% mutate(type="Amp",calling="Gain",n=0,freq=1e-36),
    hg38centro2 %>% mutate(type="Del",calling="Loss",n=0,freq=-1e-36),
    hg38centro2 %>% mutate(type="Del",calling="Homozygous deletion",n=0,freq=-1e-36)
  ) %>% mutate(calling=factor(calling,levels = cplevels)) %>% 
    arrange(chrom,start,end,type,calling) %>% select(-Centro)
  
  # %>% bind_rows(
  #   heatscale %>% mutate(type="Amp",calling="2",n=0,freq=1e-36),
  #   heatscale %>% mutate(type="Amp",calling="1",n=0,freq=1e-36),
  #   heatscale %>% mutate(type="Del",calling="-1",n=0,freq=-1e-36),
  #   heatscale %>% mutate(type="Del",calling="-2",n=0,freq=-1e-36)
  # ) 
  
  
  freqdata <- bind_rows(
    freqdata %>% rename(pos=start) %>% select(-end),
    freqdata %>% rename(pos=end) %>% select(-start)
  ) %>%
    filter(!is.na(calling)) %>% unique()
  #group_by(chrom,pos,type,calling) %>% 
  #arrange(desc(freq)) %>% 
  #slice(1) %>%
  #ungroup()
  
  
  freqdata <- freqdata %>% 
    pivot_wider(id_cols = -c(type,n),names_from = "calling",values_from = "freq") %>% 
    pivot_longer(cols = -c(chrom,pos),names_to = "calling",values_to = "freq") %>%
    mutate(sig=if_else(calling %in% cplevels[3:4],-1e-36,1e-36),freq=if_else(is.na(freq),sig,freq)) %>% 
    select(-sig) %>%
    arrange(chrom,pos,calling)
  
  
  # overlap freqdata and cnvdata
  tmp <- cnvdata %>% group_by(chr) %>% summarise(start=min(startpos,na.rm = TRUE),end=max(endpos,na.rm = TRUE)) %>% select(chrom=chr,start,end)
  freqdata <- freqdata %>% left_join(tmp) %>% filter(pos>=start,pos<=end) %>% select(-start,-end)
  freqdata_loh <- freqdata_loh %>% left_join(tmp) %>% filter(pos>=start,pos<=end) %>% select(-start,-end) %>% mutate(chrom=factor(chrom,levels = chrlevels))
  
  
  freqdata_extra1 <- 
    freqdata %>% group_by(chrom) %>% summarise(pos=min(pos)) %>% 
    left_join(tmp %>% select(chrom,start)) %>% 
    filter(pos!=start) %>% left_join(freqdata) %>% 
    mutate(pos=start) %>% select(-start)
  
  freqdata_extra2 <-  
    freqdata %>% group_by(chrom) %>% summarise(pos=max(pos)) %>% 
    left_join(tmp %>% select(chrom,end)) %>% 
    filter(pos!=end) %>% left_join(freqdata) %>% 
    mutate(pos=end) %>% select(-end)
  
  freqdata <- bind_rows(freqdata,freqdata_extra1,freqdata_extra2) %>% arrange(chrom,pos,calling)
  
  
  # freqdata <- 
  #   freqdata %>% mutate(pos=(start+end)/2) %>%
  #   select(-start,-end) %>%
  #   arrange(chrom,pos,type,calling) 
  
  #freqdata <- freqdata %>% complete(nesting(chrom,pos),calling,fill=list(n=0,freq=0.0001))
  
  #freqdata2 <- freqdata %>% filter(chrom=="1") %>% arrange(chrom,pos,freq,calling)
  #%>% slice(421:1500)
  
  cnvcolor <- rev(c('#2166ac','#92c5de','#f4a582','#b2182b'))
  #names(cnvcolor) <- c('-2','-1','1','2')
  names(cnvcolor) <- cplevels
  
  hg38centro_data <- hg38centro %>% select(chrom,pos=end) %>% group_by(chrom) %>% arrange(pos) %>% slice(1) %>% ungroup() %>% mutate(n=0,freq=0,type="Del",chrom=factor(chrom,levels = chrlevels))
  
  hg38centro_data <- hg38centro_data %>% left_join(tmp) %>% filter(pos>=start,pos<=end) %>% select(-start,-end) %>% mutate(chrom=factor(chrom,levels = chrlevels))
  
  #heatscale <- cnvdata %>% select(chrom=chr,startpos,endpos) %>% group_by(chrom) %>% summarise(start=min(startpos),end=max(endpos)) %>% ungroup() 
  #scales_ <- list(heatscale$chrom=scale_x_continuous(limits = c(heatscale$start,heatscale$end)))
  
  p_freq <- freqdata %>% 
    mutate(chrom=factor(chrom,levels = chrlevels)) %>% 
    ggplot(aes(x = pos,y=freq,fill=factor(calling,levels = cplevels[c(1,2,4,3)])))+
    geom_area(position = "stack",linetype=2)+
    #geom_bar(data=freqdata2 %>% filter(type=="Del"),stat="identity",position = "stack")+
    scale_fill_manual(values = cnvcolor)+
    labs(fill="SCNA Calling")+
    geom_hline(yintercept = 0,color="gray40")+
    #facet_wrap(~chr,nrow = 1,scales = 'free_x',strip.position = "bottom")+
    facet_grid(~chrom,scales = 'free_x',space = 'free')+
    scale_x_continuous(breaks = pretty_breaks(),expand = c(0,0),sec.axis = dup_axis())+
    scale_y_continuous(expand = expand_scale(add = c(0, 0)), sec.axis = dup_axis(),limits = c(-0.8,0.8),breaks = c(-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8),labels = function(x) paste0(x*100, "%"))+
    labs(x="",y="% CNV gain/loss, copy neutral LOH")+
    theme(text = element_text(family = "Roboto Condensed"),panel.spacing = unit(0, "lines"),panel.grid.major.y = element_line(colour = 'gray50',linetype = 5),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.placement = "inside",strip.switch.pad.wrap = unit(0, "cm"),strip.text = element_blank(),strip.background = element_blank(),panel.background = element_rect(fill = 'white'),axis.line= element_line(colour = 'black'),axis.title.y.right = element_blank(),axis.text.y.right = element_blank(),axis.ticks.y.right = element_blank(),axis.title.x.top = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank() )+
    coord_cartesian(clip="off")+panel_border(colour = 'gray80',size = 0.2)+
    geom_line(data=freqdata_loh,aes(pos,freq),colour="black",size=0.4)+
    geom_point(data=hg38centro_data,aes(pos,freq),size=1.5,pch=21,fill="black")
  
  p_freq <- flush_ticks(gg = p_freq)+theme(axis.text.x = element_blank())
  
  p_freq_legend <- get_legend(
    # create some space to the left of the legend
    p_freq + theme(legend.box.margin = margin(8, 0, 5, 0))
  )
  
  p_freq <- p_freq+theme(legend.position = "none")
  
  
  # Genome_Coverage ---------------------------------------------------------
  # cnv_freq <- 
  #   cnvdata %>% 
  #   select(Tumor_Barcode,chr,startpos,endpos,relative_copy) %>% 
  #   filter(relative_copy!=0) %>% 
  #   mutate(size=endpos-startpos) %>% 
  #   group_by(Tumor_Barcode) %>% 
  #   summarise(total=sum(size)) %>% 
  #   mutate(coverage=total/2881033286) %>% 
  #   right_join(cnvclust) %>% 
  #   mutate(coverage=if_else(is.na(coverage),0,coverage)) %>% 
  #   left_join(sample_order) 
  # cnv_freq %>% 
  #   ggplot(aes(Seq,coverage,fill=CNV_Clust))+geom_col()+coord_flip()
  
  #cnv_coverage <- cnv_freq %>% select(Tumor_Barcode,CNV_Coverage=coverage)
  #save(cnv_coverage,file='cnv_coverage.RData')
  
  
  # Combined figures --------------------------------------------------------
  p_dend <- p_dend+theme(plot.margin=margin(r=0,unit="cm"))
  p_heatmap <- p_heatmap+theme(plot.margin=margin(l=-1,unit="cm"))+panel_border(colour = 'gray80',size = 0.2)+theme(axis.text.y = element_blank())
  p_freq <- p_freq+theme(plot.margin=margin(b=-0.4,t=0.3,l=-1,unit="cm"))
  
  p_com <- align_plots(p_freq,
                       p_heatmap,
                       align = 'v',
                       axis = 'lr'
  )
  
  #p_combined_plot <- plot_grid(p_dend,p_com[[2]] ,p_a1,p_a2,p_a3,p_a4,p_a5,align = 'h',axis = "tb",rel_widths = c(0.5,1.5,0.02,0.02,0.02,0.02,0.02),nrow = 1)
  p_combined_plot <- plot_grid(p_dend,p_com[[2]] ,align = 'h',axis = "tb",rel_widths = c(0.5,1.5,0.02,0.02,0.02,0.02,0.02),nrow = 1)
  
  p_combined_plot2 <- plot_grid(NULL,p_com[[1]],NULL,NULL,NULL,NULL,NULL,align = 'h',rel_widths = c(0.5,1.5,0.02,0.02,0.02,0.02,0.02),nrow = 1)
  
  #p_combined_legend <- plot_grid(NULL,NULL,NULL,p_freq_legend,p_heatmap_legend,p_a1_legend,p_a2_legend,p_a3_legend,p_a4_legend,p_a5_legend,NULL,NULL,nrow = 1,align = 'h')
  p_combined_legend <- plot_grid(NULL,NULL,NULL,p_freq_legend,p_heatmap_legend,NULL,NULL,nrow = 1,align = 'h')
  
  p_combined <- plot_grid(p_combined_plot2,p_combined_plot,p_combined_legend,nrow = 3,align = "h",rel_heights = c(2,10,1.2))
  
  # letting users download clustering plot?
  #ggsave(file="heatmap_all_segments_raw_tmp.pdf",plot = p_combined,width = 22,height = 17,device = cairo_pdf)
  #ggsave(file="heatmap_all_clone.pdf",plot = p_combined,width = 22,height = 17,device = cairo_pdf)
  #ggsave(file="heatmap_all_subclone.pdf",plot = p_combined,width = 20,height = 17,device = cairo_pdf)
  
  return(p_combined)

}

sv_recon_plot <- function(sv_data, cn_data, barcode, genome_build, chrs, chrs_start, chrs_end) {
  
  chr_selection <- data.frame(chr=chrs, start=as.numeric(chrs_start), end=as.numeric(chrs_end))
  print(chr_selection)
  str(chr_selection)
  
  p = ReConPlot(sv_data %>% filter(Tumor_Barcode == barcode),
                cn_data %>% filter(Tumor_Barcode == barcode),
                chr_selection=chr_selection,
                size_text = 10,
                size_chr_labels = 10,
                size_title = 12,
                size_sv_line = 0.6,
                size_interchr_SV_tip = 1,
                genome_version = genome_build,
                title=barcode)
  
  #p = p + theme(plot.margin = margin(15, 5, 5, 5))
  #print(p)
  return(p)
}
  
dendro_data_k <- function(hc, k) {
  
  hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
  seg       <-  hcdata$segments
  labclust  <-  cutree(hc, k)[hc$order]
  segclust  <-  rep(0L, nrow(seg))
  heights   <-  sort(hc$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)
  
  for (i in 1:k) {
    xi      <-  hcdata$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }
  
  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  hcdata$segments$clust  <-  segclust
  hcdata$segments$line   <-  as.integer(segclust < 1L)
  hcdata$labels$clust    <-  labclust
  
  hcdata
}

plot_ggdendro <- function(hcdata,
         direction   = c("lr", "rl", "tb", "bt"),
         fan         = FALSE,
         scale.color = NULL,
         branch.size = 1,
         label.size  = 3,
         nudge.label = 0.01,
         expand.y    = 0.1) {
  
  direction <- match.arg(direction) # if fan = FALSE
  ybreaks   <- pretty(segment(hcdata)$y, n = 5)
  ymax      <- max(segment(hcdata)$y)
  
  ## branches
  p <- ggplot() +
    geom_segment(data         =  segment(hcdata),
                 aes(x        =  x,
                     y        =  y,
                     xend     =  xend,
                     yend     =  yend,
                     linetype =  factor(line),
                     colour   =  factor(clust)),
                 lineend      =  "round",
                 show.legend  =  FALSE,
                 size         =  branch.size)
  
  ## orientation
  if (fan) {
    p <- p +
      coord_polar(direction = -1) +
      scale_x_continuous(breaks = NULL,
                         limits = c(0, nrow(label(hcdata)))) +
      scale_y_reverse(breaks = ybreaks)
  } else {
    p <- p + scale_x_continuous(breaks = NULL)
    if (direction %in% c("rl", "lr")) {
      p <- p + coord_flip()
    }
    if (direction %in% c("bt", "lr")) {
      p <- p + scale_y_reverse(breaks = ybreaks,expand=c(0,0))
    } else {
      p <- p + scale_y_continuous(breaks = ybreaks)
      nudge.label <- -(nudge.label)
    }
  }
  
  # labels
  labelParams <- set_labels_params(nrow(hcdata$labels), direction, fan)
  hcdata$labels$angle <- labelParams$angle
  
  # p <- p +
  #   geom_text(data        =  label(hcdata),
  #             aes(x       =  x,
  #                 y       =  y,
  #                 label   =  label,
  #                 colour  =  factor(clust),
  #                 angle   =  angle),
  #             vjust       =  labelParams$vjust,
  #             hjust       =  labelParams$hjust,
  #             nudge_y     =  ymax * nudge.label,
  #             size        =  label.size,
  #             show.legend =  FALSE)
  
  # colors and limits
  if (!is.null(scale.color)) {
    p <- p + scale_color_manual(values = scale.color)
  }
  
  ylim <- -round(ymax * expand.y, 1)
  p    <- p + expand_limits(y = ylim)
  
  p
}

set_labels_params <- function(nbLabels,
         direction = c("tb", "bt", "lr", "rl"),
         fan       = FALSE) {
  if (fan) {
    angle       <-  360 / nbLabels * 1:nbLabels + 90
    idx         <-  angle >= 90 & angle <= 270
    angle[idx]  <-  angle[idx] + 180
    hjust       <-  rep(0, nbLabels)
    hjust[idx]  <-  1
  } else {
    angle       <-  rep(0, nbLabels)
    hjust       <-  0
    if (direction %in% c("tb", "bt")) { angle <- angle + 45 }
    if (direction %in% c("tb", "rl")) { hjust <- 1 }
  }
  list(angle = angle, hjust = hjust, vjust = 0.5)
}
