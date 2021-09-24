
# detect os to apply correct file paths
os_detect <- function(){
  os_name <- Sys.info()[['sysname']]
  return(os_name)
}

# function to read in  column names and create a dropdown menu of them
read_colnames <- function(filename, sep="\t",header=TRUE){
  file <- read.delim(file=filename)
  return(colnames(file))
}

read_in_file <- function(filename, sep="\t",header=TRUE){
  file <- read.delim(file=filename,sep="\t", header=TRUE)
  return(file)
}

change_data_type <- function(data){
  if(is.numeric(data)){
    data <- as.character(data)
  }else if(is.character(data)){
    data <- as.numeric(data)
  }
  return(data)
}

validate_vardf <- function(data, forces=NULL, nachars=c("","na","Na","NA","nan","NAN","Nan"), nacode='NA',Nmin=5, Nmin_drop=FALSE,excludes=NULL, lump=TRUE){
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
  data <- data %>% mutate(across(where(is.character), ~  if_else(is.na(.x), nacode, .x)))
  
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
  if(lump){
    data <- data %>% mutate(across(where(~ is.character(.) & (n_distinct(.) < 0.8*n())),~ fct_lump(fct_infreq(as.factor(.x)),prop = 0.2)))
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

# inspect data function
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

# output figures for ngspurity data
figure_display_ngspurity <- function(tumor_barcode=NULL,battenberg=NULL,type=NULL,project_code=NULL){
  ngspurity_qc <- read_in_file("ngspurity_qc_file.txt")
  # pull the filename from the ngspurity_qc_file table
  file_path <- ngspurity_qc$File[which(ngspurity_qc$Tumor_Barcode==tumor_barcode & ngspurity_qc$Battenberg==battenberg & ngspurity_qc$Type==type)]
  
  # absolute path
  # filename <- normalizePath(file.path(file_path))
  # relative path
  # filename <- file.path(file_path)
  # if(endsWith(filename , "pdf")){
  #   filename <- str_split(filename, "/", simplify=TRUE)[which(endsWith(str_split(filename, "/", simplify=TRUE),"pdf"))]
  # }else{
  #   #filename <- file.path(file_path)
  #   filename <- str_split(filename, "/", simplify=TRUE)[which(endsWith(str_split(filename, "/", simplify=TRUE),"png"))]
  # }
  #tags$iframe(style="height:600px; width:50%", src=filename)
  # if(endsWith(filename, "pdf") ==TRUE){
  #file.copy(from = file_path, to = paste0("../../www/",filename))
  # }else{
  #   print("This file is not a pdf.")
  # }
  
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
  
  #list(src = filename, alt = paste0(tumor_barcode, "_", battenberg, "_", type))
}

# empty_www <- function(){
#   current_path <- rstudioapi::getSourceEditorContext()$path
#   setwd(paste0(dirname(current_path),"/www"))
#   file.remove(list.files())
# }

sherlock_genome_regression <- function(data,formula=NULL) {
  
  data <- validate_vardf(data)
  
  ## for regression module
  supported_types <- c("lm", "glm")
  
  # if(!str_detect(formula,"~")){
  #   stop("Please check your formula for regression, for example, lm( mpg ~ vs + gear")
  # }
  
  str_spl_formula <- str_split(formula, "\\(") 
  formula <- paste0(str_spl_formula[[1]][1],"(formula=",str_spl_formula[[1]][2])
  
  input_formula <- paste0("mod <- data %>% ",formula)
  eval(parse(text=input_formula))
  
  p <- ggstatsplot::ggcoefstats(
    x= mod,
    point.args = list(color = "red", size = 3, shape = 15),
    exclude.intercept = TRUE,
    title = formula,
    ggtheme = hrbrthemes::theme_ipsum_rc()
  ) + # note the order in which the labels are entered
    ggplot2::labs(x = "regression coefficient", y = NULL)
  
  return(p)
}

sherlock_genome_regression_group <- function(data, formula= NULL, Group_Var, type){
  
  data <- validate_vardf(data,excludes = Group_Var)
  
  ## for regression module
  supported_types <- c("lm", "glm")
    
  # if(is.null(formula)|!str_detect(formula,"~")){
  #   stop("Please check your formula for regression, for example, lm( mpg ~ vs + gear")
  # }
    
  colnames(data)[colnames(data) == Group_Var] <- 'Group'
  
  str_spl_formula <- str_split(formula, "\\(") 
  formula <- paste0("(", str_spl_formula[[1]][2])
  
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
  
  return(result)
}


sherlock_genome_association_group <- function(data, Var1, Var2, Group_Var,filter_zero1=FALSE, filter_zero2=FALSE,log1=FALSE,log2=FALSE, type="parametric", collapse_var1=NULL, collapse_var2=NULL) {

    data <- validate_vardf(data,excludes = Group_Var)
    
    ## subset data
    data <- data %>% select(one_of(c(Group_Var,Var1,Var2)))
    colnames(data) <- c("Group","Var1","Var2")
    var1_type <- if_else(is.factor(data[["Var1"]]),"categorical", if_else(is.numeric(data[["Var1"]]),"continuous",NA_character_))
    var2_type <- if_else(is.factor(data[["Var2"]]),"categorical", if_else(is.numeric(data[["Var2"]]),"continuous",NA_character_))

    if(is.na(var1_type)|is.na(var2_type)){
      stop("Please check your data type of these two selected variables")
    }

    # process data or filtering data
    if(filter_zero1 & var1_type == 'continuous') {
      data <- data %>% filter(Var1 != 0)
    }

    if(filter_zero2 & var2_type == 'continuous') {
      data <- data %>% filter(Var2 != 0)
    }

    if(log1 & var1_type == 'continuous') {
      data <- data %>% filter(Var1>0) %>% mutate(Var1 = log2(Var1))
    }

    if(log2 & var2_type == 'continuous') {
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
      # supported types: "parametric", "nonparametric", "robust", "bayes"
      # switch the name if Var2 is categorical
      ## remove unique value
      tmp <- data %>% group_by(Group) %>% summarise(n1=n_distinct(Var1),n2=n_distinct(Var2)) %>% filter(n1==1|n2==1) %>% pull(Group)

      if(var1_type=="categorical"){
        result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::two_sample_test(data = .,x=Var1,y=Var2,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
        result$parameter1 <- Var1
        result$parameter2 <- Var2
        colnames(result)[1:3] <- c(tolower(Group_Var),"variable_name1","variable_name2")
      }else{
        result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::two_sample_test(data = .,x=Var2,y=Var1,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
        result$parameter1 <- Var2
        result$parameter2 <- Var1
        colnames(result)[1:3] <- c(tolower(Group_Var),"variable_name1","variable_name2")
      }
    }

    # transpose the result output to make it fit better and easier to read in shiny app
    result <- t(result)
    colnames(result) <- (1:ncol(result))
    
    return(result)
}

sherlock_genome_association <- function(data, Var1, Var2, filter_zero1=FALSE, filter_zero2=FALSE, log_var1=FALSE, log_var2=FALSE, type="parametric", collapse_var1=NULL, collapse_var2=NULL, xlab=Var1, ylab=Var2, output_plot, file_ext = "png",plot_height, plot_width) {
  
  data <- validate_vardf(data)
  print(Var1)
  print(Var2)
  ## subset data
  data <- data %>% select(one_of(c(Var1,Var2)))
  colnames(data) <- c("Var1","Var2")
  var1_type <- if_else(is.factor(data[[1]]),"categorical", if_else(is.numeric(data[[1]]),"continuous",NA_character_))
  var2_type <- if_else(is.factor(data[[2]]),"categorical", if_else(is.numeric(data[[2]]),"continuous",NA_character_))
  print(var1_type)
  print(var2_type)
  
  print(filter_zero1)
  print(filter_zero2)
  
  print(log_var1)
  print(log_var2)
  
  print(collapse_var1)
  print(collapse_var2)
  
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
      xlab= xlab,
      legend.title = ylab,
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
        xlab= xlab,
        ylab = ylab,
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
  
  if(is.null(output_plot)){
    return(p)
  }else{
    filename <- paste0("output_plot.", file_ext)
    ggsave(filename = filename ,plot = p,width = plot_width,height = plot_height)
    return(p)
  }
}

#----------------------------------------------------------------------------------------


server <- function(input, output, session){
  
# could not get this to work at this time, so split into two different reactives (at least for now)  
  # inspect_option <- reactive({
  #   #req(input$inspect_data_type)
  #   if(input$qc_sample_subtabs == "qc_inspect_sample_data"){
  #     req(input$inspect_data_type_qc_sample)
  #     req(input$column_name_to_inspect_qc_sample)
  #     dataframe = read_in_file("QC_sample_level.txt")
  #     type_of_inspection = input$inspect_data_type_qc_sample
  #     column_name = input$column_name_to_inspect_qc_sample
  #   }
  #   
  #   if(input$ngspurity_tabs == "ngs_inspect_data"){
  #     req(input$inspect_data_type_ngs)
  #     req(input$column_name_to_inspect_ngs)
  #     dataframe = read_in_file("all_ngspurity_output.txt")
  #     type_of_inspection = input$inspect_data_type_ngs
  #     column_name = input$column_name_to_inspect_ngs
  #   }
  #   
  #   #type_of_inspection = input$inspect_data_type
  #   #column_name = input$column_name_to_inspect
  #   #column_name = column_name
  #   return(inspect_data_function(dataframe, type_of_inspection, column_name))
  # })
  
  inspect_option_qc_sample <- reactive({

      req(input$inspect_data_type_qc_sample)
      req(input$column_name_to_inspect_qc_sample)
      dataframe = read_in_file("QC_sample_level.txt")
      type_of_inspection = input$inspect_data_type_qc_sample
      column_name = input$column_name_to_inspect_qc_sample

      return(inspect_data_function(dataframe, type_of_inspection, column_name))
  })

  inspect_option_qc_subject <- reactive({
    
    req(input$inspect_data_type_qc_subject)
    req(input$column_name_to_inspect_qc_subject)
    dataframe = read_in_file("QC_subject_level.txt")
    type_of_inspection = input$inspect_data_type_qc_subject
    column_name = input$column_name_to_inspect_qc_subject
    
    return(inspect_data_function(dataframe, type_of_inspection, column_name))
  })
  
  inspect_option_ngs <- reactive({
    
    req(input$inspect_data_type_ngs)
    req(input$column_name_to_inspect_ngs)
    dataframe = read_in_file("all_ngspurity_output.txt")
    type_of_inspection = input$inspect_data_type_ngs
    column_name = input$column_name_to_inspect_ngs
  
    return(inspect_data_function(dataframe, type_of_inspection, column_name))
  })
  
  # reactive to call figure_display() function above
  figure_output <- reactive({
    
    #ngspurity_qc <- read_in_file("ngspurity_qc_file.txt")
    req(input$tumor_barcode_to_inspect)
    req(input$battenberg_to_inspect)
    req(input$type_to_inspect)
    tumor_barcode = input$tumor_barcode_to_inspect
    battenberg = input$battenberg_to_inspect
    type =  input$type_to_inspect
    project_code = input$project_code
    if(os_detect() %in% c("Linux","Darwin")){
      return(figure_display_ngspurity(tumor_barcode, battenberg,type,project_code))
    }else{
      return(src=figure_display_ngspurity(tumor_barcode, battenberg,type,project_code))
    }
  })
  
  # reactive to call sherlock_genome_regression() function
  regression_inputs <- reactive({
    data = read_in_file("all_ngspurity_output.txt")
    # regression formula
    req(input$ngspurity_regression_formula)
    formula = input$ngspurity_regression_formula
    return(sherlock_genome_regression(data, formula))
  })   
  
  regression_group_inputs <- reactive({
    data = read_in_file("all_ngspurity_output.txt")
    # regression formula
    req(input$ngspurity_regression_formula_group)
    req(input$group_var_input)
    formula = input$ngspurity_regression_formula_group
    Group_Var = input$group_var_input
    type = input$regression_model_group_input #lm or glm
    return(sherlock_genome_regression_group(data, formula, Group_Var, type))
  })
  
  # reactive to call sherlock_genome_association()function
  bivariable_inputs <- reactive({
    # need way to take in data regardless of tab clicked- input$sidebarmenu == menuItem?? input$menuItem?
    #req(input$variable_choices_1)
    #req(input$variable_chocies_2)
    data = read_in_file("all_ngspurity_output.txt")
    Var1 = input$variable_choices_1
    Var2 = input$variable_choices_2
    filter_zero1 = input$filter_ngs_var1
    filter_zero2 = input$filter_ngs_var2
    log_var1 = input$log2_ngs_var1
    log_var2 = input$log2_ngs_var2
    type = input$ngspurity_types_list
    collapse_var1 = input$collapse_ngs_var1
    
    if(input$collapse_ngs_var1==""){
      collapse_var1 = NULL
    }else{
      collapse_var1 = input$collapse_ngs_var1
    }
    
    if(input$collapse_ngs_var2==""){
      collapse_var2 = NULL
    }else{
      collapse_var2 = input$collapse_ngs_var2
    }
    
    xlab = input$variable_choices_1
    ylab= input$variable_choices_2
    output_plot = input$ngs_output_plot
    file_ext = input$file_ext
    plot_width = 8
    plot_height = 8
    
    return(sherlock_genome_association(data, Var1, Var2, filter_zero1, filter_zero2, log_var1, log_var2, type, collapse_var1, collapse_var2, xlab, ylab, output_plot, file_ext, plot_height, plot_width))
    
  })
  
  bivariable_group_inputs <- reactive({
    #req(input$variable_choices_1_group)
    #req(input$variable_chocies_2_group)
    req(input$group_var_input_bivariable)
    data = read_in_file("all_ngspurity_output.txt")
    Var1 = input$variable_choices_1_group
    Var2 = input$variable_choices_2_group
    Group_Var = input$group_var_input_bivariable
    filter_zero1 = input$filter_ngs_var1_group
    filter_zero2 = input$filter_ngs_var2_group
    log_var1 = input$log2_ngs_var1_group
    log_var2 = input$log2_ngs_var2_group
    type = input$ngspurity_types_list_group
    
    if(input$collapse_ngs_var1_group==""){
      collapse_var1 = NULL
    }else{
      collapse_var1 = input$collapse_ngs_var1_group
    }
    
    if(input$collapse_ngs_var2_group==""){
      collapse_var2 = NULL
    }else{
      collapse_var2 = input$collapse_ngs_var2_group
    }
    #collapse_var1 = input$collapse_ngs_var1_group
    #xlab = input$variable_choices_1
    #ylab= input$variable_choices_2
    #output_plot = input$ngs_output_plot
    #file_ext = input$file_ext
    #plot_width = 8
    #plot_height = 8
    
    return(sherlock_genome_association_group(data, Var1, Var2, Group_Var, filter_zero1, filter_zero2, log_var1, log_var2, type, collapse_var1, collapse_var2))
    
  })
  
  output$module_info_table <- DT::renderDataTable({datatable(read_in_file("Module_Info.txt"), options=list(paging=FALSE))}) 
  
  shinyjs::disable("reset_project")
  
  observeEvent(input$select,{
    os_name <- os_detect()
    # if(os_name == "Darwin" | os_name == "Linux"){ #/
    #   setwd(paste0("www/Genomic Data/", input$project_code))
    # }
    # if(os_name == "Windows"){ #\
    #   setwd(paste0("www/Genomic Data/", input$project_code))
    # }
    setwd(paste0("www/Genomic Data/", input$project_code))
    
    shinyjs::enable("reset_project")
    shinyjs::disable("project_code")
    shinyjs::disable("select")
    #output$file_prompt <- renderUI({paste0("Below are the files available for the ", input$project_code, " project:")})
    #output$file_prompt <- renderUI({getwd()})
    #output$files <- renderPrint({list.files(path=paste0("./Genomic Data/",input$project_code))})
    shinyjs::show("choose_files_all")
    shinyjs::show("choose_files_ind")
    shinyjs::show("submit")
    shinyjs::show("reset")
    files <- list.files()
    # only show file list of the files needed to load into modules; add additional file names for other projects as development continues
    #NGSpurity <- c("all_ngspurity_output.txt")
    #sample_qc <- c("QC_sample_level.txt", "QC_subject_level.txt")
    #study_overview <- c("Study_Overview.Rmd")
    module_files <- c("all_ngspurity_output.txt","QC_sample_level.txt", "QC_subject_level.txt","Study_Overview.Rmd","genomePlot_list.txt")
    #module_files <- c(NGSpurity, sample_qc,study_overview)
    file_list_show <- c()
    i <- 1
    for(each in module_files){
      if(module_files[i] %in% files){
        #file_list_show[i] <- append(file_list_show, module_files[i])
        file_list_show[i] <- module_files[i]
        
      }
      i <- i + 1
      file_list_show <- file_list_show[!is.na(file_list_show)]
    }
    output$file_list_select <- renderUI({disabled(checkboxGroupInput("file_list", label=paste0("Below are the files available for the ", input$project_code, " project:"), choices= file_list_show))})
  })
  
  observeEvent(input$reset_project,{
    shinyjs::enable("project_code")
    shinyjs::enable("select")
    shinyjs::disable("reset_project")
    os_name <- os_detect()
    if(os_name == "Darwin" | os_name == "Linux"){ #/
      setwd("../../..")
    }
    if(os_name == "Windows"){ #\
      setwd("../../..")
    }
    shinyjs::hide("file_list")
    shinyjs::hide("choose_files_all")
    shinyjs::hide("choose_files_ind")
    output$file_list_select <- renderUI({updateCheckboxGroupInput(session,"file_list", label=NULL, choices= NULL)})
    shinyjs::hide("submit")
    shinyjs::hide("reset")
    
  })
  
  observeEvent(input$choose_files_all,{
    files <- list.files()
    #NGSpurity <- c("all_ngspurity_output.txt")
    #sample_qc <- c("QC_sample_level.txt", "QC_subject_level.txt")
    #study_overview <- c("Study_Overview.Rmd")
    module_files <- c("all_ngspurity_output.txt","QC_sample_level.txt", "QC_subject_level.txt","Study_Overview.Rmd","genomePlot_list.txt")
    #module_files <- c(NGSpurity, sample_qc,study_overview)
    file_list_show <- c()
    i <- 1
    for(each in module_files){
      if(module_files[i] %in% files){
        #file_list_show[i] <- append(file_list_show, module_files[i])
        file_list_show[i] <- module_files[i]
        
      }
      i <- i + 1
      file_list_show <- file_list_show[!is.na(file_list_show)]
    }
    updateCheckboxGroupInput(session,"file_list", label=paste0("Below are the files available for the ", input$project_code, " project:"), choices= file_list_show, selected= file_list_show)
    delay(3,disable("file_list"))
  })
  
  observeEvent(input$choose_files_ind, {
    shinyjs::enable("file_list")
  })
  
  # if selecting individual files, check that at least one file is selected before enabling the "Submit" button
  
  observe({
    x <- input$file_list
    if(is.null(x)){
      shinyjs::disable("submit")
    }else{
      shinyjs::enable("submit")
    }
    
  })
  
  # reset
  observeEvent(input$reset,{
    files <- list.files()
    #NGSpurity <- c("all_ngspurity_output.txt")
    #sample_qc <- c("QC_sample_level.txt", "QC_subject_level.txt")
    #study_overview <- c("Study_Overview.Rmd")
    module_files <- c("all_ngspurity_output.txt","QC_sample_level.txt", "QC_subject_level.txt","Study_Overview.Rmd","genomePlot_list.txt")
    #module_files <- c(NGSpurity, sample_qc,study_overview)
    file_list_show <- c()
    i <- 1
    for(each in module_files){
      if(module_files[i] %in% files){
        #file_list_show[i] <- append(file_list_show, module_files[i])
        file_list_show[i] <- module_files[i]
        
      }
      i <- i + 1
      file_list_show <- file_list_show[!is.na(file_list_show)]
    }
    updateCheckboxGroupInput(session,"file_list", label=paste0("Below are the files available for the ", input$project_code, " project:"), choices= file_list_show, selected=character(0))
    delay(3,disable("file_list"))
  }
  
  )
  
  # load files that were selected (right now works for both individual and all files)
  observeEvent(input$submit, {
    # right now in Sherlock: all_ngspurity_output.txt, Study_Overview.Rmd, test.txt (alphabetical order)
    if(input$choose_files_ind || input$choose_files_all){
      output$files_selected <- renderText({input$file_list})}
    if("Study_Overview.Rmd" %in% input$file_list){
      study_overview <- includeMarkdown("Study_Overview.Rmd")
      #study_overview <- includeMarkdown(paste0(getwd(), "/", "Study_Overview.Rmd"))
      output$study_overview <- renderUI({study_overview})}
    if("QC_sample_level.txt" %in% input$file_list){
      data_qc_sample <- read_in_file("QC_sample_level.txt")
      #ngs_purity_table <- read_in_file("all_ngspurity_output.txt")
      qc_sample_header <- read_colnames("QC_sample_level.txt") 
      # possibly move to the ui side (depends on if it should be available from the start)
      # View Sample Data Tab
      output$qc_sample_header <- renderUI({dropdownButton(inputId="qc_sample_dropdown",label="Select Columns",circle=FALSE,checkboxGroupInput("qcsample_header",label="Column Names",choices=c("All columns",qc_sample_header),selected=c("Study","Subject","Tumor_Barcode","Normal_Barcode","Gender","Tumor_File")))})
      output$qc_sample_table <- DT::renderDataTable({datatable(data_qc_sample[ ,input$qcsample_header,drop=FALSE], options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
    
      # Inspect Sample Data Tab
      output$inspect_df_qc_sample <- renderUI({pickerInput("inspect_data_type_qc_sample","Select one of the options below to inspect the data:" ,choices=c("cat","cat_levels","na","num","types"), multiple=FALSE, selected="cat")})
      output$inspect_sample_qc_select_column <- renderUI({selectInput("column_name_to_inspect_qc_sample","Select one column name to inspect, or all columns:", choices= c("All columns",colnames(data_qc_sample)), multiple= FALSE, selected="All columns")})
      output$inspect_sample_plot <- renderPlot({inspect_option_qc_sample()})
      }
    
    if("QC_subject_level.txt" %in% input$file_list){
      data_qc_subject <- read_in_file("QC_subject_level.txt")
      #ngs_purity_table <- read_in_file("all_ngspurity_output.txt")
      qc_subject_header <- read_colnames("QC_subject_level.txt") 
      # possibly move to the ui side (depends on if it should be available from the start)
      # View Sample Data Tab
      output$qc_subject_header <- renderUI({dropdownButton(inputId="qc_subject_dropdown",label="Select Columns",circle=FALSE,checkboxGroupInput("qcsubject_header",label="Column Names",choices=c("All columns",qc_subject_header),selected=c("Study","Subject","Barcode","Gender","Sample_ID","File","Wave")))})
      output$qc_subject_table <- DT::renderDataTable({datatable(data_qc_subject[ ,input$qcsubject_header,drop=FALSE], options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
      
      #Inspect Subject Data Tab
      output$inspect_df_qc_subject <- renderUI({pickerInput("inspect_data_type_qc_subject","Select one of the options below to inspect the data:" ,choices=c("cat","cat_levels","na","num","types"), multiple=FALSE, selected="cat")})
      output$inspect_subject_qc_select_column <- renderUI({selectInput("column_name_to_inspect_qc_subject","Select one column name to inspect, or all columns:", choices= c("All columns",colnames(data_qc_subject)), multiple= FALSE, selected="All columns")})
      output$inspect_subject_plot <- renderPlot({inspect_option_qc_subject()})
    }
    if("all_ngspurity_output.txt" %in% input$file_list){
        #observeEvent(input$sidebarmenu, {
      data_ngs <- read_in_file("all_ngspurity_output.txt")
      #ngs_purity_table <- read_in_file("all_ngspurity_output.txt")
      ngs_purity_header <- read_colnames("all_ngspurity_output.txt") 
      # possibly move to the ui side (depends on if it should be available from the start)
      # View Sample Data Tab
      output$ngs_purity_header <- renderUI({dropdownButton(inputId="ngspurity_dropdown",label="Select Columns",circle=FALSE,checkboxGroupInput("ngspurity_header",label="Column Names",choices=c("All columns",ngs_purity_header),selected=c("Tumor_Barcode","PGA","PGA_Subclonal","PGA_TETRA","PGA_LOH","PGA_Haploid_LOH")))})
      output$ngs_purity_table <- DT::renderDataTable({datatable(data_ngs[ ,input$ngspurity_header,drop=FALSE], options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
      # Inspect Data Tab
      #output$inspect_df_test <- renderUI({selectInput("inspect_data_type","Select one of the options below to inspect the data:" ,choices=c("cat","cat_levels","na","num","types"), multiple=FALSE, selected="cat")})
      output$inspect_df_test <- renderUI({pickerInput("inspect_data_type_ngs","Select one of the options below to inspect the data:" ,choices=c("cat","cat_levels","na","num","types"), multiple=FALSE, selected="cat")})
      output$ngs_purity_header_inspect_tab <- renderUI({selectInput("column_name_to_inspect_ngs","Select one column name to inspect, or all columns:", choices= c("All columns",colnames(data_ngs)), multiple= FALSE, selected="All columns")})
      output$inspect_ngs_plot <- renderPlot({inspect_option_ngs()})
      #}  
      #})
      
      # View Figures Tab
      os_name <- os_detect()
      if(os_name == "Windows"){
        ngspurity_qc <- read_in_file("ngspurity_qc_file.txt")
        #ngspurity_qc$File<- gsub("/","\\",ngspurity_qc$File)
        i <- 1 
        for(each in ngspurity_qc$File){
          ngspurity_qc$File[i] <- paste0("C:/", ngspurity_qc$File[i])
          i <- i + 1
        }
      }else{
        ngspurity_qc <- read_in_file("ngspurity_qc_file.txt")
      }
      output$ngspurity_barcode <- renderUI({selectInput("tumor_barcode_to_inspect","Select one Tumor Barcode to inspect:", choices= unique(ngspurity_qc$Tumor_Barcode), multiple= FALSE)})
      output$ngspurity_battenberg <- renderUI({selectInput("battenberg_to_inspect","Select one Battenberg to inspect:", choices= unique(ngspurity_qc$Battenberg[ngspurity_qc$Tumor_Barcode==input$tumor_barcode_to_inspect]), multiple= FALSE)})
      output$ngspurity_type <- renderUI({selectInput("type_to_inspect","Select one Type to inspect:", choices= unique(ngspurity_qc$Type[ngspurity_qc$Tumor_Barcode==input$tumor_barcode_to_inspect & ngspurity_qc$Battenberg == input$battenberg_to_inspect]), multiple= FALSE)})
      
      # observe({
      #   req(input$ngspurity_tabs)
      #   if(input$ngspurity_tabs== "ngs_view_figures"){
      #     setwd("./NGSpurity")
      #   }
      # })
      #output$ngspurity_qc_battenberg <- renderUI({selectInput("battenberg_to_inspect","Select one Battenberg to inspect:", choices= unique(ngspurity_qc$Battenberg), multiple= FALSE)})
      #output$ngspurity_qc_type <- renderUI({selectInput("type_to_inspect","Select one Type to inspect:", choices= unique(ngspurity_qc$Type), multiple= FALSE)})
      # output$figure <- renderImage({figure_output()},deleteFile=FALSE)
      if(os_detect() %in% c("Linux","Darwin")){
        output$figure_pdf <-  renderImage({figure_output()},deleteFile=FALSE)
      }else{
        output$figure_pdf <- renderUI({ tags$iframe(style="height:1000px; width:60%", src= figure_output())})
      }
      
      # Regression/Bivariable Association Tab
      # not sure this ngspurity_assoc variable is needed; need to remove those that end with "_Detail" from choices
      ngspurity_assoc <- read_in_file("all_ngspurity_output.txt") 
      ngspurity_assoc2 <- colnames(ngspurity_assoc %>% select(!ends_with("_Detail"),-"Tumor_Barcode"))
      ngspurity_assoc3 <- ngspurity_assoc %>% validate_vardf() %>% inspectdf::inspect_types()
      factor_var <- data.frame(ngspurity_assoc3$col_name[["factor"]],rep_len("factor",length.out= length(ngspurity_assoc3$col_name[["factor"]]))) 
      colnames(factor_var) <- c("Variable Name", "Variable Type")
      numeric_var <- data.frame(ngspurity_assoc3$col_name[["numeric"]],rep_len("numeric",length.out= length(ngspurity_assoc3$col_name[["numeric"]])))   
      colnames(numeric_var) <- c("Variable Name", "Variable Type")
      character_var <- data.frame(ngspurity_assoc3$col_name[["character"]],rep_len("character",length.out= length(ngspurity_assoc3$col_name[["character"]])))
      colnames(character_var) <- c("Variable Name", "Variable Type")
      variable_type_table <- rbind(factor_var,character_var,numeric_var)
      output$ngs_variable_table <- DT::renderDataTable({variable_type_table}) #options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE)})
      output$ngspurity_variable_list_1 <- renderUI({pickerInput("variable_choices_1", "Variable One", choices= ngspurity_assoc2, multiple=FALSE)})
      output$ngspurity_variable_list_1_group <- renderUI({pickerInput("variable_choices_1_group", "Variable One", choices= ngspurity_assoc2, multiple=FALSE)})
      output$ngspurity_variable_list_2 <- renderUI({pickerInput("variable_choices_2", "Variable Two",choices= ngspurity_assoc2, multiple=FALSE)})
      output$ngspurity_variable_list_2_group <- renderUI({pickerInput("variable_choices_2_group", "Variable Two",choices= ngspurity_assoc2, multiple=FALSE)})
      #output$ylab <- renderUI({input$variable_choices_2})
      
      output$group_var_text_box_bivariable <- renderUI({textInput(inputId= "group_var_input_bivariable", label= NULL,placeholder= "Group Variable Name")})
      output$bivariable_grp_result <- DT::renderDataTable({
        req(input$calculate_ngs_association_group)
        isolate(bivariable_group_inputs())})
      
      output$ngspurity_bivariable_plot <- renderPlot({
        req(input$calculate_ngs_association)
        isolate(bivariable_inputs())})
      
    # Regression Group

      output$group_var_text_box <- renderUI({textInput(inputId= "group_var_input", label= NULL,placeholder= "Group Variable Name" )})
      output$regression_model_group <- renderUI({checkboxGroupInput(inputId= "regression_model_group_input", label= "Select a supported regression model type", choices= c("lm", "glm"))})
      output$formula_input_group <- renderUI({textInput("ngspurity_regression_formula_group", "Regression Formula", placeholder="lm( mpg ~ vs + gear)")})
      output$regression_grp_result <- DT::renderDataTable({
        input$calculate_ngs_regression_group
        isolate(regression_group_inputs())})
        #output$test <- renderText({"test"})
      # }else{
        # shinyjs::hide("reset_ngs_regression")
      output$formula_input <- renderUI({textInput("ngspurity_regression_formula", "Regression Formula", placeholder="lm( mpg ~ vs + gear)")})
        #output$calculate_regression <- renderUI({actionButton("calculate_ngs_regression","Calculate")})
        # shinyjs::hide("group_var_input")
        # shinyjs::hide("regression_model_group_input")
      output$ngspurity_regression_plot <- renderPlot({
        input$calculate_ngs_regression
        isolate(regression_inputs())})
    
   
    if("genomePlot_list.txt" %in% input$file_list){
      data_genomePlot <- read_in_file("genomePlot_list.txt")    
      output$genomePlot_Barcodes <- renderUI({pickerInput("tumor_barcode_to_inspect_genomePlot","Select one Tumor Barcode to inspect:", choices= unique(data_genomePlot$Tumor_Barcode), multiple= FALSE, options= pickerOptions(liveSearch=TRUE))})
      tumor_barcode_to_inspect_genomePlot_reactive <- reactive({
        req(input$tumor_barcode_to_inspect_genomePlot)
        tumor_barcode_to_inspect_genomePlot = input$tumor_barcode_to_inspect_genomePlot
        genomePlot_figurePath <- data_genomePlot$genomePlot[which(data_genomePlot["Tumor_Barcode"] == tumor_barcode_to_inspect_genomePlot)]
      })
      #genomePlot_figurePath <- data_genomePlot$genomePlot[which(data_genomePlot["Tumor_Barcode"] == tumor_barcode_to_inspect_genomePlot())]
      output$genomePlot_figure <- renderImage({list(src=tumor_barcode_to_inspect_genomePlot_reactive(),width="800px", height = "auto", alt=paste0("Tumor Barcode: ",input$tumor_barcode_to_inspect_genomePlot))}, deleteFile = FALSE)
    }
      
  }
    
})
  
  # edit the option for displaying all columns, not finished yet
  #Warning in if (input$ngspurity_header == "All columns") { :
  # the condition has length > 1 and only the first element will be used
  #Warning: Error in [.data.frame: undefined columns selected
  observeEvent(input$ngspurity_header,{
    if("All columns" %in% input$ngspurity_header){
      ngs_purity_header <- read_colnames("all_ngspurity_output.txt")
      updateCheckboxGroupInput(session,"ngspurity_header", label="Column names", choices= c("All columns",ngs_purity_header), selected= c("All_columns", ngs_purity_header))
    }
  })

  # # update inspect df dropdown to update depending on which column type is chosen
   observeEvent(input$column_name_to_inspect_qc_sample,{
     data <- read_in_file("QC_sample_level.txt")
     if(input$column_name_to_inspect_qc_sample != "All columns"){
       #indiv_col <- as.data.frame(dataframe[,input$column_name_to_inspect])
       if(typeof(data[,input$column_name_to_inspect_qc_sample]) =="character"){
         #output$test2 <- renderPrint({typeof(dataframe[,input$column_name_to_inspect])})
         updatePickerInput(session,"inspect_data_type_qc_sample", "Select one of the options below to inspect the data:",choices=c("cat","cat_levels","na","num","types"),choicesOpt = list(disabled=c("cat","cat_levels","na","num","types") %in% c("num")),options=pickerOptions(hideDisabled=TRUE))
       }
       if(typeof(data[,input$column_name_to_inspect_qc_sample]) == "numeric" || typeof(data[,input$column_name_to_inspect_qc_sample]) == "integer" || typeof(data[,input$column_name_to_inspect_qc_sample]) == "double"){
         #output$test2 <- renderPrint({typeof(dataframe[,input$column_name_to_inspect])})
         updatePickerInput(session,"inspect_data_type_qc_sample", "Select one of the options below to inspect the data:",choices=c("cat","cat_levels","na","num","types"),choicesOpt = list(disabled=c("cat","cat_levels","na","num","types") %in% c("cat","cat_levels")),options=pickerOptions(hideDisabled=TRUE))
       }
     }

    })

   observeEvent(input$column_name_to_inspect_qc_subject,{
     data <- read_in_file("QC_subject_level.txt")
     if(input$column_name_to_inspect_qc_subject != "All columns"){
       #indiv_col <- as.data.frame(dataframe[,input$column_name_to_inspect])
       if(typeof(data[,input$column_name_to_inspect_qc_subject]) =="character"){
         #output$test2 <- renderPrint({typeof(dataframe[,input$column_name_to_inspect])})
         updatePickerInput(session,"inspect_data_type_qc_subject", "Select one of the options below to inspect the data:",choices=c("cat","cat_levels","na","num","types"),choicesOpt = list(disabled=c("cat","cat_levels","na","num","types") %in% c("num")),options=pickerOptions(hideDisabled=TRUE))
       }
       if(typeof(data[,input$column_name_to_inspect_qc_subject]) == "numeric" || typeof(data[,input$column_name_to_inspect_qc_subject]) == "integer" || typeof(data[,input$column_name_to_inspect_qc_subject]) == "double"){
         #output$test2 <- renderPrint({typeof(dataframe[,input$column_name_to_inspect])})
         updatePickerInput(session,"inspect_data_type_qc_subject", "Select one of the options below to inspect the data:",choices=c("cat","cat_levels","na","num","types"),choicesOpt = list(disabled=c("cat","cat_levels","na","num","types") %in% c("cat","cat_levels")),options=pickerOptions(hideDisabled=TRUE))
       }
     }
     
   })
   
   observeEvent(input$column_name_to_inspect_ngs,{
     data <- read_in_file("all_ngspurity_output.txt")
     if(input$column_name_to_inspect_ngs != "All columns"){
       #indiv_col <- as.data.frame(dataframe[,input$column_name_to_inspect])
       if(typeof(data[,input$column_name_to_inspect_ngs]) =="character"){
         #output$test2 <- renderPrint({typeof(dataframe[,input$column_name_to_inspect])})
         updatePickerInput(session,"inspect_data_type_ngs", "Select one of the options below to inspect the data:",choices=c("cat","cat_levels","na","num","types"),choicesOpt = list(disabled=c("cat","cat_levels","na","num","types") %in% c("num")), options=pickerOptions(hideDisabled=TRUE))
       }
       if(typeof(data[,input$column_name_to_inspect_ngs]) == "numeric" || typeof(data[,input$column_name_to_inspect_ngs]) == "integer" || typeof(data[,input$column_name_to_inspect_ngs]) == "double"){
         #output$test2 <- renderPrint({typeof(dataframe[,input$column_name_to_inspect])})
         updatePickerInput(session,"inspect_data_type_ngs", "Select one of the options below to inspect the data:",choices=c("cat","cat_levels","na","num","types"),choicesOpt = list(disabled=c("cat","cat_levels","na","num","types") %in% c("cat","cat_levels")), options=pickerOptions(hideDisabled=TRUE))
       }
     }
     
   })
  
  # remove first variable chosen from the second variable list
  observeEvent(input$variable_choices_1,{
    ngspurity_assoc <- read_in_file("all_ngspurity_output.txt")
    updatePickerInput(session, "variable_choices_2","Variable Two",choices=colnames(ngspurity_assoc)[colnames(ngspurity_assoc) != "Tumor_Barcode"], choicesOpt = list(disabled=colnames(ngspurity_assoc)[colnames(ngspurity_assoc) != "Tumor_Barcode"] %in% input$variable_choices_1), options=pickerOptions(hideDisabled=TRUE))
  })
  
  observeEvent(input$variable_choices_1_group,{
    ngspurity_assoc <- read_in_file("all_ngspurity_output.txt")
    updatePickerInput(session, "variable_choices_2_group","Variable Two",choices=colnames(ngspurity_assoc)[colnames(ngspurity_assoc) != "Tumor_Barcode"], choicesOpt = list(disabled=colnames(ngspurity_assoc)[colnames(ngspurity_assoc) != "Tumor_Barcode"] %in% input$variable_choices_1_group), options=pickerOptions(hideDisabled=TRUE))
  })
  
  # Other Association Tests - figure out how to access labels on the variables (continuous and categorical, change to these labels in if statements)
  observe( {
    req(input$variable_choices_1, input$variable_choices_2)
    #print("testing 1 2 3")
    data <- read_in_file(filename="all_ngspurity_output.txt")
    data <- validate_vardf(data)
    if(is.numeric(data[[input$variable_choices_1]]) && is.numeric(data[[input$variable_choices_2]])) {
      updatePickerInput(session,"ngspurity_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes", "skit"), options=pickerOptions(hideDisabled=TRUE))
    }
    if(is.factor(data[[input$variable_choices_1]]) && is.factor(data[[input$variable_choices_2]])) {
      updatePickerInput(session,"ngspurity_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes", "fisher"), options=pickerOptions(hideDisabled=TRUE))
    }
    if(is.numeric(data[[input$variable_choices_1]]) && is.factor(data[[input$variable_choices_2]])) {
      updatePickerInput(session,"ngspurity_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes"), options=pickerOptions(hideDisabled=TRUE))
    }
    if(is.factor(data[[input$variable_choices_1]]) && is.numeric(data[[input$variable_choices_2]])) {
      updatePickerInput(session,"ngspurity_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes"), options=pickerOptions(hideDisabled=TRUE))
    }
  })
  
  observe( {
    req(input$variable_choices_1_group, input$variable_choices_2_group)
    data <- read_in_file(filename="all_ngspurity_output.txt")
    data <- validate_vardf(data)
    if(is.numeric(data[[input$variable_choices_1_group]]) && is.numeric(data[[input$variable_choices_2_group]])) {
      updatePickerInput(session,"ngspurity_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes", "skit"), options=pickerOptions(hideDisabled=TRUE))
    }
    if(is.factor(data[[input$variable_choices_1_group]]) && is.factor(data[[input$variable_choices_2_group]])) {
      updatePickerInput(session,"ngspurity_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes", "fisher"), options=pickerOptions(hideDisabled=TRUE))
    }
    if(is.numeric(data[[input$variable_choices_1_group]]) && is.factor(data[[input$variable_choices_2_group]])) {
      updatePickerInput(session,"ngspurity_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes"), options=pickerOptions(hideDisabled=TRUE))
    }
    if(is.factor(data[[input$variable_choices_1_group]]) && is.numeric(data[[input$variable_choices_2_group]])) {
      updatePickerInput(session,"ngspurity_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes"), options=pickerOptions(hideDisabled=TRUE))
    }
  })
  
  observe({
    req(input$variable_choices_1, input$variable_choices_2)
    data <- read_in_file(filename="all_ngspurity_output.txt")
    data <- validate_vardf(data)
    if(is.numeric(data[[input$variable_choices_1]])){
      shinyjs::disable("collapse_ngs_var1")
    }else{
      shinyjs::enable("collapse_ngs_var1")
    }
    if(is.numeric(data[[input$variable_choices_2]])){
      shinyjs::disable("collapse_ngs_var2")
    }else{
      shinyjs::enable("collapse_ngs_var2")
    }
  })
  
  observe({
    req(input$variable_choices_1_group, input$variable_choices_2_group)
    data <- read_in_file(filename="all_ngspurity_output.txt")
    data <- validate_vardf(data)
    if(is.numeric(data[[input$variable_choices_1_group]])){
      shinyjs::disable("collapse_ngs_var1_group")
    }else{
      shinyjs::enable("collapse_ngs_var1_group")
    }
    if(is.numeric(data[[input$variable_choices_2_group]])){
      shinyjs::disable("collapse_ngs_var2_group")
    }else{
      shinyjs::enable("collapse_ngs_var2_group")
    }
  })
  
  
  observe({
    if(input$ngs_output_plot ==TRUE){
      shinyjs:: show("file_ext")
    }else{
      shinyjs:: hide("file_ext")
    }
  })
  
  observeEvent(input$reset_ngs_regression, {
    reset("ngspurity_regression_formula")
  })
  
  observeEvent(input$reset_ngs_regression_group,{
    reset("ngspurity_regression_formula_group")
    reset("group_var_input")
    updateCheckboxGroupInput(inputId= "regression_model_group_input", label= "Select a supported regression model type", choices= c("lm", "glm"))
  })
  
  observeEvent(input$reset_ngs_association,{
    ngspurity_assoc <- read_in_file("all_ngspurity_output.txt")
    ngspurity_assoc2 <- colnames(ngspurity_assoc %>% select(!ends_with("_Detail"),-"Tumor_Barcode"))
    updatePickerInput(session, inputId="variable_choices_1","Variable One", choices= ngspurity_assoc2)
    updatePickerInput(session, inputId= "variable_choices_2", "Variable Two", choices = ngspurity_assoc2)
    updatePickerInput(session, "variable_choices_2","Variable Two",choices=colnames(ngspurity_assoc)[colnames(ngspurity_assoc) != "Tumor_Barcode"], choicesOpt = list(disabled=colnames(ngspurity_assoc)[colnames(ngspurity_assoc) != "Tumor_Barcode"] %in% input$variable_choices_1), options=pickerOptions(hideDisabled=TRUE))
    updateCheckboxInput(session, "filter_ngs_var1","Filter (>0)",value=FALSE)
    updateCheckboxInput(session,"log2_ngs_var1", "log2", value=FALSE)
    updateTextInput(session, "collapse_ngs_var1", "Collapse Level",value=NULL)
    updateCheckboxInput(session, "filter_ngs_var2","Filter (>0)",value=FALSE)
    updateCheckboxInput(session,"log2_ngs_var2", "log2", value=FALSE)
    updateTextInput(session,"collapse_ngs_var2", "Collapse Level",value=NULL)
    updatePickerInput(session,"ngspurity_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes","fisher","skit"))
    updateCheckboxInput(session,"ngs_output_plot", "Save plot to computer",value = FALSE)
    updateTextInput(session,"file_ext", "File Extension (png, svg, or jpg)", value="png")
   })
  
  observeEvent(input$reset_ngs_association_group,{
    ngspurity_assoc <- read_in_file("all_ngspurity_output.txt")
    ngspurity_assoc2 <- colnames(ngspurity_assoc %>% select(!ends_with("_Detail"),-"Tumor_Barcode"))
    updatePickerInput(session, inputId="variable_choices_1_group","Variable One", choices= ngspurity_assoc2)
    updatePickerInput(session, inputId= "variable_choices_2_group", "Variable Two", choices = ngspurity_assoc2)
    updatePickerInput(session, "variable_choices_2_group","Variable Two",choices=colnames(ngspurity_assoc)[colnames(ngspurity_assoc) != "Tumor_Barcode"], choicesOpt = list(disabled=colnames(ngspurity_assoc)[colnames(ngspurity_assoc) != "Tumor_Barcode"] %in% input$variable_choices_1_group), options=pickerOptions(hideDisabled=TRUE))
    updateCheckboxInput(session, "filter_ngs_var1_group","Filter (>0)",value=FALSE)
    updateCheckboxInput(session,"log2_ngs_var1_group", "log2", value=FALSE)
    updateTextInput(session, "collapse_ngs_var1_group", "Collapse Level",value=NULL)
    updateCheckboxInput(session, "filter_ngs_var2_group","Filter (>0)",value=FALSE)
    updateCheckboxInput(session,"log2_ngs_var2_group", "log2", value=FALSE)
    updateTextInput(session,"collapse_ngs_var2_group", "Collapse Level",value=NULL)
    updatePickerInput(session,"ngspurity_types_list_group", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes","fisher","skit"))
    #updateTextInput(session,inputId= "group_var_input_bivariable", label= NULL,placeholder= "Group Variable Name",value= NULL)
    reset("group_var_input_bivariable")
  })
   
  
}

# output$oncoplot <- renderPlot({})
