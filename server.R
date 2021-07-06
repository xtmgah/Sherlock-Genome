# function to read in  column names and create a dropdown menu of them
read_colnames <- function(filename, sep="\t",header=TRUE){
  file <- read.delim(file=filename)
  colnames(file)
}

read_in_file <- function(filename, sep="\t",header=TRUE){
  file <- read.delim(file=filename,sep="\t", header=TRUE)
}

# should we reassign the column types? (numeric and character)
# example: ngs_purity_table$DPClust_Mutations <- as.numeric(ngs_purity_table$DPClust_Mutations)
# typeof(ngs_purity_table$DPClust_Mutations) == typeof(ngs_purity_table$Tumor_Barcode) FALSE (int and chr)
# typeof(ngs_purity_table$WGD_MP2_Ratio) == typeof(ngs_purity_table$Tumor_Barcode ) FALSE (num and chr)
# typeof(ngs_purity_table$WGD_MP2_Ratio) == typeof(ngs_purity_table$nWGD_CCF_MCluster )  FALSE (int and num)
#str(ngs_purity_table)
#'data.frame':	1255 obs. of  31 variables:
#$ Tumor_Barcode                  : chr  "NALC-0001-T01" "NALC-0002-T01" "NALC-0003-T01" "NALC-0004-T01" ...
#$ PGA                            : num  0.17326 0.00351 0.96262 0.00349 0.89282 ...
#$ PGA_Subclonal                  : num  0.03725 0.00213 0 0.00281 0.03377 ...
#$ PGA_TETRA                      : num  0.0349 0 0.7454 0 0.184 ...
#$ PGA_LOH                        : num  0.02447 0.00209 0.08505 0.00325 0.1894 ...
#$ PGA_Haploid_LOH                : num  0.00 1.79e-03 8.12e-02 1.14e-05 1.50e-01 ...
#$ MCN                            : num  0.139976 0.003186 0.962026 0.000263 0.884416 ...
#$ MCN_WGD                        : chr  "nWGD" "nWGD" "WGD" "nWGD" ...
#$ MCN_CHR                        : int  2 0 22 0 21 19 1 19 19 1 ...
#$ DPClust_Mutations              : int  61909 428 15437 221 98709 17368 22470 17102 17761 1651 ...
#$ ASCAT_Purity                   : num  0.23 1 0.13 1 0.13 0.11 0.33 0.27 0.18 0.16 ...
#$ ASCAT_Ploidy                   : num  2.26 2 3.89 2 3.3 ...
#$ BB_Purity                      : num  0.227 0.987 0.13 0.961 0.124 ...
#$ BB_Ploidy                      : num  2.14 2 3.95 2 3.43 ...
#$ CCUBE_Purity                   : num  0.319 0.327 NA 0.82 0.184 ...
#$ GOF                            : num  0.686 0.013 0.9977 0.0627 0.9537 ...
#$ High_Purity                    : chr  "" "BB_Purity|ASCAT_Purity" "" "BB_Purity|ASCAT_Purity" ...
#$ WGD_MP2_Number                 : int  2314 0 0 0 0 352 0 549 0 0 ...
#$ WGD_MP2_Ratio                  : num  0.339 NA NA NA NA ...
#$ WGD_EVEN_Ratio                 : num  4.05e-02 1.88e-03 8.55e-01 2.33e-07 4.12e-01 ...
#$ nWGD_CCF_MCluster              : int  0 0 0 0 0 0 1 0 1 0 ...
#$ nWGD_CCF_MCluster_Detail       : chr  "" "" "" "" ...
#$ nWGD_Substantial_Segment       : int  0 0 0 0 3 2 0 3 1 7 ...
#$ nWGD_Substantial_Segment_Detail: chr  "" "" "" "" ...
#$ Super_Cluster_1                : int  0 0 1 0 0 0 0 0 1 0 ...
#$ Super_Cluster_5                : int  0 0 1 0 0 0 0 0 0 0 ...
#$ Super_Cluster_10               : int  0 0 1 0 0 0 0 0 0 0 ...
#$ Homozoygous_Deletion           : int  0 0 0 0 0 0 0 0 0 0 ...
#$ Homozoygous_Deletion2          : chr  "N" "N" "N" "N" ...
#$ Homozoygous_Deletion_Detail    : chr  "" "" "" "" ...
#$ DPCLUST_Detail                 : chr  "C3|1|61732" "C3|0.11|413" "C1|1.3|15387" "C3|0.402|23; C4|0.089|164" ...
# only pulls data type of "chr" for inspect_cat() and inspect_imb()

# inspect data functions
inspect_data_function <- function(dataframe, type_of_inspection=NULL,column_name=NULL ){
  if(column_name != "All columns"){
    dataframe <- as.data.frame(dataframe[,column_name])
  }
  # inspect_na()
  if(type_of_inspection =="na"){
    x <- inspect_na(dataframe)
    x <- show_plot(x)
  }
  # inspect_cat()
  if(type_of_inspection == "cat"){
    x <- inspect_cat(dataframe)
    x <- show_plot(x)
  } 
  # inspect_imb()
  if(type_of_inspection == "cat_levels"){
    x <- inspect_imb(dataframe)
    x <- show_plot(x)
  }
  # inspect_num()
  if(type_of_inspection == "num"){
    #dataframe <- as.data.frame(as.numeric(unlist(dataframe)))
    x <- inspect_num(dataframe)
    x <- show_plot(x)
  }
  # inspect_types()
  if(type_of_inspection == "types"){
    x <- inspect_types(dataframe)
    x <- show_plot(x)
  }
  return(x)
}

# output figures for ngspurity data
figure_display_ngspurity <- function(tumor_barcode=NULL,battenberg=NULL,type=NULL){
  ngspurity_qc <- read_in_file("ngspurity_qc_file.txt")
  # pull the filename from the ngspurity_qc_file table
  file_path <- ngspurity_qc$File[which(ngspurity_qc$Tumor_Barcode==tumor_barcode & ngspurity_qc$Battenberg==battenberg & ngspurity_qc$Type==type)]
  # absolute path
  # filename <- normalizePath(file.path(file_path))
  # relative path
  filename <- file.path(file_path)
  
  list(src = filename, alt = paste0(tumor_barcode, "_", battenberg, "_", type))
}

# variable check when a dataframe for any module is loaded
# Look into forcats more
# fct_lump(), factor() - more than 3 levels, make least common levels “other”?
# Need to make Tumor_Barcode a factor? All categorical variables?

variable_check <- function(dataframe= NULL){
  # Each variable input will be checked to see what R is assigning the variable type as
  # Then, depending on the number of unique value for each variable, the type of variable could potentially be changed
  # For example, if there are 23 unique values and the values are 0-22, and R is assigning the variable type of 'continuous', it will be changed to categorical
  # ngspurity example:not sure if this is a good way to do this- for example, Homozygous Deletion has 9 unique values, but not sure if it is the count
  # or chromosome number
  #dataframe <- dataframe
  i <- 1
  for(each in colnames(dataframe)){
    uniq_len <- length(unique(dataframe[,colnames(dataframe)[i]]))
    #print(colnames(dataframe)[i])
    #print(uniq_len)
    #print(i)
    # will need to be adjusted depending on number of unique values there can be
    if(uniq_len <5 & typeof(dataframe[,colnames(dataframe)[i]]) %in% c("numeric","integer","double")){
      #print("test")
      #dataframe %>% mutate(dataframe[,colnames(dataframe)[i]] = as.character(dataframe[,colnames(dataframe)[i]]))
      dataframe[,colnames(dataframe)[i]] <- as.character(dataframe[,colnames(dataframe)[i]])
      #data1 <- dataframe
      #print(typeof(dataframe[,colnames(dataframe)[i]]))
    }
    i <- i + 1
  }
 return(dataframe)

}
# association tests (starting with ngspurity data)
association_test <- function(dataframe=NULL, variable_one=NULL, variable_two=NULL) {
# check if the length of unique values in continuous variable is less than 23 - will probably need to be changed in the future  
  #if(typeof(dataframe[,variable_one]) %in% c("integer","double")){ #continuous
  #  if(length(unique(dataframe[,variable_one])) > 23){
      #print(paste0(variable_one, "is a continuous variable."))
  #  }
   # if(length(unique(dataframe[,variable_one])) <= 23){
   #   as.integer(dataframe[,variable_one])
   #   typeof(variable_one)
        #print(typeof(variable_one))
   # } 
  #}

# add a check to see what type each variable is in order to know which loop to enter (and therefore test to run)
  var_one_type <- typeof(dataframe[,variable_one])
  var_two_type <- typeof(dataframe[,variable_two])
  type_table <- table(as.vector(c(var_one_type,var_two_type)))
  num_cat_vars <- sum(type_table[which(rownames(type_table) =="character")])
  num_contin_vars <-sum(type_table[which(rownames(type_table) %in% c("numeric", "integer","double"))])
  
  if(num_cat_vars == 0 & num_contin_vars == 2){ #for two continuous variables
  
    shapiro_one <- shapiro.test(dataframe[,variable_one])
    shapiro_two <- shapiro.test(dataframe[,variable_two])
  
  # make p-values for each variable into a table
    shapiro_table <- data.table(x=c(variable_one,variable_two), y=c(shapiro_one$p.value,shapiro_two$p.value))
    colnames(shapiro_table) <- c("variable","shapiro p-value")
  
  # Histogram to plot distribution
    hist_one <- ggplot(data= dataframe, aes(x=dataframe[,variable_one])) + geom_histogram() 
  #hist_one <- hist(dataframe[,variable_one])
    hist_two <- ggplot(data= dataframe, aes(x=dataframe[,variable_two])) + geom_histogram()
  
  # QQ-Plot to observe distribution
    qqplot_one <- ggqqplot(dataframe[,variable_one])
    qqplot_two <- ggqqplot(dataframe[,variable_two])
  
    # Association Test for 2 continuous variables
    # condition to determine the correct type (parametric (p) or nonparametric (n)) depending on p-value from shapiro test
    if(shapiro_one$p.value || shapiro_two$p.value < 0.05){
      type <- "n"
    }else{
      type <- "p"
    }

    data <- data.frame(x=dataframe[,variable_one], y=dataframe[,variable_two])
    colnames(data) <- c("variable_one", "variable_two")
  
    # add in arguments for generating the plot
    #ggscatterstats(data=data, x= variable_one, y=variable_two, type="n")
    continuous_scatter <- ggscatterstats(data=data, x="variable_one",y= "variable_two", type=type, xlab=colnames(dataframe)[which(colnames(dataframe) == variable_one)],ylab=colnames(dataframe)[which(colnames(dataframe) == variable_two)])
    #continuous_scatter <- ggscatterstats(data=data, x="variable_one",y= "variable_two", type=type, xlab=variable_one, ylab = variable_two)
  
    # consider this as an option for different combinations of multiple continuous variables
    continuous_corr <- ggcorrmat(data=data, matrix.type="full",type=type)
  
    combined_hist_qq <- plot_grid(hist_one,hist_two,qqplot_one,qqplot_two)
    final_plot <- plot_grid(continuous_scatter,continuous_corr)
  }
  # Association Test for one continuous variable and one categorical variable
  # if statement to set the categorical as variable_one and the continuous variable as variable_two
  
  
  # continuous and categorical variable  - test: MCN_WGD (cat) and PGA (cont) - will need to identify which variable is categorical and which is continuous
  # check distribution of the continuous variable
  if(num_cat_vars ==1 & num_contin_vars ==1){
    
    if(typeof(dataframe[,variable_one]) %in% c("integer","double")){
      tmp <- variable_two
      variable_two <- variable_one
      #var_two_type <- typeof(dataframe[,variable_two])
      variable_one <- tmp
      #var_one_type <- typeof(dataframe[,variable_one])
    }else{
      variable_one <- variable_one
      variable_two <- variable_two
    }
    
    data <- data.frame(x=dataframe[,variable_one], y=dataframe[,variable_two])
    colnames(data) <- c("variable_one", "variable_two")  
    
    shapiro_two <- shapiro.test(dataframe[,variable_two])
    
    #make p-values for each variable into a table
    shapiro_table <- data.table(x=variable_two, y=shapiro_two$p.value)
    colnames(shapiro_table) <- c("variable","shapiro p-value")
    
    hist_var_two <- ggplot(data= dataframe, aes(x=dataframe[,variable_two])) + geom_histogram() 
    qqplot_var_two <- ggqqplot(dataframe[,variable_two])
    
    if(shapiro_two$p.value < 0.05){
      type <- "n"
    }else{
      type <- "p"
    }
    
    final_plot <- ggbetweenstats(data=data,x=variable_one,y=variable_two,type=type,xlab=colnames(dataframe)[which(colnames(dataframe) == variable_one)],ylab=colnames(dataframe)[which(colnames(dataframe) == variable_two)])
  }
  
  # two categorical variables
  if(num_cat_vars ==2 & num_contin_vars ==0){
    data <- data.frame(x=dataframe[,variable_one], y=dataframe[,variable_two])
    colnames(data) <- c("variable_one", "variable_two")
    # did not see a different when changing the type
    final_plot <- ggbarstats(data=data,x=variable_one,y=variable_two,xlab=colnames(dataframe)[which(colnames(dataframe) == variable_two)],ylab=colnames(dataframe)[which(colnames(dataframe) == variable_one)], legend.title =colnames(dataframe)[which(colnames(dataframe) == variable_one)] )
  }
  
  #return(list(shapiro_one$p.value,shapiro_two$p.value,combined_hist_qq))
  #return(list(shapiro_table,final_plot))
  return(final_plot)
  #return(list(paste0("Shapiro-Wilk results: ", variable_one, " p-value: ",shapiro_one$p.value, ".") , paste0("Shapiro-Wilk results: ", variable_two, " p-value: ", shapiro_two$p.value, "."),combined_hist_qq))
}


server <- function(input, output, session){
  
# going to probably only want to fun variable_check one time and not for multiple modules/tabs within modules  
  # will need to be use for all dataframes and not only all_ngspurity_output
   inspect_option <- reactive({
        req(input$inspect_data_type)
        req(input$column_name_to_inspect)
        # in future will need to fix as to how to allow for different filename
        dataframe = read_in_file("all_ngspurity_output.txt") %>% variable_check() # somehow change to input$dataframe somewhere?
        type_of_inspection = input$inspect_data_type
        column_name = input$column_name_to_inspect
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
     return(figure_display_ngspurity(tumor_barcode, battenberg,type))
   })
   
  # reactive to call association_test() function above
  association_input <- reactive({
    req(input$variable_choices_1)
    req(input$variable_choices_2)
    dataframe = read_in_file(filename="all_ngspurity_output.txt") %>% variable_check() # somehow change to be th file that corresponds to the current module
    variable_one = input$variable_choices_1
    variable_two = input$variable_choices_2
    return(association_test(dataframe, variable_one, variable_two))
  }) 
   
  observeEvent(input$select,{
    setwd(paste0("./Genomic Data/", input$project_code))
    #output$file_prompt <- renderUI({paste0("Below are the files available for the ", input$project_code, " project:")})
    #output$file_prompt <- renderUI({getwd()})
    #output$files <- renderPrint({list.files(path=paste0("./Genomic Data/",input$project_code))})
    files <- list.files()
    output$file_list_select <- renderUI({disabled(checkboxGroupInput("file_list", label=paste0("Below are the files available for the ", input$project_code, " project:"), choices= files))})
  })
  
  #eventReactive(input$select,{
  #  setwd("../")
  #  setwd(paste0("./Genomic Data/",input$project_code))
  #})
  
  observeEvent(input$choose_files_all,{
    files <- list.files()
    updateCheckboxGroupInput(session,"file_list", label=paste0("Below are the files available for the ", input$project_code, " project:"), choices= files, selected= files)
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
      updateCheckboxGroupInput(session,"file_list", label=paste0("Below are the files available for the ", input$project_code, " project:"), choices= files, selected=character(0))
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
      if("all_ngspurity_output.txt" %in% input$file_list){
        ngs_purity_table <- read_in_file("all_ngspurity_output.txt")
        ngs_purity_header <- read_colnames("all_ngspurity_output.txt")
        # possibly move to the ui side (depends on if it should be available from the start)
        # View Sample Data Tab
        output$ngs_purity_header <- renderUI({dropdownButton(inputId="ngspurity_dropdown",label="Select Columns",circle=FALSE,checkboxGroupInput("ngspurity_header",label="Column Names",choices=c("All columns",ngs_purity_header),selected=c("Tumor_Barcode","PGA","PGA_Subclonal","PGA_TETRA","PGA_LOH","PGA_Haploid_LOH")))})
        output$ngs_purity_header_b <- DT::renderDataTable({datatable(ngs_purity_table[ ,input$ngspurity_header,drop=FALSE], options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
        # Inspect Data Tab
        #output$inspect_df_test <- renderUI({selectInput("inspect_data_type","Select one of the options below to inspect the data:" ,choices=c("cat","cat_levels","na","num","types"), multiple=FALSE, selected="cat")})
        output$inspect_df_test <- renderUI({pickerInput("inspect_data_type","Select one of the options below to inspect the data:" ,choices=c("cat","cat_levels","na","num","types"), multiple=FALSE, selected="cat")})
        output$ngs_purity_header_inspect_tab <- renderUI({selectInput("column_name_to_inspect","Select one column name to inspect, or all columns:", choices= c("All columns",colnames(ngs_purity_table)), multiple= FALSE, selected="All columns")})
        output$testprint <- renderPlot({inspect_option()})
        # View Figures Tab
        ngspurity_qc <- read_in_file("ngspurity_qc_file.txt")
        output$ngspurity_barcode <- renderUI({selectInput("tumor_barcode_to_inspect","Select one Tumour Barcode to inspect:", choices= unique(ngspurity_qc$Tumor_Barcode), multiple= FALSE)})
        output$ngspurity_battenberg <- renderUI({selectInput("battenberg_to_inspect","Select one Battenberg to inspect:", choices= sort(unique(ngspurity_qc$Battenberg[ngspurity_qc$Tumor_Barcode==input$tumor_barcode_to_inspect])), multiple= FALSE)})
        output$ngspurity_type <- renderUI({selectInput("type_to_inspect","Select one Type to inspect:", choices= sort(unique(ngspurity_qc$Type[ngspurity_qc$Tumor_Barcode==input$tumor_barcode_to_inspect & ngspurity_qc$Battenberg == input$battenberg_to_inspect])), multiple= FALSE)})
        #output$ngspurity_qc_battenberg <- renderUI({selectInput("battenberg_to_inspect","Select one Battenberg to inspect:", choices= unique(ngspurity_qc$Battenberg), multiple= FALSE)})
        #output$ngspurity_qc_type <- renderUI({selectInput("type_to_inspect","Select one Type to inspect:", choices= unique(ngspurity_qc$Type), multiple= FALSE)})
        output$figure <- renderImage({figure_output()},deleteFile=FALSE)
        # Association Tab
        # not sure this ngspurity_assoc variable is needed; need to remove those that end with "_Detail" from choices
        ngspurity_assoc <- read_in_file("all_ngspurity_output.txt")
        ngspurity_assoc2 <- colnames(ngspurity_assoc %>% select(!ends_with("_Detail"),-"Tumor_Barcode"))
        output$ngspurity_variable_list_1 <- renderUI({pickerInput("variable_choices_1", "Variable One", choices= ngspurity_assoc2, multiple=FALSE)})
        output$ngspurity_variable_list_2 <- renderUI({pickerInput("variable_choices_2", "Variable Two",choices= ngspurity_assoc2, multiple=FALSE)})
        observeEvent(input$submit1, {
          #association_output <- association_input()
          output$shapiro_test <- renderPrint({association_input()})
          output$scatter_corr <- renderPlot({association_input()})
          
        })
        
      }
    
      if("test.txt" %in% input$file_list){
        test <- readLines("test.txt")
        output$test <- renderUI({test})
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
  
  # update inspect df dropdown to update depending on which column type is chosen
  observeEvent(input$column_name_to_inspect,{
    dataframe <- read_in_file(filename="all_ngspurity_output.txt")
    if(input$column_name_to_inspect != "All columns"){
      #indiv_col <- as.data.frame(dataframe[,input$column_name_to_inspect])
      if(typeof(dataframe[,input$column_name_to_inspect]) =="character"){
        #output$test2 <- renderPrint({typeof(dataframe[,input$column_name_to_inspect])})
        updatePickerInput(session,"inspect_data_type", "Select one of the options below to inspect the data:",choices=c("cat","cat_levels","na","num","types"),choicesOpt = list(disabled=c("cat","cat_levels","na","num","types") %in% c("num")))
      }
      if(typeof(dataframe[,input$column_name_to_inspect]) == "numeric" || typeof(dataframe[,input$column_name_to_inspect]) == "integer" || typeof(dataframe[,input$column_name_to_inspect]) == "double"){
        #output$test2 <- renderPrint({typeof(dataframe[,input$column_name_to_inspect])})
        updatePickerInput(session,"inspect_data_type", "Select one of the options below to inspect the data:",choices=c("cat","cat_levels","na","num","types"),choicesOpt = list(disabled=c("cat","cat_levels","na","num","types") %in% c("cat","cat_levels")))
      }
    }
    
  })
  
  # remove first variable chosen from the second variable list
  observeEvent(input$variable_choices_1,{
    ngspurity_assoc <- read_in_file("all_ngspurity_output.txt")
    updatePickerInput(session, "variable_choices_2","Variable Two",choices=colnames(ngspurity_assoc)[colnames(ngspurity_assoc) != "Tumor_Barcode"], choicesOpt = list(disabled=colnames(ngspurity_assoc)[colnames(ngspurity_assoc) != "Tumor_Barcode"] %in% input$variable_choices_1), options=pickerOptions(hideDisabled=TRUE))
  })
}

# output$oncoplot <- renderPlot({})
  