# function to read in  column names and create a dropdown menu of them
read_colnames <- function(filename, sep="\t",header=TRUE){
  file <- read.delim(file=filename)
  colnames(file)
}

read_in_file <- function(filename, sep="\t",header=TRUE){
  file <- read.delim(file=filename,sep="\t", header=TRUE)
}

# call inspectdf function (na, cat, etc)
inspect_data <- function(filename,sep="\t",header=TRUE){
 dataframe <- read.delim(file=filename,sep=sep,header=header)
 na_result <- inspect_na(dataframe)
 #cat_result <- inspect_cat(dataframe)
 #final_result <- cat(na_result, cat_result)
 return(na_result=na_result)
}

# should we reassign the column types? (numeric and character)
# example: ngs_purity_table$DPClust_Mutations <- as.numeric(ngs_purity_table$DPClust_Mutations)
# typeof(ngs_purity_table$DPClust_Mutations) == typeof(ngs_purity_table$Tumor_Barcode) FALSE (int and chr)
# typeof(ngs_purity_table$WGD_MP2_Ratio) == typeof(ngs_purity_table$Tumor_Barcode ) FALSE (num and chr)
# typeof(ngs_purity_table$WGD_MP2_Ratio) == typeof(ngs_purity_table$nWGD_CCF_MCluster )  FALSE (int and num)
#str(ngs_purity_table)
#'data.frame':	1255 obs. of  31 variables:
#  $ Tumor_Barcode                  : chr  "NALC-0001-T01" "NALC-0002-T01" "NALC-0003-T01" "NALC-0004-T01" ...
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

# inspect data functions
inspect_data_function <- function(dataframe, type_of_inspection=NULL,column_name=NULL ){
  #dataframe <- read_in_file(filename="all_ngspurity_output.txt")
  #ngs_purity_table <- read_in_file(filename="all_ngspurity_output.txt")
  # possible arguments for type_of_inspection: na, cat, cat_levels, num, types; default is NULL
  # if statements depending on the option chosen by the user (eventually for any dataframe- test is ngspurity)
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

server <- function(input, output, session){
  
   inspect_option <- reactive({
        dataframe = read_in_file(filename="all_ngspurity_output.txt")
        type_of_inspection = input$inspect_data_type
        column_name = input$column_name_to_inspect
        return(inspect_data_function(dataframe, type_of_inspection))
  })
  
  observeEvent(input$select,{
    setwd(paste0("./Genomic Data/",input$project_code))
    #output$file_prompt <- renderUI({paste0("Below are the files available for the ", input$project_code, " project:")})
    #output$file_prompt <- renderUI({getwd()})
    #output$files <- renderPrint({list.files(path=paste0("./Genomic Data/",input$project_code))})
    files <- list.files()
    output$file_list_select <- renderUI({disabled(checkboxGroupInput("file_list", label=paste0("Below are the files available for the ", input$project_code, " project:"), choices= files))})
  })
  
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
      #files <- list.files()
      updateCheckboxGroupInput(session,"file_list", label=paste0("Below are the files available for the ", input$project_code, " project:"), choices= files, selected=character(0))
      delay(3,disable("file_list"))
    }
    
  )
  
  # load files depending on if user selects "Select All Files" or "Select Individual Files"
 # observeEvent(input$submit, {
 #    if(input$choose_files_all){
#      files <- list.files()
#      i <- 1
#     for(each in 1:length(files))
  #paste0("./Genomic Data/",input$project_code)
 #   }else{
      
  #  }
#  })
  

  # load files that were selected (right now works for both individual and all files)
  # dropdown for column headers? and dropdown for the different inspect options?
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
        #ngs_purity_table <- read.delim(file="all_ngspurity_output.txt",sep="\t", header=TRUE)
        #output$ngs_purity_inspect <- renderTable({inspect_na(ngs_purity_table)})
        # add in selection for range of rows to display
        #output$na <- renderTable({head(inspect_data)})
        ngs_purity_header <- read_colnames("all_ngspurity_output.txt")
        # possibly move to the ui side (depends on if it should be available from the start)
        output$ngs_purity_header <- renderUI({selectInput("ngspurity_header","Select at least one column name:", choices= ngs_purity_header, multiple= TRUE)})
        #output$ngs_rows_to_display <- renderUI({textInput("ngs_rows_to_display","Select the number or rows you would like to view (must be between 1 and 100):", value=10)})
        #output$ngs_rows_to_display_validate <- validate_row_input(input$ngs_rows_to_display)
        output$ngs_purity_header_b <- DT::renderDataTable({datatable(ngs_purity_table[ ,input$ngspurity_header,drop=FALSE], options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
        output$inspect_df_test <- renderUI({selectInput("inspect_data_type","Select one of the options below to inspect the data:" ,choices=c("cat","cat_levels","na","num","types"), multiple=FALSE, selected="na")})
        #output$inspect_df_test <- renderPrint({inspect_data(filename="all_ngspurity_output.txt")})
        output$testprint <- renderPlot(inspect_option())
      }
      if("test.txt" %in% input$file_list){
        test <- readLines("test.txt")
        output$test <- renderUI({test})
      }
  })
  
  #observeEvent(input$inspect_data_type,{
  #  updateSelectInput(session,"inspect_data_type", "Select one of the options below to inspect the data:",choices=c("cat","cat_levels","na","num","types"))
  #})

  
#output$ngs_purity_header_b <- renderTable()
 
 #ngs_purity_table <- read.delim(file="all_ngspurity_output.txt",sep="\t", header=TRUE)
 #observeEvent(input$ngs_rows_to_display,{
   #r <- seq(1:100)
    #ngs_purity_table <- read.delim(file="all_ngspurity_output.txt",sep="\t", header=TRUE)
    #if(input$ngs_rows_to_display %in% r){
    #  output$ngs_purity_header_b <- renderTable({ngs_purity_table[1:input$ngs_rows_to_display,inspect_data_colname(),drop=FALSE]})
    #}else{
    # output$range_error <- renderPrint({"The number of rows entered must be between 1 and 100."})
    #}
  #})

  #observe({
  #  x <- input$ngspurity_header
  #  if(is.null(x)){
  #    shinyjs::disable("create_table")
   # }else{
   #   shinyjs::enable("create_table")
   # }
  


  #observeEvent(input$ngs_rows_to_display,{
    # select rows to display (default 10, maximum 100?)
   # ngs_purity_table <- read.delim(file="all_ngspurity_output.txt",sep="\t", header=TRUE)
  #  if(input$ngs_rows_to_display >100 || input$ngs_rows_to_display <=0){
   #   output$range_error <- renderPrint({"The number of rows entered must be between 1 and 100."})
  #  }else{
  #      output$ngs_purity_header_b <- renderTable({ngs_purity_table[1:input$ngs_rows_to_display ,input$ngspurity_header,drop=FALSE]})
  #      }
  #  })
 
  


}

# output$oncoplot <- renderPlot({})
  