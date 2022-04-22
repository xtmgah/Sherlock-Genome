
source('Sherlock_Genome_Functions.R')

server <- function(input, output, session){
  
  roboto_font_import()
  # disable reset until user clicks "Select Project"
  shinyjs::disable("reset_project")
  
  # setwd based on project code selected
  observeEvent(input$select,{
    os_name <- os_detect()

    setwd(paste0("www/Genomic Data/", input$project_code))
    
  # disable, enable, and show different buttons/selection options
    shinyjs::enable("reset_project")
    shinyjs::disable("project_code")
    shinyjs::disable("select")
    shinyjs::show("choose_files_all")
    shinyjs::show("choose_files_ind")
    shinyjs::show("submit")
    shinyjs::show("reset")
    
    # list of files available for the project selected, separated by module
    output$file_list_select <- renderUI({pickerInput("file_list", 
                                                              label=paste0("Below are the files available for the ", input$project_code, " project:"), choices= list("Study Overview"=files_list()[[9]],"Sample QC"= files_list()[[7]],"NGSpurity"= files_list()[[6]], "Mutations"= files_list()[[5]], "SCNA"= files_list()[[8]], "SV"= files_list()[[11]], "Mutational Signatures"= files_list()[[4]],"Genomic Landscape"= files_list()[[2]], "Clonal Evolution"= files_list()[[1]], "Survival Analysis"= files_list()[[10]], "Integrative Analysis"=files_list()[[3]]), 
                                                              multiple=TRUE, options = list(`actions-box`= TRUE))})
    
  })
  
  # reset button to select a different project
  observeEvent(input$reset_project,{
    shinyjs::enable("project_code")
    shinyjs::enable("select")
    shinyjs::disable("reset_project")
    os_name <- os_detect()
    # if(os_name == "Darwin" | os_name == "Linux"){ #/
    #   setwd("../../..")
    #   # setwd()
    # }
    # if(os_name == "Windows"){ #\
    #   setwd("../../..")
    # }
    shinyjs::hide("file_list")
    shinyjs::hide("choose_files_all")
    shinyjs::hide("choose_files_ind")
    output$file_list_select <- renderUI({updatePickerInput(session, "file_list", label=NULL, choices=NULL)})
    shinyjs::hide("submit")
    shinyjs::hide("reset")
    
  })
  
  # only enable submit button when the user has selected at least one file
  observe({
    x <- input$file_list
    if(is.null(x)){
      shinyjs::disable("submit")
    }else{
      shinyjs::enable("submit")
    }
    
  })
  
  observeEvent(input$submit,{
    x <- input$file_list
    x <- paste (x,sep="", collapse=", ")
    #print(x)
    output$file_load_message <- renderText(paste0("The following files were loaded into their corresponding modules: ", x[1:length(x)], "."))
  })

######### Study Overview ###############################################
  
  # load the study overview for the selected project
  observeEvent(input$submit, {
    if("Study_Overview.Rmd" %in% input$file_list){
      study_overview <- includeMarkdown("Study_Overview/Study_Overview.Rmd")
      output$study_overview_rmd <- renderUI({includeMarkdown("Study_Overview/Study_Overview.Rmd")})}
  })  
#######################################################################
  
  # reactive for dataframe filter (used in sample qc for sample and subject data, and ngspurity)
  df_filter_reactive <- reactive({
    if(input$sbmenu == "sample_qc"){
      if(input$qc_tabs == "qc_sample"){
        req(input$user_filter_input_sample)
        data = read_in_file("Sample_QC/QC_sample_level.txt")
        conditions = input$user_filter_input_sample
        column_selection = input$qcsample_header
      }
      if(input$qc_tabs == "qc_subject"){
        req(input$user_filter_input_subject)
        data = read_in_file("Sample_QC/QC_subject_level.txt")
        conditions = input$user_filter_input_subject
        column_selection = input$qcsubject_header
      }
    }
    
    if(input$sbmenu == "NGSpurity"){
      req(input$user_filter_input)
      data = read_in_file("NGSpurity/all_ngspurity_output.txt")
      conditions = input$user_filter_input
      column_selection = input$ngspurity_header
    }
    
    return(sherlock_genome_filter(data,conditions, column_selection))
    # return(sherlock_genome_filter(data,conditions) # to try in the near future so columns selected when filtering update automatically (right now only update after clicking filter)
  })
  
  # reactive for dataframe inspect (used in sample qc for sample and subject data, and ngspurity)
  inspect_df_option_qc_ngs <- reactive({
    
    if(input$sbmenu == "sample_qc"){
      if(input$qc_tabs == "qc_sample"){
        
        req(input$inspect_data_type_qc_sample)
        req(input$column_name_to_inspect_qc_sample)
        
        dataframe = read_in_file("Sample_QC/QC_sample_level.txt")
        
        type_of_inspection = input$inspect_data_type_qc_sample
        column_name = input$column_name_to_inspect_qc_sample
      }
      
      if(input$qc_tabs == "qc_subject"){ 
        
        req(input$inspect_data_type_qc_subject)
        req(input$column_name_to_inspect_qc_subject)
        
        dataframe = read_in_file("Sample_QC/QC_subject_level.txt")
        
        
        type_of_inspection = input$inspect_data_type_qc_subject
        column_name = input$column_name_to_inspect_qc_subject
      }
    }
    
    if(input$sbmenu == "NGSpurity"){
      req(input$inspect_data_type_ngs)
      req(input$column_name_to_inspect_ngs)
      
      dataframe = read_in_file("NGSpurity/all_ngspurity_output.txt")
      
      type_of_inspection = input$inspect_data_type_ngs
      column_name = input$column_name_to_inspect_ngs
    }
    
    return(inspect_data_function(dataframe, type_of_inspection, column_name))
  })
  
  df_columns_subset <- reactive({
    #req(input$qcsample_header)
    if(input$sbmenu == "sample_qc"){
      if(input$qc_tabs == "qc_sample"){
        dataframe= read_in_file("Sample_QC/QC_sample_level.txt")
        col_names = input$qcsample_header
      }
      if(input$qc_tabs == "qc_subject"){ 
        dataframe= read_in_file("Sample_QC/QC_subject_level.txt")
        col_names = input$qcsubject_header
      }
    }
    if(input$sbmenu == "NGSpurity"){
      dataframe= read_in_file("NGSpurity/all_ngspurity_output.txt")
      col_names = input$ngspurity_header
    }
      
    return(df_columns_subset_function(dataframe, col_names))
    
  })
  
######### Sample QC ###############################################
  
  observeEvent(input$submit,{
    if("QC_sample_level.txt" %in% input$file_list){
      data_qc_sample <- read_in_file("Sample_QC/QC_sample_level.txt")
      qc_sample_header <- read_colnames("Sample_QC/QC_sample_level.txt") 
      # qc_sample_assoc <- data_qc_sample %>% validate_vardf() %>% inspectdf::inspect_types()
      
      ##### View Sample Data #####
      output$qc_sample_header_output <- renderUI({dropdownButton(inputId="qc_sample_dropdown",label="Select Columns",circle=FALSE,checkboxGroupInput("qcsample_header",label="Column Names",choices=c("All columns",colnames(data_qc_sample)),selected=c("Study","Subject","Tumor_Barcode","Normal_Barcode","Gender","Tumor_File")))})
      #output$qc_sample_header_output <- renderUI({dropdownButton(inputId="qc_sample_dropdown",label="Select Columns",circle=FALSE,checkboxGroupInput("qcsample_header",label="Column Names",choices=c("All columns",qc_sample_header),selected=c(qc_sample_header[1:length(qc_sample_header)])))})
      output$filter_input_sample <- renderUI({textInput(inputId="user_filter_input_sample",label="Filter Dataframe",placeholder = "Wave == 'W3'",value="")})
      #output$qc_sample_table <- DT::renderDataTable({datatable(data_qc_sample[ ,input$qcsample_header,drop=FALSE], options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
      output$qc_sample_table <- DT::renderDataTable({datatable(df_columns_subset(), options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
      output$qc_sample_table2 <- DT::renderDataTable({
        req(input$filter_df_sample)
        datatable(isolate(df_filter_reactive()),options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
      
      ##### Inspect Sample Data #####
      output$inspect_df_qc_sample <- renderUI({pickerInput("inspect_data_type_qc_sample","Select one of the options below to inspect the data:" ,choices=c("cat","cat_levels","na","num","types"), multiple=FALSE, selected="cat")})
      output$inspect_sample_qc_select_column <- renderUI({selectInput("column_name_to_inspect_qc_sample","Select one column name to inspect, or all columns:", choices= c("All columns",colnames(data_qc_sample)), multiple= FALSE, selected="All columns")})
      output$inspect_sample_plot <- renderPlot({inspect_df_option_qc_ngs()})
    }
    
    if("QC_subject_level.txt" %in% input$file_list){
      data_qc_subject <- read_in_file("Sample_QC/QC_subject_level.txt")
      qc_subject_header <- read_colnames("Sample_QC/QC_subject_level.txt") 
      # qc_subject_assoc <- data_qc_subject %>% validate_vardf() %>% inspectdf::inspect_types()

      ##### View Subject Data #####
      output$qc_subject_header <- renderUI({dropdownButton(inputId="qc_subject_dropdown",label="Select Columns",circle=FALSE,checkboxGroupInput("qcsubject_header",label="Column Names",choices=c("All columns",qc_subject_header),selected=c("Study","Subject","Barcode","Gender","Sample_ID","File","Wave")))})
      output$filter_input_subject<- renderUI({textInput(inputId="user_filter_input_subject",label="Filter Dataframe",placeholder = "source_material == 'Lung'",value="")})
      
      output$qc_subject_table <- DT::renderDataTable({datatable(df_columns_subset(), options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
      #output$qc_subject_table <- DT::renderDataTable({datatable(data_qc_subject[ ,input$qcsubject_header,drop=FALSE], options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
      
      output$qc_subject_table2 <- DT::renderDataTable({
        req(input$filter_df_subject)
        datatable(isolate(df_filter_reactive()),options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
      
      ##### Inspect Subject Data #####
      output$inspect_df_qc_subject <- renderUI({pickerInput("inspect_data_type_qc_subject","Select one of the options below to inspect the data:" ,choices=c("cat","cat_levels","na","num","types"), multiple=FALSE, selected="cat")})
      output$inspect_subject_qc_select_column <- renderUI({selectInput("column_name_to_inspect_qc_subject","Select one column name to inspect, or all columns:", choices= c("All columns",colnames(data_qc_subject)), multiple= FALSE, selected="All columns")})
      output$inspect_subject_plot <- renderPlot({inspect_df_option_qc_ngs()})
    }
  })

  observeEvent(input$qcsample_header, {
    if("All columns" %in% input$qcsample_header){
      data_qc_sample <- read_in_file("Sample_QC/QC_sample_level.txt")
      qc_sample_header <- read_colnames("Sample_QC/QC_sample_level.txt") 
      updateCheckboxGroupInput(session, "qcsample_header", label="Column Names",choices=c("All columns",colnames(data_qc_sample)),selected=c("All columns",qc_sample_header[1:length(qc_sample_header)]))
    }
  })
  
  # update inspect types depending on variable types (sample data)
  observeEvent(input$column_name_to_inspect_qc_sample,{
    data <- read_in_file("Sample_QC/QC_sample_level.txt")
    if(input$column_name_to_inspect_qc_sample == "All columns"){
      updatePickerInput(session,"inspect_data_type_qc_sample", "Select one of the options below to inspect the data:",choices=c("cat","cat_levels","na","num","types"))
    }
    if(input$column_name_to_inspect_qc_sample != "All columns"){
      if(typeof(data[,input$column_name_to_inspect_qc_sample]) =="character"){
        updatePickerInput(session,"inspect_data_type_qc_sample", "Select one of the options below to inspect the data:",choices=c("cat","cat_levels","na","num","types"),choicesOpt = list(disabled=c("cat","cat_levels","na","num","types") %in% c("num")), options=pickerOptions(hideDisabled=TRUE))
      }
      if(typeof(data[,input$column_name_to_inspect_qc_sample]) == "numeric" || typeof(data[,input$column_name_to_inspect_qc_sample]) == "integer" || typeof(data[,input$column_name_to_inspect_qc_sample]) == "double"){
        updatePickerInput(session,"inspect_data_type_qc_sample", "Select one of the options below to inspect the data:",choices=c("cat","cat_levels","na","num","types"),choicesOpt = list(disabled=c("cat","cat_levels","na","num","types") %in% c("cat","cat_levels")), options=pickerOptions(hideDisabled=TRUE))
      }
    }
    
  })
  
  observeEvent(input$qcsubject_header, {
    if("All columns" %in% input$qcsubject_header){
      data_qc_subject <- read_in_file("Sample_QC/QC_subject_level.txt")
      qc_subject_header <- read_colnames("Sample_QC/QC_subject_level.txt") 
      updateCheckboxGroupInput(session, "qcsubject_header", label="Column Names",choices=c("All columns",colnames(data_qc_subject)),selected=c("All columns",qc_subject_header[1:length(qc_subject_header)]))
    }
  })
  
  # update inspect types depending on variable types (subject data)
  observeEvent(input$column_name_to_inspect_qc_subject,{
    data <- read_in_file("Sample_QC/QC_subject_level.txt")
    if(input$column_name_to_inspect_qc_subject == "All columns"){
      updatePickerInput(session,"inspect_data_type_qc_subject", "Select one of the options below to inspect the data:",choices=c("cat","cat_levels","na","num","types"))
    }
    if(input$column_name_to_inspect_qc_subject != "All columns"){
      if(typeof(data[,input$column_name_to_inspect_qc_subject]) =="character"){
        updatePickerInput(session,"inspect_data_type_qc_subject", "Select one of the options below to inspect the data:",choices=c("cat","cat_levels","na","num","types"),choicesOpt = list(disabled=c("cat","cat_levels","na","num","types") %in% c("num")), options=pickerOptions(hideDisabled=TRUE))
      }
      if(typeof(data[,input$column_name_to_inspect_qc_subject]) == "numeric" || typeof(data[,input$column_name_to_inspect_qc_subject]) == "integer" || typeof(data[,input$column_name_to_inspect_qc_subject]) == "double"){
        updatePickerInput(session,"inspect_data_type_qc_subject", "Select one of the options below to inspect the data:",choices=c("cat","cat_levels","na","num","types"),choicesOpt = list(disabled=c("cat","cat_levels","na","num","types") %in% c("cat","cat_levels")), options=pickerOptions(hideDisabled=TRUE))
      }
    }
    
  })

######### NGSpurity ###############################################
  
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
  
  observeEvent(input$submit,{
    if("all_ngspurity_output.txt" %in% input$file_list){
      data_ngs <- read_in_file("NGSpurity/all_ngspurity_output.txt")
      ngs_purity_header <- read_colnames("NGSpurity/all_ngspurity_output.txt") 
    }
      ##### View Sample Data #####
      output$ngs_purity_header <- renderUI({dropdownButton(inputId="ngspurity_dropdown",label="Select Columns",circle=FALSE,checkboxGroupInput("ngspurity_header",label="Column Names",choices=c("All columns",ngs_purity_header),selected=c("Tumor_Barcode","PGA","PGA_Subclonal","PGA_TETRA","PGA_LOH","PGA_Haploid_LOH")))})

      #output$ngs_purity_table1 <- DT::renderDataTable({datatable(data_ngs[ ,input$ngspurity_header,drop=FALSE], options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
      # output$ngs_purity_table1 <- DT::renderDataTable({
        # datatable(isolate(ngs_column_select()), options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
      output$ngs_purity_table1 <- DT::renderDataTable({datatable(df_columns_subset(), options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
      
      output$filter_input <- renderUI({textInput(inputId="user_filter_input",label="Filter Dataframe",placeholder = "MCN_WGD == 'nWGD'",value="")})
      output$ngs_purity_table2 <- DT::renderDataTable({
        req(input$filter_df)
        datatable(isolate(df_filter_reactive()),options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
      
      ###### Inspect Data #####
      output$inspect_df_test <- renderUI({pickerInput("inspect_data_type_ngs","Select one of the options below to inspect the data:" ,choices=c("cat","cat_levels","na","num","types"), multiple=FALSE, selected="cat")})
      output$ngs_purity_header_inspect_tab <- renderUI({selectInput("column_name_to_inspect_ngs","Select one column name to inspect, or all columns:", choices= c("All columns",colnames(data_ngs)), multiple= FALSE, selected="All columns")})
      output$inspect_ngs_plot <- renderPlot({inspect_df_option_qc_ngs()})
      
      ##### View Figures #####
      os_name <- os_detect()
      if(os_name == "Windows"){
        ngspurity_qc <- read_in_file("NGSpurity/ngspurity_qc_file.txt")
        i <- 1 
        for(each in ngspurity_qc$File){
          ngspurity_qc$File[i] <- paste0("C:/", ngspurity_qc$File[i])
          i <- i + 1
        }
      }else{
        ngspurity_qc <- read_in_file("NGSpurity/ngspurity_qc_file.txt")
      }
      output$ngspurity_barcode <- renderUI({selectInput("tumor_barcode_to_inspect","Select one Tumor Barcode to inspect:", choices= unique(ngspurity_qc$Tumor_Barcode), multiple= FALSE)})
      output$ngspurity_battenberg <- renderUI({selectInput("battenberg_to_inspect","Select one Battenberg to inspect:", choices= unique(ngspurity_qc$Battenberg[ngspurity_qc$Tumor_Barcode==input$tumor_barcode_to_inspect]), multiple= FALSE)})
      output$ngspurity_type <- renderUI({selectInput("type_to_inspect","Select one Type to inspect:", choices= unique(ngspurity_qc$Type[ngspurity_qc$Tumor_Barcode==input$tumor_barcode_to_inspect & ngspurity_qc$Battenberg == input$battenberg_to_inspect]), multiple= FALSE)})
      
      if(os_detect() %in% c("Linux","Darwin")){
        output$figure_pdf <-  renderImage({figure_output()},deleteFile=FALSE)
      }else{
        output$figure_pdf <- renderUI({ tags$iframe(style="height:1000px; width:60%", src= figure_output())})
      }
      
  })
  
  observeEvent(input$ngspurity_header, {
    if("All columns" %in% input$ngspurity_header){
      data_ngs <- read_in_file("NGSpurity/all_ngspurity_output.txt")
      ngs_purity_header <- read_colnames("NGSpurity/all_ngspurity_output.txt") 
      updateCheckboxGroupInput(session, "ngspurity_header", label="Column Names",choices=c("All columns",colnames(data_ngs)),selected=c("All columns",ngs_purity_header[1:length(ngs_purity_header)]))
    }
  })
  
  # update inspect types depending on variable types (nsgpurity data)
   observeEvent(input$column_name_to_inspect_ngs,{
     data <- read_in_file("NGSpurity/all_ngspurity_output.txt")
     
     if(input$column_name_to_inspect_ngs != "All columns"){
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

######### Mutations ###############################################
   
   lolliplot_setup_reactive <- reactive({
    gene = input$lolli_gene_name
    group = input$lolli_group
    #minN = input$lolli_minN
    
    return(lolliplot_setup(gene, group))
     
   })
   
   lolliplot_plot_reactive <- reactive({
     #return(list(tdata0, tslist, samplelist))
     #print(lolliplot_setup_reactive()[1])
     gene = input$lolli_gene_name
     tdata0 =  lolliplot_setup_reactive()[[1]]
     tslist = lolliplot_setup_reactive()[[2]] %>% unlist()
     print(input$tslist_menu)
     tslist_input = input$tslist_menu
     samplelist = lolliplot_setup_reactive()[[3]]
     minN = input$lolli_minN
     
     
     return(lolliplot_plot(gene, tdata0, tslist, tslist_input,samplelist, minN))
   })
   
   observeEvent(input$submit,{
     if(input$project_code == "Sherlock"){
       updateSelectInput(session,"lolli_group",choices= sort(unique(sherlock_overall$SP_Group)))
       #updateSelectInput(session, "lolli_transcript", choices=c("test1", "test2")) # will need to update
       shinyjs::disable("get_ts_info")
     }
     
   })
   
   observeEvent(input$lolli_gene_name,{
     if(input$lolli_gene_name == ""){
       shinyjs::disable("get_ts_info")
     }
     if(input$lolli_gene_name != ""){
       shinyjs::enable("get_ts_info")
     }
   })
   
   observeEvent(input$get_ts_info,{
     #print(lolliplot_setup_reactive()[2])
     choices <- lolliplot_setup_reactive()[[2]] %>% unlist()
     #return(list(tdata0, tslist, samplelist))
     output$tslist_input <- renderUI({selectInput("tslist_menu", "Select a Transcript from the list below:", choices= choices, selected=NULL, multiple= FALSE)})
   })
   
  observeEvent(input$lolliplot_calculate,{
    output$lolliplot <- renderPlot(lolliplot_plot_reactive()[[1]])
    #output$lolliplot <- renderPlot(lolliplot_plot_reactive())
    output$lolliplot_table <-DT::renderDataTable({datatable(lolliplot_plot_reactive()[[2]], options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
    #output$lolliplot_table <- renderDataTable(lolliplot_plot_reactive()[[2]])
  })
  
######### SCNA ###############################################
   
######### SV ###############################################
   
######### Mutational Signatures ###############################################
   
######### Genomic Landscape ###############################################

  observeEvent(input$submit,{
    if("genomePlot_list.txt" %in% input$file_list){
      data_genomePlot <- read_in_file("Genomic_Landscape/genomePlot/genomePlot_list.txt") 
      output$genomePlot_Barcodes <- renderUI({pickerInput("tumor_barcode_to_inspect_genomePlot","Select one Tumor Barcode to inspect:", choices= unique(data_genomePlot$Tumor_Barcode), multiple= FALSE, options= pickerOptions(liveSearch=TRUE, dropupAuto = FALSE))})
      tumor_barcode_to_inspect_genomePlot_reactive <- reactive({
        req(input$tumor_barcode_to_inspect_genomePlot)
        tumor_barcode_to_inspect_genomePlot = input$tumor_barcode_to_inspect_genomePlot
        genomePlot_figurePath <- data_genomePlot$genomePlot[which(data_genomePlot["Tumor_Barcode"] == tumor_barcode_to_inspect_genomePlot)]
      })
      
      output$genomePlot_figure <- renderImage({list(src=paste0("Genomic_Landscape/",tumor_barcode_to_inspect_genomePlot_reactive()),width="800px", height = "auto", alt=paste0("Tumor Barcode: ",input$tumor_barcode_to_inspect_genomePlot))}, deleteFile = FALSE)
    }
  })
    
######### Clonal Evolution ###############################################
   
   figure_output_mutationTime <- reactive({
     req(input$tumor_barcode_mutation_time)
     tumor_barcode = input$tumor_barcode_mutation_time
     if(os_detect() %in% c("Linux","Darwin")){
       return(figure_display_mutationTime(tumor_barcode))
     }else{
       return(src=figure_display_mutationTime(tumor_barcode))
     }
     
   })
   
   observeEvent(input$submit,{
     if("MutationTime_Proportion.txt" %in% input$file_list){
       data_mutation_time <- read_in_file("Clonal_Evolution/MutationTime/MutationTime_Proportion.txt") 
       output$mutation_time_barcode <- renderUI({pickerInput("tumor_barcode_mutation_time","Select one Tumor Barcode to inspect:", choices= unique(data_mutation_time$Tumor_Barcode), multiple= FALSE, options= pickerOptions(liveSearch=TRUE, dropupAuto = FALSE))})
     }
   })
  
   # mutationTime_figurePath <- data_mutation_time$Tumor_Barcode[which(data_mutation_time["Tumor_Barcode"] == tumor_barcode_to_inspect_mutation_time)]
   # mutationTime_figurePath <- paste0("/MutationTime/", mutationTime_figurePath,"_MTime.pdf")
   if(os_detect() %in% c("Linux","Darwin")){
     output$figure_pdf_mutation_time <-  renderImage({figure_output_mutationTime()},deleteFile=FALSE)
   }else{
     output$figure_pdf_mutation_time <- renderUI({ tags$iframe(style="height:1000px; width:60%", src= figure_output_mutationTime())})
   }
  
######### Survival Analysis ###############################################
  
   survival_analysis <- reactive({
     vartmp = input$vartmp_options_select_survival
     sp_group = input$sp_group_selected_survival
     reference = input$reference_survival
     
     # if(input$keyname_checkbox_survival == FALSE){
     #   keyname = NULL
     # }else{
     #   keyname = input$keyname_select_survival
     # }
     keyname= ""
     
     filename = ""
     # if(input$filename_survival == ""){
     #   filename = NULL
     # }else{
     #   filename = input$filename_survival
     # }
     print(paste0("Keyname:",keyname))
     print(vartmp)
     print(sp_group)
     print(reference)
     print(filename)
     return(Survgroup(vartmp=vartmp, sp_group=sp_group, reference=reference, filename=filename, keyname=keyname))

   })
   
   observeEvent(input$submit, {
     if(input$project_code =="Sherlock"){
       updateSelectizeInput(session, 'vartmp_options_select_survival', choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% colnames(), server = TRUE)
       output$sp_group_choices_survival <- renderUI({selectInput("sp_group_selected_survival", "Select one SP_Group to run the fisher test:", choices= sort(unique(sherlock_overall$SP_Group)), multiple= FALSE)})
       updateSelectizeInput(session, 'keyname_select_survival', choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% colnames(), server = TRUE)
     }
   }) 
  
  observeEvent(input$calculate_survival, {
    output$survival_value <- renderTable({isolate(survival_analysis()[[1]])})
    output$survival_plot <- renderPlot(isolate(survival_analysis()[[2]]))
    output$download_survival_plot <- renderUI(actionButton("download_surv_plot", "Download Plot"))
  })
  
  observeEvent(input$download_surv_plot,{
    # pdfhr2() # run when opening app
    surv_filename <- paste0(str_replace_all(input$vartmp_options_select_survival,'[/ \\|]','_'),'_Survival.pdf') #add date for when plot was generated?
    ggsave(surv_filename)
    
  })
######### Integrative Analysis ###############################################
  
  ##### Enrichment Analysis (Fisher Test) #####
  fisher_test <- reactive({
    #vartmp = c('Tumor_Barcode',input$vartmp_options_select)
    vartmp = input$vartmp_options_select
    sp_group = input$sp_group_selected
    if(input$fisher_samplelist_checkbox == FALSE){
      samplelist <- NULL
    }else{
      samplelist = input$fisher_samplelist
      #samplelist <- read_delim("barcode_fisher_testfile.txt", delim = "\t", col_names = FALSE) %>% pull(1) %>% c()
      #samplelist <- read_delim(samplelist, delim = "\t", col_names = FALSE)
      samplelist <- read_delim(samplelist$datapath, delim = "\t", col_names = FALSE) %>% pull(1) %>% c()
    }
   
    if(input$var2name_checkbox == FALSE){
      var2name = NULL
    }else{
      var2name = input$vartmp_options_select_var2
    }
    
    excludes = input$fisher_excludes
    excludes_cat = input$fisher_excludes_cat
    keeps = input$fisher_keeps
    keeps_cat = input$fisher_keeps_cat
    minfreq = input$fisher_min_freq
    freq_column = input$fisher_freq_colnames
    method = input$fisher_test_type
    glm_formula = input$fisher_glm_input
    fdrcutoff = input$fisher_fdr_cutoff
    # clarification on covdata/covdata0
    covdata =  covdata0
    
    print(paste0("sp_group:",sp_group))
    print(paste0("Var2name:",var2name))
    print(paste0("excludes:",excludes))
    print(paste0("excludes_cat:",excludes_cat))
    print(paste0("keeps:",keeps))
    print(paste0("keeps_cat",keeps_cat))
    print(paste0("minfreq:",minfreq))
    return(fishergroup(vartmp,sp_group, samplelist, var2name,excludes,excludes_cat, keeps,keeps_cat, minfreq=0.03, freq_column, method, glm_formula,covdata,fdrcutoff))
    
    # excludes category (Overall Features error)
    # Warning: Error in : Aesthetics must be either length 1 or the same as the data (1): x, y and fill
  })
  
  observeEvent(input$submit, {
    if(input$project_code =="Sherlock"){
      updateSelectizeInput(session, 'vartmp_options_select', choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% colnames(), server = TRUE)
      updateSelectInput(session, "fisher_freq_colnames", choices=colnames(sherlock_freq)[which(str_detect(colnames(sherlock_freq), regex("freq", ignore_case = T)))])
      output$fisher_glm_message <- renderText({
        var_included <- paste0(colnames(covdata0 %>% select(-"Tumor_Barcode")), collapse= ", ")
        paste0("The following variables have already been included in the glm function: ", var_included, ".")})
    }
  }) 
  
  observe({
    # print(input$vartmp_options_select)
    if(input$vartmp_options_select!=""){
      updateSelectizeInput(session, 'vartmp_options_select_var2', choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% select(-input$vartmp_options_select) %>% colnames(), server = TRUE, selected =NULL)
      
      updateSelectizeInput(session, 'fisher_excludes', choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% select(-input$vartmp_options_select) %>% colnames(), server = TRUE, selected =NULL)
      updateSelectizeInput(session, 'fisher_excludes_cat', choices =  unique(sherlock_data_full$Type), server= TRUE)
      
      updateSelectizeInput(session, 'fisher_keeps', choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% select(-input$vartmp_options_select) %>% colnames(), server = TRUE, selected =NULL)
      updateSelectizeInput(session, 'fisher_keeps_cat', choices =  unique(sherlock_data_full$Type), server= TRUE)
      
      # updateSelectizeInput(session, 'fisher_excludes_cat', choices =  mdata0 %>% select(c(2:length(colnames(mdata0)))), server= TRUE)
      # updateSelectizeInput(session, 'fisher_excludes_cat', choices = mdata0 %>% pivot_longer(cols= 2:length(colnames(mdata0))) %>% separate(col = name, into =c('A','B'), sep='\\|') %>% select(B) %>% unique(), server = TRUE, selected =NULL)
      
    }
    
  })
  
  # observeEvent(input$vartmp_options_select,{
  #   
  # })
  observe({
    if(input$vartmp_options_select!="" & input$vartmp_options_select_var2!=""){
      updateSelectizeInput(session, 'fisher_excludes', choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% select(c(-input$vartmp_options_select, -input$vartmp_options_select_var2)) %>% colnames(), server = TRUE)

    }
    
  })
  
  output$sp_group_choices <- renderUI({selectInput("sp_group_selected", "Select one SP_Group to run the fisher test:", choices= sort(unique(sherlock_overall$SP_Group)), multiple= FALSE )})
  
  # observeEvent(input$calculate_fisher,{
  #   
  #   result_table <- isolate(fisher_test())
  #   output$fisher_output_table <- DT::renderDataTable({datatable(result_table[[1]], options=list(searchHighlight=TRUE, order=c(3,'asc')),filter=list(position="top",clear=TRUE,plain=FALSE))})
  # 
  #   output$fisher_output_plots1 <- renderPlot(result_table[[2]])
  #   if(length(result_table) == 3){
  #     output$fisher_output_plots3 <- renderPlot(result_table[[3]])
  #   }
  #   
  #   output$download_fisher_volcano <- renderUI(actionButton("download_fisher_volc_plot", "Download Plot"))
  #   output$download_fisher_barplot <- renderUI(actionButton("download_fisher_bar_plot", "Download Plot"))
  # 
  # })
  
  
  fisher_output_reactive <- eventReactive(input$calculate_fisher, {
    
    # result_table <- isolate(fisher_test())
    isolate(fisher_test())
    
  })
  
  output$fisher_output_table <- DT::renderDataTable({datatable(fisher_output_reactive()[[1]], options=list(searchHighlight=TRUE, order=c(3,'asc')),filter=list(position="top",clear=TRUE,plain=FALSE))}) # %>% bindEvent(input$calculate_fisher)
  
  output$fisher_output_plots1 <- renderPlot(fisher_output_reactive()[[2]])
  # if(length(fisher_output_reactive()) == 3){
  output$fisher_output_plots3 <- renderPlot(fisher_output_reactive()[[3]])
  # }
  
  output$download_fisher_volcano <- renderUI(actionButton("download_fisher_volc_plot", "Download Plot"))
  output$download_fisher_barplot <- renderUI(actionButton("download_fisher_bar_plot", "Download Plot"))
# })
  observeEvent(input$download_fisher_volc_plot,{
    # pdfhr2() # run when opening app
    fisher_volc_filename <- paste0(str_replace_all(input$vartmp_options_select,'[/ \\|]','_'),'_fisher_enrichment_volcano.pdf') #add date for when plot was generated?
    ggsave(fisher_volc_filename,fisher_output_reactive()[[2]])

  })

  observeEvent(input$download_fisher_bar_plot,{
    # pdfhr2() # run when opening app
    # add second variable name to filename?
    fisher_bar_filename <- paste0(str_replace_all(input$vartmp_options_select,'[/ \\|]','_'),'_fisher_enrichment_barplot.pdf') #add date for when plot was generated?
    ggsave(fisher_bar_filename, fisher_output_reactive()[[3]])

  })
  
  observe({
    if(input$var2name_checkbox ==TRUE){
      showTab(inputId = "fisher_results", target = "Figure 2: Bar plot")
    }else{
      hideTab(inputId = "fisher_results", target = "Figure 2: Bar plot")
    }
    
  })
  
  ##### Fisher Bar Plot #####
  
  fisher_test_bar <- reactive({
    #vartmp = c('Tumor_Barcode',input$vartmp_options_select)
    vartmp = input$vartmp_options_select_bar
    sp_group = input$sp_group_selected_bar
    samplelist = NULL
    var2name = input$vartmp_options_select_var2_bar
  
    print("fisherbarplot inputs")
    print(vartmp)
    print(sp_group)
    print(samplelist)
    print(var2name)
    
    # return(fishergroup(vartmp, sp_group, var2name, excludes))
    return(fisherbarplot(vartmp, sp_group, samplelist, var2name))
    
    # excludes category (Overall Features error)
    # Warning: Error in : Aesthetics must be either length 1 or the same as the data (1): x, y and fill
  })
  
  observeEvent(input$submit, {
    if(input$project_code =="Sherlock"){
      updateSelectizeInput(session, 'vartmp_options_select_bar', choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% colnames(), server = TRUE)

    }
  }) 
  
  output$sp_group_choices_bar <- renderUI({selectInput("sp_group_selected_bar", "Select one SP_Group to run the fisher test:", choices= sort(unique(sherlock_overall$SP_Group)), multiple= FALSE )})
  
  
  observe({
    # print(input$vartmp_options_select)
    if(input$vartmp_options_select_bar!=""){
      updateSelectizeInput(session, 'vartmp_options_select_var2_bar', choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% select(-input$vartmp_options_select_bar) %>% colnames(), server = TRUE, selected =NULL)
    }
  })
    
  fisher_output_bar_reactive <- eventReactive(input$generate_fisher_barplot, {
    
    # result_table <- isolate(fisher_test())
    isolate(fisher_test_bar())
    
  })
  
  output$fisher_output_bar <- renderPlot(fisher_output_bar_reactive(), width= 500, height= 550)
  
  output$download_fisher_barplot_separate <- renderUI(actionButton("download_fisher_bar_plot_separate", "Download Plot"))


  ##### Association Testing #####
  # Data select input for association testing
  
  association_inputs <- reactive({
  
    data = data_input_assoc()[[1]]

    plot_width = 12
    plot_height = 8
    
    # multivariable- regression
    if(input$association_testing == "multivar_analysis"){
      if(input$group_var_regression == FALSE){
        # print("regression with no group_var")
        regression =TRUE
        formula = input$regression_formula
        Var1 = NULL
        Var2 = NULL
        filter_zero1 = NULL
        filter_zero2 = NULL
        log_var1 = FALSE
        log_var2 = FALSE
        type = "parametric"
        collapse_var1 = NULL
        collapse_var2 = NULL
        xlab = NULL
        ylab= NULL
        output_plot = input$regression_output_plot
        file_ext = input$file_ext_regression
        
        return(sherlock_genome_association(data, Var1, Var2, regression, formula, filter_zero1, filter_zero2, log_var1, log_var2, type, collapse_var1, collapse_var2, xlab, ylab, output_plot, file_ext, plot_height, plot_width))
        
      }else{
        print("regression WITH group_var")
        regression = TRUE
        formula = input$regression_formula_group
        Group_Var = input$group_var_input
        Var1 = input$variable_choices_1_group
        Var2 = input$variable_choices_2_group
        
        filter_zero1 = NULL
        filter_zero2 = NULL
        log_var1 = FALSE
        log_var2 = FALSE
        collapse_var1 = NULL
        collapse_var2 = NULL
        type = "parametric"
        
        return(sherlock_genome_association_group(data, Var1, Var2, Group_Var, regression, formula, filter_zero1, filter_zero2,log_var1,log_var2,type,collapse_var1,collapse_var2))
        
      }
    }
    
    # bivariable
    if(input$association_testing == "bivar_analysis"){
      if(input$group_var_bivariable ==FALSE){
        regression = FALSE
        formula = NULL
        Var1 = input$variable_choices_1
        Var2 = input$variable_choices_2
        filter_zero1 = input$filter_assoc_var1
        filter_zero2 = input$filter_assoc_var2
        log_var1 = input$log2_assoc_var1
        log_var2 = input$log2_assoc_var2
        type = input$assoc_types_list
        xlab = input$variable_choices_1
        ylab= input$variable_choices_2
        output_plot = input$association_output_plot
        file_ext = input$file_ext_assoc
        
        if(input$collapse_assoc_var1==""){
          collapse_var1 = NULL
        }else{
          collapse_var1 = input$collapse_assoc_var1
        }
        
        if(input$collapse_assoc_var2==""){
          collapse_var2 = NULL
        }else{
          collapse_var2 = input$collapse_assoc_var2
        }
        return(sherlock_genome_association(data, Var1, Var2, regression, formula, filter_zero1, filter_zero2, log_var1, log_var2, type, collapse_var1, collapse_var2, xlab, ylab, output_plot, file_ext, plot_height, plot_width))
        
      }else{
        
        Var1 = input$variable_choices_1_group
        Var2 = input$variable_choices_2_group
        Group_Var = input$group_var_input_bivariable
        filter_zero1 = input$filter_assoc_var1_group
        filter_zero2 = input$filter_assoc_var2_group
        log_var1 = input$log2_assoc_var1_group
        log_var2 = input$log2_assoc_var2_group
        type = input$assoc_types_list_group
        regression = FALSE
        formula =NULL
        
        if(input$collapse_assoc_var1_group==""){
          collapse_var1 = NULL
        }else{
          collapse_var1 = input$collapse_assoc_var1_group
        }
        
        if(input$collapse_assoc_var2_group==""){
          collapse_var2 = NULL
        }else{
          collapse_var2 = input$collapse_assoc_var2_group
        }
        
        return(sherlock_genome_association_group(data, Var1, Var2, Group_Var, regression, formula, filter_zero1, filter_zero2,log_var1,log_var2,type, collapse_var1, collapse_var2))
        
      } 
    }
  })
  
  data_input_assoc <- reactive({
    
    # browser()
    df_list <- list()
    if(length(input$data_list_association_selected) > 1){
      no_tum_barcode <- list()
      if("Sample Level" %in% input$data_list_association_selected){
        sample_data <- read_in_file("Sample_QC/QC_sample_level.txt")
        if("Tumor_Barcode" %in% colnames(sample_data)){
          df_list <- append(df_list, list(sample_data), after = length(df_list))
          print("Sample Level added")
        }else{
          print("Sample level could not be joined with the selected datasets because it does contain the common column to be joined by.")
          no_tum_barcode <- append(no_tum_barcode, "sample level", after = length(df_list))
        }
        
      }
      if("Subject Level" %in% input$data_list_association_selected){
        subject_data <- read_in_file("Sample_QC/QC_subject_level.txt")
        if("Tumor_Barcode" %in% colnames(subject_data)){
          df_list <- append(df_list, list(subject_data), after = length(df_list))
          print("Subject Level added")
        }else{
          print("Subject level could not be joined with the selected datasets because it does contain the common column to be joined by.")
          no_tum_barcode <- append(no_tum_barcode, "subject level", after = length(df_list))
        }
        
      }
      if("NGSpurity" %in% input$data_list_association_selected){
        ngspurity_data <- read_in_file("NGSpurity/all_ngspurity_output.txt")
        if("Tumor_Barcode" %in% colnames(ngspurity_data)){
          df_list <- append(df_list, list(ngspurity_data), after = length(df_list))
          print("NGSpurity added")
        }else{
          print("NGSpurity could not be joined with the selected datasets because it does contain the common column to be joined by.")
          no_tum_barcode <- append(no_tum_barcode, "ngspurity", after = length(df_list))
        }
        
      if("Survival Data" %in% input$data_list_association_selected){
        #suvdata
        if("Tumor_Barcode" %in% colnames(suvdata)){
          df_list <- append(df_list, list(suvdata), after = length(df_list))
          print("Survival Data added")
        }else{
          print("Survival Data could not be joined with the selected datasets because it does contain the common column to be joined by.")
          no_tum_barcode <- append(no_tum_barcode, "survival data", after = length(df_list))
        }
      }
        
        assoc_join <- join_assoc_data(df_list)
        assoc_colnames <- colnames(assoc_join %>% select(!ends_with("_Detail"),-"Tumor_Barcode"))
        assoc_join_by <- print("Joined by: Tumor Barcode")
        if(length(no_tum_barcode >=1)) {
          assoc_no_tum_barcode <- print(paste0("The following datasets were not able to be joined because they lack the Tumor Barcode column:", no_tum_barcode, "."))
          return(list(assoc_join, assoc_colnames, assoc_join_by, assoc_no_tum_barcode))
        }else{
          return(list(assoc_join, assoc_colnames, assoc_join_by))
        }
        
      }
      
    }else{
      if("Sample Level" %in% input$data_list_association_selected){
        sample_data <- read_in_file("Sample_QC/QC_sample_level.txt")
        df_list <- append(df_list, list(sample_data), after = length(df_list))
      }
      if("Subject Level" %in% input$data_list_association_selected){
        subject_data <- read_in_file("Sample_QC/QC_subject_level.txt")
        df_list <- append(df_list, list(subject_data), after = length(df_list))
      }
      if("NGSpurity" %in% input$data_list_association_selected){
        ngspurity_data <- read_in_file("NGSpurity/all_ngspurity_output.txt")
        df_list <- append(df_list, list(ngspurity_data), after = length(df_list))
      }
      
      if("Survival Data" %in% input$data_list_association_selected){
        #suvdata
        df_list <- append(df_list, list(suvdata), after = length(df_list))
      }
      
      assoc_join <- data.frame(df_list[1])
      if("Tumor_Barcode" %in% colnames(assoc_join)){
        print("option 1")
        assoc_colnames <- colnames(assoc_join %>% select(!ends_with("_Detail"),-"Tumor_Barcode"))
      }else{
        print("option 2")
        assoc_colnames <- colnames(assoc_join %>% select(!ends_with("_Detail")))
      }
      
      return(list(assoc_join, assoc_colnames))
    }
    
      # print(typeof(df_list[1]))
      # if(length(df_list) > 1){
      #   # assoc_join <- join_assoc_data(df_list)
      #   assoc_join <- join_assoc_data(df_list) #list of results
      #   print(assoc_join)
      #   assoc_colnames <- colnames(assoc_join %>% select(!ends_with("_Detail"),-"Tumor_Barcode"))
      #   assoc_join_by <- print("Joined by: Tumor Barcode")
      # 
      # }#else{
      #   assoc_join <- data.frame(df_list[1])
      #   #currently does not work for "Subject Data" because it does not have a "Tumor_Barcode" column
      #   assoc_colnames <- colnames(assoc_join %>% select(!ends_with("_Detail"),-"Tumor_Barcode"))
      # }
      
      # return(list(assoc_join, assoc_colnames, assoc_join_by))
    })

  observeEvent(input$submit, {
    if(input$project_code == "Sherlock"){
      output$data_list_association <- renderUI({pickerInput("data_list_association_selected", 
                                                            label=paste0("Select one or more datasets for association testing:"), choices= list("Sample Level", "Subject Level", "NGSpurity", "Survival Data"), 
                                                            multiple=TRUE, options = list(`actions-box`= TRUE))})
    }
    
    output$assoc_variable_list_1 <- renderUI({pickerInput("variable_choices_1", "Variable One", choices= NULL, multiple=FALSE)})
    output$assoc_variable_list_2 <- renderUI({pickerInput("variable_choices_2", "Variable Two",choices= NULL, multiple=FALSE)})
    
    output$assoc_variable_list_1_group <- renderUI({pickerInput("variable_choices_1_group", "Variable One", choices= NULL, multiple=FALSE)})
    output$assoc_variable_list_2_group <- renderUI({pickerInput("variable_choices_2_group", "Variable One", choices= NULL, multiple=FALSE)})
  })
  
  observeEvent(input$load_datasets,{
    
    print(length(data_input_assoc()))
    data_final <- data_input_assoc()[[1]]
    data_final_dim <- data_input_assoc()[[1]] %>% dim()
    print(data_final_dim)
    data_colnames <- data_input_assoc()[[2]]

    if(length(data_input_assoc()) >2 ){
      if(length(data_input_assoc()) ==3){
        data_join_by <- data_input_assoc()[[3]]
        output$assoc_data_join_by <- renderText(data_join_by)
        
      }
      if(length(data_input_assoc()) ==4){
        data_join_by <- data_input_assoc()[[3]]
        data_no_tum_barcode <- data_input_assoc()[[4]]
        
        output$assoc_data_join_by <- renderText(data_join_by)
        output$assoc_no_tum_barcode <- renderText(data_no_tum_barcode)
      }
      
      
    }
    # if(data_input_assoc()[[3]] %in% data_input_assoc()){
    #   data_join_by <- data_input_assoc()[[3]]
    #   print(data_join_by)
    #   output$assoc_data_join_by <- renderText(data_join_by)
    # }
    # if(data_input_assoc()[[4]] %in% data_input_assoc()){
    #   data_no_tum_barcode <- data_input_assoc()[[4]]
    #   print(data_no_tum_barcode)
    #   output$assoc_no_tum_barcode <- renderText(data_no_tum_barcode)
    # }


    # select columns dropdown
    output$assoc_dropdown <- renderUI({dropdownButton(inputId="assoc_data_dropdown",label="Select Columns",circle=FALSE,checkboxGroupInput("assoc_header",label="Column Names",choices=c("All columns",data_colnames),selected=data_colnames[1:6]))})
    output$assoc_datatable_dim <- renderTable({tibble(rows=data_final_dim[[1]], columns=data_final_dim[[2]])})
    # if(length(data_input_assoc() == 4)){
    #   output$assoc_data_join_by <- renderText(data_join_by)
    #   output$assoc_no_tum_barcode <- renderText(data_no_tum_barcode)
    # }
    # if(data_join_by ==TRUE){
      # output$assoc_data_join_by <- renderText(data_join_by)
    # }
    # if(data_no_tum_barcode ==TRUE){
      # output$assoc_no_tum_barcode <- renderText(data_no_tum_barcode)
    # }

    # output$assoc_datatable <- renderDT({datatable(data_final,filter=list(position="top",clear=TRUE,plain=FALSE))})
    output$assoc_datatable <- DT::renderDataTable({datatable(data_final[ ,input$assoc_header,drop=FALSE], options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
    
    updatePickerInput(session, "variable_choices_1", "Variable One", choices= data_colnames)
    updatePickerInput(session, "variable_choices_2", "Variable Two", choices= data_input_assoc()[[2]])
    
    updatePickerInput(session, "variable_choices_1_group", "Variable One", choices= data_colnames)
    updatePickerInput(session, "variable_choices_2_group", "Variable Two", choices= data_input_assoc()[[2]])
    # updatePickerInput(session, "variable_choices_2_group","Variable Two",choices=data_colnames, choicesOpt = list(disabled=data_colnames %in% input$variable_choices_1), options=pickerOptions(hideDisabled=TRUE))
    
  })
  
########### Multivariable (Regression) and Regression Group ##############
  
  ##### Regression #####
  
  # Input
  output$formula_input <- renderUI({textInput("regression_formula", "Regression Formula", placeholder="lm( mpg ~ vs + gear)")})
  
  # Results
  output$regression_plot <- renderPlot({
    req(input$calculate_regression)
    isolate(association_inputs())})
  
  observe({
    if(input$regression_output_plot ==TRUE){
      shinyjs:: show("file_ext_regression")
    }else{
      shinyjs:: hide("file_ext_regression")
    }
  })
  
  # Reset
  observeEvent(input$reset_regression, {
    reset("regression_formula")
  })
  
  
  ##### Regression Group #####
  
  # Input
  output$group_var_text_box <- renderUI({textInput(inputId= "group_var_input", label= NULL,placeholder= "Group Variable Name" )})
  output$formula_input_group <- renderUI({textInput("regression_formula_group", "Regression Formula", placeholder="lm( mpg ~ vs + gear)")})
  
  # Results
  output$regression_grp_result <- DT::renderDataTable({
    req(input$calculate_regression_group)
    isolate(association_inputs())})
  
  # Reset
  observeEvent(input$reset_regression_group,{
    reset("regression_formula_group")
    reset("group_var_input")
    updateCheckboxGroupInput(inputId= "regression_model_group_input", label= "Select a supported regression model type", choices= c("lm", "glm"))
  })
  

########### Bivariable and Bivariable Group ##############
  
  ##### Bivariable #####
  
  # remove first variable from the second variable list
  observeEvent(input$variable_choices_1,{
    data_colnames <- data_input_assoc()[[2]]
    updatePickerInput(session, "variable_choices_2","Variable Two",choices=data_colnames, choicesOpt = list(disabled=data_colnames %in% input$variable_choices_1), options=pickerOptions(hideDisabled=TRUE))
  })
  
  # update the types of association tests permitted depending on the variable one and variable two types
  observe( {
    req(input$variable_choices_1, input$variable_choices_2)
    data <- data_input_assoc()[[1]]
    data <- validate_vardf(data)
    if(is.numeric(data[[input$variable_choices_1]]) && is.numeric(data[[input$variable_choices_2]])) {
      updatePickerInput(session,"assoc_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes", "skit"), options=pickerOptions(hideDisabled=TRUE))
    }
    if(is.factor(data[[input$variable_choices_1]]) && is.factor(data[[input$variable_choices_2]])) {
      updatePickerInput(session,"assoc_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes", "fisher"), options=pickerOptions(hideDisabled=TRUE))
    }
    if(is.numeric(data[[input$variable_choices_1]]) && is.factor(data[[input$variable_choices_2]])) {
      updatePickerInput(session,"assoc_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes"), options=pickerOptions(hideDisabled=TRUE))
    }
    if(is.factor(data[[input$variable_choices_1]]) && is.numeric(data[[input$variable_choices_2]])) {
      updatePickerInput(session,"assoc_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes"), options=pickerOptions(hideDisabled=TRUE))
    }
  })
  
  observe({
    req(input$variable_choices_1, input$variable_choices_2)
    data <- data_input_assoc()[[1]]
    data <- validate_vardf(data)
    if(is.numeric(data[[input$variable_choices_1]])){
      shinyjs::disable("collapse_assoc_var1")
    }else{
      shinyjs::enable("collapse_assoc_var1")
    }
    if(is.numeric(data[[input$variable_choices_2]])){
      shinyjs::disable("collapse_assoc_var2")
    }else{
      shinyjs::enable("collapse_assoc_var2")
    }
  })
  
  
  observe({
    if(input$association_output_plot ==TRUE){
      shinyjs:: show("file_ext_assoc")
    }else{
      shinyjs:: hide("file_ext_assoc")
    }
  })
  
  output$bivariable_plot <- renderPlot({
    req(input$calculate_association)
    isolate(association_inputs())})
  
  # reset association inputs
  observeEvent(input$reset_association,{
    data <- data_input_assoc()[[1]]
    data_colnames <- data_input_assoc()[[2]]
    
    updatePickerInput(session, inputId="variable_choices_1","Variable One", choices= data_colnames)
    updatePickerInput(session, inputId= "variable_choices_2", "Variable Two", choices = data_colnames)
    updatePickerInput(session, "variable_choices_2","Variable Two",choices=data_colnames, choicesOpt = list(disabled=data_colnames %in% input$variable_choices_1), options=pickerOptions(hideDisabled=TRUE))
    updateCheckboxInput(session, "filter_assoc_var1","Filter (>0)",value=FALSE)
    updateCheckboxInput(session,"log2_assoc_var1", "log2", value=FALSE)
    updateTextInput(session, "collapse_assoc_var1", "Collapse Level",value=NULL)
    updateCheckboxInput(session, "filter_assoc_var2","Filter (>0)",value=FALSE)
    updateCheckboxInput(session,"log2_assoc_var2", "log2", value=FALSE)
    updateTextInput(session,"collapse_assoc_var2", "Collapse Level",value=NULL)
    updatePickerInput(session,"assoc_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes","fisher","skit"))
    updateCheckboxInput(session,"assoc_output_plot", "Save plot to computer",value = FALSE)
    updateTextInput(session,"file_ext", "File Extension (png, svg, or jpg)", value="png")
  })
  
  ##### Bivariable Group #####
  
  # remove first variable from the second variable list in group function
  observeEvent(input$variable_choices_1_group,{
    data_colnames <- data_input_assoc()[[2]]
    updatePickerInput(session, "variable_choices_2_group","Variable Two",choices=data_colnames, choicesOpt = list(disabled=data_colnames %in% input$variable_choices_1_group), options=pickerOptions(hideDisabled=TRUE))
  })

  # update the types of association tests permitted depending on the variable one and variable two types in group function
  observe( {
    req(input$variable_choices_1_group, input$variable_choices_2_group)
    data <- data_input_assoc()[[1]]
    data <- validate_vardf(data)
    if(is.numeric(data[[input$variable_choices_1_group]]) && is.numeric(data[[input$variable_choices_2_group]])) {
      updatePickerInput(session,"assoc_types_list_group", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes", "skit"), options=pickerOptions(hideDisabled=TRUE))
    }
    if(is.factor(data[[input$variable_choices_1_group]]) && is.factor(data[[input$variable_choices_2_group]])) {
      updatePickerInput(session,"assoc_types_list_group", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes", "fisher"), options=pickerOptions(hideDisabled=TRUE))
    }
    if(is.numeric(data[[input$variable_choices_1_group]]) && is.factor(data[[input$variable_choices_2_group]])) {
      updatePickerInput(session,"assoc_types_list_group", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes"), options=pickerOptions(hideDisabled=TRUE))
    }
    if(is.factor(data[[input$variable_choices_1_group]]) && is.numeric(data[[input$variable_choices_2_group]])) {
      updatePickerInput(session,"assoc_types_list_group", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes"), options=pickerOptions(hideDisabled=TRUE))
    }
  })

  # update the collapse level option depending on the variable one and variable two types in group function
  observe({
    req(input$variable_choices_1_group, input$variable_choices_2_group)
    data <- data_input_assoc()[[1]]
    data <- validate_vardf(data)
    if(is.numeric(data[[input$variable_choices_1_group]])){
      shinyjs::disable("collapse_assoc_var1_group")
    }else{
      shinyjs::enable("collapse_assoc_var1_group")
    }
    if(is.numeric(data[[input$variable_choices_2_group]])){
      shinyjs::disable("collapse_assoc_var2_group")
    }else{
      shinyjs::enable("collapse_assoc_var2_group")
    }
  })
  
  output$group_var_text_box_bivariable <- renderUI({textInput(inputId= "group_var_input_bivariable", label= NULL,placeholder= "Group Variable Name")})
  
  output$bivariable_grp_result <- DT::renderDataTable({
    req(input$calculate_association_group)
    isolate(association_inputs())})
  
  # reset association group inputs
  observeEvent(input$reset_association_group,{
    data <- data_input_assoc()[[1]]
    data_colnames <- data_input_assoc()[[2]]

    updatePickerInput(session, inputId="variable_choices_1_group","Variable One", choices= data_colnames)
    updatePickerInput(session, inputId= "variable_choices_2_group", "Variable Two", choices = data_colnames)
    updatePickerInput(session, "variable_choices_2_group","Variable Two",choices= data_colnames, choicesOpt = list(disabled=data_colnames %in% input$variable_choices_1_group), options=pickerOptions(hideDisabled=TRUE))
    updateCheckboxInput(session, "filter_assoc_var1_group","Filter (>0)",value=FALSE)
    updateCheckboxInput(session,"log2_assoc_var1_group", "log2", value=FALSE)
    updateTextInput(session, "collapse_assoc_var1_group", "Collapse Level",value=NULL)
    updateCheckboxInput(session, "filter_assoc_var2_group","Filter (>0)",value=FALSE)
    updateCheckboxInput(session,"log2_assoc_var2_group", "log2", value=FALSE)
    updateTextInput(session,"collapse_assoc_var2_group", "Collapse Level",value=NULL)
    updatePickerInput(session,"assoc_types_list_group", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes","fisher","skit"))
    reset("group_var_input_bivariable")
  })
  
  ##### Oncoplot #####
  
  observeEvent(input$submit, {
    if(input$project_code =="Sherlock"){
      type_list <- sherlock_data_full %>% select(Type) %>% unique()
      # print(type_list)
      updateSelectizeInput(session, "oncoplot_genalts_select", choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% colnames(), server = TRUE)
      updateSelectizeInput(session, "oncoplot_categories_select", choices =  unique(sherlock_data_full$Type), server= TRUE)
      updateSelectizeInput(session, "oncoplot_cat_freq_select", choices =  unique(sherlock_data_full$Type), server= TRUE)
      
    }
  })
  
  oncoplot_reactive <- reactive({
    
    if(input$submit){
      if(input$project_code == "Sherlock"){
        # will need to edit these later
        sample_level0 <- sherlock_data_full %>% pull(Tumor_Barcode) %>% unique()
        # sample_new_level <- sherlock_data_full %>% pull(Tumor_Barcode) %>% unique()
        freq_table <- sherlock_freq
      }
    }
    
    # 1
    if(input$oncoplot_input_selection == 'oncoplot_genalts'){
      genomic_alts = input$oncoplot_genalts_select 
      order_by_input = input$order_by_input1
      opt_four_freq = NULL
      print(paste0("1:",genomic_alts))
    }
    
    # 2
    if(input$oncoplot_input_selection == 'oncoplot_text_input_genalts'){
      genomic_alts <- input$oncoplot_user_input_genalts
      order_by_input = input$order_by_input2
      print(genomic_alts)
      opt_four_freq = NULL
      # print(paste0("2:",genomic_alts))
      # print(genomic_alts)
    }
    
    # 3 
    if(input$oncoplot_input_selection == 'oncoplot_categories'){
      genomic_alts = input$oncoplot_categories_select
      opt_four_freq = NULL
      print(paste0("3:",genomic_alts))
    }
    
    # 4
    if(input$oncoplot_input_selection == 'oncoplot_cat_freq_inputs'){
      opt_four_freq <- c()
      genomic_alts_final <- c()
      genomic_alts = input$oncoplot_cat_freq_select
      genomic_alts <- str_split(genomic_alts, "\\n") %>% unlist()
      i <- 1
      for(each in genomic_alts){
        genomic_alts[i] <- str_split(genomic_alts[i], ",")
        genomic_alts_final <- append(genomic_alts_final, trimws(genomic_alts[[i]][1]))
        opt_four_freq <- append(opt_four_freq, trimws(genomic_alts[[i]][2]))
        i <- i + 1
        # print(genomic_alts[i])
      }
      # print(paste0("4:",genomic_alts))
      genomic_alts <- genomic_alts_final
    }
    

    # print(genomic_alt_four)
    frequency = input$oncoplot_min_freq
    # 
    # print(genomic_alt_one, genomic_alt_two, genomic_alt_three, genomic_alt_four)
    # data <- oncoplot_data_prep(genomic_alt)
    data <- oncoplot_data_prep(genomic_alts, opt_four_freq, freq_table)
    # data <- oncoplot_data_prep(genomic_alt_one, genomic_alt_two, genomic_alt_three, genomic_alt_four)
    # print(data)
    
    landscape_colors <- get_vcColors()
    if(data[[1]] == "source_cat"){
      reps <- rep(paste0("data[[2]][[",1:length(data[[2]]),"]]"))
      reps <- paste0(reps, collapse= " , ")
      # rbind_list <- paste0("data_for_colors <- rbind(", reps, ")")
      rbind_list <- paste0("data_all_cats <- rbind(", reps, ")")
      eval(parse(text=rbind_list))
      # # data_for_colors <- rbind(data[[1]], data[[2]], data[[3]])
      # # print(data_for_colors)
      # # should the color be taken out of the options for the next data_alts
      data_for_colors <- data_all_cats %>% arrange(Gene, Type)
      data_alts <- data_for_colors %>% pull(Alteration) %>% unique()
    }else{
      data_alts <- data[[2]] %>% pull(Alteration) %>% unique()
    }
    
    # print(data_alts)
    # data_alts <- data %>% pull(Alteration) %>% unique()
    color_list <- read.delim("../landscape_colors.csv", header = TRUE, sep = ",")
    i <- 1
    new_colors <- c()
    for(each in data_alts){
      if(!(data_alts[i] %in% names(landscape_colors))){
        data_alt <- data_alts[i]
        # print(data_alt)
        # assign color to each alteration without a color already set 
        # landscape_colors[data_alts[i]] <- sample(color_list$Value, 1)
        new_colors[data_alt] <- sample(color_list$Value, 1)
        # print(new_colors)
        # print(paste0(data_alts[i], landscape_colors[data_alts[i]]))
        # names(new_colors)[i] <- data_alt
      }
        i <- i + 1
        # print(landscape_colors)
    }
    
    # print(names(new_colors))
    landscape_colors <- c(landscape_colors, new_colors)
    # print(landscape_colors)
    
    oncolist_result <- list()
    oncolist_result_by_freq <- list()
    if(data[[1]] == "ind_genalt"){
      data <- data[[2]]
      # data_set_sample_level <- data %>% filter(Gene == data$Gene[1])
      # plot_color_vec <- landscape_colors[which(names(landscape_colors) %in% data_set_sample_level$Alteration)] %>% c()
      plot_color_vec <- landscape_colors[which(names(landscape_colors) %in% data$Alteration)] %>% c()
      
      # print(plot_color_vec)
      # gene_level <- c("Kataegis", "WGD_Status")
      gene_type_level <- data %>% select(Gene,Type) %>% unique()
      print("gene_type_level")
      print(gene_type_level)
      gene_level <- c(data %>% pull(Gene) %>% unique())
      # res <- oncoplot(data=data_set_sample_level,gene_level=NULL,sample_level=NULL, sample_level0 = sample_level0,landscape_colors = plot_color_vec, p2_axis_hidden = TRUE,p2_hidden = TRUE, GeneSortOnly=FALSE, frequency= frequency, gen_alt_second = TRUE)
      #res <- oncoplot(data=data,gene_level=NULL,sample_level=NULL, sample_level0 = sample_level0,landscape_colors = plot_color_vec, p2_axis_hidden = TRUE,p2_hidden = TRUE, GeneSortOnly=FALSE, frequency= frequency, gen_alt_second = TRUE, order_by_input= order_by_input)
      res <- oncoplot(data=data,gene_level=NULL,sample_level=NULL, sample_level0 = sample_level0,landscape_colors = plot_color_vec, p2_axis_hidden = TRUE,p2_hidden = TRUE, GeneSortOnly=FALSE, frequency= frequency, gen_alt_second = TRUE, order_by_input= order_by_input)
      
      tmp_sample_level <- res$sample_level %>% pull(Tumor_Barcode)
      tmp_gene_level <- res$gene_level
      print("tmp_gene_level")
      print(tmp_gene_level)
      tmp_gene_freq <- res$gene_freq %>% arrange(desc(Freq)) %>% pull(Gene)
      print("tmp_gene_freq")
      print(tmp_gene_freq)
      
      res <- oncoplot(data=data,gene_level=tmp_gene_level,sample_level=tmp_sample_level, sample_level0 = sample_level0,landscape_colors = plot_color_vec, p2_axis_hidden = TRUE,p2_hidden = TRUE, GeneSortOnly=FALSE, frequency= frequency, gen_alt_second = TRUE)
      
      if(order_by_input){
        i <- 1
        for(each in unique(data$Gene)){
          data_ind_gen_alt <- data %>% filter(Gene ==unique(data$Gene)[i])
          print(data$Gene[i])
          plot_color_vec <- landscape_colors[which(names(landscape_colors) %in% data_ind_gen_alt$Alteration)] %>% c()
          res <- oncoplot(data=data_ind_gen_alt,gene_level= NULL,sample_level=tmp_sample_level, sample_level0 = sample_level0,landscape_colors = plot_color_vec, p2_axis_hidden = TRUE,p2_hidden = TRUE, GeneSortOnly=FALSE, frequency= frequency, gen_alt_second= TRUE, order_by_input= TRUE)
          oncolist_result <- list.append(oncolist_result, res)
          i <- i + 1
          print(paste0("the value of i is:", i))
        }
        
      }
      # if(order_by_frequency){
      #   res <- oncoplot(data=data,gene_level=tmp_gene_level,sample_level=tmp_sample_level, sample_level0 = sample_level0,landscape_colors = plot_color_vec, p2_axis_hidden = TRUE,p2_hidden = TRUE, GeneSortOnly=FALSE, frequency= frequency)
      #   oncolist_result_by_freq <- list.append(oncolist_result_by_freq, res)
      #   # pull Type from original genomic alts?
      #   data2 <- tibble(Gene=tmp_gene_freq, Type=gene_type_level$Type[which(gene_type_level$Gene==tmp_gene_freq)])
      #   print("data2")
      #   print(data2)
      # }else{
      #   data2 <- data %>% select(Gene,Type) %>% unique()
      # }
      # 
      # # data2 <- subset(data2,!(Gene %in% c(tmp_gene_level)))
      # print(data2)
      # i <- 1
      # for(i in 1:length(data2$Gene)){
      # data_ind_alt <- data %>% filter(Type== data2$Type[i]) %>% filter(Gene== data2$Gene[i])
      # print(data_ind_alt)
      # print(paste0(data2$Type[i], data2$Gene[i]))
      # plot_color_vec <- landscape_colors[which(names(landscape_colors) %in% data_ind_alt$Alteration)] %>% c()
      # print(plot_color_vec)
      # res <- oncoplot(data=data_ind_alt,gene_level= NULL,sample_level=tmp_sample_level, sample_level0 = sample_level0,landscape_colors = plot_color_vec, p2_axis_hidden = TRUE,p2_hidden = TRUE, GeneSortOnly=FALSE, frequency= frequency, gen_alt_second= TRUE)
      #     # plot(res$oncoplot_legend)
      #     # print(res$gene_level)
      # oncolist_result <- list.append(oncolist_result, res)
      # print(res$sample_level)
      # i <- i + 1
      # }

      # sample level
      # oncolist_result <- list.append(oncolist_result, res)
      
      
        # print(i)
      
      # res <- oncoplot(data=data,gene_level= tmp_gene_level,sample_level=tmp_sample_level, sample_level0 = sample_level0,landscape_colors = plot_color_vec, p2_axis_hidden = TRUE,p2_hidden = TRUE, GeneSortOnly=FALSE, frequency= frequency, gen_alt_second= TRUE)
      # p_legends <- plot_grid(p_study_legend,p_purity_legend, p_cosine_legend, p_proportion_legend, align = 'v',axis = 'l',ncol = 1,rel_heights = c(1,1,1,1.5))
      # p_main <- plot_grid(p_cluster,p_study,p_purity,p_mutation,p_cosine,p_proportion,align = 'v',axis = 'lr',ncol = 1,rel_heights = c(2,h_study,h_purity,4,0.05,7))
      # p_all <- plot_grid(p_main,p_legends+theme(plot.margin = margin(l = -0.5,unit = 'cm')),ncol = 2,rel_widths = c(10,2))
      # if(order_by_frequency){
      #   reps <- rep(paste0("oncolist_result_by_freq[[",1:length(oncolist_result_by_freq),"]]$oncoplot[[1]]"))
      #   reps <- paste0(reps, collapse= " , ")
      #   print(reps)
      #   onco_all <- paste0("oncoplots <- plot_grid(", reps, ", align= 'v', ncol =1)")
      #   #onco_all <- paste0("oncoplots <- plot_grid(", reps, ", align = 'v',axis = 'lr',ncol = 1,rel_heights = c(2,4,0.05,7)")
      #   # print(onco_all)
      #   eval(parse(text=onco_all))
      #   # print(oncoplots)
      # }else{
      if(order_by_input){
        reps <- rep(paste0("oncolist_result[[",1:length(oncolist_result),"]]$oncoplot[[1]]"))
        reps <- paste0(reps, collapse= " , ")
        print(reps)
        onco_all <- paste0("oncoplots <- plot_grid(", reps, ", align= 'v', ncol =1)")
        #onco_all <- paste0("oncoplots <- plot_grid(", reps, ", align = 'v',axis = 'lr',ncol = 1,rel_heights = c(2,4,1,1,0.05,7)")
        # print(onco_all)
        eval(parse(text=onco_all))
      }else{
        oncoplots <- res$oncoplot[[1]]
      }
        
        # print(oncoplots)
      #}

      # 
      # reps <- rep(paste0("oncolist_result[[",1:length(oncolist_result),"]]$oncoplot_legend"))
      # reps <- paste0(reps, collapse= " , ")
      # print(reps)
      # num_plots <- length(oncolist_result)
      # print(num_plots)
      # legs_all <- paste0("plegends <- plot_grid(", reps, ",align = 'h',rel_widths = c(rep(1,num_plots),1),nrow = 1)")
      # #legs_all <- paste0("plegends <- plot_grid(", reps, ",align = 'v',axis = 'l',ncol = 1,rel_heights = c(1,1,1,1.5)")
      # print(legs_all)
      # eval(parse(text=legs_all))
      # print(plegends)
      
      if(order_by_input){
        reps <- rep(paste0("oncolist_result[[",1:length(oncolist_result),"]]$oncoplot_legend[[1]]"))
        #reps <- rep(paste0("res$oncoplot_legend[[",1:length(res$oncoplot_legend),"]]"))
        reps <- paste0(reps, collapse= " , ")
        ncol <- length(oncolist_result)
        #legs_all <- paste0("plegends <- plot_grid(", reps, ",align = 'h',axis = 'l',ncol = 1,rel_heights = c(1,1,1,1.5,1))")
        #legs_all <- paste0("plegends <- plot_grid(", reps, ",align = 'h',ncol = ncol,rel_heights = c(1,1,1,1.5,1))")
        legs_all <- paste0("plegends <- plot_grid(", reps, ",align = 'h',ncol = ncol, rel_widths=0.2)")
        eval(parse(text=legs_all))
      }else{
        reps <- rep(paste0("res$oncoplot_legend[[",1:length(res$oncoplot_legend),"]]"))
        reps <- paste0(reps, collapse= " , ")
        ncol <- length(res$oncoplot_legend)
        #legs_all <- paste0("plegends <- plot_grid(", reps, ",align = 'h',axis = 'l',ncol = 1,rel_heights = c(1,1,1,1.5,1))")
        #legs_all <- paste0("plegends <- plot_grid(", reps, ",align = 'h',ncol = ncol,rel_heights = c(1,1,1,1.5,1))")
        legs_all <- paste0("plegends <- plot_grid(", reps, ",align = 'h',ncol = ncol, rel_widths=0.2)")
        eval(parse(text=legs_all))
      }
      
      #plot_combined <- plot_grid(res$oncoplot[[1]], plegends, ncol = 1,rel_heights = c(1,0.2),align = 'v',axis = 'l')
      plot_combined <- plot_grid(oncoplots, plegends, ncol = 1,rel_heights = c(1,0.2),align = 'v',axis = 'l')
      # plot_combined <- plot_grid(oncoplots, plegends+theme(plot.margin = margin(l = -0.5,unit = 'cm')),ncol = 2,rel_widths = c(10,2))
      
      return(plot_combined)
      
    }else{
      # data <- data_all_cats
      data <- data[[2]]
      # i <- 1
      # sample_list <- list()
      # for(each in data_all_cats){
        # data_by_type <- data_[i]
        # data_by_type <- as.data.frame(data_by_type)
        name <- paste("test_result", i, sep = "")
        plot_color_vec <- landscape_colors[which(names(landscape_colors) %in% data_all_cats$Alteration)] %>% c()
        print(plot_color_vec)
        # run oncoplot function the first time to determine the gene level and sample level for each input?
        # second run of oncoplot funcion in sherlock_landscape_v2.R was mostly to add in the correct sample level (included all samples from other oncoplot results)
        res <- assign(name,oncoplot(data=data_all_cats,gene_level=NULL,sample_level=NULL, sample_level0 = sample_level0,landscape_colors = plot_color_vec, p2_axis_hidden = TRUE,p2_hidden = TRUE, GeneSortOnly=FALSE, frequency= frequency))
        # print(res$sample_level)
        tmp_sample_level <- res$sample_level %>% pull(Tumor_Barcode)

        i <- i + 1
        # print(oncolist_result)
      # }
      
      # tmp_sample_level <- sample_list %>% reduce(left_join, by= "Tumor_Barcode") %>% pull(Tumor_Barcode)

      # run again with the correct sample level across all inputs
      i <- 1
      # sample_check <- list()
      oncolist_result <- list()
      # legends <- c()
      # print(tmp_sample_level)
      for(each in data){
        data_by_type <- data[i]
        data_by_type <- as.data.frame(data_by_type)
        plot_color_vec <- landscape_colors[which(names(landscape_colors) %in% data_by_type$Alteration)] %>% c()
        plot_color_vec <- plot_color_vec[order(match(names(plot_color_vec),data_alts))]
        # print(plot_color_vec)
        name <- paste("test_result", i, sep = "")
        res <- assign(name,oncoplot(data=data_by_type,gene_level=NULL,sample_level=tmp_sample_level, sample_level0 = sample_level0,landscape_colors = plot_color_vec, p2_axis_hidden = FALSE,p2_hidden = FALSE, GeneSortOnly=FALSE, frequency= frequency))
        print(res$sample_level)
        # legends <- append(legends, res$oncoplot_legend)
        # sample_check <- append(sample_check, list(res$sample_level))
        oncolist_result <- list.append(oncolist_result, res)
        i <- i + 1
        # print(res$oncoplot_legend)
        # print(legends)
      }
      
      
      reps <- rep(paste0("oncolist_result[[",1:length(oncolist_result),"]]$oncoplot[[1]]"))
      reps <- paste0(reps, collapse= " , ")
      #nrow to adjust plot size 
      #gene_level_length <- rep(paste0("oncolist_result[[",1:length(oncolist_result),"]]$gene_level"))
      gene_level_length <- paste0("oncolist_result[[",1:length(oncolist_result),"]]$gene_level")
      i <- 1
      total_length <- 0
      for(each in gene_level_length){
        gene_level <- paste0("gene_levelb <-oncolist_result[[i]]$gene_level")
        eval(parse(text=gene_level))
        total_length <- total_length + length(gene_levelb)
        i <- i + 1
        
      }
      #gene_level_length <- paste0(gene_level_length, collapse= " , ")
      print(reps)
      onco_all <- paste0("oncoplots <- plot_grid(", reps, ",align = 'v', ncol= 1)")
      # print(onco_all)
      eval(parse(text=onco_all))
      # print(oncoplots)
      #nrow <- nrow(length(gene_level))

      reps <- rep(paste0("oncolist_result[[",1:length(oncolist_result),"]]$oncoplot_legend"))
      reps <- paste0(reps, collapse= " , ")
      print(reps)
      ncol <- length(genomic_alts)
      h2 <- length(data_alts)/25
      if(h2 < 0.2){h2 <- 0.2}
      num_plots <- length(oncolist_result)
      #legs_all <- paste0("plegends <- plot_grid(", reps, ",align = 'h',rel_widths = c(rep(1,num_plots),1),nrow = 1)")
      legs_all <- paste0("plegends <- plot_grid(", reps, ",align = 'h',ncol= ncol, rel_heights=h2,rel_widths = c(1,1,1,1),nrow = 1)")
      print(legs_all)
      eval(parse(text=legs_all))
      # print(plegends)
      
      #.2=5
      #.4=10
      #.6=25
      h2 <- length(data_alts)/25
      if(h2 < 0.2){h2 <- 0.2}
      #h2 <- if_else(length(data_alts) <=5, 0.2, if_else(length(data_alts) >5 & <10, 0.4), 0.6 )
      #plot_combined <- plot_grid(oncoplots, plegends, ncol = 1,rel_heights = c(1,h2),axis = 'l',align = 'v')
      plot_combined <- plot_grid(oncoplots, plegends, ncol = 1,align = 'v')
      return(plot_combined)
      # combined_output <- oncoplot_combined(oncolist=oncolist_result)
      # return(combined_output)
    }
    
    # plot_combined <- plot_grid(oncoplots, plegends, ncol = 1)
    
    # check if the sample names are in the same order
    # print(all(sample_check[[1]] == sample_check[[2]]))
    # print(all(sample_check[[1]][,1] == sample_check[[2]][,1]))
    # print(all(sample_check[[1]][,2] == sample_check[[2]][,2]))

  })
  
 
  # observeEvent(input$generate_oncoplot, {
    # output$data_output <- renderTable(oncoplot_reactive())
    output$data_output_plot <- renderPlot({
      # final_plot <- plot_grid(oncoplot_reactive()[[1]], oncoplot_reactive()[[2]])
      # final_plot <- plot_grid(oncoplot_reactive())
      final_plot <- oncoplot_reactive()
      # pdfhr2()
      # print(oncoplot_reactive())
      # test_result <- oncoplot(data=oncoplot_reactive(),gene_level=NULL,sample_level=sample_new_level, sample_level0 = sample_level0)
      # plot(oncoplot_reactive()$oncoplot[[1]])
      plot(final_plot)
    },bg = 'transparent') %>%  bindEvent(input$generate_oncoplot)

  # })
  
  ######### Documentation ###############################################
  
  output$module_info_table <- DT::renderDataTable({datatable(read_in_file("~/Sherlock-Genome/Module_Info.txt"), options=list(paging=FALSE))}) 
   


}

