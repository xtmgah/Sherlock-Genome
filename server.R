
source('Sherlock_Genome_Functions.R')

server <- function(input, output, session){
  
  observe({
    Sys.sleep(3)
    remove_start_up(timeout = 200)
  })
  
  updateTabItems(session, 'sbmenu', 'data_load')
  
  onStop(function() {
    setwd(parent)
    print(getwd())
    unlink("www/Genomic Data/User_project", recursive = TRUE)
    
    
    
  })
  
  session$onSessionEnded(function() {
    cat("Session Ended\n")
    project_code <- isolate(input$project_code)
    tmb_cohort_name <- isolate(input$tmb_cohort_name)
    print(file.exists(paste0('www/Genomic Data/', project_code, '/', project_code, '_', tmb_cohort_name,'_tmb.pdf')))
    
    if(file.exists(paste0('www/Genomic Data/', project_code, '/', project_code, '_', tmb_cohort_name,'_tmb.pdf'))){
      unlink(paste0('www/Genomic Data/', project_code, '/', project_code, '_', tmb_cohort_name,'_tmb.pdf'))
      unlink(paste0('www/Genomic Data/', project_code, '/Rplots.pdf'))
    }

  })
  
  # set up colors for plots
  get_vcColors()

  # disable reset until user clicks "Select Project"
  shinyjs::disable("reset_project")
  shinyjs::disable("reset_filter_df_sample")
  
  observeEvent(input$int_ext, {
    disable("int_ext")
    
    output$ui_clear_int_ext <- renderUI(actionButton("clear_int_ext", "Clear Selection", class = 'action-buttons-app'))
    
    print(getwd())
      
    if(input$int_ext == "Upload Project Data"){
      print("testVolumes")
      volumes <- c(Home = fs::path_home(), getVolumes()())
      print(volumes)
      shinyDirChoose(input, "user_choose_directory", roots = volumes, session = session, restrictions = system.file(package = "base"), allowDirCreate = FALSE)
      output$user_directory_path <- renderText({
        if (is.integer(input$user_choose_directory)) {
          cat("No directory has been selected (shinyDirChoose)")
        } else {
          paste0("User directory selected: ", parseDirPath(volumes, input$user_choose_directory))
        }
        
      })
      shinyjs::show("project_code")
      shinyjs::show("user_choose_directory")
      shinyjs::show("user_directory_path")
      shinyjs::show("choose_method_outside")
      shinyjs::show("all_mods_outside")
      shinyjs::show("ind_mods_outside")
      shinyjs::show("file_list_external")
      shinyjs::show("file_load_message_external")
      
      
    }

    if(input$int_ext == "Project in App"){
    
      shinyjs::show("project_code")
      shinyjs::show("select")
      shinyjs::disable('select')
      shinyjs::show("reset_project")
      
    }
    
  })
  
  observe({
    if(!is.null(input$project_code)){
      shinyjs::enable('select')
    }
  })
  
    observeEvent(input$clear_int_ext, {
      refresh()
      print(getwd())
      setwd('../../../')
      print(getwd())

    })
    
    observeEvent(input$user_choose_directory, {
      print(input$user_choose_directory)
      if (!is.integer(input$user_choose_directory)) {
        print(input$user_choose_directory)
        volumes <- c(Home = fs::path_home(), getVolumes()())
        output$choose_method_outside <- renderText("Choose one of the two methods below to load data into the corresponding modules:")
        output$all_mods_outside <- renderText("1. Select all modules available for the selected project. This will automatically populate each module with its corresponding data.")
        output$ind_mods_outside <- renderText("2. Select individual modules from those available for the selected project. This will automatically populate each module with its corresponding data.")
        # copy directory over to www/Genomic Data directory (current working directory)
        R.utils::copyDirectory(parseDirPath(volumes, input$user_choose_directory), "./www/Genomic Data/User_project")
        
        # check to see what is in the user folder
        setwd('www/Genomic Data/User_project/')
        dirs_in_user_folder <- list.files() # overall directories
        print(dirs_in_user_folder)
        dirs_in_user_folder <- gsub("_", " ",dirs_in_user_folder) # look into type here
        print(dirs_in_user_folder)
        print(typeof(dirs_in_user_folder))
        print(length(dirs_in_user_folder))
        choices_ext <- c("Study Overview","Manifest Information","NGSpurity","Mutations","SCNA","SV","Mutational Signatures","Genomic Landscape","Clonal Evolution","Survival Analysis","Integrative Analysis")
        print(choices_ext)
        choices_disabled <- !choices_ext %in% dirs_in_user_folder
        output$file_list_select_external <- renderUI({pickerInput("file_list_external",
                                                                  #label="Below are the modules available for the user-selected project. Select all or some of the modules available to load the corresponding data.",
                                                                  label=NULL,
                                                                  choices= choices_ext,
                                                                  multiple=TRUE, choicesOpt = list(disabled=choices_disabled,  style = ifelse(choices_disabled, yes = "color: rgba(119, 119, 119, 0.5);",
                                                                                                                                                                      no = "")), options = list(style = 'picker-inputs-app', `actions-box`= TRUE, title = 'Click here to select modules'))})#, choicesOpt = list(disabled=choices_ext[which(!(choices_ext %in% gsub("_", " ",dirs_in_user_folder)))]))})
        }
    })

  # setwd based on project code selected
  observeEvent(input$select,{ #change id to be specific to internal selection- internal_select
    os_name <- os_detect()

    print("setting project dir now")
    setwd(paste0("www/Genomic Data/",input$project_code))
    print(getwd())
    
  # disable, enable, and show different buttons/selection options
    shinyjs::disable("project_code")
    shinyjs::disable("select")
    shinyjs::show("file_list_select_internal")
    shinyjs::show("file_list_select_external")
    shinyjs::show("submit_int_ext_data_out")
  
    output$choose_method_in_app <- renderText("Choose one of the two methods below to load data into the corresponding modules:")
    output$all_mods_in_app <- renderText("1. Select all modules available for the selected project. This will automatically populate each module with its corresponding data.")
    output$ind_mods_in_app <- renderText("2. Select individual modules from those available for the selected project. This will automatically populate each module with its corresponding data.")
    
    shinyjs::show("choose_method_in_app")
    shinyjs::show("all_mods_in_app")
    shinyjs::show("ind_mods_in_app")
    
    dirs_in_internal_folder <- list.files() # overall directories
    print(dirs_in_internal_folder)
    dirs_in_internal_folder <- gsub("_", " ",dirs_in_internal_folder) # look into type here
    print(dirs_in_internal_folder)
    print(typeof(dirs_in_internal_folder))
    print(length(dirs_in_internal_folder))
    choices_int <- c("Study Overview","Manifest Information","NGSpurity","Mutations","SCNA","SV","Mutational Signatures","Genomic Landscape","Clonal Evolution","Survival Analysis","Integrative Analysis")
    print(choices_int)
    choices_disabled_int <- !choices_int %in% dirs_in_internal_folder
    output$file_list_select_internal <- renderUI({pickerInput("file_list_internal", 
                                                              #label="Below are the modules available for the user-selected project. Select all or some of the modules available to load the corresponding data.",
                                                              label=NULL,
                                                              choices= choices_int,
                                                              multiple=TRUE, choicesOpt = list(disabled=choices_disabled_int, style = ifelse(choices_disabled_int, yes = "color: rgba(119, 119, 119, 0.5);",
                                                                                                                                          no = "")),options = list(style = 'picker-inputs-app', `actions-box`= TRUE, title = 'Click here to select modules'))})
    print(paste0('file_list_internal',input$file_list_internal))
   
  
  })
  
  # only enable submit button when the user has selected at least one module #file
  observe({
    req(input$int_ext)
    if(input$int_ext=='Project in App'){
      x <- input$file_list_internal
      print(paste('x:',x))
    }
    if(input$int_ext == 'Upload Project Data'){
      x <- input$file_list_external
    }

    if(!is.null(x)){
      output$submit_int_ext_data_out <- renderUI(actionButton("submit_int_ext_data", "Submit",class = 'action-buttons-app'))

    }else{
      shinyjs::hide("submit_int_ext_data")
    }
  })
  
  observeEvent(input$submit_int_ext_data,{
    if(input$int_ext == "Project in App"){
    #if(input$sel_project_in_app){
      x <- input$file_list_internal
      x <- paste (x,sep="", collapse=", ")
      output$file_load_message_internal <- renderText(paste0("The following modules were loaded with the corresponding data: ", x[1:length(x)], "."))
      shinyjs::disable("file_list_internal")

    }
    if(input$int_ext == "Upload Project Data"){
      x <- input$file_list_external
      x <- paste (x,sep="", collapse=", ")
      output$file_load_message_external <- renderText(paste0("The following modules were loaded with the corresponding data: ", x[1:length(x)], "."))
      shinyjs::disable("file_list_external")
    }
    
    shinyjs::disable('submit_int_ext_data')
    
   
  })

  # reactive file list for both internal or external selection
  file_list_reactive <- reactive({
    req(input$submit_int_ext_data)
    if(input$int_ext == "Project in App"){
      x <- input$file_list_internal
     
    }
    if(input$int_ext == "Upload Project Data"){
      x <- input$file_list_external

    }
    
    print(x)
  })
  
######### Study Overview ###############################################
  
  # load the study overview for the selected project
  observeEvent(input$submit_int_ext_data, {
    if("Study Overview" %in% file_list_reactive()){
      study_overview <- includeMarkdown("Study_Overview/Study_Overview.Rmd")
      output$study_overview_rmd <- renderUI({includeMarkdown("Study_Overview/Study_Overview.Rmd")})}
  })

########################################################################

  df_columns_subset <- reactive({
    if(input$sbmenu == "sample_qc"){
      if(input$qc_tabs == "qc_sample"){
        dataframe= read_delim("Manifest_Information/QC_sample_level.txt")
        col_names = input$qcsample_header
      }
      if(input$qc_tabs == "qc_subject"){
        dataframe= read_delim("Manifest_Information/QC_subject_level.txt")
        col_names = input$qcsubject_header
      }
    }
    if(input$sbmenu == "NGSpurity"){
      dataframe= read_delim("NGSpurity/all_ngspurity_output.txt")
      col_names = input$ngspurity_header
    }

    return(df_columns_subset_function(dataframe, col_names))

  })

  # reactive to check input before running through filtering function if text input is provided
  check_input_reactive <- reactive({

    if(input$sbmenu == "sample_qc"){
      if(input$qc_tabs == "qc_sample"){
        data = read_delim("Manifest_Information/QC_sample_level.txt")
        inputs = input$user_filter_input_sample
      }
      if(input$qc_tabs == "qc_subject"){
        data = read_delim("Manifest_Information/QC_subject_level.txt")
        inputs = input$user_filter_input_subject
      }
    }
    if(input$sbmenu == "NGSpurity"){
      data = read_delim("NGSpurity/all_ngspurity_output.txt")
      inputs = input$user_filter_input_ngspurity
    }

    if(input$sbmenu == 'scna'){
      print('scna check check')
      if(is.null(input$project_code)){
        data = read_delim("SCNA/scna_data2.txt")
      }else{
        if(input$project_code == 'Sherlock_TCGA'){
          load("SCNA/BBprofile.RData")
          
          data = BBprofile
        }
      }
      

      inputs = input$user_filter_input_scna

    }

    return(check_inputs(data, inputs))
  })

  # reactive for dataframe filter (used in sample qc for sample data, subject data, and ngspurity)
  df_filter_reactive <- reactive({
    if(input$sbmenu == "sample_qc"){
      if(input$qc_tabs == "qc_sample"){
        data = read_delim("Manifest_Information/QC_sample_level.txt")
        conditions = input$user_filter_input_sample
        col_names = input$qcsample_header
      }
      if(input$qc_tabs == "qc_subject"){
        print(input$qc_tabs)
        data = read_delim("Manifest_Information/QC_subject_level.txt")
        conditions = input$user_filter_input_subject
        col_names = input$qcsubject_header
      }
    }

    if(input$sbmenu == "NGSpurity"){
      data = read_delim("NGSpurity/all_ngspurity_output.txt")
      conditions = input$user_filter_input_ngspurity
      col_names = input$ngspurity_header
    }

    if(input$sbmenu == 'scna'){
      if(is.null(input$project_code)){
        data <- read_delim('SCNA/scna_data2.txt')
        
      }else{
        if(input$project_code == 'Sherlock_TCGA'){
          load("SCNA/BBprofile.RData")
          data = BBprofile
        }
      }
        
      conditions = input$user_filter_input_scna
      col_names = input$scna_cols
    }

    if(input$sbmenu == "documentation"){
      if(input$documentation_tabs == "data_req_info"){
        
        #data_req_paths <- read_delim("~/Documents/Sherlock_Genome/Data_requirements.txt", delim="\t")
        data_req_paths <- read_delim(paste0(parent, '/Data_requirements.txt'), delim="\t")
        
        module_selected <- module_chosen()
        print(paste0('module_selected',module_selected))
        
        data_req_paths <- data_req_paths %>% filter(Module == module_selected)
        
        submodule_selected <- data_req_paths$Submodule[which(data_req_paths["Module"]== module_selected)]
        print(paste0('submodule_selected',submodule_selected))
        
        data_req_paths <- data_req_paths %>% filter(Submodule == submodule_selected)
        
        files_in_submodule <- data_req_paths$File[which(data_req_paths["Submodule"] == input$submodule_user_select)]
        file_paths_in_submodule <- data_req_paths$File_path[which(data_req_paths["File"] == files_in_submodule)]
        
        x <- input$file_user_select
        x_path_in_file <- data_req_paths$File_path[which(data_req_paths["File"] == x)]
        print(paste0('x_path_in_file:', x_path_in_file))
        x_path_all <- paste0(parent, '/www/Genomic Data/Data_Requirements/', x_path_in_file)
        x_path_all <- x_path_all[1]
        print(paste0('x_path_all:', x_path_all))
        x_file <- read_delim(x_path_all, delim="\t")
        
        data = x_file
        conditions = ""
        col_names = input$example_file_col_input
      }
    }

    if(conditions == ""){
      conditions = NULL
    }else{
      conditions = conditions
    }

    if("All columns" %in% col_names){
      col_names = col_names[col_names != "All columns"]
    }else{
      col_names = col_names
    }

    return(sherlock_genome_filter(data,conditions, col_names))
  })

  # reactive for dataframe inspect (used in sample qc for sample and subject data, and ngspurity)
  inspect_df_option_qc_ngs <- reactive({

    if(input$sbmenu == "sample_qc"){
      if(input$qc_tabs == "qc_sample"){

        req(input$inspect_data_type_qc_sample)
        req(input$column_name_to_inspect_qc_sample)

        dataframe = read_delim("Manifest_Information/QC_sample_level.txt")
        type_of_inspection = input$inspect_data_type_qc_sample
        column_name = input$column_name_to_inspect_qc_sample
      }

      if(input$qc_tabs == "qc_subject"){

        req(input$inspect_data_type_qc_subject)
        req(input$column_name_to_inspect_qc_subject)

        dataframe = read_delim("Manifest_Information/QC_subject_level.txt")
        type_of_inspection = input$inspect_data_type_qc_subject
        column_name = input$column_name_to_inspect_qc_subject
      }
    }

    if(input$sbmenu == "NGSpurity"){
      req(input$inspect_data_type_ngs)
      req(input$column_name_to_inspect_ngs)

      dataframe = read_delim("NGSpurity/all_ngspurity_output.txt")
      type_of_inspection = input$inspect_data_type_ngs
      column_name = input$column_name_to_inspect_ngs
    }

    return(inspect_data_function(dataframe, type_of_inspection, column_name))
  })

########### Manifest Information ###############################################

  observeEvent(input$submit_int_ext_data,{
    if("Manifest Information" %in% file_list_reactive()){
      print("Manifest Information" %in% file_list_reactive())
      data_qc_sample <- read_delim("Manifest_Information/QC_sample_level.txt")
      qc_sample_header <- colnames(data_qc_sample)
      qc_sample_header_length <- length(qc_sample_header)
      shinyjs::disable('check_input_sample')
      
      ### Sample Data ###
      
      ##### View Sample Data #####
      output$qc_sample_header_output <- renderUI({
        
        if(qc_sample_header_length < 6){
          n_col <- qc_sample_header_length
        }else{
          n_col <- 6
        }
        dropdownButton(inputId="qc_sample_dropdown",label="Select Columns",circle=FALSE, status = 'dropdownbutton', checkboxGroupInput("qcsample_header",label="Column Names",choices=c("All columns",qc_sample_header),selected=qc_sample_header[1:n_col]))})

            output$qc_sample_table_orig <- DT::renderDataTable({
        if(qc_sample_header_length < 6){
          n_col <- qc_sample_header_length
        }else{
          n_col <- 6
        }
        datatable(data_qc_sample %>% select(all_of((1:n_col))),options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})

      output$qc_sample_table<- DT::renderDataTable({datatable(df_filter_reactive(), options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))}) %>% bindEvent(input$filter_df_sample)
      
      ##### Inspect Sample Data #####
      output$inspect_sample_qc_select_column <- renderUI({selectInput("column_name_to_inspect_qc_sample","Select one column name to inspect, or all columns:", choices= c("All columns",colnames(data_qc_sample)), multiple= FALSE, selected="All columns")})
      output$inspect_sample_plot <- renderPlot({inspect_df_option_qc_ngs()}, height= 1000, width = 1000)

      ### Subject Data ###
      data_qc_subject <- read_delim("Manifest_Information/QC_subject_level.txt")
      qc_subject_header <- colnames(data_qc_subject)
      qc_subject_header_length <- length(qc_subject_header)
      shinyjs::disable('check_input_subject')

      ##### View Subject Data #####
      output$qc_subject_header <- renderUI({
        if(qc_subject_header_length < 6){
          n_col <- qc_subject_header_length
        }else{
          n_col <- 6
        }
        dropdownButton(inputId="qc_subject_dropdown",label="Select Columns",circle=FALSE,status = 'dropdownbutton',checkboxGroupInput("qcsubject_header",label="Column Names",choices=c("All columns",qc_subject_header),selected=qc_subject_header[1:n_col]))})

      output$qc_subject_table_orig <- DT::renderDataTable({
        if(qc_subject_header_length < 6){
          n_col <- qc_subject_header_length
        }else{
          n_col <- 6
        }
        datatable(data_qc_subject %>% select(all_of((1:n_col))),options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})

      output$qc_subject_table <- DT::renderDataTable({datatable(df_filter_reactive(),options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))}) %>% bindEvent(input$filter_df_subject)

      ##### Inspect Subject Data #####
      output$inspect_subject_qc_select_column <- renderUI({selectInput("column_name_to_inspect_qc_subject","Select one column name to inspect, or all columns:", choices= c("All columns",colnames(data_qc_subject)), multiple= FALSE, selected="All columns")})
      output$inspect_subject_plot <- renderPlot({inspect_df_option_qc_ngs()}, height= 1000, width = 1000)
    }
  })

  observeEvent(input$user_filter_input_sample, {
    if(input$user_filter_input_sample == ''){
      shinyjs::disable('check_input_sample')
      shinyjs::enable('filter_df_sample')
      output$input_check_result_sample <- renderText({''})
    }else{
      shinyjs::enable('check_input_sample')
      shinyjs::disable('filter_df_sample')
    }

  })

 observeEvent(input$check_input_sample, {

  output$input_check_result_sample <- renderText({isolate(check_input_reactive())}) #%>% bindEvent(input$check_input)
  print(paste0('check_input_reactive:',check_input_reactive()))
  if(check_input_reactive() == 'All variables are found in the data.'){
    shinyjs::enable('filter_df_sample')
  } else{
    shinyjs::disable('filter_df_sample')
  }

  })

  observeEvent(input$qcsample_header, {
    if("All columns" %in% input$qcsample_header){
      data_qc_sample <- read_delim("Manifest_Information/QC_sample_level.txt")
      qc_sample_header <-colnames(data_qc_sample)
      updateCheckboxGroupInput(session, "qcsample_header", label="Column Names",choices=c("All columns",colnames(data_qc_sample)),selected=c("All columns",qc_sample_header[1:length(qc_sample_header)]))
    }
  })

  observeEvent(input$filter_df_sample,{
    output$input_check_result_sample <- renderText({''})

  })

  observeEvent(input$user_filter_input_subject, {
    if(input$user_filter_input_subject == ''){
      shinyjs::disable('check_input_subject')
      shinyjs::enable('filter_df_subject')
      output$input_check_result_subject <- renderText({''})
    }else{
      shinyjs::enable('check_input_subject')
      shinyjs::disable('filter_df_subject')
    }

  })

  observeEvent(input$check_input_subject, {

    output$input_check_result_subject <- renderText({isolate(check_input_reactive())}) #%>% bindEvent(input$check_input)
    print(paste0('check_input_reactive:',check_input_reactive()))
    if(check_input_reactive() == 'All variables are found in the data.'){
      shinyjs::enable('filter_df_subject')
    } else{
      shinyjs::disable('filter_df_subject')
    }

  })

  observeEvent(input$qcsubject_header, {
    if("All columns" %in% input$qcsubject_header){
      data_qc_subject <- read_delim("Manifest_Information/QC_subject_level.txt")
      qc_subject_header <- colnames(data_qc_subject)
      updateCheckboxGroupInput(session, "qcsubject_header", label="Column Names",choices=c("All columns",colnames(data_qc_subject)),selected=c("All columns",qc_subject_header[1:length(qc_subject_header)]))
    }
  })

  observeEvent(input$filter_df_subject,{
    output$input_check_result_subject <- renderText({''})

  })

  data_to_inspect <- reactive({

    if(input$sbmenu == "sample_qc"){
      if(input$qc_tabs == "qc_sample"){
        data <- read_delim("Manifest_Information/QC_sample_level.txt")
      }
      if(input$qc_tabs == "qc_subject"){
        data <- read_delim("Manifest_Information/QC_subject_level.txt")
      }
    }
    if(input$sbmenu== "NGSpurity"){
      data <- read_delim("NGSpurity/all_ngspurity_output.txt")
    }

    data <- data %>% validate_vardf()

  })

  col_to_inspect <- reactive({
    if(input$sbmenu == "sample_qc"){
      if(input$qc_tabs == "qc_sample"){
        req(input$column_name_to_inspect_qc_sample)
        column_name_tmp <- input$column_name_to_inspect_qc_sample
      }
      if(input$qc_tabs == "qc_subject"){
        req(input$column_name_to_inspect_qc_subject)
        column_name_tmp <- input$column_name_to_inspect_qc_subject
      }
    }
    if(input$sbmenu== "NGSpurity"){
      req(input$column_name_to_inspect_ngs)
      column_name_tmp <- input$column_name_to_inspect_ngs
    }

    tmp <- data_to_inspect()

    if(column_name_tmp == 'All columns'){
      choices <- c('cat','cat_levels','na','num','types')
    }else{
      tmp2 <- tmp %>% select(all_of(column_name_tmp))

      if("character" %in% sapply(tmp2, class)[1]){
        choices <- c('cat','cat_levels','na','types')
      }
      if("numeric" %in% sapply(tmp2, class)[1]){
        choices <- c('na','num','types')
      }
    }
    print(choices)
  })

  observeEvent(input$column_name_to_inspect_qc_sample,{

    updatePickerInput(session,"inspect_data_type_qc_sample", "Select one of the options below to inspect the data:",choices=col_to_inspect())

  })

  observeEvent(input$column_name_to_inspect_qc_subject,{

    updatePickerInput(session,"inspect_data_type_qc_subject", "Select one of the options below to inspect the data:",choices=col_to_inspect())

  })

########### NGSpurity ###############################################

  # reactive to call figure_display() function above
  # figure_output <- reactive({
  #   req(input$tumor_barcode_to_inspect)
  #   req(input$battenberg_to_inspect)
  #   req(input$type_to_inspect)
  #   tumor_barcode = input$tumor_barcode_to_inspect
  #   battenberg = input$battenberg_to_inspect
  #   type = input$type_to_inspect
  #   project_code = input$project_code
  #   if(os_detect() %in% c("Linux","Darwin")){
  #     return(figure_display_ngspurity(tumor_barcode, battenberg,type,project_code))
  #   }else{
  #     return(src=figure_display_ngspurity(tumor_barcode, battenberg,type,project_code))
  #   }
  # })
  
  figure_output_ngspurity <- reactive({
    ngspurity_qc <- read_delim("NGSpurity/ngspurity_qc_file.txt")
    tumor_barcode = input$tumor_barcode_to_inspect
    battenberg = input$battenberg_to_inspect
    type = input$type_to_inspect

    ngspurity_qc_figure <- ngspurity_qc %>% filter(Tumor_Barcode == input$tumor_barcode_to_inspect) %>% filter(Battenberg == input$battenberg_to_inspect) %>% filter(Type == input$type_to_inspect) %>% pull(File)
    print(ngspurity_qc_figure)
    
    filename_ngs <- sub(".","NGSpurity", ngspurity_qc_figure)
    # print(filename)
    
  })

  observeEvent(input$submit_int_ext_data,{
    if("NGSpurity" %in% file_list_reactive()){
      data_ngs <- read_delim("NGSpurity/all_ngspurity_output.txt")
      ngs_purity_header <- colnames(data_ngs)
      ngs_purity_header_length <- length(ngs_purity_header)
      shinyjs::disable('check_input_ngspurity')

      ##### View Data QC #####

      output$ngs_purity_header <- renderUI({
        if(ngs_purity_header_length < 6){
          n_col <- ngs_purity_header_length
        }else{
          n_col <- 6
        }
        dropdownButton(inputId="ngspurity_dropdown",label="Select Columns",circle=FALSE,status = 'dropdownbutton', checkboxGroupInput("ngspurity_header",label="Column Names",choices=c("All columns",ngs_purity_header),selected=ngs_purity_header[1:n_col]))})
      
      output$qc_ngspurity_table_orig <- DT::renderDataTable({
        
        if(ngs_purity_header_length < 6){
          n_col <- ngs_purity_header_length
        }else{
          n_col <- 6
        }
        datatable(data_ngs %>% select(all_of((1:n_col))),options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})
      
      output$qc_ngspurity_table <- DT::renderDataTable({datatable(df_filter_reactive(),options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))}) %>% bindEvent(input$filter_df_ngspurity)

      ###### Inspect Data #####
      output$ngs_purity_header_inspect_tab <- renderUI({selectInput("column_name_to_inspect_ngs","Select one column name to inspect, or all columns:", choices= c("All columns",ngs_purity_header), multiple= FALSE, selected="All columns")})
      output$inspect_ngs_plot <- renderPlot({inspect_df_option_qc_ngs()}, height= 1000, width = 1000)

      ##### View Figures #####
      os_name <- os_detect()
      if(os_name == "Windows"){
        ngspurity_qc <- read_delim("NGSpurity/ngspurity_qc_file.txt")
        i <- 1
        for(each in ngspurity_qc$File){
          ngspurity_qc$File[i] <- paste0("C:/", ngspurity_qc$File[i])
          i <- i + 1
        }
      }else{
        ngspurity_qc <- read_delim("NGSpurity/ngspurity_qc_file.txt")
      }

      output$ngspurity_type <- renderUI({selectInput("type_to_inspect","Select one Type to inspect:", choices= unique(ngspurity_qc$Type[ngspurity_qc$Tumor_Barcode==input$tumor_barcode_to_inspect & ngspurity_qc$Battenberg == input$battenberg_to_inspect]), multiple= FALSE)})

      updateSelectInput(session, 'tumor_barcode_to_inspect', label = "Select one Tumor Barcode to inspect:", choices= unique(ngspurity_qc$Tumor_Barcode))

      
      observe({
        x <- input$tumor_barcode_to_inspect

        updateSelectInput(session, "battenberg_to_inspect","Select one Battenberg to inspect:", choices= unique(ngspurity_qc$Battenberg[ngspurity_qc$Tumor_Barcode==x]))

      })

      observe({

        x <- input$battenberg_to_inspect
        updateSelectInput(session, 'type_to_inspect',"Select one Type to inspect:", choices= unique(ngspurity_qc$Type[ngspurity_qc$Tumor_Barcode==input$tumor_barcode_to_inspect & ngspurity_qc$Battenberg == x]))

      })

      # if(os_detect() %in% c("Linux","Darwin")){
      #   output$figure_pdf <-  renderImage({figure_output()},deleteFile=FALSE)
      # }else{
      #   output$figure_pdf <- renderUI({ tags$iframe(style="height:1000px; width:100%", src= figure_output())})
      # }
      
      output$figure_ngspurity <- renderImage({list(src=figure_output_ngspurity(),width="800px", height = "auto", alt=figure_output_ngspurity())}, deleteFile = FALSE)


    }
  })

  observeEvent(input$user_filter_input_ngspurity, {
    if(input$user_filter_input_ngspurity == ''){
      shinyjs::disable('check_input_ngspurity')
      shinyjs::enable('filter_df_ngspurity')
      output$input_check_result_ngspurity <- renderText({''})
    }else{
      shinyjs::enable('check_input_ngspurity')
      shinyjs::disable('filter_df_ngspurity')
    }

  })

  observeEvent(input$check_input_ngspurity, {

    print('trying to check input for ngspurity')
    output$input_check_result_ngspurity <- renderText({isolate(check_input_reactive())}) #%>% bindEvent(input$check_input)
    print(paste0('check_input_reactive:',check_input_reactive()))
    if(check_input_reactive() == 'All variables are found in the data.'){
      shinyjs::enable('filter_df_ngspurity')
    } else{
      shinyjs::disable('filter_df_ngspurity')
    }

  })

  observeEvent(input$ngspurity_header, {
    if("All columns" %in% input$ngspurity_header){
      data_ngs <- read_delim("NGSpurity/all_ngspurity_output.txt")
      ngs_purity_header <- colnames(data_ngs)
      updateCheckboxGroupInput(session, "ngspurity_header", label="Column Names",choices=c("All columns",ngs_purity_header),selected=c("All columns",ngs_purity_header[1:length(ngs_purity_header)]))
    }
  })


  observeEvent(input$filter_df_ngspurity,{
    output$input_check_result_ngspurity <- renderText({''}) #%>% bindEvent(input$check_input)

  })

  observeEvent(input$column_name_to_inspect_ngs,{

    updatePickerInput(session,"inspect_data_type_ngs", "Select one of the options below to inspect the data:",choices=col_to_inspect())

  })

########## Mutations ###############################################

#### Mutation Summary ####
  mut_summary_reactive <- reactive({
    if(is.null(input$project_code)){
      data_manifest <- read_delim('Mutations/mutations_manifest.txt')
      data_maf <- read_delim('Mutations/mut_summary_maf_data.txt')

      if(is.null(input$mut_summary)){
        sample_group = NULL
        data_maf = data_maf
      }else{
        sample_group = input$mut_summary
        data_maf = data_maf %>% left_join(data_manifest, by = 'Tumor_Barcode') %>% filter(SP_Group == sample_group)
      }
    }else{
      if(input$project_code == 'Sherlock_TCGA'){
        load('Mutations/sherlock_maf.RData')
        if(is.null(input$mut_summary)){
          sample_group = NULL
          data_maf = sherlock_maf
        }else{
          sample_group = input$mut_summary
          data_maf = sherlock_maf %>% left_join(sherlock_overall, by = 'Tumor_Barcode') %>% filter(SP_Group == sample_group)
        }
      }
    }


    return(mutation_summary(data_maf, sample_group))
  })


  observeEvent(input$submit_int_ext_data,{



    if("Mutations" %in% file_list_reactive()){
      if(is.null(input$project_code)){
        data_manifest <- read_delim('Mutations/mutations_manifest.txt')
        print(unique(data_manifest$SP_Group))
        updateSelectInput(session,"mut_summary",choices= sort(unique(data_manifest$SP_Group)))
      }else{
        if(input$project_code == "Sherlock_TCGA"){
          updateSelectInput(session,"mut_summary",choices= sort(unique(sherlock_overall$SP_Group)))
        }
      }

    }
  })

  output$mut_summary_plot <- renderPlot({
    withProgress(message = 'Tabulating mutations',{

      mut_summary_reactive()
    })
    }, height = 500) %>% bindEvent(input$mut_summary_generate)
  
    
  output$download_mut_summary = downloadHandler(
    
    filename = function() {
      if(is.null(input$project_code)){
        project_code <- 'User_Project'
      }else{
        project_code <- input$project_code
      }
      if(!is.null(input$mut_summary)){
        group_list <- paste0(input$mut_summary, collapse = '_')
        filename = paste0(project_code,'_', group_list, '_MutationSummary.pdf')
      }else{
        filename = paste0(project_code, '_MutationSummary.pdf')
      }
      
      return(filename)
    },
    
    content = function(file) {
        ggsave(file, plot = mut_summary_reactive(), width = 15, height = 10, device = cairo_pdf)
      
    })

#### TMB Plot ####

  tmb_plot_reactive <- reactive({
    # need to work out where to store files needed in modules regardless of project
    tcga_cohorts_file <- paste0(parent, '/www/Genomic Data/Sherlock_TCGA/Mutations/tcga.cohort.RData')
    
    if(is.null(input$project_code)){
      data <- read_delim('Mutations/tmb_data.txt')
      project_code <- 'User_Project'
    }else{
      project_code <- input$project_code
      if(input$project_code == "Sherlock_TCGA"){
        load('Mutations/tmb_data.RData')
        data <- tmb_data
      }
    }
   
    if(!input$tmb_sp_group){
      maf.mutload1 <- data %>% dplyr::select(Tumor_Sample_Barcode=Tumor_Barcode,total=TMB)
    }else{
      maf.mutload1 <- data %>% filter(SP_Group == input$tmb_group) %>% select(Tumor_Sample_Barcode=Tumor_Barcode,total=TMB)
    }
    maf.mutload <- list(maf.mutload1)
    cohortsize = c(dim(maf.mutload1)[1])
    print(input$tmb_cohort_name)
    if(input$tmb_cohort_name == ""){
      cohortName = 'Cohort'
    }else{
      cohortName = input$tmb_cohort_name
    }
    
    return(tcgaCompareWGS2(project_code, maf.mutload = maf.mutload,cohortName=cohortName,col = c('gray70','darkblue'),capture_size = 1,bg_col = c("white","gray80"),tcga_cohorts_file = tcga_cohorts_file ,cohortsize=cohortsize))
  
  })

  figure_output_tmb <- reactive({
    if(is.null(input$project_code)){
      project_code = 'User_Project'
    }else{
      project_code = input$project_code
    }
    cohortName = input$tmb_cohort_name
    filename = paste0(project_code, '_', cohortName,'_tmb.pdf')
    if(os_detect() %in% c("Linux","Darwin")){
      return(figure_display_tmb(filename))
    }else{
      return(src=figure_display_tmb(filename))
    }
  })

  observeEvent(input$submit_int_ext_data,{
    if("Mutations" %in% file_list_reactive()){
      if(is.null(input$project_code)){
        data_manifest <- read_delim('Mutations/mutations_manifest.txt')
        print(unique(data_manifest$SP_Group))
        updateSelectInput(session,"tmb_group",choices= sort(unique(data_manifest$SP_Group)))
      }else{
        if(input$project_code == "Sherlock_TCGA"){
          updateSelectInput(session,"tmb_group",choices= sort(unique(sherlock_overall$SP_Group)))
        }
      }
      
    }

  })

  observeEvent(input$tmb_plot_generate,{
    print('here now')
    
    tmb_plot_reactive()
    if(os_detect() %in% c("Linux","Darwin")){
      output$figure_pdf_tmb <-  renderImage({figure_output_tmb()},deleteFile=FALSE)
    }else{
      output$figure_pdf_tmb <- renderUI({ tags$iframe(style="height:1000px; width:60%", src= figure_output_tmb())})
    }

  })
  
  output$download_tmb_plot = downloadHandler(

    filename = function() {
      if(is.null(input$project_code)){
        project_code <- 'User_Project'
      }else{
        project_code <- input$project_code
      }
      cohortName = input$tmb_cohort_name

      filename = paste0(project_code, '_', cohortName,'_tmb.pdf')

      return(filename)
    },

    content = function(theFile) {
      tmb_plot_reactive()
      unlink("./Rplots.pdf", recursive = FALSE)
      if(is.null(input$project_code)){
        project_code <- 'User_Project'
      }else{
        project_code <- input$project_code
      }
      
      cohortName = input$tmb_cohort_name

      filename = paste0(project_code, '_', cohortName,'_tmb.pdf')

      file.copy(from =filename, to = theFile)
      file.remove(paste0('./',filename))
      
    })



#### Lollipop Plot ####
  
  check_gene_reactive <- reactive({
    if(is.null(input$project_code)){
      data = read_delim('Mutations/lollipop_input_data.txt')
      gene = input$lolli_gene_name
    }else{
      if(input$project_code == 'Sherlock_TCGA'){
        load('Mutations/sherlock_maf.RData')
        data = sherlock_maf
        gene = input$lolli_gene_name
      }
    }
   
    
    return(check_gene(data, gene))
    
  })

  check_gene_message_output_reactive <- reactive({
    
    
    if(input$lolli_gene_name != ''){
      out_message <- check_gene_reactive()
      if(check_gene_reactive() == 'Gene entered exists in the data.'){
        shinyjs::enable("get_ts_info")
        #shinyjs::enable("lolliplot_calculate")
      }else{
        shinyjs::disable("get_ts_info")
      }
    }
    if(input$lolli_gene_name == ''){
      out_message <- ''
      shinyjs::disable("get_ts_info")
      shinyjs::disable("lolliplot_calculate")
    }
    print(out_message)
  })
  
  output$check_gene_input <- renderText({check_gene_message_output_reactive()})
  
  
   lolliplot_setup_reactive <- reactive({
    gene = input$lolli_gene_name
    if(input$lolli_sp_group){
      group = input$lolli_group
    }else{
      group = NULL
    }
    
    if(is.null(input$project_code)){
      data_overall <- read_delim('Mutations/mutations_manifest.txt')
      data_maf <- read_delim('Mutations/lollipop_input_data.txt')
    }else{
      if(input$project_code == "Sherlock_TCGA"){
        data_overall = sherlock_overall
        load('Mutations/sherlock_maf.RData')
        data_maf = sherlock_maf
      }
    }
    

    return(lolliplot_setup(gene, group, data_overall, data_maf))

   })

   lolliplot_plot_reactive <- reactive({
     gene = input$lolli_gene_name
     tdata0 =  lolliplot_setup_reactive()[[1]]
     tslist = lolliplot_setup_reactive()[[2]] %>% unlist()
     print(input$tslist_menu)
     tslist_input = input$tslist_menu
     samplelist = lolliplot_setup_reactive()[[3]]
     minN = as.numeric(input$lolli_minN)
     domain_annotation = input$lolli_domain_annot


     return(lolliplot_plot(gene, tdata0, tslist, tslist_input,samplelist, minN, domain_annotation))
   })

   observeEvent(input$submit_int_ext_data,{
     if("Mutations" %in% file_list_reactive()){
       if(is.null(input$project_code)){
         data_manifest <- read_delim('Mutations/mutations_manifest.txt')
         updateSelectInput(session,"lolli_group",choices= sort(unique(data_manifest$SP_Group)))
       }else{
         if(input$project_code == "Sherlock_TCGA"){
           updateSelectInput(session,"lolli_group",choices= sort(unique(sherlock_overall$SP_Group)))
         }
       }
       shinyjs::disable("get_ts_info")
       shinyjs::disable("lolliplot_calculate")
       shinyjs::disable("lolliplot_reset")
       
     }

   })
   
   observeEvent(input$get_ts_info,{
     choices <- lolliplot_setup_reactive()[[2]] %>% unlist()
     output$tslist_input <- renderUI({selectInput("tslist_menu", "Select a transcript from the list below:", choices= choices, selected=NULL, multiple= FALSE)})
     shinyjs::enable("lolliplot_calculate")
   })

  output$lolliplot <- renderPlot(lolliplot_plot_reactive()[[1]]) %>% bindEvent(input$lolliplot_calculate)
  output$lolliplot_table <-DT::renderDataTable({datatable(lolliplot_plot_reactive()[[2]], options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))}) %>% bindEvent(input$lolliplot_calculate)


  output$download_lolliplot = downloadHandler(
    
    filename = function() {
      if(is.null(input$project_code)){
        project_code <- 'User_Project'
      }else{
        project_code <- input$project_code
      }
      if(input$lolli_sp_group){
        group = input$lolli_group
        filename = paste0(project_code, '_', group, '_', input$lolli_gene_name, '_',input$tslist_menu,'_lollipop.pdf')
      }else{
        group = NULL
        filename = paste0(project_code, '_', input$lolli_gene_name, '_',input$tslist_menu,'_lollipop.pdf')
        
      }
      
      print(filename)
      return(filename)
    },
    
    content = function(file) {
      ggsave(file,plot = lolliplot_plot_reactive()[[1]],width = 20,height = 10,device = cairo_pdf)
      
      
    })
  
  observeEvent(input$lolliplot_calculate, {
    shinyjs::enable("lolliplot_reset")
    shinyjs::disable('lolli_gene_name')
    shinyjs::disable("lolli_sp_group")
    shinyjs::disable("lolli_group")
    shinyjs::disable('get_ts_info')
    shinyjs::disable('tslist_menu')
    shinyjs::disable('lolli_minN')
    
  })

  observeEvent(input$lolliplot_reset,{
    reset('lolli_gene_name')
    shinyjs::enable('lolli_gene_name')
    reset('lolli_sp_group')
    shinyjs::enable("lolli_sp_group")
    reset('lolli_group')
    shinyjs::enable("lolli_group")
    reset('get_ts_info')
    shinyjs::hide('tslist_menu')
    shinyjs::enable('tslist_menu')
    shinyjs::enable('lolli_minN')
    
  })

######### SCNA ###############################################

  scna_clustering_reactive_part_one <- reactive({
    #req(input$scna_cluster_calculate)
    sample_input = input$scna_samples
    if(is.null(input$project_code)){
      data1 <- read_delim('SCNA/scna_data1.txt')
      data2 <- read_delim('SCNA/scna_data2.txt')
      
      data1 <- data1 %>% filter(Tumor_Barcode %in% sample_input)
      data2 <- data2 %>% filter(Tumor_Barcode %in% sample_input)
      
   
    }else{
      if(input$project_code == 'Sherlock_TCGA'){

        load(paste0(parent,'/www/Genomic Data/Sherlock_TCGA/SCNA/bb_heatmap_inputs.RData'))

        data1 <- BBsolution4 %>% filter(Tumor_Barcode %in% sample_input)

        data2 <- BBprofile %>% filter(Tumor_Barcode %in% data1$Tumor_Barcode)

      }
    }
   


    return(scna_clustering_part_one(data1 = data1, data2 = data2))

  })

  scna_clustering_reactive_part_two <- reactive({

    cnvdata = scna_clustering_reactive_part_one()[[1]]
    p_dend = scna_clustering_reactive_part_one()[[2]]
    sample_order = scna_clustering_reactive_part_one()[[3]]
    totalsamples = scna_clustering_reactive_part_one()[[5]]
    cnseg = scna_clustering_reactive_part_one()[[6]]

    load(paste0(parent,'/hg38centro_cyto_info.RData'))
    print(head(cyto))
    print(head(hg38centro))
    return(scna_clustering_part_two(cnvdata, p_dend, sample_order, totalsamples,cnseg, hg38centro, cyto))
  })


  figure_display_scna_gistic_output <- reactive({

    if(input$gistic_out_select == "Amplification Plot"){
      amp_or_del = 'amp'
    }else{
      amp_or_del = 'del'
    }


    #return(figure_display_scna_gistic(amp_or_del))

  })


  observeEvent(input$submit_int_ext_data,{
    if("SCNA" %in% file_list_reactive()){
      
      if(is.null(input$project_code)){
        scna_data <- read_delim('SCNA/scna_data2.txt')
        updateSelectizeInput(session, 'scna_samples',choices = scna_data %>% filter(startsWith(Tumor_Barcode, 'TCGA')) %>% pull(Tumor_Barcode) %>% unique(), selected = NULL)
        
      }else{
        if(input$project_code == 'Sherlock_TCGA'){


          load(paste0(parent,'/www/Genomic Data/Sherlock_TCGA/SCNA/bb_heatmap_inputs.RData'))

          scna_data <- BBprofile
          updateSelectizeInput(session, 'scna_samples',choices = scna_data %>% filter(startsWith(Tumor_Barcode, 'TCGA')) %>% pull(Tumor_Barcode) %>% unique(), selected = NULL)
          
        }
      }
      
        scna_data <- validate_vardf(scna_data)
      

        output$scna_header_output <- renderUI({dropdownButton(inputId="scna_dropdown",label="Select Columns",circle=FALSE, status = 'dropdownbutton', checkboxGroupInput("scna_cols",label="Column Names",choices=c("All columns",colnames(scna_data)), selected = colnames(scna_data)[1:6]))})

        output$scna_table_orig <- DT::renderDataTable({datatable(scna_data %>% select(1:6),options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})

        output$scna_table<- DT::renderDataTable({datatable(df_filter_reactive(), options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))}) %>% bindEvent(input$filter_scna_table)

        # gistic output
        # if(os_detect() %in% c("Linux","Darwin")){
        #   output$scna_gistic_plot <-  renderImage({figure_display_scna_gistic_output()},deleteFile=FALSE)
        # }else{
        #   output$scna_gistic_plot <- renderUI({ tags$iframe(style="height:1000px; width:60%", src= figure_display_scna_gistic_output())})
        # }
        
        output$scna_gistic_plot <- renderImage({list(src=paste0("SCNA/",figure_display_scna_gistic_output(), "_qplot.png"),width="auto", height = "1000px", alt=paste0(figure_display_scna_gistic_output(), '_qplot.png'))}, deleteFile = FALSE)
        
      }
    
  })

  scna_out <- eventReactive(input$scna_cluster_calculate,{
    table_out <- scna_clustering_reactive_part_one()[[4]]
    plot_out <- scna_clustering_reactive_part_two()
    
    return(list(table_out, plot_out))
  })

    output$scna_clustering_table <- DT::renderDataTable({
      withProgress(message = 'Running cluster analysis', {
      
        datatable(scna_out()[[1]], options = list(searchHighlight = TRUE), filter = list(position = 'top', clear=TRUE, plain=FALSE))
    
    })})
    
    # cluster heatmap
    output$scna_clustering_heatmap <- renderPlot({

        scna_out()[[2]]
    })

  
  output$download_scna_clustering_plot = downloadHandler(
    
    filename = function() {
      if(is.null(input$project_code)){
        project_code <- 'User_Project'
      }else{
        project_code <- input$project_code
      }
      
      filename = paste0(project_code, '_scna_clustering.pdf')
      
      return(filename)
    },
    
    content = function(file) {
      ggsave(file, plot = scna_clustering_reactive_part_two(), width = 15, height = 20, device = cairo_pdf)
      
    })
  
  observeEvent(input$scna_cols, {
    if("All columns" %in% input$scna_cols){

      load("./www/Genomic Data/Sherlock_TCGA/SCNA/BBprofile.RData")
      scna_data <- BBprofile
      updateCheckboxGroupInput(session, "scna_cols", label="Column Names",choices=c("All columns",colnames(scna_data)),selected=c("All columns",colnames(scna_data)))
    }
  })

  observeEvent(input$user_filter_input_scna, {
    if(input$user_filter_input_scna == ''){
      shinyjs::disable('check_input_scna')
      shinyjs::enable('filter_scna_table')
      output$input_check_result_scna <- renderText({''})
    }else{
      shinyjs::enable('check_input_scna')
      shinyjs::disable('filter_scna_table')
    }

  })

  observeEvent(input$check_input_scna, {

    print('trying to check input for scna')
    output$input_check_result_scna <- renderText({isolate(check_input_reactive())}) #%>% bindEvent(input$check_input)
    print(paste0('check_input_reactive:',check_input_reactive()))
    if(check_input_reactive() == 'All variables are found in the data.'){
      shinyjs::enable('filter_scna_table')
    } else{
      shinyjs::disable('filter_scna_table')
    }

  })


########## SV #################################################################
  
  recon_plot_reactive <- reactive({
    
    sv_data <- read_delim('SV/sv_data.txt')
    cn_data <- read_delim('SV/cn_data.txt')
    
    sv_data <- sv_data %>% filter(chr1 %in% c(1:22, 'X','Y'))
    sv_data <- sv_data %>% filter(chr2 %in% c(1:22, 'X', 'Y'))
    
    if(!str_detect(sv_data$chr1[1],'chr')){
      sv_data <- sv_data %>% mutate(chr1=paste0('chr',chr1),chr2=paste0('chr',chr2))
    }
    
    
    if(!str_detect(cn_data$chr[1],'chr')){
      cn_data <- cn_data %>% mutate(chr=paste0('chr',chr))
    }
    
    barcode = input$recon_tumor_barcode
    genome_build = input$recon_genome_build
    chrs = input$recon_chr_select
    print(paste0('chrs:', chrs))
    chrs_start = input$recon_chr_start
    chrs_end = input$recon_chr_end
    
    return(sv_recon_plot(sv_data = sv_data, cn_data = cn_data, barcode, genome_build, chrs, chrs_start, chrs_end))
  })
  
  observeEvent(input$submit_int_ext_data,{
    if("SV" %in% file_list_reactive()){
      sv_data <- read_delim('SV/sv_data.txt')
      cn_data <- read_delim('SV/cn_data.txt')
      
      # tumor barcode list
      tumor_barcode_sv <- unique(sv_data$Tumor_Barcode)
      tumor_barcode_cn <- unique(cn_data$Tumor_Barcode)
      
      tumor_barcode_input <- unique(c(tumor_barcode_sv, tumor_barcode_cn))
      updatePickerInput(session, "recon_tumor_barcode", label = "Select one Tumor Barcode:", choices= tumor_barcode_input)
      
      # chr list
      sv_data <- sv_data %>% filter(chr1 %in% c(1:22, 'X','Y'))
      sv_data <- sv_data %>% filter(chr2 %in% c(1:22, 'X', 'Y'))
      
      if(!str_detect(sv_data$chr1[1],'chr')){
        sv_data <- sv_data %>% mutate(chr1=paste0('chr',chr1),chr2=paste0('chr',chr2))
      }
      
      
      if(!str_detect(cn_data$chr[1],'chr')){
        cn_data <- cn_data %>% mutate(chr=paste0('chr',chr))
      }
      
      chr_sv <- unique(sv_data$chr1)
      chr_cn <- unique(cn_data$chr)
      
      chr_input <- unique(c(chr_sv, chr_cn))
      updatePickerInput(session, "recon_chr_select", label = "Select one chromosome:", choices= chr_input)
      
    }
  })
  
  observeEvent(input$generate_recon_plot,{
    
    output$recon_plot <- renderPlot({isolate(recon_plot_reactive())})
    
  })
  
  output$download_sv_recon_plot = downloadHandler(
    
    filename = function() {
      if(is.null(input$project_code)){
        project_code <- 'User_Project'
      }else{
        project_code <- input$project_code
      }
    
      filename = paste0(project_code, '_', input$recon_tumor_barcode, '_', input$recon_genome_build, '_',input$recon_chr_select,'_', input$recon_chr_start, '_', input$recon_chr_end, '_reconPlot.pdf')
      print(filename)
      return(filename)
    },
    
    content = function(file) {
      ggsave(file,plot = recon_plot_reactive(),width = 25,height = 20,device = cairo_pdf)
      
      
    })
  
######### Mutational Signatures ###############################################

  observeEvent(input$submit_int_ext_data,{
    if("Mutational Signatures" %in% file_list_reactive()){
      if(is.null(input$project_code)){
        print(paste0('mut_sigs_files:',list.files('Mutational_Signatures/Clustered_Mutations')))
        file_list <- list.files('Mutational_Signatures/Clustered_Mutations')
        print(file_list)
        barcodes <- c()
        for(each in file_list){
          print(each)
          barcode <- str_split(each, "_")[[1]][1]
          barcodes <- append(barcodes, barcode)
        }
        
        print(paste0('barcodes:',barcodes))
        updateSelectizeInput(session, 'clustered_mut_barcodes', choices = barcodes, server = TRUE)
      }else{
        if(input$project_code == "Sherlock_TCGA"){
          updateSelectizeInput(session, 'clustered_mut_barcodes', choices = unique(sherlock_data_full$Tumor_Barcode), server = TRUE)
        }
      }
      
      output$figure_pdf_clustered_mut <-  renderImage({list(src=paste0("Mutational_Signatures/Clustered_Mutations/",input$clustered_mut_barcodes, "_Mutation_Clustering.png"),width="1500px", height = "auto", alt=paste0("Tumor Barcode: ",input$tumor_barcode_to_inspect_genomePlot))}, deleteFile = FALSE)
      #output$genomePlot_figure <- renderImage({list(src=paste0("Genomic_Landscape/",tumor_barcode_to_inspect_genomePlot_reactive()),width="800px", height = "auto", alt=paste0("Tumor Barcode: ",input$tumor_barcode_to_inspect_genomePlot))}, deleteFile = FALSE)
    }
  })
  
  # figure_output_clustered_mut <- reactive({
  #   req(input$clustered_mut_barcodes)
  #   tumor_barcode = input$clustered_mut_barcodes
  #   if(!is.null(input$project_code)){ # project included in the app
  #     project_code <- input$project_code
  #   }else{ # user project
  #     project_code <- 'User_project'
  #   }
  #   if(os_detect() %in% c("Linux","Darwin")){
  #     return(figure_display_clustered_mut(tumor_barcode, project_code))
  #   }else{
  #     return(src=figure_display_clustered_mut(tumor_barcode, project_code))
  #   }
  # 
  # })

  # observeEvent(input$submit_int_ext_data,{
  #   if("Mutational Signatures" %in% file_list_reactive()){
  #     if(is.null(input$project_code)){
  #         print(paste0('mut_sigs_files:',list.files('Mutational_Signatures/Clustered_Mutations')))
  #         file_list <- list.files('Mutational_Signatures/Clustered_Mutations')
  #         print(file_list)
  #         barcodes <- c()
  #         for(each in file_list){
  #           print(each)
  #           barcode <- str_split(each, "_")[[1]][1]
  #           barcodes <- append(barcodes, barcode)
  #         }
  # 
  #         print(paste0('barcodes:',barcodes))
  #         updateSelectizeInput(session, 'clustered_mut_barcodes', choices = barcodes, server = TRUE)
  #     }else{
  #       if(input$project_code == "Sherlock_TCGA"){
  #         updateSelectizeInput(session, 'clustered_mut_barcodes', choices = unique(sherlock_data_full$Tumor_Barcode), server = TRUE)
  #       }
  #     }
  # 
  #    
  #   }
  #   
  # })
  # 
  # if(os_detect() %in% c("Linux","Darwin")){
  #   output$figure_pdf_clustered_mut <-  renderImage({figure_output_clustered_mut()},deleteFile=FALSE)
  # }else{
  #   output$figure_pdf_clustered_mut <- renderUI({ tags$iframe(style="height:1000px; width:60%", src= figure_output_clustered_mut())})
  # }
  
######### Genomic Landscape ###############################################

  observeEvent(input$submit_int_ext_data,{
    if("Genomic Landscape" %in% file_list_reactive()){
      data_genomePlot <- read_delim("Genomic_Landscape/genomePlot/genomePlot_list.txt")
      updatePickerInput(session, 'tumor_barcode_to_inspect_genomePlot', "Select one Tumor Barcode:", choices = data_genomePlot %>% pull(1) %>% unique())

      tumor_barcode_to_inspect_genomePlot_reactive <- reactive({
        req(input$tumor_barcode_to_inspect_genomePlot)
        tumor_barcode_to_inspect_genomePlot = input$tumor_barcode_to_inspect_genomePlot
        genomePlot_figurePath <- data_genomePlot$genomePlot[which(data_genomePlot["Tumor_Barcode"] == tumor_barcode_to_inspect_genomePlot)]
      })

      output$genomePlot_figure <- renderImage({list(src=paste0("Genomic_Landscape/",tumor_barcode_to_inspect_genomePlot_reactive()),width="800px", height = "auto", alt=paste0("Tumor Barcode: ",input$tumor_barcode_to_inspect_genomePlot))}, deleteFile = FALSE)
    }
  })

######### Clonal Evolution ###############################################

  observeEvent(input$submit_int_ext_data,{
    if("Clonal Evolution" %in% file_list_reactive()){
      data_mutation_time <- read_delim("Clonal_Evolution/MutationTime_Proportion.txt")
      updatePickerInput(session, "tumor_barcode_mutation_time","Select one Tumor Barcode:", choices= data_mutation_time %>% pull(1) %>% unique())
      
      tumor_barcode_to_inspect_clonal_evolution_reactive <- reactive({
        req(input$tumor_barcode_mutation_time)
        print(input$tumor_barcode_mutation_time)
        #tumor_barcode_to_inspect_clonal_evolution = input$tumor_barcode_mutation_time
        data_mutation_time$Tumor_Barcode[which(data_mutation_time$Tumor_Barcode==input$tumor_barcode_mutation_time)]
        
      })
      
      output$figure_pdf_mutation_time <- renderImage({list(src=paste0("Clonal_Evolution/MutationTime/",tumor_barcode_to_inspect_clonal_evolution_reactive(), "_MTime.png"),width="auto", height = "1000px", alt=paste0("Tumor Barcode: ",tumor_barcode_to_inspect_clonal_evolution_reactive()))}, deleteFile = FALSE)
    }
  })
  
   # figure_output_mutationTime <- reactive({
   #   req(input$tumor_barcode_mutation_time)
   #   tumor_barcode = input$tumor_barcode_mutation_time
   #   if(!is.null(input$project_code)){ # project included in the app
   #     project_code <- input$project_code
   #   }else{ # user project
   #     project_code <- 'User_project'
   #   }
   # 
   #   # if(os_detect() %in% c("Linux","Darwin")){
   #   #   return(figure_display_mutationTime(tumor_barcode, project_code))
   #   # }else{
   #   #   return(src=figure_display_mutationTime(tumor_barcode, project_code))
   #   # }
   # 
   # })

   # observeEvent(input$submit_int_ext_data,{
   #   if("Clonal Evolution" %in% file_list_reactive()){
   #     data_mutation_time <- read_delim("Clonal_Evolution/MutationTime_Proportion.txt")
   #     updatePickerInput(session, "tumor_barcode_mutation_time","Select one Tumor Barcode:", choices= data_mutation_time %>% pull(1) %>% unique())
   # 
   #   }
   # })

   # if(os_detect() %in% c("Linux","Darwin")){
   #   output$figure_pdf_mutation_time <-  renderImage({figure_output_mutationTime()},deleteFile=FALSE)
   # }else{
   #   output$figure_pdf_mutation_time <- renderUI({ tags$iframe(style="height:1000px; width:60%", src= figure_output_mutationTime())})
   # }
   
   #output$figure_pdf_mutation_time <- renderImage({list(src=paste0("Genomic_Landscape/",tumor_barcode_to_inspect_genomePlot_reactive()),width="800px", height = "auto", alt=paste0("Tumor Barcode: ",input$tumor_barcode_to_inspect_genomePlot))}, deleteFile = FALSE)
   #output$figure_pdf_mutation_time <- renderImage({figure_output_mutationTime()},deleteFile=FALSE)

######### Survival Analysis ###############################################

   survprep_reactive <- reactive({
     vartmp = input$vartmp_options_select_survival
     if(input$sp_group_checkbox_survival){
       sp_group = input$sp_group_selected_survival
     }else{
       sp_group = NULL
     }
     
     print(sp_group)
     
     if(is.null(input$project_code)){ 
       user_manifest <- read_delim('Survival_Analysis/survival_manifest.txt')
       data <- read_delim('Survival_Analysis/survival_data_full.txt')
       load('Survival_Analysis/suvdata.RData') # generates suvdata

       print('setting up group for user project for survival analysis')
       if(is.null(sp_group)){
         slists <- user_manifest %>% pull(Tumor_Barcode)
       }else{
         slists <- user_manifest %>% filter(SP_Group==sp_group) %>% pull(Tumor_Barcode)
       }
     }else{
       if(input$project_code == 'Sherlock_TCGA'){
         data = sherlock_data_full
         load('Survival_Analysis/suvdata.RData') # generates suvdata
         
         if(is.null(sp_group)){
           slists <- sherlock_overall %>% pull(Tumor_Barcode)
         }else{
           slists <- sherlock_overall %>% filter(SP_Group==sp_group) %>% pull(Tumor_Barcode)
         }
       }
     }
     
     reference = NULL
     
     return(survprep(data=data, suvdata=suvdata, vartmp=vartmp, sp_group=sp_group, reference=reference, slists=slists))
   })
   
   check_reference_input_survival <- reactive({
     
     data = survprep_reactive()[[1]]
     
     variable = data$Key
     reference_level = input$enter_reference_level
     
     return(check_factor_levels(data, variable, reference_level))
   })
   
   survival_analysis <- reactive({
     
     data <- survprep_reactive()[[1]]
     olevel <- survprep_reactive()[[2]]
     
     print(head(data))
     print(olevel)
     vartmp = input$vartmp_options_select_survival
     
     reference = input$enter_reference_level

     keyname= NULL

     filename = NULL

     return(Survgroup(data = data, vartmp=vartmp, reference=reference, keyname=keyname, filename=filename))
   })

   observeEvent(input$submit_int_ext_data, {
     if("Survival Analysis" %in% file_list_reactive()){
       #shinyjs::disable('survival_reset')
       if(is.null(input$project_code)){
         print('user project survival analysis part 1')
         user_manifest <- read_delim('Survival_Analysis/survival_manifest.txt')
         mdata0 <- read_delim('Survival_Analysis/survival_genomic_alterations.txt')
         updateSelectizeInput(session, 'vartmp_options_select_survival', choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% colnames(), server = TRUE)
         output$sp_group_choices_survival <- renderUI({selectInput("sp_group_selected_survival", "Select one SP_Group to run the fisher test:", choices= sort(unique(user_manifest$SP_Group)), multiple= FALSE)})
       }else{
         if(input$project_code =="Sherlock_TCGA"){
           mdata0 <- sherlock_data_full %>% 
             mutate(Gene=paste0(Gene,"|",Type)) %>% 
             select(Tumor_Barcode,Gene,Alteration) %>% 
             pivot_wider(names_from = "Gene",values_from = "Alteration")
           
           updateSelectizeInput(session, 'vartmp_options_select_survival', choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% colnames(), server = TRUE)
           output$sp_group_choices_survival <- renderUI({selectInput("sp_group_selected_survival", "Select one SP_Group to run the fisher test:", choices= sort(unique(sherlock_overall$SP_Group)), multiple= FALSE)})
         }
       }
     }
     #shinyjs::disable('enter_reference_level')
   })
  
  survival_reference <- eventReactive(input$get_reference_levels,{
    
    ref_levels <- survprep_reactive()[[2]]
    ref_levels <- paste(ref_levels,sep="", collapse=", ")
    shinyjs::enable('enter_reference_level')
    return(ref_levels)
    
  })
  
  observeEvent(input$get_reference_levels,{
    shinyjs::show('reference_levels_out')
    shinyjs::enable('enter_reference_level')
  })
  
  output$reference_levels_out <- renderText({
    #req(input$get_reference_levels)
    survival_reference()})
  
  reference_level_output_reactive <- reactive({
    
    if(input$enter_reference_level != ''){
      out_message <- check_reference_input_survival()
      
    }
    
    if(input$enter_reference_level == ''){
      #print('input box for collapse 1 is blank')
      out_message <- ''
      shinyjs::disable('calculate_survival')
    }
    
    return(out_message)
  })
  
  output$reference_level_check <- renderText({reference_level_output_reactive()})
  
  observeEvent(input$enter_reference_level, {
    
    if(reference_level_output_reactive() == 'Reference level exists in data.'){
      shinyjs::enable('calculate_survival')
    } else{
      shinyjs::disable('calculate_survival')
      
    }
    
    #shinyjs::enable('survival_reset')
  })
  
  output$survival_value <- renderTable({survival_test()[[1]]})
  output$survival_plot <- renderPlot({survival_test()[[2]]})

  survival_test <- eventReactive(input$calculate_survival,{
    tryCatch({
      p <- survival_analysis()
      return(p)
    },
      error = function(cond){
        hide('download_survival_plot')
        showNotification('Error occurred in survival analysis. Please adjust inputs and try again.')
        return()
      }
    )
  })
  
  # observeEvent(input$calculate_survival, {
  #   shinyjs::disable('get_reference_levels')
  # })
  
  observeEvent(input$survival_reset, {
    
    shinyjs::reset('sp_group_checkbox_survival')
    shinyjs::reset('sp_group_selected_survival')
    #shinyjs::reset('get_reference_levels')
    shinyjs::hide('reference_levels_out')
    #output$reference_levels_out <- renderText({''})
    shinyjs::disable('enter_reference_level')
    shinyjs::reset('enter_reference_level')
    shinyjs::enable('calculate_survival')
    #shinyjs::disable('survival_reset')
    
  })
  
  # observeEvent(input$survival_reset, {
  # 
  #   if(input$project_code =="Sherlock_TCGA"){
  #     mdata0 <- sherlock_data_full %>% 
  #       mutate(Gene=paste0(Gene,"|",Type)) %>% 
  #       select(Tumor_Barcode,Gene,Alteration) %>% 
  #       pivot_wider(names_from = "Gene",values_from = "Alteration")
  #   }
  #   #shinyjs::reset('vartmp_options_select_survival')
  #   print(paste0('input$vartmp_options_select_survival:', input$vartmp_options_select_survival))
  #   updateSelectizeInput(session, 'vartmp_options_select_survival', choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% colnames(), server = TRUE, selected = (mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% colnames())[1])
  #   shinyjs::reset('sp_group_checkbox_survival')
  #   shinyjs::reset('sp_group_selected_survival')
  #   shinyjs::enable('get_reference_levels')
  #   output$reference_levels_out <- renderText({''})
  #   shinyjs::reset('enter_reference_level')
  #   shinyjs::reset('calculate_survival')
  # 
  # })

  # observeEvent(input$vartmp_options_select_survival,{
  # 
  #   if(!is.null(input$vartmp_options_select_survival)){
  #     output$reference_levels_out <- renderText({''})
  #     shinyjs::enable('get_reference_levels')
  #     shinyjs::reset('enter_reference_level')
  #   }
  # })

  # observe({
  #   #print(input$vartmp_options_select_survival)
  #   req(reference_reset_reactive())
  #   output$reference_levels_out <- renderText({''})
  #   #shinyjs::enable('get_reference_levels')
  #   shinyjs::reset('enter_reference_level')
  #   
  # })
  
  #reference_reset_reactive <- reactive(input$vartmp_options_select_survival)
  
  
  observeEvent(input$download_surv_plot,{
    
    surv_filename <- paste0(str_replace_all(input$vartmp_options_select_survival,'[/ \\|]','_'),'_Survival.pdf') #add date for when plot was generated?

    ggsave_workaround <- function(g){survminer:::.build_ggsurvplot(x = g,
                                                                   surv.plot.height = NULL,
                                                                   risk.table.height = NULL,
                                                                   ncensor.plot.height = NULL)}

    g_to_save <- ggsave_workaround(survival_analysis()[[2]])

    ggsave(filename = surv_filename, plot = g_to_save,
           width = 10, height = 7, device = cairo_pdf)
    unlink("./Rplots.pdf", recursive = FALSE, force = FALSE, expand = FALSE)

  })
  
  output$download_survival_plot = downloadHandler(
    
    filename = function() {
      if(is.null(input$project_code)){
        project_code <- 'User_Project'
      }else{
        project_code <- input$project_code
      }
      if(input$sp_group_checkbox_survival){
        sp_group = input$sp_group_selected_survival
      }else{
        sp_group = NULL
      }
      
      if(input$sp_group_checkbox_survival){
        group = input$sp_group_selected_survival
        filename = paste0(project_code , '_', str_replace_all(input$vartmp_options_select_survival,'[/ \\|]','_'), '_', sp_group, '_Survival.pdf')
      }else{
        group = NULL
        filename = paste0(project_code , '_', str_replace_all(input$vartmp_options_select_survival,'[/ \\|]','_'), '_Survival.pdf')
        
      }

      return(filename)
    },
    
    content = function(file) {
      ggsave_workaround <- function(g){survminer:::.build_ggsurvplot(x = g,
                                                                     surv.plot.height = NULL,
                                                                     risk.table.height = NULL,
                                                                     ncensor.plot.height = NULL)}
      
      g_to_save <- ggsave_workaround(survival_analysis()[[2]])
      
      ggsave(file, plot = g_to_save,
             width = 10, height = 7, device = cairo_pdf)
      unlink("./Rplots.pdf", recursive = FALSE, force = FALSE, expand = FALSE)
      
    })
  
######### Integrative Analysis ###############################################

  ##### Enrichment Analysis (Fisher Test) #####
  fisher_test <- reactive({
    vartmp = input$vartmp_options_select
    if(input$sp_group_checkbox){
      sp_group = input$sp_group_selected
    }else{
      sp_group = NULL
    }

    if(input$fisher_samplelist_checkbox == FALSE){
      samplelist = NULL
    }else{
      samplelist = input$fisher_samplelist
      samplelist <- read_delim(samplelist$datapath, delim = "\t", col_names = FALSE) %>% pull(1) %>% c()
    }

    if(input$var2name_checkbox == FALSE){
      var2name = NULL
    }else{
      var2name = input$vartmp_options_select_var2
    }

    if(input$excludes_checkbox){
      excludes = input$fisher_excludes
    }else{
      excludes = NULL
    }

    if(input$excludes_cat_checkbox){
      excludes_cat = input$fisher_excludes_cat
    }else{
      excludes_cat = NULL
    }

    if(input$keeps_checkbox){
      keeps = input$fisher_keeps
    }else{
      keeps = NULL
    }

    if(input$keeps_cat_checkbox){
      keeps_cat = input$fisher_keeps_cat
    }else{
      keeps_cat = NULL
    }

    minfreq = as.numeric(input$fisher_min_freq)
    freq_column = input$fisher_freq_colnames
    method = input$fisher_test_type
    glm_formula = input$fisher_glm_input
    fdrcutoff = as.numeric(input$fisher_fdr_cutoff)
    if(is.null(input$project_code)){
      mdata0 <- read_delim('Integrative_Analysis/integrative_mdata0.txt')
      freq_data <- read_delim('Integrative_Analysis/integrative_freq.txt')
      group_data <- read_delim('Integrative_Analysis/integrative_group_data.txt')
      covdata <- read_delim('Integrative_Analysis/integrative_covdata.txt')
    }else{
      if(input$project_code == 'Sherlock_TCGA'){
        covdata0 <- read_delim('Integrative_Analysis/integrative_covdata.txt')

        covdata =  covdata0
        mdata0 <- sherlock_data_full %>%
          mutate(Gene=paste0(Gene,"|",Type)) %>%
          select(Tumor_Barcode,Gene,Alteration) %>%
          pivot_wider(names_from = "Gene",values_from = "Alteration")
        freq_data <- sherlock_freq
        group_data <- sherlock_overall

      }
    }
   
    return(fishergroup(mdata0,freq_data, group_data, vartmp,sp_group, samplelist, var2name,excludes,excludes_cat, keeps,keeps_cat, minfreq=0.03, freq_column, method, glm_formula,covdata,fdrcutoff))

  })

  selected_input_vartmp1 <- reactiveVal()
  
  observeEvent(input$vartmp_options_select,{
    selected_input_vartmp1(input$vartmp_options_select)
  })
  
  selected_input_vartmp2 <- reactiveVal()
  
  observeEvent(input$vartmp_options_select_var2,{
    selected_input_vartmp2(input$vartmp_options_select_var2)
  })
  
  observeEvent(input$submit_int_ext_data, {
    if("Integrative Analysis" %in% file_list_reactive()){
      if(is.null(input$project_code)){
        mdata0 <- read_delim('Integrative_Analysis/integrative_mdata0.txt')
        freq <- read_delim('Integrative_Analysis/integrative_freq.txt')
        updateSelectizeInput(session, 'vartmp_options_select', choices = mdata0 %>% select(-Tumor_Barcode, c(colnames(mdata0))) %>% colnames(), server = TRUE)
        updateSelectInput(session, "fisher_freq_colnames", choices=colnames(freq)[which(str_detect(colnames(freq), regex("freq", ignore_case = T)))])
      }else{
        if(input$project_code =="Sherlock_TCGA"){
          mdata0 <- sherlock_data_full %>%
            mutate(Gene=paste0(Gene,"|",Type)) %>%
            select(Tumor_Barcode,Gene,Alteration) %>%
            pivot_wider(names_from = "Gene",values_from = "Alteration")

          updateSelectizeInput(session, 'vartmp_options_select', choices = mdata0 %>% select(-Tumor_Barcode, c(colnames(mdata0))) %>% colnames(), server = TRUE)
          updateSelectInput(session, "fisher_freq_colnames", choices=colnames(sherlock_freq)[which(str_detect(colnames(sherlock_freq), regex("freq", ignore_case = T)))])
          print(input$vartmp_options_select)
          print(selected_input_vartmp1())
        
        }
      }
      
      output$fisher_glm_message <- renderText({
        if(is.null(input$project_code)){
          covdata0 <- read_delim('Integrative_Analysis/integrative_covdata.txt')
        }else{
          if(input$project_code == 'Sherlock_TCGA'){
            covdata0 <- read_delim('Integrative_Analysis/integrative_covdata.txt')
          }
        }

        var_included <- paste0(colnames(covdata0 %>% select(-"Tumor_Barcode")), collapse= ", ")
        paste0("The following variables can be included in the glm function: ", var_included,
               ". Use Var1 and Var2 in the glm formula. Var1 denotes the first genomic alteration selected, and Var2 denotes the value of each Tumor Barcode for all other genomic alterations.")})
     }
  })
  
  observeEvent(input$vartmp_options_select,{
    if(input$vartmp_options_select !=''){
        if(is.null(input$project_code)){
          mdata0 <- read_delim('Integrative_Analysis/integrative_mdata0.txt')
        }else{
          if(input$project_code == 'Sherlock_TCGA'){
            mdata0 <- sherlock_data_full %>%
          mutate(Gene=paste0(Gene,"|",Type)) %>%
          select(Tumor_Barcode,Gene,Alteration) %>%
          pivot_wider(names_from = "Gene",values_from = "Alteration")
          }
        }
      choices <- mdata0 %>% select(-c(Tumor_Barcode, input$vartmp_options_select)) %>% colnames()
      updateSelectizeInput(session, 'vartmp_options_select_var2', choices = choices, server = TRUE, selected =NULL)
      
    }
    
  })

  output$sp_group_choices <- renderUI({selectInput("sp_group_selected", "Select one SP_Group to run the fisher test:", choices= sort(unique(sherlock_overall$SP_Group)), multiple= FALSE )})

  output$fisher_output_table <- DT::renderDataTable({datatable(fisher_output_reactive()[[1]], options=list(searchHighlight=TRUE, order=c(3,'asc')),filter=list(position="top",clear=TRUE,plain=FALSE))}) # %>% bindEvent(input$calculate_fisher)

  output$fisher_output_plots1 <- renderPlot(fisher_output_reactive()[[2]])
  output$fisher_output_plots3 <- renderPlot(fisher_output_reactive()[[3]])

  fisher_output_reactive <- eventReactive(input$calculate_fisher,{
    tryCatch({
      p <- fisher_test()
      return(p)
    },
    error = function(cond){
      hide('download_fisher_volc_plot')
      hide('download_fisher_bar_plot')
      showNotification('Error occurred in enrichment test. Please adjust inputs and try again.', duration = 10, type = 'error')
      return()
    })
  })

  output$download_fisher_volc_plot = downloadHandler(
    
    filename = function() {
      if(is.null(input$project_code)){
        project_code <- 'User_Project'
      }else{
        project_code <- input$project_code
      }
      if(input$sp_group_checkbox){
        group = input$sp_group_selected
        filename <- paste0(project_code, '_', str_replace_all(input$vartmp_options_select,'[/ \\|]','_'), '_', group, '_fisher_enrichment_volcano_plot.pdf') #add date for when plot was generated?
      }else{
        group = NULL
        filename <- paste0(project_code, '_', str_replace_all(input$vartmp_options_select,'[/ \\|]','_'),'_fisher_enrichment_volcano_plot.pdf') #add date for when plot was generated?
        
      }
      
      print(filename)
      return(filename)
    },
    
    content = function(file) {
      ggsave(file,fisher_output_reactive()[[2]], height = 10, width = 20, device = cairo_pdf)
      
      
    })
  
  output$download_fisher_bar_plot = downloadHandler(
    
    filename = function() {
      if(is.null(input$project_code)){
        project_code <- 'User_Project'
      }else{
        project_code <- input$project_code
      }
      if(input$sp_group_checkbox){
        group = input$sp_group_selected
        filename <- paste0(project_code, '_', str_replace_all(input$vartmp_options_select,'[/ \\|]','_'), '_', group, '_fisher_enrichment_barplot.pdf') #add date for when plot was generated?
      }else{
        group = NULL
        filename <- paste0(project_code, '_', str_replace_all(input$vartmp_options_select,'[/ \\|]','_'),'_fisher_enrichment_barplot.pdf') #add date for when plot was generated?
        
      }
      
      return(filename)
    },
    
    content = function(file) {
      ggsave(file, fisher_output_reactive()[[3]], device = cairo_pdf)
      
      
    })

  observeEvent(input$var2name_checkbox, {
    if(input$var2name_checkbox ==TRUE){

      showTab(inputId = "fisher_results", target = "Figure 2: Bar Plot")
 
    }else{
      hideTab(inputId = "fisher_results", target = "Figure 2: Bar Plot")

    }
  })

  ##### Fisher Bar Plot #####

  fisher_test_bar <- reactive({
    vartmp = input$vartmp_options_select_bar
    if(input$sp_group_checkbox_bar){
      sp_group = input$sp_group_selected_bar
    }else{
      sp_group = NULL
    }

    if(input$fisher_bar_samplelist_checkbox == FALSE){
      samplelist = NULL
    }else{
      samplelist = input$fisher_bar_samplelist
      samplelist <- read_delim(samplelist$datapath, delim = "\t", col_names = FALSE) %>% pull(1) %>% c()
    }

    var2name = input$vartmp_options_select_var2_bar

    if(is.null(input$project_code)){
      mdata0 <- read_delim('Integrative_Analysis/integrative_mdata0.txt')
      group_data <- read_delim('Integrative_Analysis/integrative_group_data.txt')
      
      
    }else{
      if(input$project_code == 'Sherlock_TCGA'){
        mdata0 <- sherlock_data_full %>%
          mutate(Gene=paste0(Gene,"|",Type)) %>%
          select(Tumor_Barcode,Gene,Alteration) %>%
          pivot_wider(names_from = "Gene",values_from = "Alteration")
        group_data <- sherlock_overall
      }
    }
    
    return(fisherbarplot(mdata0,group_data, vartmp, sp_group, samplelist, var2name))

  })

  observeEvent(input$submit_int_ext_data, {
    if(is.null(input$project_code)){
      mdata0 <- read_delim('Integrative_Analysis/integrative_mdata0.txt')
      group_data <- read_delim('Integrative_Analysis/integrative_group_data.txt')
      
      updateSelectizeInput(session, 'vartmp_options_select_bar', choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% colnames(), server = TRUE)
      output$sp_group_choices_bar <- renderUI({selectInput("sp_group_selected_bar", "Select one Group to run the fisher test:", choices= sort(unique(group_data$SP_Group)), multiple= FALSE )})
      
    }else{
      if(input$project_code =="Sherlock_TCGA"){
        mdata0 <- sherlock_data_full %>%
          mutate(Gene=paste0(Gene,"|",Type)) %>%
          select(Tumor_Barcode,Gene,Alteration) %>%
          pivot_wider(names_from = "Gene",values_from = "Alteration")
        updateSelectizeInput(session, 'vartmp_options_select_bar', choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% colnames(), server = TRUE)
        output$sp_group_choices_bar <- renderUI({selectInput("sp_group_selected_bar", "Select one Group to run the fisher test:", choices= sort(unique(sherlock_overall$SP_Group)), multiple= FALSE )})
        
      }
    }
    
  })



  observe({
    req(input$submit_int_ext_data)
    if("Integrative Analysis" %in% file_list_reactive()){
      if(is.null(input$project_code)){
        mdata0 <- read_delim('Integrative_Analysis/integrative_mdata0.txt')
       
      }else{
        if(input$project_code == 'Sherlock_TCGA'){
          mdata0 <- sherlock_data_full %>%
            mutate(Gene=paste0(Gene,"|",Type)) %>%
            select(Tumor_Barcode,Gene,Alteration) %>%
            pivot_wider(names_from = "Gene",values_from = "Alteration")
         
        }
      }
      
      if(input$vartmp_options_select_bar!=""){
        updateSelectizeInput(session, 'vartmp_options_select_var2_bar', choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% select(-input$vartmp_options_select_bar) %>% colnames(), server = TRUE, selected =NULL)
      }
     }
    })

  fisher_output_bar_reactive <- eventReactive(input$generate_fisher_barplot, {

    isolate(fisher_test_bar())

  })

  output$fisher_output_bar <- renderPlot(fisher_output_bar_reactive(), width = 500)
  
  output$download_fisher_bar_plot_only = downloadHandler(
    
    filename = function() {
      if(is.null(input$project_code)){
        project_code <- 'User_Project'
      }else{
        project_code <- input$project_code
      }
      if(input$sp_group_checkbox_bar){
        group <- input$sp_group_selected_bar
        filename <- paste0(project_code, str_replace_all(input$vartmp_options_select_bar,'[/ \\|]','_'),'_',str_replace_all(input$vartmp_options_select_var2_bar,'[/ \\|]','_'), '_', group, '_Group_enriched.pdf')
        
      }else{
        filename <- paste0(project_code, str_replace_all(input$vartmp_options_select_bar,'[/ \\|]','_'),'_',str_replace_all(input$vartmp_options_select_var2_bar,'[/ \\|]','_'),'_enriched.pdf')
        
      }
      
      print(filename)
      return(filename)
    },
    
    content = function(file) {
      ggsave(file, fisher_output_bar_reactive(), device = cairo_pdf)
      
      
    })


##### Association Testing #####
  # Data select input for association testing

  # check levels for inputs in Collapse Level (variable one and/or variable two)
  check_levels_reactive_collapse_var1 <- reactive({

    data = data_input_assoc()[[1]]

    variable = input$variable_choices_1
    collapse_level = input$collapse_assoc_var1

    return(check_levels(data, variable, collapse_level))
  })

  check_levels_reactive_collapse_var2 <- reactive({

    data = data_input_assoc()[[1]]

    variable = input$variable_choices_2
    collapse_level = input$collapse_assoc_var2

    return(check_levels(data, variable, collapse_level))
  })

  check_input_group_var_bivar_reactive <- reactive({

    data = data_input_assoc()[[1]]

    inputs = input$group_var_input_bivariable

    return(check_inputs(data, inputs))
  })

  check_input_formula_multivar_reactive <- reactive({

    data = data_input_assoc()[[1]]

    formula = input$regression_formula

    return(check_inputs(data=data, formula=formula, inputs=NULL))
  })

  check_input_group_var_multivar_reactive <- reactive({

    data = data_input_assoc()[[1]]

    inputs = input$group_var_input

    return(check_inputs(data, inputs))
  })

  association_inputs <- reactive({

    data = data_input_assoc()[[1]]

    plot_width = 12
    plot_height = 8

    if(input$association_testing == "multivar_analysis"){
      regression = TRUE
      formula = input$regression_formula
      if(input$group_var_input == ""){
        Group_Var = NULL
      }else{
        Group_Var = input$group_var_input
      }

      Var1 = NULL
      Var2 = NULL

      filter_zero1 = NULL
      filter_zero2 = NULL
      log_var1 = FALSE
      log_var2 = FALSE
      collapse_var1 = NULL
      collapse_var2 = NULL
      type = "parametric"
    }

    if(input$association_testing == "bivar_analysis"){
      regression = FALSE
      formula = NULL
      if(input$group_var_input_bivariable == ""){
        Group_Var = NULL
      }else{
        Group_Var = input$group_var_input_bivariable
      }
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

      if(input$filter_assoc_var1==FALSE | is.character(data[[input$variable_choices_1]])){
        filter_zero1 = FALSE
      }else{
        filter_zero1 = input$filter_assoc_var1
      }

      if(input$filter_assoc_var2==FALSE | is.character(data[[input$variable_choices_2]])){
        filter_zero2 = FALSE
      }else{
        filter_zero2 = input$filter_assoc_var2
      }
      print(paste('filter_zero2:', filter_zero2))

      if(input$log2_assoc_var1==FALSE | is.character(data[[input$variable_choices_1]])){
        log_var1 = FALSE
      }else{
        log_var1 = input$log2_assoc_var1
      }

      if(input$log2_assoc_var2==FALSE | is.character(data[[input$variable_choices_2]])){
        log_var2 = FALSE
      }else{
        log_var2 = input$log2_assoc_var2
      }

      if(input$collapse_assoc_var1=="" | is.numeric(data[[input$variable_choices_1]])){
        collapse_var1 = NULL
      }else{
        collapse_var1 = input$collapse_assoc_var1
      }

      if(input$collapse_assoc_var2=="" | is.numeric(data[[input$variable_choices_2]])){
        collapse_var2 = NULL
      }else{
        collapse_var2 = input$collapse_assoc_var2
      }

    }

    return(sherlock_genome_association_test(data, Var1, Var2, Group_Var, regression, formula, filter_zero1, filter_zero2,log_var1,log_var2,type,collapse_var1,collapse_var2,xlab,ylab))

  })

  data_input_assoc <- reactive({

    df_list <- list()
    if(length(input$data_list_association_selected) > 1){
      no_tum_barcode <- list()
      if("Sample Level" %in% input$data_list_association_selected){
        sample_data <- read_delim("Manifest_Information/QC_sample_level.txt")
        if("Tumor_Barcode" %in% colnames(sample_data)){
          df_list <- append(df_list, list(sample_data), after = length(df_list))
          print("Sample Level added")
        }else{
          print("Sample level could not be joined with the selected datasets because it does contain the necessary column to be joined by.")
          no_tum_barcode <- append(no_tum_barcode, "sample level", after = length(df_list))
        }

      }
      if("Subject Level" %in% input$data_list_association_selected){
        print('trying to join with subject level data now')
        subject_data <- read_delim("Manifest_Information/QC_subject_level.txt")
        if("Tumor_Barcode" %in% colnames(subject_data)){
          df_list <- append(df_list, list(subject_data), after = length(df_list))
          print("Subject Level added")
        }else{
          print("Subject level could not be joined with the selected datasets because it does not contain the necessary column to be joined by.")
          no_tum_barcode <- append(no_tum_barcode, "subject level", after = length(df_list))
          print(paste0('no_tum_barcode: ', no_tum_barcode))
        }

      }
      if("NGSpurity" %in% input$data_list_association_selected){
        ngspurity_data <- read_delim("NGSpurity/all_ngspurity_output.txt")
        if("Tumor_Barcode" %in% colnames(ngspurity_data)){
          df_list <- append(df_list, list(ngspurity_data), after = length(df_list))
          print("NGSpurity added")
        }else{
          print("NGSpurity could not be joined with the selected datasets because it does contain the necessary column to be joined by.")
          no_tum_barcode <- append(no_tum_barcode, "ngspurity", after = length(df_list))
        }
      }
      if("Survival Data" %in% input$data_list_association_selected){
        #suvdata
        if(is.null(input$project_code)){
          suvdata <- read_delim("Survival_Analysis/survival_metadata.txt")
        }else{
          if(input$project_code == 'Sherlock_TCGA'){
            load('Survival_Analysis/suvdata.RData')
          }
         
        }
        if("Tumor_Barcode" %in% colnames(suvdata)){
          df_list <- append(df_list, list(suvdata), after = length(df_list))
          print("Survival Data added")
        }else{
          print("Survival Data could not be joined with the selected datasets because it does contain the necessary column to be joined by.")
          no_tum_barcode <- append(no_tum_barcode, "survival data", after = length(df_list))
        }
      }

        assoc_join <- join_assoc_data(df_list)[[1]]
        
        assoc_colnames <- colnames(assoc_join %>% select(!ends_with("_Detail")))
        assoc_join_by <- print(paste0("Joined by: ", paste0(join_assoc_data(df_list)[[2]], collapse= ", ")))
        if(length(no_tum_barcode >=1)) {
          assoc_no_tum_barcode <- print(paste0("The following datasets were not able to be joined because they lack the Tumor Barcode column:", no_tum_barcode, "."))
          return(list(assoc_join, assoc_colnames, assoc_join_by, assoc_no_tum_barcode))
        }else{
          assoc_no_tum_barcode <- print('')
          return(list(assoc_join, assoc_colnames, assoc_join_by, assoc_no_tum_barcode))
        }

    }else{
      if("Sample Level" %in% input$data_list_association_selected){
        sample_data <- read_delim("Manifest_Information/QC_sample_level.txt")
        df_list <- append(df_list, list(sample_data), after = length(df_list))
      }
      if("Subject Level" %in% input$data_list_association_selected){
        subject_data <- read_delim("Manifest_Information/QC_subject_level.txt")
        df_list <- append(df_list, list(subject_data), after = length(df_list))
      }
      if("NGSpurity" %in% input$data_list_association_selected){
        ngspurity_data <- read_delim("NGSpurity/all_ngspurity_output.txt")
        df_list <- append(df_list, list(ngspurity_data), after = length(df_list))
      }

      if("Survival Data" %in% input$data_list_association_selected){ # add in user project survival data file

        if(is.null(input$project_code)){
          suvdata <- read_delim("Survival_Analysis/survival_metadata.txt")
        }else{
          if(input$project_code == 'Sherlock_TCGA'){
          load('Survival_Analysis/suvdata.RData')
          }
        }
        df_list <- append(df_list, list(suvdata), after = length(df_list))
      }

      assoc_join <- data.frame(df_list[1])
      assoc_colnames <- colnames(assoc_join %>% select(!ends_with("_Detail")))
      assoc_join_by <- print('')
      assoc_no_tum_barcode <- print('')
      
      return(list(assoc_join, assoc_colnames, assoc_join_by, assoc_no_tum_barcode))
    }

    })

  observeEvent(input$submit_int_ext_data, {
    
    shinyjs::disable('load_datasets')
    if(is.null(input$project_code)){ 
      possible_data_choices <- c("Sample Level", "Subject Level", "NGSpurity", "Survival Data")
      possible_data_choices <- data.frame(data_type = possible_data_choices, file=c('Manifest_Information/QC_sample_level.txt', 'Manifest_Information/QC_subject_level.txt', 'NGSpurity/all_ngspurity_output.txt','Survival_Analysis/survival_metadata.txt'))
      possible_data_choices_data_types <- possible_data_choices$data_type
      # check to see which files exist to be used in association testing (follow order sample level, subject level, ngspurity, survival data)
      possible_data_choices_exist <- possible_data_choices$data_type[which(file.exists('Manifest_Information/QC_sample_level.txt', 'Manifest_Information/QC_subject_level.txt', 'NGSpurity/all_ngspurity_output.txt','Survival_Analysis/survival_metadata.txt'))]
      possible_data_choices_do_not_exist <- !possible_data_choices$data_type %in% possible_data_choices_exist

      output$data_list_association <- renderUI({pickerInput("data_list_association_selected",
                                                            label=paste0("Select one or more datasets from the project selected:"), choices= possible_data_choices_data_types, multiple = TRUE,
                                                            choicesOpt = list(disabled=possible_data_choices_do_not_exist, 
                                                            style = ifelse(possible_data_choices_do_not_exist, yes = "color: rgba(119, 119, 119, 0.5);",no = "")), options = list(style = 'picker-inputs-app', `max-options` = 2))})

    }else{
      if(input$project_code == "Sherlock_TCGA"){
        output$data_list_association <- renderUI({pickerInput("data_list_association_selected",
                                                              label=paste0("Select one or more datasets from the project selected:"), choices= list("Sample Level", "Subject Level", "NGSpurity", "Survival Data"),
                                                              multiple=TRUE, options = list(style = 'picker-inputs-app', `max-options`= 2))})
      }
    }
    
    data_to_load <- reactiveVal('')
    
    observe({
      new_data_to_load <- c(data_to_load(), input$data_list_association_selected)
      print(new_data_to_load)
      if(length(new_data_to_load) > 1){
        shinyjs::enable('load_datasets')
      }else{
        shinyjs::disable('load_datasets')
      }
    })
    
    output$assoc_variable_list_1 <- renderUI({pickerInput("variable_choices_1", "Variable One", choices= NULL, multiple=FALSE, options = list(style = 'picker-inputs-app'))})
    output$assoc_variable_list_2 <- renderUI({pickerInput("variable_choices_2", "Variable Two",choices= NULL, multiple=FALSE, options = list(style = 'picker-inputs-app'))})

    shinyjs::disable('reset_association')
    shinyjs::disable('reset_regression')
    shinyjs::disable('calculate_regression')
    
  })

  observeEvent(input$load_datasets,{

    print(length(data_input_assoc()))
    data_final <- data_input_assoc()[[1]]
    data_final_dim <- data_input_assoc()[[1]] %>% dim()
    print(data_final_dim)
    data_colnames <- data_input_assoc()[[2]]

    print(paste0('length data assoc:', length(data_input_assoc())))
    
    # if(length(data_input_assoc()) >2 ){
    #   if(length(data_input_assoc()) ==3){
    #     data_join_by <- data_input_assoc()[[3]]
    #     output$assoc_data_join_by <- renderText(data_join_by)
    # 
    #   }
    #   if(length(data_input_assoc()) ==4){
    #     data_join_by <- data_input_assoc()[[3]]
    #     data_no_tum_barcode <- data_input_assoc()[[4]]
    #     print(data_no_tum_barcode)
    # 
    #     output$assoc_data_join_by <- renderText(data_join_by)
    #     output$assoc_no_tum_barcode <- renderText(data_no_tum_barcode)
    #   }
    # 
    # 
    # }
    # if(length(data_input_assoc()) ==3 ){
    #   data_join_by <- data_input_assoc()[[3]]
    #   output$assoc_data_join_by <- renderText(data_join_by)
    # }
   # if(length(data_input_assoc()) ==4){
      data_join_by <- data_input_assoc()[[3]]
      data_no_tum_barcode <- data_input_assoc()[[4]]
      print(data_no_tum_barcode)
      
      output$assoc_data_join_by <- renderText(data_join_by)
      output$assoc_no_tum_barcode <- renderText(data_no_tum_barcode)
    #}

    # select columns dropdown
    output$assoc_dropdown <- renderUI({dropdownButton(inputId="assoc_data_dropdown",label="Select Columns",circle=FALSE,status = 'dropdownbutton',checkboxGroupInput("assoc_header",label="Column Names",choices=c("All columns",data_colnames),selected=data_colnames[1:6]))})
    output$assoc_datatable_dim <- renderTable({tibble(rows=data_final_dim[[1]], columns=data_final_dim[[2]])})

    output$assoc_datatable <- DT::renderDataTable({datatable(data_final[ ,input$assoc_header,drop=FALSE], options=list(searchHighlight=TRUE),filter=list(position="top",clear=TRUE,plain=FALSE))})

    updatePickerInput(session, "variable_choices_1", "Variable One", choices= data_colnames)
    updatePickerInput(session, "variable_choices_2", "Variable Two", choices= data_colnames)

  })

  observeEvent(input$group_var_input_bivariable, {
    if(input$group_var_input_bivariable == ''){
      shinyjs::disable('check_group_var_bivar')
    }else{
      shinyjs::enable('check_group_var_bivar')
    }

  })

  observeEvent(input$variable_choices_1,{
    shinyjs::reset('collapse_assoc_var1')
  })
  
  observeEvent(input$variable_choices_2,{
    shinyjs::reset('collapse_assoc_var2')
  })
    check_collapse_var1_output_reactive <- reactive({

      if(input$collapse_assoc_var1 != ''){
        out_message <- check_levels_reactive_collapse_var1()
        
      }

      if(input$collapse_assoc_var1 == ''){
        #print('input box for collapse 1 is blank')
        out_message <- ''
      }
      
      return(out_message)
  })

    
    output$check_collapse_var1_output <- renderText({check_collapse_var1_output_reactive()})
    
    observeEvent(input$collapse_assoc_var1, {
      
      if(check_collapse_var1_output_reactive() == "Categorical value exists in data." | input$collapse_assoc_var1 == ''){
        shinyjs::enable('calculate_association')
      } else{
        shinyjs::disable('calculate_association')
      }
      
    })

    check_collapse_var2_output_reactive <- reactive({

      if(input$collapse_assoc_var2 != ''){

        out_message <- check_levels_reactive_collapse_var2()

      }

      if(input$collapse_assoc_var2 == ''){
        print('input box for collapse 2 is blank')
        out_message <- ''
      }
      #print(out_message)
      return(out_message)
      
    })

    output$check_collapse_var2_output <- renderText({check_collapse_var2_output_reactive()})

    observeEvent(input$collapse_assoc_var2, {
      
      if(check_collapse_var2_output_reactive() == "Categorical value exists in data." | input$collapse_assoc_var2 == ''){
        shinyjs::enable('calculate_association')
      } else{
        shinyjs::disable('calculate_association')
      }
      
    })
    
  observeEvent(input$group_var_input_bivariable,{
    if(input$group_var_input_bivariable == ''){
      shinyjs::disable('check_group_var_bivar')
      shinyjs::enable('calculate_association')
      output$check_group_var_bivar_output <- renderText({''})
    }else{
      shinyjs::enable('check_group_var_bivar')
      shinyjs::disable('calculate_association')
    }
  })


  observeEvent(input$check_group_var_bivar, {

    output$check_group_var_bivar_output <- renderText({isolate(check_input_group_var_bivar_reactive())}) #%>% bindEvent(input$check_input)
    print(paste0('check_group_var_input:',check_input_group_var_bivar_reactive()))
    if(check_input_group_var_bivar_reactive() == 'All variables are found in the data.'){
      shinyjs::enable('calculate_association')
    } else{
      shinyjs::disable('calculate_association')
    }

  })


########### Multivariable (Regression) ##############

  observeEvent(input$regression_formula,{
    if(input$regression_formula == ''){
      shinyjs::disable('check_multivar_formula')
      shinyjs::enable('calculate_regression')
      output$check_multivar_formula_output <- renderText({''})
    }else{
      shinyjs::enable('check_multivar_formula')
      shinyjs::disable('calculate_regression')
    }
  })

  observeEvent(input$check_multivar_formula, {

    output$check_multivar_formula_output <- renderText({isolate(check_input_formula_multivar_reactive())}) #%>% bindEvent(input$check_input)
    print(paste0('check_formula_multivar_input:',check_input_formula_multivar_reactive()))
    if(check_input_formula_multivar_reactive() == 'All variables are found in the data.'){
      shinyjs::enable('calculate_regression')
    } else{
      shinyjs::disable('calculate_regression')
    }

  })


  observeEvent(input$group_var_input,{
    if(input$group_var_input == ''){
      shinyjs::disable('check_group_var_multivar')
      shinyjs::enable('calculate_regression')
      output$check_group_var_multivar_output <- renderText({''})
    }else{
      shinyjs::enable('check_group_var_multivar')
      shinyjs::disable('calculate_regression')
    }
  })


  observeEvent(input$check_group_var_multivar, {

    output$check_group_var_multivar_output <- renderText({isolate(check_input_group_var_multivar_reactive())}) #%>% bindEvent(input$check_input)
    print(paste0('check_group_var_multivar_input:',check_input_group_var_multivar_reactive()))
    if(check_input_group_var_multivar_reactive() == 'All variables are found in the data.' & check_input_formula_multivar_reactive() == 'All variables are found in the data.'){
      shinyjs::enable('calculate_regression')
    } else{
      shinyjs::disable('calculate_regression')
    }

  })

  output$regression_table <-
    renderDataTable({multivar_test()})

  output$regression_plot <-
    renderPlot({multivar_test()})

 multivar_test <- eventReactive(input$calculate_regression,{
    tryCatch({
      p <- association_inputs()
      return(p)
    },
    error = function(cond){
      hide('download_multivar_results')
      showNotification('Error occurred in association test. Please adjust inputs and try again.', duration = 10, type = 'error')
      return()
    })
 })

  observeEvent(input$calculate_regression, {
    shinyjs::disable("calculate_regression")
    shinyjs::enable("reset_regression")
    shinyjs::show("regression_plot")
    shinyjs::show("regression_table")
    shinyjs::show("download_multivar_result")
    
  })
  
  output$download_multivar_result = downloadHandler(
    
    filename = function() {
      if(is.null(input$project_code)){
        project_code <- 'User_Project'
      }else{
        project_code <- input$project_code
      }
      str_spl_formula <- str_split(input$regression_formula, "\\(")
      formula <- paste0("(", str_spl_formula[[1]][2])
      formula <- gsub("\\(", "", formula)
      formula <- gsub("\\)", "", formula)
      formula <- gsub("~","_",formula)
      if(input$group_var_input == ''){
        filename <- paste0(project_code, '_', formula,"_multivar_analysis.pdf")
      }else{
        filename <- paste0(project_code, '_', formula, "_", input$group_var_input,"_multivar_analysis.tsv")
      }
      
      print(filename)
      return(filename)
    },
    
    content = function(file) {
      if(input$group_var_input_bivariable == ''){
        ggsave(file, association_inputs(), device = cairo_pdf, width = 20, height = 7)
      }else{
        table_output <- association_inputs() %>% data.frame()
        write_delim(x=table_output, file=file,delim= '\t')
      }      
      
    })

  # Reset
  observeEvent(input$reset_regression, {
    shinyjs::enable("calculate_regression")
    shinyjs::disable("reset_regression")
    updateTextInput(session, "regression_formula", "Regression Formula", value = "")
    updateTextInput(session, "group_var_input", label = NULL, value = "")

    shinyjs::hide("regression_plot")
    shinyjs::hide("regression_table")
    shinyjs::hide("download_multivar_result")

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
    if(is.character(data[[input$variable_choices_1]]) && is.character(data[[input$variable_choices_2]])) {
      updatePickerInput(session,"assoc_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes", "fisher"), options=pickerOptions(hideDisabled=TRUE))
    }
    if(is.numeric(data[[input$variable_choices_1]]) && is.character(data[[input$variable_choices_2]])) {
      updatePickerInput(session,"assoc_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes"), options=pickerOptions(hideDisabled=TRUE))
    }
    if(is.character(data[[input$variable_choices_1]]) && is.numeric(data[[input$variable_choices_2]])) {
      updatePickerInput(session,"assoc_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes"), options=pickerOptions(hideDisabled=TRUE))
    }
  })

  observe({
    req(input$variable_choices_1, input$variable_choices_2)
    data <- data_input_assoc()[[1]]
    data <- validate_vardf(data)
    if(is.numeric(data[[input$variable_choices_1]])){
      shinyjs::disable("collapse_assoc_var1")
      updateTextInput(session = session, "collapse_assoc_var1", value= "")
    }else{
      shinyjs::enable("collapse_assoc_var1")
    }
    if(is.numeric(data[[input$variable_choices_2]])){
      shinyjs::disable("collapse_assoc_var2")
      updateTextInput(session = session, "collapse_assoc_var2", value= "")
    }else{
      shinyjs::enable("collapse_assoc_var2")
    }
  })

  observe({
    req(input$variable_choices_1, input$variable_choices_2)
    data <- data_input_assoc()[[1]]
    data <- validate_vardf(data)
    if(!is.numeric(data[[input$variable_choices_1]])){
      shinyjs::disable("filter_assoc_var1")
      updateCheckboxInput(session = session,'filter_assoc_var1', value=FALSE)
      shinyjs::disable("log2_assoc_var1")
      updateCheckboxInput(session = session,'log2_assoc_var1', value=FALSE)
    }else{
      shinyjs::enable("filter_assoc_var1")
      shinyjs::enable("log2_assoc_var1")
    }
    if(!is.numeric(data[[input$variable_choices_2]])){
      shinyjs::disable("filter_assoc_var2")
      updateCheckboxInput(session = session,'filter_assoc_var2', value=FALSE)
      shinyjs::disable("log2_assoc_var2")
      updateCheckboxInput(session = session,'log_assoc_var2', value=FALSE)
    }else{
      shinyjs::enable("filter_assoc_var2")
      shinyjs::enable("log2_assoc_var2")
    }
  })


  output$bivariable_table <-
    renderDataTable({bivar_test()})

  output$bivariable_plot <-
    renderPlot({bivar_test()})

  bivar_test <- eventReactive(input$calculate_association,{
    tryCatch({
      p <- association_inputs()
      return(p)
    },
    error = function(cond){

      hide('download_bivar_result')
      showNotification('Error occurred in association test. Please adjust inputs and try again.', duration = 10, type = 'error')
      return()
    })
  })

  observeEvent(input$calculate_association,{
    shinyjs::enable('reset_association')
    shinyjs::disable('calculate_association')
    shinyjs::show("bivariable_plot")
    shinyjs::show("bivariable_table")

    shinyjs::disable("assoc_variable_list_1")
    shinyjs::disable("filter_assoc_var1")
    shinyjs::disable("log2_assoc_var1")
    shinyjs::disable("collapse_assoc_var1")

    shinyjs::disable("assoc_variable_list_2")
    shinyjs::disable("filter_assoc_var2")
    shinyjs::disable("log2_assoc_var2")
    shinyjs::disable("collapse_assoc_var2")

    shinyjs::disable("assoc_types_list")
    shinyjs::disable("group_var_input_bivariable")
    shinyjs::show('download_bivar_result')

  })

  output$download_bivar_result = downloadHandler(
    
    filename = function() {
      if(is.null(input$project_code)){
        project_code <- 'User_Project'
      }else{
        project_code <- input$project_code
      }
      if(input$group_var_input_bivariable == ''){
        filename <- paste0(input$variable_choices_1, "_",input$variable_choices_2,"_bivariable_analysis.pdf")
      }else{
        filename <- paste0(input$variable_choices_1, "_",input$variable_choices_2,"_",input$group_var_input_bivariable,"_bivariable_analysis.tsv")

      }
      
      print(filename)
      return(filename)
    },
    
    content = function(file) {
      if(input$group_var_input_bivariable == ''){
        ggsave(file, association_inputs(), device = cairo_pdf, width = 20, height = 7)
      }else{
        table_output <- association_inputs()
        write_delim(x=table_output, file=file,delim= '\t')
      }      
      
    })
  
  # reset association inputs
  observeEvent(input$reset_association,{
    shinyjs::enable('calculate_association')
    shinyjs::disable('reset_association')
    data <- data_input_assoc()[[1]]
    data_colnames <- data_input_assoc()[[2]]

    shinyjs::enable("assoc_variable_list_1")
    shinyjs::enable("filter_assoc_var1")
    shinyjs::enable("log2_assoc_var1")
    shinyjs::enable("collapse_assoc_var1")

    shinyjs::enable("variable_choices_2")
    shinyjs::enable("filter_assoc_var2")
    shinyjs::enable("log2_assoc_var2")
    shinyjs::enable("collapse_assoc_var2")

    shinyjs::enable("assoc_types_list")
    shinyjs::enable("group_var_input_bivariable")

    updatePickerInput(session, inputId="variable_choices_1","Variable One", choices= data_colnames)
    updatePickerInput(session, inputId= "variable_choices_2", "Variable Two", choices = data_colnames)
    updatePickerInput(session, "variable_choices_2","Variable Two",choices=data_colnames, choicesOpt = list(disabled=data_colnames %in% input$variable_choices_1), options=pickerOptions(hideDisabled=TRUE))
    updateCheckboxInput(session, "filter_assoc_var1","Filter (>0)",value=FALSE)
    updateCheckboxInput(session,"log2_assoc_var1", "log2", value=FALSE)
    updateTextInput(session, "collapse_assoc_var1", "Collapse Level",value = "")
    updateCheckboxInput(session, "filter_assoc_var2","Filter (>0)",value=FALSE)
    updateCheckboxInput(session,"log2_assoc_var2", "log2", value=FALSE)
    updateTextInput(session,"collapse_assoc_var2", "Collapse Level",value = "")
    updatePickerInput(session,"assoc_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes","fisher","skit"))
    updateTextInput(session, "group_var_input_bivariable", value = "")
    shinyjs::hide("download_bivar_result")
    shinyjs::hide("bivariable_plot")
    shinyjs::hide("bivariable_table")

  })

##### Oncoplot #####

  observeEvent(input$submit_int_ext_data, {
    if(is.null(input$project_code)){
      mdata0 <- read_delim('Integrative_Analysis/integrative_mdata0.txt')
      example_data_full <- read_delim('Integrative_Analysis/integrative_data_full.txt')
      updateSelectizeInput(session, "oncoplot_genalts_select", choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% colnames(), server = TRUE)
      updateSelectizeInput(session, "oncoplot_categories_select", choices =  unique(example_data_full$Type), server= TRUE)
      updateSelectizeInput(session, "oncoplot_cat_freq_select", choices =  unique(example_data_full$Type), server= TRUE)
    }else{
      if(input$project_code =="Sherlock_TCGA"){
        mdata0 <- sherlock_data_full %>%
          mutate(Gene=paste0(Gene,"|",Type)) %>%
          select(Tumor_Barcode,Gene,Alteration) %>%
          pivot_wider(names_from = "Gene",values_from = "Alteration")
        type_list <- sherlock_data_full %>% select(Type) %>% unique()
        
        updateSelectizeInput(session, "oncoplot_genalts_select", choices = mdata0 %>% select(c(2:length(colnames(mdata0)))) %>% colnames(), server = TRUE)
        updateSelectizeInput(session, "oncoplot_categories_select", choices =  unique(sherlock_data_full$Type), server= TRUE)
        updateSelectizeInput(session, "oncoplot_cat_freq_select", choices =  unique(sherlock_data_full$Type), server= TRUE)
        
      }
    }

  })

  oncoplot_reactive <- reactive({
    print('entering reactive')
    if(is.null(input$project_code)){
      example_data_full <- read_delim('Integrative_Analysis/integrative_data_full.txt')
      sample_level0 <- example_data_full %>% pull(Tumor_Barcode) %>% unique()
      freq_table <- read_delim('Integrative_Analysis/integrative_freq.txt')
    }else{

        if(input$project_code == "Sherlock_TCGA"){

          sample_level0 <- sherlock_data_full %>% pull(Tumor_Barcode) %>% unique()

          freq_table <- sherlock_freq
        }

    }
 

    # 1
    if(input$oncoplot_input_selection == 'oncoplot_genalts'){
      genomic_alts = input$oncoplot_genalts_select
      order_by_input = input$order_by_input1
      opt_four_freq = NULL
      print(paste0("1:",genomic_alts))
      input_opt = 1
    }

    # 2
    if(input$oncoplot_input_selection == 'oncoplot_text_input_genalts'){
      genomic_alts <- input$oncoplot_user_input_genalts
      order_by_input = input$order_by_input2
      print(genomic_alts)
      opt_four_freq = NULL
      input_opt = 2

    }

    # 3
    if(input$oncoplot_input_selection == 'oncoplot_categories'){
      genomic_alts = input$oncoplot_categories_select
      opt_four_freq = NULL
      input_opt = 3
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

      }
      genomic_alts <- genomic_alts_final
      input_opt = 4
    }



    frequency = as.numeric(input$oncoplot_min_freq)

    data <- oncoplot_data_prep(genomic_alts, opt_four_freq, freq_table, input_opt)

    print('made it past data prep step')


    if(data[[1]] == "source_cat"){
      reps <- rep(paste0("data[[2]][[",1:length(data[[2]]),"]]"))
      reps <- paste0(reps, collapse= " , ")
      rbind_list <- paste0("data_all_cats <- rbind(", reps, ")")
      eval(parse(text=rbind_list))

      data_for_colors <- data_all_cats %>% arrange(Gene, Type)
      data_alts <- data_for_colors %>% pull(Alteration) %>% unique()
    }else{
      data_alts <- data[[2]] %>% pull(Alteration) %>% unique()
    }

    color_list <- read_delim(paste0(parent, "/landscape_colors_all.csv"), delim = ",")
    landscape_colors <- landscape_colors_fcn()

    i <- 1
    new_colors <- c()
    for(each in data_alts){
      if(!(data_alts[i] %in% names(landscape_colors))){
        data_alt <- data_alts[i]

        new_colors[data_alt] <- sample(color_list$Color, 1)

      }
        i <- i + 1

    }

    landscape_colors <- c(landscape_colors, new_colors)

    oncolist_result <- list()
    oncolist_result_by_freq <- list()
    if(data[[1]] == "ind_genalt"){
      data <- data[[2]]
     
      plot_color_vec <- landscape_colors[which(names(landscape_colors) %in% data$Alteration)] %>% c()

      gene_type_level <- data %>% select(Gene,Type) %>% unique()
      print("gene_type_level")
      print(gene_type_level)
      gene_level <- c(data %>% pull(Gene) %>% unique())
     
      print('trying first oncoplot generation')
      res <- oncoplot(data=data,gene_level=NULL,sample_level=NULL, sample_level0 = sample_level0,landscape_colors = plot_color_vec, p2_axis_hidden = TRUE,p2_hidden = TRUE, GeneSortOnly=FALSE, frequency= frequency, gen_alt_second = TRUE, order_by_input= order_by_input)
      print('finished first oncoplot generation')
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

      if(order_by_input){
        reps <- rep(paste0("oncolist_result[[",1:length(oncolist_result),"]]$oncoplot[[1]]"))
        reps <- paste0(reps, collapse= " , ")
        print(reps)
        onco_all <- paste0("oncoplots <- plot_grid(", reps, ", align= 'v', ncol =1)")

        eval(parse(text=onco_all))
      }else{
        oncoplots <- res$oncoplot[[1]]
      }

      if(order_by_input){
        reps <- rep(paste0("oncolist_result[[",1:length(oncolist_result),"]]$oncoplot_legend[[1]]"))
        reps <- paste0(reps, collapse= " , ")
        ncol <- length(oncolist_result)
       
        legs_all <- paste0("plegends <- plot_grid(", reps, ",align = 'h',ncol = ncol, rel_widths=0.2)")
        eval(parse(text=legs_all))
      }else{
        reps <- rep(paste0("res$oncoplot_legend[[",1:length(res$oncoplot_legend),"]]"))
        reps <- paste0(reps, collapse= " , ")
        ncol <- length(res$oncoplot_legend)
        
        legs_all <- paste0("plegends <- plot_grid(", reps, ",align = 'h',ncol = ncol, rel_widths=0.2)")
        eval(parse(text=legs_all))
      }

      plot_combined <- plot_grid(oncoplots, plegends, ncol = 1,rel_heights = c(1,0.2),align = 'v',axis = 'l')

      return(plot_combined)

    }else{
      data <- data[[2]]
    
        name <- paste("test_result", i, sep = "")
        plot_color_vec <- landscape_colors[which(names(landscape_colors) %in% data_all_cats$Alteration)] %>% c()
        # run oncoplot function the first time to determine the gene level and sample level for each input?
        # second run of oncoplot funcion in sherlock_landscape_v2.R was mostly to add in the correct sample level (included all samples from other oncoplot results)
        res <- assign(name,oncoplot(data=data_all_cats,gene_level=NULL,sample_level=NULL, sample_level0 = sample_level0,landscape_colors = plot_color_vec, p2_axis_hidden = TRUE,p2_hidden = TRUE, GeneSortOnly=FALSE, frequency= frequency))
        tmp_sample_level <- res$sample_level %>% pull(Tumor_Barcode)

        i <- i + 1

      # run again with the correct sample level across all inputs
      i <- 1

      oncolist_result <- list()

      for(each in data){
        data_by_type <- data[i]
        data_by_type <- as.data.frame(data_by_type)
        plot_color_vec <- landscape_colors[which(names(landscape_colors) %in% data_by_type$Alteration)] %>% c()
        plot_color_vec <- plot_color_vec[order(match(names(plot_color_vec),data_alts))]

        name <- paste("test_result", i, sep = "")
        res <- assign(name,oncoplot(data=data_by_type,gene_level=NULL,sample_level=tmp_sample_level, sample_level0 = sample_level0,landscape_colors = plot_color_vec, p2_axis_hidden = FALSE,p2_hidden = FALSE, GeneSortOnly=FALSE, frequency= frequency))
        print(res$sample_level)

        oncolist_result <- list.append(oncolist_result, res)
        i <- i + 1

      }


      reps <- rep(paste0("oncolist_result[[",1:length(oncolist_result),"]]$oncoplot[[1]]"))
      reps <- paste0(reps, collapse= " , ")
     
      gene_level_length <- paste0("oncolist_result[[",1:length(oncolist_result),"]]$gene_level")
      i <- 1
      total_length <- 0
      for(each in gene_level_length){
        gene_level <- paste0("gene_levelb <-oncolist_result[[i]]$gene_level")
        eval(parse(text=gene_level))
        total_length <- total_length + length(gene_levelb)
        i <- i + 1

      }
      print(reps)
      onco_all <- paste0("oncoplots <- plot_grid(", reps, ",align = 'v', ncol= 1)")

      eval(parse(text=onco_all))

      reps <- rep(paste0("oncolist_result[[",1:length(oncolist_result),"]]$oncoplot_legend"))
      reps <- paste0(reps, collapse= " , ")
      print(reps)
      ncol <- length(genomic_alts)
      h2 <- length(data_alts)/25
      if(h2 < 0.2){h2 <- 0.2}
      num_plots <- length(oncolist_result)
      legs_all <- paste0("plegends <- plot_grid(", reps, ",align = 'h',ncol= ncol, rel_heights=h2,rel_widths = c(1,1,1,1),nrow = 1)")
      print(legs_all)
      eval(parse(text=legs_all))
     
      plot_combined <- plot_grid(oncoplots, plegends, ncol = 1,align = 'v')
      return(plot_combined)

    }

  })

  # tests to see if the bivariable association test results in an error somewhere in the function
  # if error, prints error message to the ui
  oncoplot_create <- eventReactive(input$generate_oncoplot,{
    tryCatch({
      final_plot <- oncoplot_reactive()
      return(final_plot)
    },
    error = function(cond){
      hide('download_oncoplot')
      showNotification('Error occurred in generating the oncoplot. Please adjust inputs and try again.', duration = 10, type = 'error')
      return()
    })
  })
  
  output$data_output_plot <- renderPlot({oncoplot_create()}, bg = 'transparent')
    
    output$download_oncoplot = downloadHandler(
      
      filename = function() {
        if(is.null(input$project_code)){
          project_code <- 'User_Project'
        }else{
          project_code <- input$project_code
        }
        
        filename <- paste0(project_code, '_oncoplot.pdf')
        print(filename)
        return(filename)
      },
      
      content = function(file) {
        ggsave(file,oncoplot_reactive(), height = 10, width = 25, device = cairo_pdf)
        
      })
    
##### Genomic Features Wordcloud #####

    wordcloud_reactive <- reactive({
      if(is.null(input$project_code)){
        print('here in wordcloud reactive now')
        example_user_data <- read_delim('Integrative_Analysis/integrative_wordcloud_user_data.txt')
        freq_new <- read_delim('Integrative_Analysis/integrative_freq.txt') 
        freq_new <- freq_new %>% separate(col=name, into=c("Gene","Type"), sep= "\\|", remove= FALSE)
        user_wordcloud <- left_join(example_user_data, freq_new, by= c("Gene", "Type"))
        user_wordcloud <- user_wordcloud %>% mutate(Wordcloud_cat= if_else(Type %in% c("Mutation_Driver", "Fusion","SCNA_Focal_Cytoband"), "Driver_Events","Others"))
        word_alts <- user_wordcloud %>% pull(Alteration) %>% unique()
        wordcloud <- user_wordcloud
        
      }else{
        if(input$project_code == "Sherlock_TCGA"){

          sherlock_freq_new  <- sherlock_freq %>% separate(col=name, into=c("Gene","Type"), sep= "\\|", remove= FALSE)
          sherlock_wordcloud <- left_join(sherlock_data, sherlock_freq_new, by= c("Gene", "Type"))
          sherlock_wordcloud <- sherlock_wordcloud %>% mutate(Wordcloud_cat= if_else(Type %in% c("Mutation_Driver", "Fusion","SCNA_Focal_Cytoband"), "Driver_Events","Others"))
          word_alts <- sherlock_wordcloud %>% pull(Alteration) %>% unique()
          wordcloud <- sherlock_wordcloud
          
        }
      }
      
      word_alts = word_alts
      tumor_barcode = input$wordcloud_tumor_barcode
      

      return(genomic_features_wordcloud(data=wordcloud, word_alts, tumor_barcode))
    })

    observeEvent(input$submit_int_ext_data, {
      if(is.null(input$project_code)){
        example_data_full <- read_delim('Integrative_Analysis/integrative_data_full.txt')
        updateSelectizeInput(session, 'wordcloud_tumor_barcode', choices = example_data_full %>% pull(Tumor_Barcode), server = TRUE)
        
      }else{
        if(input$project_code =="Sherlock_TCGA"){
          updateSelectizeInput(session, 'wordcloud_tumor_barcode', choices = sherlock_data %>% pull(Tumor_Barcode), server = TRUE)
        }
      }


      output$sherlock_wordcloud_plot <- renderPlot({
        wordcloud_reactive()}, height ='auto', width='auto') %>%
        bindEvent(input$generate_wordcloud)
    })
    
    output$download_wordcloud = downloadHandler(
      
      filename = function() {
        if(is.null(input$project_code)){
          project_code <- 'User_Project'
        }else{
          project_code <- input$project_code
        }
        
        filename <- paste0(project_code, '_wordcloud.pdf')
        print(filename)
        return(filename)
      },
      
      content = function(file) {
        ggsave(file,wordcloud_reactive(), height = 10, width = 15, device = cairo_pdf)
        
      })


########## Documentation ###############################################

  # Module Info
 observe(
    if(input$sbmenu == 'documentation'){
      if(input$documentation_tabs == 'module_info'){
        
        file_path <- paste0(parent,'/Module_Info.txt')
        output$module_info_table <- DT::renderDataTable({datatable(read_delim(file_path, delim = "\t"), options=list(paging=FALSE))})

      }
    }
)
  
  # Data Requirement Info
  module_chosen <- reactive({
    req(input$module_list)
    data_req_paths <- paste0(parent, '/Data_requirements.txt')
    data_req_paths <- read_delim(data_req_paths, delim = '\t')

    module_list_data_req = input$module_list

  })

  observeEvent(input$module_list, {
    
    data_req_paths <- paste0(parent, '/Data_requirements.txt')
    data_req_paths <- read_delim(data_req_paths, delim = '\t')

    module_list_submodule <- data_req_paths$Submodule[which(data_req_paths["Module"]== module_chosen())]

    updatePickerInput(session,"submodule_user_select",choices=unique(module_list_submodule))

  })

  files_for_submodule <- reactive({

    data_req_paths <- paste0(parent, '/Data_requirements.txt')
    data_req_paths <- read_delim(data_req_paths, delim = '\t')
    
    module_selected <- module_chosen()
    print(paste0('module_selected',module_selected))
    
    data_req_paths <- data_req_paths %>% filter(Module == module_selected)
    
    submodule_selected <- data_req_paths$Submodule[which(data_req_paths["Module"]== module_selected)]
    print(paste0('submodule_selected',submodule_selected))
    
    data_req_paths <- data_req_paths %>% filter(Submodule == submodule_selected)
    
    files_in_submodule <- data_req_paths$File[which(data_req_paths["Submodule"] == input$submodule_user_select)]

    return(files_in_submodule)
  })
  
  output$file_list <- renderUI({pickerInput("file_user_select", label="File List",choices=files_for_submodule(), options = list(style = 'picker-inputs-app'))})

  shinyjs::disable('clear_selected_example_file')
  
  observeEvent(input$select_example_file,{
    
    shinyjs::disable('select_example_file')
    shinyjs::enable('clear_selected_example_file')
    shinyjs::show('example_file_col_dropdown')
    shinyjs::disable('module_list')
    shinyjs::disable('submodule_user_select')
    shinyjs::disable('file_user_select')
    
    data_req_paths <- paste0(parent, '/Data_requirements.txt')
    data_req_paths <- read_delim(data_req_paths, delim = '\t')
    
    module_selected <- module_chosen()
    print(paste0('module_selected',module_selected))
    
    data_req_paths <- data_req_paths %>% filter(Module == module_selected)
    
    submodule_selected <- data_req_paths$Submodule[which(data_req_paths["Module"]== module_selected)]
    print(paste0('submodule_selected',submodule_selected))
    
    data_req_paths <- data_req_paths %>% filter(Submodule == submodule_selected)
    
    files_in_submodule <- data_req_paths$File[which(data_req_paths["Submodule"] == input$submodule_user_select)]
    
    x <- input$file_user_select
    x_path_in_file <- data_req_paths$File_path[which(data_req_paths["File"] == x)]
    print(paste0('x_path_in_file:', x_path_in_file))
    x_path_all <- paste0(parent, '/www/Genomic Data/Data_Requirements/', x_path_in_file)
    x_path_all <- x_path_all[1]
    print(paste0('x_path_all:', x_path_all))
    x_file <- read_delim(x_path_all, delim="\t")
    print(head(x_file))
    x_col_total <- length(colnames(x_file))
    if(x_col_total < 6){
      x_col_total <- x_col_total
    }else{
      x_col_total <- 6
    }
    output$example_file_column_names <- renderUI({dropdownButton(inputId="example_file_col_dropdown",label="Select Columns",circle=FALSE,status = 'dropdownbutton',checkboxGroupInput("example_file_col_input",label="Column Names",choices=c("All columns",colnames(x_file)),selected=colnames(x_file[1:x_col_total])))})
    output$example_file <- DT::renderDataTable({datatable(head(df_filter_reactive()), options=list(paging = FALSE, searching = FALSE))})
    
    # required columns output
    submodule_select = input$submodule_user_select
    req_col <- data_req_paths$Required_columns[which(data_req_paths["File_path"] == x_path_in_file)[1]]
    
    print(paste0('submodule_select ', submodule_select))
    print(paste0('req_col:',req_col))
    output$req_col_out <- renderText({paste0('The required columns for this file are: ', req_col)})

  })
  
  observeEvent(input$example_file_col_input, {

    data_req_paths <- paste0(parent, '/Data_requirements.txt')
    data_req_paths <- read_delim(data_req_paths, delim = '\t')
    
    module_selected <- module_chosen()
    print(paste0('module_selected',module_selected))

    data_req_paths <- data_req_paths %>% filter(Module == module_selected)

    submodule_selected <- data_req_paths$Submodule[which(data_req_paths["Module"]== module_selected)]
    print(paste0('submodule_selected',submodule_selected))

    data_req_paths <- data_req_paths %>% filter(Submodule == submodule_selected)

    files_in_submodule <- data_req_paths$File[which(data_req_paths["Submodule"] == input$submodule_user_select)]
    file_paths_in_submodule <- data_req_paths$File_path[which(data_req_paths["File"] == files_in_submodule)]

    x <- input$file_user_select
    x_path_in_file <- data_req_paths$File_path[which(data_req_paths["File"] == x)]
    print(paste0('x_path_in_file:', x_path_in_file))
    x_path_all <- paste0(parent, '/www/Genomic Data/Data_Requirements/', x_path_in_file)
    x_path_all <- x_path_all[1]
    print(paste0('x_path_all:', x_path_all))
    x_file <- read_delim(x_path_all, delim="\t")
    print(head(x_file))

    x_col_total <- length(colnames(x_file))
    if(x_col_total < 6){
      x_col_total <- x_col_total
    }else{
      x_col_total <- 6
    }
    if("All columns" %in% input$example_file_col_input){
      updateCheckboxGroupInput(session, "example_file_col_input", label="Column Names",choices=c("All columns",colnames(x_file)),selected=c("All columns",colnames(x_file)[1:length(colnames(x_file))]))
      
    }
  })
  
  # try adding in clear/reset functionality so there isn't the error in the console for column names
  observeEvent(input$clear_selected_example_file, {
    shinyjs::enable('select_example_file')
    shinyjs::disable('clear_selected_example_file')
    
    shinyjs::hide('example_file_col_dropdown')
    output$example_file <- DT::renderDataTable(NULL)
    output$req_col_out <- renderText(NULL)
    shinyjs::enable('module_list')
    shinyjs::enable('submodule_user_select')
    shinyjs::enable('file_user_select')
    
  })
  
  
}

