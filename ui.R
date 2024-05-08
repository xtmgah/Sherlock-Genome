
#Sherlock-Genome: A R Shiny App for Genomic Analysis and Visualization

# Warning messages:
#   1: In fileDependencies.R(file) :
#   Failed to parse /private/var/folders/z0/03n38lh1409btl6jtzdq9x49w00fmd/T/RtmpheTtUV/file606e34acd684/Sherlock_Genome_Functions.R ; dependencies in this file will not be discovered.
# 2: In sprintf(gettext(fmt, domain = domain), ...) :
#   one argument not used by format 'invalid uid value replaced by that for user 'nobody''
# 3: invalid uid value replaced by that for user 'nobody' 
# 4: In sprintf(gettext(fmt, domain = domain), ...) :
#   one argument not used by format 'invalid gid value replaced by that for user 'nobody''
# 5: invalid gid value replaced by that for user 'nobody'

# keep adding to this vector when additional packages are required to run the app 
# packages_req <- c("shiny","shinydashboard", "hrbrthemes", "rlist","tidyverse", "tidyr", "markdown","shinyjs", "ggside","tibble","inspectdf","DT","dplyr","shinyWidgets","stringr","PMCMRplus","ggsci", "ggstatsplot","ggplot2", "ggpubr","cowplot","data.table","forcats","broom","rstudioapi","fs","purrr","scales", "survival", "survminer", "extrafontdb", "purrr") #ggpubr- ggqqplot()
packages_req <- c("shinybusy", "devtools","BiocManager", "randomcoloR","shinyBS","shinycssloaders","ggrepel","spsComps","bslib", "htmltools", "shiny","shinydashboard", "shinyFiles","hrbrthemes", "rlist","tidyverse", "markdown","shinyjs", "ggside","inspectdf","DT","shinyWidgets","PMCMRplus","ggsci", "ggstatsplot", "ggpubr","cowplot","data.table","broom","rstudioapi","fs","scales", "survival", "survminer", "extrafontdb","ggwordcloud","graphics",'SKIT',"shinydashboardPlus","valr","ggdendro", "circlize", "ggplot2") #BiocManager? ggpubr- ggqqplot()
packages_req_bioconductor <- c("CNTools", "maftools")
packages_req_from_github <- c("ReConPlot")

# check for required packages and install those not installed; load all packages required
pkg_check_fcn <-lapply(packages_req,function(x){
  #if a package in the packages_req vector is not installed, install it
  if(!(x %in% installed.packages())){
    install.packages(x, dependencies = TRUE, repos=structure(c(CRAN="http://cloud.r-project.org/")))
    #load the package after installing it
    library(x, character.only = TRUE)
  #else statement for when the packages are already installed
  }else{
  #load the package already installed into the workspace
      library(x, character.only = TRUE)
      #print(paste0(x,":" , "Package already installed"))
  }

})

pkg_check_fcn_bioconductor <- lapply(packages_req_bioconductor,function(x){
  if(!(x %in% installed.packages())){
    BiocManager::install(x)
    #load the package after installing it
    library(x, character.only = TRUE)
    #else statement for when the packages are already installed
  }else{
    #load the package already installed into the workspace
    library(x, character.only = TRUE)
    #print(paste0(x,":" , "Package already installed"))
  }
})

pkg_check_fcn_from_github <- lapply(packages_req_from_github,function(x){
  if(!(x %in% installed.packages())){
    if(x == 'ReConPlot'){
      devtools::install_github(paste0('cortes-ciriano-lab/',x))
    }else{
      devtools::install_github(x)
    }
    #load the package after installing it
    library(x, character.only = TRUE)
    #else statement for when the packages are already installed
  }else{
    #load the package already installed into the workspace
    library(x, character.only = TRUE)
    #print(paste0(x,":" , "Package already installed"))
  }
})
  

extrafont::loadfonts()

parent <- getwd()


# detect OS (Windows vs. Mac)
os_detect <- function(){
  os_name <- Sys.info()[['sysname']]
  return(os_name)
}

inline <-  function (x) {
  tags$div(style="display:inline-block;", x)
}

# onStop(function() {
#   setwd("..")
#   unlink("User_project", recursive = TRUE)
# })


# ui setup 
ui <- fluidPage(
  busy_start_up(
    loader = spin_kit(
      spin = "cube-grid",
      color = "#2379ba",
      style = "width:50px; height:50px;"
    ),
    text = "Loading...",
    mode = "manual",
    color = "#2379ba",
    background = "#ffffff"
  ),
#ui <- page_fluid(  
  
  # all action buttons
  tags$head(tags$style(HTML(
                            '.action-buttons-app {background-color: #3983C0; border-color:#3983C0; border-radius: 28px; color: white;}',
                            '.picker-inputs-app {background-color: white; border-color:gray; border-radius: 28px}',
                            '.selectize-input {background-color: white; border-color:gray; border-radius: 28px}',
                            '.btn-dropdownbutton{background-color: #3983C0; border-color:#3983C0; border-radius: 28px; color: white;}',
                            '.flipbox {
                              background-color: #ffffff;
                              width: 300px;
                              height: 200px;
                              border: 1px solid #000000;
                              }',
                            '.btn-custom-class {
                              color: white;
                              background-color: #2379ba;
                              border: 1px solid #ffffff;
                            }'
                            
                            #'.box {border-radius: 28px;}'
                            ))),
  
  # all pickerInputs
  #tags$head(tags$style(HTML('.picker-inputs-app {background-color: white; border-color:gray; border-radius: 28px}'))),
  # all selectInputs
  #tags$head(tags$style(HTML('.selectize-input {background-color: white; border-color:gray; border-radius: 28px}'))),
  
  # tags$head(tags$style(HTML(
  # '
  # .btn-dropdownbutton{
  #                     background-color: #3983C0;border-color:#3983C0; border-radius: 28px; color: white;}'))),
  # title at the top of the page; windowTitle = what the browser tab reads

  useShinyjs(),
  titlePanel(title = "Sherlock-Genome: A R Shiny App for Genomic Analysis and Visualization",
             windowTitle = "Sherlock-Genome"),

  # dashboard setup with left sidebar
  dashboardPage(skin = "blue",
    dashboardHeader(title = "Sherlock-Genome"),
                                        #tags$li(strong('About Sherlock-Genome'), class = 'dropdown')),
    
# each of the following will be a tab on the left sidebar
dashboardSidebar(sidebarMenu(id = "sbmenu",
      menuItem("About",tabName = "introduction", icon = icon('house')),
      menuItem("Load Data", tabName = "data_load", selected = TRUE),
      menuItem("Study Overview", tabName = "study_overview"),
      menuItem("Manifest Information", tabName = "sample_qc"),
      menuItem("NGSpurity", tabName = "NGSpurity"),
      menuItem("Mutations", tabName = "mutations"),
      menuItem("SCNA", tabName = "scna"),
      menuItem("SV", tabName = "sv"),
      menuItem("Mutational Signatures", tabName = "mutational_signature"),
      menuItem("Genomic Landscape", tabName = "genomic_landscape"),
      menuItem("Clonal Evolution", tabName = "clonal_evolution"),
      menuItem("Survival Analysis", tabName = "survival_analysis"),
      menuItem("Integrative Analysis", tabName = "integrative_analysis"),
      menuItem("Documentation", tabName= "documentation"),
      menuItem("Frequently Asked Questions",tabName = "faq")
    )),
# each of the following are what the tabs assigned page will display when the tab is clicked
    dashboardBody(
      
      tags$head(tags$style(HTML('
      .content-wrapper {
        background-color: #fff;
        overflow: auto;
      }
    '))),
      tags$script(HTML("
        var openTab = function(tabName){
          $('a', $('.sidebar')).each(function() {
            if(this.getAttribute('data-value') == tabName) {
              this.click()
            };
          });
        }
      ")),
      
      
      tabItems(
        tabItem(tabName = "introduction", h2("Welcome to Sherlock-Genome!"),
                #h4("Sherlock-Genome is an application for..... Please expand any of the boxes below to learn more about the different modules included in the application. Otherwise, select the 'Load Data' module in the sidebar panel to begin using Sherlock-Genome."),
                h4("Sherlock-Genome is an application for various types of genomic analysis and data visualization. Several modules are included in the app to allow for such analysis and visualization. Please click on any of the boxes below to learn more about the different modules included in the application. Otherwise, select the 'Load Data' module in the sidebar panel to begin using Sherlock-Genome."),
                #fluidRow(
                #box(title = 'Load Data', status = 'primary', solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, markdown('Load the data for a selected project into the corresponding modules.')),
                #uiOutput('active_side_load_data'),
                #img(src = "intro_test.png"),
                #flipBox(id = 'load_data', front = div(h4('Load Data')), back = div(h4('Load the data for a selected project into the corresponding modules.'))), 
                #infoBox('Load Data', color = 'blue', icon = icon('upload'), markdown('Load the data for a selected project into the corresponding modules.')), 
                 box(title = 'Study Overview', status = 'primary', solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, markdown('Brief summary about the project selected.')),
                 #infoBox('Study Overview', color = 'blue',icon = icon('folder-open'),markdown('Brief summary about the project selected.')),
                 box(title = 'Manifest Information', status = 'primary', solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, markdown('View or inspect sample and subject level data')),
                 #infoBox(title = 'Manifest Information', color = 'blue', icon = icon('pen-to-square'), markdown('View or inspect sample and subject level data')),
                 box(title = 'NGSpurity', status = 'primary', solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, markdown('NGSpurity pipeline output to investigate purity, ploidy, and clonality of samples.'))
                 #infoBox('NGSpurity', color = 'blue',icon = icon(lib = 'glyphicon','tasks'), markdown('NGSpurity pipeline output to investigate purity, ploidy, and clonality of samples.'))
                  #card(full_screen = FALSE, card_header("Load Data"), card_body(markdown('Load the data for a selected project into the corresponding modules.')))
                #),
                ,
                #fluidRow(
                  box(title = 'Mutations', status = 'primary', solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, markdown('Exploration gene mutations across transcripts. Visualization includes a lollipop plot.')),
                  #infoBox('Mutations', color = 'blue', icon = icon(lib = 'glyphicon','stats'), markdown('Exploration gene mutations across transcripts. Visualization includes a lollipop plot.')),
                  box(title = 'SCNA', status = 'primary', solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, markdown('Inspect SCNA data, perform clustering analysis, and/or visualize GISTIC ouput.')),
                  #infoBox('SCNA', color = 'blue', icon = icon(lib = 'glyphicon','stats'), markdown('need to complete')),
                  box(title = 'SV', status = 'primary', solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, markdown('Generate a ReCon plot to visualize SVs.')),
                  #infoBox('SV', color = 'blue', icon = icon(lib = 'glyphicon','stats'), markdown('need to complete')),
                  box(title = 'Mutational Signatures', status = 'primary', solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, markdown('Exploration of data related to mutational signatures. This includes mutation clustering.'))
                  #infoBox('Mutational Signatures', color = 'blue',icon = icon(lib = 'glyphicon','stats'), markdown('Exploration of data related to mutational signatures. This includes mutation clustering.'))
                #),
                ,
                #fluidRow(
                  box(title = 'Genomic Landscape', status = 'primary', solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, markdown('Exploration of the genomic landscape of samples from the selected project.')),
                  #infoBox('Genomic Landscape', color = 'blue', icon = icon(lib = 'glyphicon','equalizer'), markdown('Exploration of the genomic landscape of samples from the selected project.')), 
                  box(title = 'Clonal Evolution', status = 'primary', solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, markdown('Explore measures of mutation time in each sample.')),
                  #infoBox('Clonal Evolution', color = 'blue', )
                  box(title = 'Survival Analysis', status = 'primary', solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, markdown('Conduct a survival analysis with all of the genomic alteration data available for the selected project.')),
                  box(title = 'Integrative Analysis', status = 'primary', solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, markdown('Multiple options for analyses and testing of the data available from the selected project. These include fisher enrichment testing, fisher test barplot visualization, association testing, oncoplot generation, and word cloud generation.')),
                #))
                # tags$iframe(style="height:1000px; width:60%", src=paste0(parent, '/www/sherlock_genome_app_introduction_tab'))
                img(src= 'sherlock_genome_app_introduction_tab.png', width = 1500, height = 700)
                ),
        tabItem(tabName ="data_load",h2("Load Data"),hr(),h3("Welcome to the Sherlock-Genome data loading page. If the project you wish to explore is already included in Sherlock-Genome, select 'Project in App'. Otherwise, select 'Upload Project Data'. You can click 'Clear Selection' at any time to reset your selections."),
                       h4("Click", a("here", onclick = "openTab('documentation')", href="#"), "for documentation regarding data input requirements for user data, projects included, and modules available in the Sherlock-Genome app."),
                       #prettyRadioButtons(inputId = "int_ext", label = NULL, choices = c("Project in App", "Upload Project Data"), selected = character(0), inline = TRUE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                       radioGroupButtons(
                          inputId = "int_ext",
                          label = NULL, 
                          choices = c("Project in App", "Upload Project Data"),
                          status = "primary",
                          selected = character(0)
                        ),
                       #radioButtons("int_ext", label= NULL, choices = c("Project in App", "Upload Project Data"), selected = character(0)),
                       #actionButton("sel_project_in_app", "Project in App"),actionButton("sel_project_ext", "Upload Project Data"),
                       #conditionalPanel(condition= "input.int_ext", actionButton("clear_int_ext", "Clear Selection")),
                       # internal project
                       conditionalPanel(condition="input.int_ext == 'Project in App'",
                       #conditionalPanel(condition="input.sel_project_in_app",
                       # uiOutput("ui_project_code"),
                       # uiOutput("ui_select_project"),
                       # uiOutput("ui_reset_project"),
                       #selectInput("project_code","Select a project:", c("BRCA_HK", "mWGS","Sherlock")),
                       selectizeInput("project_code","Select a project:", choices = c("Sherlock_TCGA"), selected = NULL, multiple = TRUE, 
                                      options = list(maxItems = 1)),
                       #pickerInput("project_code","Select a project:", c("BRCA_HK", "mWGS","Sherlock"), selected= "Sherlock"),
                       #actionButton("select", "Select Project", style = "background-color: #3983C0; border-color:#3983C0; border-radius: 28px; color: white;")),
                       actionButton("select", "Select Project", class = 'action-buttons-app')),
                       #actionBttn("select", "Select Project", style = "jelly", color = "primary", size = "sm")),
                       #actionButton("reset_project", "Select Different Project")),
                       # external project
                       conditionalPanel(condition="input.int_ext == 'Upload Project Data'", shinyDirButton("user_choose_directory", "Select Project Folder", "Please select a folder", class = 'action-buttons-app')),
                       conditionalPanel(condition="input.int_ext == 'Upload Project Data' && input.user_choose_directory", h5(strong(textOutput("user_directory_path")))),

                       #conditionalPanel(condition = "input.select", h4(textOutput("choose_method_in_app")), h5(textOutput("all_mods_in_app")), h5(textOutput("ind_mods_in_app"))),
                       h4(textOutput("choose_method_in_app")), h4(textOutput("all_mods_in_app")), h4(textOutput("ind_mods_in_app")),
                       # conditionalPanel(condition = "input.select", h4("Choose one of the two methods below to load data into the corresponding modules in the panel to the left:"),
                       # h5("1. Select a project, and all modules available for that project. This will automatically populate each module with its corresponding data."),
                       # h5("2. Select a project, and individual modules for that project. This will populate each module selected with its corresponding data.")),
                       # either internal or external project
                       #conditionalPanel(condition="input.select",uiOutput("file_list_select_internal"),actionButton(inputId="submit_internal", label="Submit"),textOutput("file_load_message_internal")),
                       #conditionalPanel(condition="input.user_choose_directory", uiOutput("file_list_select_external"), actionButton(inputId="submit_external", "Submit"), textOutput("file_load_message_external")))
                       conditionalPanel(condition="input.select",uiOutput("file_list_select_internal")),
                       #conditionalPanel(condition="input.user_choose_directory", uiOutput("file_list_select_external")),
                       conditionalPanel(condition = "input.int_ext == 'Upload Project Data'", h4(textOutput("choose_method_outside")), h4(textOutput("all_mods_outside")), h4(textOutput("ind_mods_outside"))),
                       uiOutput("file_list_select_external"),
                       uiOutput("submit_int_ext_data_out"),
                       conditionalPanel(condition="input.select",h4(textOutput("file_load_message_internal"))),
                       conditionalPanel(condition="input.user_choose_directory",h4(textOutput("file_load_message_external"))),
                       hr(style = "border-top: 1px solid #ffffff;"),uiOutput("ui_clear_int_ext"),
                       ),
               
               tabItem(tabName = "study_overview", h2("Study Overview"),hr(),uiOutput("study_overview_rmd")),

               tabItem(tabName = "sample_qc", h2("Manifest Information"),hr(), tabsetPanel(id="qc_tabs",tabPanel(value="qc_sample", "Sample Level", tabsetPanel(id="qc_sample_subtabs",
                                                                                                                                              #tabPanel(value="qc_view_sample_data","View Sample Level",h4("Select the columns you would like to view from the dropdown below:"),conditionalPanel("!input.filter_df_sample",uiOutput("qc_sample_header_output")),conditionalPanel("input.filter_df_sample",uiOutput("qc_sample_header_output2")),uiOutput("filter_input_sample"), actionButton("filter_df_sample","Apply Filter(s)"), actionButton("reset_filter_df_sample", "Clear Filters"), conditionalPanel(condition= "!input.filter_df_sample",DT::dataTableOutput("qc_sample_table")), conditionalPanel("input.filter_df_sample",DT::dataTableOutput("qc_sample_table2"))),
                                                                                                                                              tabPanel(value="qc_view_sample_data","View Sample Level",hr(style = "border-top: 1px solid #ffffff;"), box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, h4("Select the columns you would like to view from the dropdown below."),uiOutput("qc_sample_header_output"),h4("Enter any filters you would like to apply using the notation in the example. Click 'Apply' anytime you change the inputs."), 
                                                                                                                                                                                                                                                         #div(style="display:inline-block;vertical-align:top;",
                                                                                                                                                                                                                                                         tags$style(HTML('#check_input_sample {margin-top: 25px}')),
                                                                                                                                                                                                                                                             splitLayout(cellWidths = c("60%","40%"),
                                                                                                                                                                                                                                                                  #uiOutput("filter_input_sample"), 
                                                                                                                                                                                                                                                                  textInput(inputId="user_filter_input_sample",label="Filter Dataframe",placeholder = "Wave == 'W3'",value=''),
                                                                                                                                                                                                                                                                  actionButton('check_input_sample', 'Check Input', class = 'action-buttons-app')), 
                                                                                                                                                                                                                                                         textOutput('input_check_result_sample'),
                                                                                                                                                                                                                                                         actionButton("filter_df_sample","Apply", class = 'action-buttons-app')),
                                                                                                                                                       conditionalPanel(condition= "!input.filter_df_sample",DT::dataTableOutput("qc_sample_table_orig")), conditionalPanel("input.filter_df_sample",DT::dataTableOutput("qc_sample_table"))),
                                                                                                                                               #tabPanel(value = "qc_view_sample_data", "View Sample Level", h4("Select the columns you would like to view from the dropdown below:"), panel(uiOutput("qc_sample_header_output"), uiOutput("qc_sample_new_filter"), status = "primary"), DT::dataTableOutput(outputId = "qc_sample_table")),
                                                                                                                                              tabPanel(value="qc_inspect_sample_data","Inspect Sample Level", hr(style = "border-top: 1px solid #ffffff;"), box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, h4("Inspect the data based on the following options:"),
                                                                                                                                                       h4(tags$li("cat: summary and comparison of categorical columns")),
                                                                                                                                                       h4(tags$li("cat_levels: summary and comparison of the most common levels in categorical columns")),
                                                                                                                                                       h4(tags$li("na: summary and comparison of the rate of missingness across columns")),
                                                                                                                                                       h4(tags$li("num: summary and comparison of numeric columns")),
                                                                                                                                                       h4(tags$li("types: summary and comparison of column types")),
                                                                                                                                                       uiOutput("inspect_sample_qc_select_column"),
                                                                                                                                                       #selectInput("column_name_to_inspect_qc_sample","Select one column name to inspect, or all columns:", choices= NULL, multiple= FALSE),
                                                                                                                                                       pickerInput("inspect_data_type_qc_sample","Select one of the options below to inspect the data:",choices=NULL, multiple=FALSE, selected=NULL, options = list(style = 'picker-inputs-app'))),
                                                                                                                                                       #uiOutput("inspect_df_qc_sample")),
                                                                                                                                                       plotOutput("inspect_sample_plot")))),
                                                                                      tabPanel(value="qc_subject","Subject Level", tabsetPanel(id="qc_subject_subtabs",
                                                                                                                                              tabPanel(value="qc_view_subject_data", "View Subject Level",hr(style = "border-top: 1px solid #ffffff;"), box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, h4("Select the columns you would like to view from the dropdown below."),uiOutput("qc_subject_header"), h4("Enter any filters you would like to apply using the notation in the example. Click 'Apply' anytime you change the inputs."),
                                                                                                                                                                                                                                                            tags$style(HTML('#check_input_subject {margin-top: 25px}')),
                                                                                                                                                                                                                                                            splitLayout(cellWidths = c("60%","40%"),
                                                                                                                                                                                                                                                                        #uiOutput("filter_input_subject"),
                                                                                                                                                                                                                                                                        textInput(inputId="user_filter_input_subject",label="Filter Dataframe",placeholder = "source_material == 'Lung'",value=""),
                                                                                                                                                                                                                                                                        actionButton('check_input_subject', 'Check Input', class = 'action-buttons-app')), 
                                                                                                                                                                                                                                                            textOutput('input_check_result_subject'),
                                                                                                                                                                                                                                                            actionButton("filter_df_subject","Apply", class = 'action-buttons-app')),conditionalPanel(condition= "!input.filter_df_subject",DT::dataTableOutput("qc_subject_table_orig")), conditionalPanel("input.filter_df_subject",DT::dataTableOutput("qc_subject_table"))),
                                                                                                                                              tabPanel(value="qc_inspect_subject_data","Inspect Subject Level", hr(style = "border-top: 1px solid #ffffff;"), box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, h4("Inspect the data based on the following options:"),
                                                                                                                                                       h4(tags$li("cat: summary and comparison of categorical columns")),
                                                                                                                                                       h4(tags$li("cat_levels: summary and comparison of the most common levels in categorical columns")),
                                                                                                                                                       h4(tags$li("na: summary and comparison of the rate of missingness across columns")),
                                                                                                                                                       h4(tags$li("num: summary and comparison of numeric columns")),
                                                                                                                                                       h4(tags$li("types: summary and comparison of column types")),
                                                                                                                                                       uiOutput("inspect_subject_qc_select_column"),
                                                                                                                                                       #uiOutput("inspect_df_qc_subject")),
                                                                                                                                                       pickerInput("inspect_data_type_qc_subject","Select one of the options below to inspect the data:",choices=NULL, multiple=FALSE, selected=NULL, options = list(style = 'picker-inputs-app'))),
                                                                                                                                                       plotOutput("inspect_subject_plot")))))),

               tabItem(tabName = "NGSpurity", h2("NGSpurity"), hr(),tabsetPanel(id="ngspurity_tabs",
                        tabPanel(value="ngs_view_data_qc","View Data QC",hr(style = "border-top: 1px solid #ffffff;"), box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, h4("Select the columns you would like to view from the dropdown below."), uiOutput("ngs_purity_header"), h4("Enter any filters you would like to apply using the notation in the example. Click 'Apply' anytime you change the inputs."),
                          tags$style(HTML('#check_input_ngspurity {margin-top: 25px}')),
                          splitLayout(cellWidths = c("60%","40%"),
                                      #uiOutput("filter_input_ngspurity"), 
                                      textInput(inputId="user_filter_input_ngspurity",label="Filter Dataframe",placeholder = "MCN_WGD == 'nWGD'",value=""),
                                      actionButton('check_input_ngspurity', 'Check Input', class = 'action-buttons-app')),
                          textOutput('input_check_result_ngspurity'),
                          actionButton("filter_df_ngspurity", "Apply", class = 'action-buttons-app')),
                          conditionalPanel(condition= "!input.filter_df_ngspurity",DT::dataTableOutput("qc_ngspurity_table_orig")), conditionalPanel("input.filter_df_ngspurity",DT::dataTableOutput("qc_ngspurity_table"))),
                        tabPanel(value="ngs_inspect_data","Inspect Data", hr(style = "border-top: 1px solid #ffffff;"), box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, h4("Inspect the data based on the following options. The meaning of each is as follows:"),
                          h4(tags$li("cat: summary and comparison of categorical columns")),
                          h4(tags$li("cat_levels: summary and comparison of the most common levels in categorical columns")),
                          h4(tags$li("na: summary and comparison of the rate of missingness across columns")),
                          h4(tags$li("num: summary and comparison of numeric columns")),
                          h4(tags$li("types: summary and comparison of column types")),
                          uiOutput("ngs_purity_header_inspect_tab"),
                          #uiOutput("inspect_df_test")),
                          pickerInput("inspect_data_type_ngs","Select one of the options below to inspect the data:" ,choices=NULL, multiple=FALSE, selected=NULL, options = list(style = 'picker-inputs-app'))),
                          tags$style(type="text/css",  ".shiny-output-error { visibility: hidden; }", ".shiny-output-error:before { visibility: hidden; }"),
                          plotOutput("inspect_ngs_plot", height="800px")),

                                  #if(os_detect() %in% c("Linux","Darwin")){
                                    tabPanel(value="ngs_view_figures","View Figures", hr(style = "border-top: 1px solid #ffffff;"), box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE,h4("Using the dropdown menus below, select the figure you would like to view. The figure will be updated as choices from the dropdown menus are changed."),
                                                                                                                                        #uiOutput("ngspurity_barcode"),
                                                                                                                                        selectInput("tumor_barcode_to_inspect","Select one Tumor Barcode to inspect:", choices= NULL, multiple= FALSE),
                                                                                                                                        #uiOutput("ngspurity_battenberg"),
                                                                                                                                        selectInput("battenberg_to_inspect","Select one Battenberg to inspect:", choices= NULL, multiple= FALSE),
                                                                                                                                        #uiOutput("ngspurity_type")),
                                                                                                                                        selectInput("type_to_inspect","Select one Type to inspect:", choices= NULL, multiple= FALSE)),
                                                                                                                                        #tags$head(tags$style(type="text/css","#figure img {max-width:100%; width:100%; height:auto}")),
                                                                                                                                        #uiOutput("figure_pdf"))
                                                                                                                                        imageOutput("figure_ngspurity",width = "50%",height = "50%"))
                                  #}else{
                                    # tabPanel(value="ngs_view_figures","View Figures", hr(style = "border-top: 1px solid #ffffff;"), box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, h4("Using the dropdown menus below, select the figure you would like to view. The figure will be updated as choices from the dropdown menus are changed."),
                                    #                                                                                                     #uiOutput("ngspurity_barcode"),
                                    #                                                                                                     selectInput("tumor_barcode_to_inspect","Select one Tumor Barcode to inspect:", choices= NULL, multiple= FALSE),
                                    #                                                                                                     #uiOutput("ngspurity_battenberg"),
                                    #                                                                                                     selectInput("battenberg_to_inspect","Select one Battenberg to inspect:", choices= NULL, multiple= FALSE),
                                    #                                                                                                     #uiOutput("ngspurity_type")),
                                    #                                                                                                     selectInput("type_to_inspect","Select one Type to inspect:", choices= NULL, multiple= FALSE)),
                                    #                                                                                                     #tags$head(tags$style(type="text/css","#figure img {max-width:100%; width:100%; height:auto}")),
                                    #                                                                                                     uiOutput("figure_pdf"))
                                   #},
                                    )),

               tabItem(tabName = "mutations", h2("Mutations"),hr(),tabsetPanel(id="mutation_tabs", 
                                                                          tabPanel(value = "mutation_summary_tab", "Mutation Summary", hr(style = "border-top: 1px solid #ffffff;"), 
                                                                                   box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, h4('Generate a mutation summary of one or more groups in the selected project using the dropdown list below. If you do not select a group, the mutation summary will be generated for the entire sample population of the selected project. Find brief explanations of each of the plots below.'), selectInput("mut_summary", label= "Group:", choices= NULL, multiple=TRUE),
                                                                                       actionButton("mut_summary_generate", "Get Mutation Summary", class = 'action-buttons-app'),
                                                                                   #hr(style = "border-top: 1px solid #ffffff;"),
                                                                                   h4(tags$li('Variant Classification: Total counts across all samples in selected group(s) for each variant classification defined in the project data.')),
                                                                                   h4(tags$li('Variant Type: Counts for each variant type across all samples in selected group(s) for each type defined in the project data.')),
                                                                                   h4(tags$li('SNV Class: Total counts across all samples in selected group(s) for each of the six mutation classes.')),
                                                                                   h4(tags$li('Variants per Sample: Counts of all variants per individual sample in the selected group(s). These are color coded to match the variant classification categories in the top left plot.')),
                                                                                   h4(tags$li('Variant Classification Summary: Similar to the top left plot which showed counts per variant classification, this plot displays variant classification in a boxplot format, providing additional information as to where samples fall for each variant classification category.')),
                                                                                   h4(tags$li('Top 10 Mutated Genes: Top 10 most mutated genes when considering the percentage of samples with at least one mutation in a given gene. The bars are color coded by variant classification, and a black coloration denotes that an individual had multiple mutation types for a given gene.'))),
                                                                                   fluidRow(column(10,
                                                                                   plotOutput("mut_summary_plot")),
                                                                                   column(3, style = "margin-top: 100px;",conditionalPanel(condition = 'input.mut_summary_generate', downloadButton('download_mut_summary', 'Download Plot', class = 'action-buttons-app'))))),
                                                                          #tabPanel( 
                                                                                   if(os_detect() %in% c("Linux","Darwin")){
                                                                                     tabPanel(value="tmb_comparison_tab","TMB Comparison", hr(style = "border-top: 1px solid #ffffff;"), 
                                                                                              box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, 
                                                                                              h4('Generate a tumor mutational burden (TMB) plot of TCGA cohorts against one group or the entire sample population in the project selected. Check the Group checkbox to select a group.'),
                                                                                              #checkboxInput("tmb_sp_group", "Group", value=FALSE),
                                                                                              prettyCheckbox("tmb_sp_group", "Group",value = FALSE, inline = TRUE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                                                                                              conditionalPanel(condition="input.tmb_sp_group == true", selectInput("tmb_group", label= "Select a group:", choices= NULL, multiple=FALSE)),
                                                                                              h4("Input a name to label the project's cohort on the plot. Otherwise, it will be labeled 'Cohort'."),textInput("tmb_cohort_name", "Cohort Name", value = "Cohort", placeholder = "Project Name or Group"),
                                                                                              actionButton("tmb_plot_generate", "Generate TMB Plot", class = 'action-buttons-app')),
                                                                                              tags$head(tags$style(type="text/css","#figure img {max-width:100%; width:100%; height:auto}")),imageOutput("figure_pdf_tmb",width = "50%",height = "50%"),
                                                                                              conditionalPanel(condition = "input.tmb_plot_generate", downloadButton("download_tmb_plot","Download Plot",class = 'action-buttons-app')))
                                                                                              
                                                                                   }else{
                                                                                     tabPanel(value="tmb_comparison_tab","TMB Comparison", hr(style = "border-top: 1px solid #ffffff;"), 
                                                                                              box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, 
                                                                                                  h4('Generate a tumor mutational burden (TMB) plot of TCGA cohorts against one group or the entire sample population in the project selected. Check the Group checkbox to select a group.'),
                                                                                                  #checkboxInput("tmb_sp_group", "Group", value=FALSE),
                                                                                                  prettyCheckbox("tmb_sp_group", "Group",value = FALSE, inline = TRUE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                                                                                                  conditionalPanel(condition="input.tmb_sp_group == true", selectInput("tmb_group", label= "Select a group:", choices= NULL, multiple=FALSE)),
                                                                                                  h4("Input a name to label the project's cohort on the plot. Otherwise, it will be labeled 'Cohort'."),textInput("tmb_cohort_name", "Cohort Name", value = "Cohort", placeholder = "Project Name or Group"),
                                                                                                  actionButton("tmb_plot_generate", "Generate TMB Plot", class = 'action-buttons-app')),
                                                                                                  tags$head(tags$style(type="text/css","#figure img {max-width:100%; width:100%; height:auto}")),imageOutput("figure_pdf_tmb",width = "50%",height = "50%"),
                                                                                                  conditionalPanel(condition = "input.tmb_plot_generate", actionButton("download_tmb_plot","Download Plot", class = 'action-buttons-app')))
                                                                                              
                                                                                   },
                                                                                   #conditionalPanel('input.tmb_plot_generate', actionButton('download_tmb_plot',"Download Plot"))),
                                                                          tabPanel(value="lolliplot_tab", "Lollipop Plot",hr(style = "border-top: 1px solid #ffffff;"),
                                                                                                       #tags$ol(tags$li("Begin by entering a gene name to create the plot for.",tags$i("Note: the gene name input is case-sensitive.")),
                                                                                                             #  tags$li("Select 'Group' if desired using the checkbox. A dropdown list of group choices will appear that is specific to the project and its corresponding data."),
                                                                                                             #  tags$li("Click 'Select Transcript' to select a transcript that includes data for the gene input", tags$i("You cannot select a transcript until entering a gene name.")),
                                                                                                              # tags$li("Check the 'Domain Annotations' box if you choose to include domain annotations on the plot."),
                                                                                                              # tags$li("Enter a minN value to determine which mutations on the plot are labeled. N is a number of a given mutation type at a given position, and therefore minN establishes a threshold for text labels on the plot."), style = "font-size: 18px"),
                                                                                   
                                                                                   box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE,
                                                                                       h4("Construct a lollipop mutation plot by following the steps below."),  
                                                                                   h4("Select 'Group' if desired using the checkbox. If you do not select a group the entire data cohort will be used."),
                                                                                   #checkboxInput("lolli_sp_group", "Group", value=FALSE),
                                                                                   prettyCheckbox("lolli_sp_group", "Group", value = FALSE,inline = TRUE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                                                                                   conditionalPanel(condition="input.lolli_sp_group == true", selectInput("lolli_group", label= "Select at least one group:", choices= NULL, multiple=TRUE)),
                                                                                   h4("Enter a gene name for which to create the plot.",tags$i("Note: The gene name input is case-sensitive.")),
                                                                                   textInput("lolli_gene_name", label= "Gene name", placeholder = "TP53"),
                                                                                   textOutput('check_gene_input'),
                                                                                   h4("Click 'Select Transcript' to select a transcript available for the gene entered.", tags$i("Note: You cannot select a transcript until entering a gene name.")),
                                                                                   actionButton("get_ts_info", "Select transcript", class = 'action-buttons-app'),
                                                                                   uiOutput("tslist_input"),
                                                                                  #selectInput("lolli_transcript", label= "Select a transcript name", choices=NULL, multiple= FALSE),
                                                                                  #checkboxInput("lolli_domain_annot", "Domain Annotations", value=FALSE),
                                                                                  h4("Check the 'Domain Annotations' box if you would like to include domain annotations on the plot."),
                                                                                  prettyCheckbox("lolli_domain_annot", "Domain Annotations", value = FALSE,inline = TRUE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                                                                                  h4("Enter a minN to select the mutations observed in at least minN samples."),
                                                                                  textInput("lolli_minN", label= "minN", value= 5),
                                                                                  actionButton("lolliplot_calculate", "Generate Lolliplot", class = 'action-buttons-app'),actionButton("lolliplot_reset","Reset", class = 'action-buttons-app')),
                                                                                  tabsetPanel(
                                                                                    tabPanel(value = 'lolliplot_plot', 'Lollipop Plot', plotOutput("lolliplot"),
                                                                                    conditionalPanel(condition = "input.lolliplot_calculate", downloadButton("download_lolliplot","Download Plot", class = 'action-buttons-app'))),
                                                                                    tabPanel(value='lollipop_table', 'Mutation Table', DT::dataTableOutput("lolliplot_table"))
                                                                                 )))),

               tabItem(tabName = "scna", h2("SCNA"), hr(),
                tabsetPanel(id = 'scna_tabs',
                  tabPanel(value="scna_data_tab", "Inspect Data",hr(style = "border-top: 1px solid #ffffff;"),
                       box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE,
                           h4("Use the table below to explore SCNA events in the selected project."), 
                           h4("Select the columns you would like to view from the dropdown below."),
                           uiOutput('scna_header_output'),
                           h4("Enter any filters you would like to apply using the notation in the example."),
                           tags$style(HTML('#check_input_scna {margin-top: 25px}')),
                           splitLayout(cellWidths = c("60%","40%"),
                                       #uiOutput("filter_input_sample"), 
                                       textInput(inputId="user_filter_input_scna",label="Filter Dataframe",placeholder = "chr == '1'",value=''),
                                       actionButton('check_input_scna', 'Check Input', class = 'action-buttons-app')),
                           textOutput('input_check_result_scna'),
                           actionButton("filter_scna_table","Apply/Reset", class = 'action-buttons-app')),
                            conditionalPanel('!input.filter_scna_table', DT::dataTableOutput('scna_table_orig')),conditionalPanel('input.filter_scna_table', DT::dataTableOutput('scna_table'))),
                  tabPanel(value = 'scna_clustering_tab', "Clustering",hr(style = "border-top: 1px solid #ffffff;"),
                           box(width = 12, title = 'Instructions', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, 
                               h4('Below is SCNA clustering output of the samples in the selected project. The first tab will contain a table of which cluster each sample belongs to,
                                  and in the second tab will be a heatmap.'),
                               h5(strong('Note: Clustering is only supported for data that is based on hg38.')),
                               selectizeInput("scna_samples", "Select at least three samples for SCNA clustering:", choices= NULL, multiple= TRUE),
                               actionButton('scna_cluster_calculate','Run Clustering Analysis', class = 'action-buttons-app')),
                           tabsetPanel(id = 'scna_clustering_results',
                                       tabPanel(value = 'clustering_table', 'Clustering Data Table', DT::dataTableOutput('scna_clustering_table')),
                                       tabPanel(value = 'clustering_plot', 'Clustering Plot', plotOutput('scna_clustering_heatmap', height = 1000),
                                                conditionalPanel(condition = 'input.scna_cluster_calculate', downloadButton("download_scna_clustering_plot","Download Plot", class = 'action-buttons-app'))))),
                                       tabPanel(value="scna_gistic_output","GISTIC Output", hr(style = "border-top: 1px solid #ffffff;"), box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE,
                                                                                   h4('Below users can explore GISTIC output. Select ‘Amplification’ or ‘Deletion’ to display the corresponding plot. Documentation for GISTIC, which identifies regions of amplification or deletion across samples, can be found', a(href = 'https://broadinstitute.github.io/gistic2/', 'here.')),
                                                                                   radioGroupButtons(
                                                                                     inputId = "gistic_out_select",
                                                                                     label = NULL,
                                                                                     choices = c("Amplification Plot", "Deletion Plot"),
                                                                                     status = "primary",
                                                                                     selected = 'Amplification Plot'
                                                                                   )),
                                   tags$head(tags$style(type="text/css","#figure img {max-width:100%; width:100%; height:auto}")),imageOutput("scna_gistic_plot",width = "50%",height = "50%"))

                        # if(os_detect() %in% c("Linux","Darwin")){
                        #   tabPanel(value="scna_gistic_output","GISTIC Output", hr(style = "border-top: 1px solid #ffffff;"), box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, 
                        #                                                            h4('Below users can explore GISTIC output. Select ‘Amplification’ or ‘Deletion’ to display the corresponding plot. Documentation for GISTIC, which identifies regions of amplification or deletion across samples, can be found', a(href = 'https://broadinstitute.github.io/gistic2/', 'here.')),
                        #                                                            radioGroupButtons(
                        #                                                              inputId = "gistic_out_select",
                        #                                                              label = NULL, 
                        #                                                              choices = c("Amplification Plot", "Deletion Plot"),
                        #                                                              status = "primary",
                        #                                                              selected = 'Amplification Plot'
                        #                                                            )),
                        #            tags$head(tags$style(type="text/css","#figure img {max-width:100%; width:100%; height:auto}")),imageOutput("scna_gistic_plot",width = "50%",height = "50%"))
                        # }else{
                        #   tabPanel(value="scna_gistic_output","GISTIC Output", hr(style = "border-top: 1px solid #ffffff;"), box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, 
                        #                                                            h4('Below users can explore GISTIC output. Select ‘Amplification’ or ‘Deletion’ to display the corresponding plot. Documentation for GISTIC, which identifies regions of amplification or deletion across samples, can be found', a(href = 'https://broadinstitute.github.io/gistic2/', 'here.')), 
                        #                                                            radioGroupButtons(
                        #                                                              inputId = "gistic_out_select",
                        #                                                              label = NULL, 
                        #                                                              choices = c("Amplification Plot", "Deletion Plot Test"),
                        #                                                              status = "primary",
                        #                                                              selected = 'Amplification Plot'
                        #                                                            )),
                        #            tags$head(tags$style(type="text/css","#figure img {max-width:100%; width:100%; height:auto}")),uiOutput("scna_gistic_plot"))
                        # }
                  )),

               tabItem(tabName = "sv", h2("SV"), hr(),
                       tabsetPanel(id = 'sv_tabs',
                                   tabPanel(value = 'reconplot', 'ReConPlot', hr(style = "border-top: 1px solid #ffffff;"),
                                            box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, h4('Use this module to generate a ReConPlot (REarrangement and COpy Number Plot). Find more information about the ReConPlot package', a(href = 'https://github.com/cortes-ciriano-lab/ReConPlot', 'here.')),
                                                pickerInput("recon_tumor_barcode","Select one Tumor Barcode:", choices= NULL, multiple= FALSE, options= list(style = 'picker-inputs-app', `live-search`=TRUE, `dropup-auto` = FALSE)),
                                                selectizeInput("recon_genome_build","Select a genome build:", choices = c("hg19", "hg38", "T2T"), selected = "hg38", multiple = FALSE),
                                                pickerInput("recon_chr_select","Select one chromosome:", choices= NULL, multiple= FALSE, options= list(style = 'picker-inputs-app', `live-search`=TRUE, `dropup-auto` = FALSE)),
                                                textInput("recon_chr_start", "Enter a starting genomic coordinate for the plot"),
                                                textInput("recon_chr_end", "Enter an ending genomic coordinate for the plot"),
                                                actionButton("generate_recon_plot", 'Generate ReConPlot', class = 'action-buttons-app')),
                                                fluidRow(column(10,
                                                plotOutput('recon_plot', height = '750px')),
                                                column(3,conditionalPanel('input.generate_recon_plot', downloadButton('download_sv_recon_plot', 'Download Plot', class = 'action-buttons-app'))))))),

               tabItem(tabName = "mutational_signature", h2("Mutational Signatures"), hr(), box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, 
                                                                                                h4('This module allows for easy viewing of mutational signature analysis results.'),
                                                                                                h4('The mutational signature analysis was completed using the', a(href = 'https://github.com/AlexandrovLab/SigProfilerClusters', 'SigProfilerClusters tool'), 'developed by the Alexandrov lab. The paper for the tool can be found', a(href = 'https://academic.oup.com/bioinformatics/article/38/13/3470/6589887?login=false', 'here.')),
                                                                                                h4("Select a tumor barcode from the dropdown list below to display a mutational signature analysis plot for the selected sample."),selectizeInput("clustered_mut_barcodes", "Select one Tumor Barcode:", choices= NULL, multiple= FALSE)), imageOutput("figure_pdf_clustered_mut")),

               tabItem(tabName = "genomic_landscape", h2("Genomic Landscape"),hr(),
               tabsetPanel(id="genomic_landscape_tabs",
                           tabPanel(value="genome_plot_figures", "View genomePlot", hr(style = "border-top: 1px solid #ffffff;"), box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, h4("This module allows for exploration of the genomic landscape of samples from the selected project."),
                              h4("The circos plots displayed below were generated using the",code("genomePlot()", style = "font-size: 16px"), "function from the Signature Tools Lib R package, which is used for mutational signature analysis. The plots display somatic variants across the genome, including single nucleotide variants (SNV), small insertions and deletions (indels), copy number variants (CNV), and rearrangements."),
                              h4("Additional information about the Signature Tools Lib package and its functions can be found on the ",a("Signature Tools Lib GitHub", href="https://github.com/Nik-Zainal-Group/signature.tools.lib"), "page.
                                      The paper,", a(em("A practical framework and online tool for mutational signature analyses show intertissue variation and driver dependencies"), href="https://www.nature.com/articles/s43018-020-0027-5"), ", also contains a lot of useful information regarding use of the package."),
                              #uiOutput("genomePlot_Barcodes")
                              h4('Select a tumor barcode from the dropdown list below to display the genomic landscape plots for the selected sample.'),
                              pickerInput("tumor_barcode_to_inspect_genomePlot","Select one Tumor Barcode:", choices= NULL, multiple= FALSE, options= list(style = 'picker-inputs-app', pickerOptions(liveSearch=TRUE, dropupAuto = FALSE)))),
                              imageOutput("genomePlot_figure")))),

               tabItem(tabName = "clonal_evolution", h2("Clonal Evolution"),hr(),
                       box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, 
                                                                                            h4('This module allows for the exploration of mutation time for samples in a selected project.'),
                                                                                            h4("The plots displayed below were generated using the", a("MutationTimeR R package", href="https://github.com/Nik-Zainal-Group/signature.tools.lib"), code("plotSample()", style = "font-size: 16px"), "function."),
                                                                                            h4("Select one tumor barcode from the dropdown below to explore clonal evolution of samples in the selected project, with regard to VAF, copy number, and mutation time."),
                                                                                            #uiOutput("mutation_time_barcode"),
                                                                                            pickerInput("tumor_barcode_mutation_time","Select one Tumor Barcode:", choices= NULL, multiple= FALSE, options= list(pickerOptions(liveSearch=TRUE, dropupAuto = FALSE), style = 'picker-inputs-app'))),
                                                                                            tags$head(tags$style(type="text/css","#figure img {max-width:100%; width:100%; height:auto}")),imageOutput("figure_pdf_mutation_time",width = "50%",height = "50%")),

                       # if(os_detect() %in% c("Linux","Darwin")){
                       #   tabPanel(value="mutationTime_view_figures","View Figures", box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, 
                       #                                                                  h4('This module allows for the exploration of mutation time for samples in a selected project.'),
                       #                                                                  h4("The plots displayed below were generated using the", a("MutationTimeR R package", href="https://github.com/Nik-Zainal-Group/signature.tools.lib"), code("plotSample()", style = "font-size: 16px"), "function."),
                       #                                                                  h4("Select one tumor barcode from the dropdown below to explore clonal evolution of samples in the selected project, with regard to VAF, copy number, and mutation time."),
                       #                                                                  #uiOutput("mutation_time_barcode"),
                       #                                                                  pickerInput("tumor_barcode_mutation_time","Select one Tumor Barcode:", choices= NULL, multiple= FALSE, options= list(pickerOptions(liveSearch=TRUE, dropupAuto = FALSE), style = 'picker-inputs-app'))),
                       #                                                                  tags$head(tags$style(type="text/css","#figure img {max-width:100%; width:100%; height:auto}")),imageOutput("figure_pdf_mutation_time",width = "50%",height = "50%"))
                       # }else{
                       #   tabPanel(value="mutationTime_view_figures","View Figures", box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, h4("Select one tumor barcode from the dropdown below to explore clonal evolution of samples in the selected project, with regard to VAF, copy number, and mutation time."), 
                       #                                                                  #uiOutput("mutation_time_barcode")
                       #                                                                  pickerInput("tumor_barcode_mutation_time","Select one Tumor Barcode:", choices= NULL, multiple= FALSE, options= list(pickerOptions(liveSearch=TRUE, dropupAuto = FALSE), style = 'picker-inputs-app'))),
                       #                                                                  tags$head(tags$style(type="text/css","#figure img {max-width:100%; width:100%; height:auto}")),uiOutput("figure_pdf_mutation_time"))
                       # }
                       #),

               tabItem(tabName = "survival_analysis", h2("Survival Analysis"), hr(), box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, h4("This module allows the user to conduct a survival analysis using the genomic alteration data available for the selected project."),
                       # tags$ol(tags$li("Select one genomic alteration from the dropdown list."),
                       #         tags$li("Select a specific group (if desired) by checking the 'Group' box, and then selecting a group from the dropdown list. Note: Selecting a group is not required."),
                       #         tags$li("Finally, select the reference group, which is used in the analysis as a control for the survival comparison. A 'No' reference means that patients whose mutation status was 'No' will be used as the reference. A 'Yes' reference means that patients whose mutation status was 'Yes' will be used as a reference."), style = "font-size: 18px"),
                       selectizeInput("vartmp_options_select_survival","Select one genomic alteration:", choices= NULL, multiple= TRUE, options = list(maxItems=1, placeholder="Select a genomic alteration"), selected= NULL),
                       #checkboxInput("sp_group_checkbox_survival", "SP_Group", value=FALSE),conditionalPanel(condition= "input.sp_group_checkbox_survival ==true",uiOutput("sp_group_choices_survival")),
                       h4("Select a specific group (if desired) by checking the 'Group' box, and then selecting a group from the dropdown list."),
                       prettyCheckbox("sp_group_checkbox_survival", "Group", value = FALSE,inline = TRUE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE), conditionalPanel(condition= "input.sp_group_checkbox_survival ==true",uiOutput("sp_group_choices_survival")),
                       #radioButtons("reference_survival", "Reference", choices = c("No","Yes"), selected = "No"),
                       #prettyRadioButtons(inputId = "reference_survival", label = "Reference", choices = c("No", "Yes"), selected = "No", inline = TRUE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                       h4("Select the reference group, which is used in the analysis as a control for the survival comparison. A 'No' reference means that patients whose mutation status was 'No' will be used as the reference. A 'Yes' reference means that patients whose mutation status was 'Yes' will be used as a reference."),
                       # radioGroupButtons(
                       #   inputId = "reference_survival",
                       #   label = 'Reference', 
                       #   choices = c("No", "Yes"),
                       #   status = "primary"
                       #   #selected = "No"
                       # ),
                       actionButton('get_reference_levels', 'Get Reference Levels', class = 'action-buttons-app'),
                       textOutput('reference_levels_out'),
                       textInput('enter_reference_level', 'Enter Reference Level', value = ''),
                       textOutput('reference_level_check'),
                       #conditionalPanel(condition = 'input.get_reference_levels', textInput('enter_reference_level', 'Enter Reference Level', value = '')),
                       # checkboxInput("keyname_checkbox_survival", "Keyname", value=FALSE),
                       #conditionalPanel(condition= "input.keyname_checkbox_survival ==true",selectizeInput("keyname_select_survival","Select one type of genomic alteration to inspect:", choices= NULL, multiple= FALSE, options = list(maxItems=1, placeholder="Select a genomic alteration"), selected= NULL)),
                       # textInput("filename_survival", "Output Filename:", placeholder= "survival_plot.pdf"),
                       actionButton("calculate_survival", "Calculate", class = 'action-buttons-app'), actionButton("survival_reset","Reset", class = 'action-buttons-app')), 
                       tableOutput("survival_value"),
                       plotOutput("survival_plot", height= 500, width= 500), 
                       #uiOutput("download_survival_plot"),
                       conditionalPanel('input.calculate_survival', downloadButton('download_survival_plot', 'Download Plot', class = 'action-buttons-app'))),

               tabItem(tabName = "integrative_analysis", h2("Integrative Analysis"), hr(), tabsetPanel(id="integrative_analysis",
                                                                                                 tabPanel(value="fisher_test","Enrichment Analysis", hr(style = "border-top: 1px solid #ffffff;"), box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, 
                                                                                                 h4("This tab allows you to conduct an enrichment test (fisher or glm) between genomic alterations. After entering the desired inputs, the 'Results' tab will contain the table for the test, and the 'Volcano Plot' will help to visualize results from the table. If selecting a second variable (i.e. Variable 2 Name),
                                                                                                    an additional 'Bar Plot' tab will appear, containing a bar chart of the genomic alterations selected."),
                                                                                                    selectizeInput("vartmp_options_select","Genomic Alteration I", choices= NULL, multiple= FALSE),#, options = list(maxItems=1, placeholder="Select a genomic alteration"), selected= NULL),
                                                                                                    #conditionalPanel(condition= "input.vartmp_options_select",
                                                                                                    #checkboxInput("sp_group_checkbox", "SP_Group", value=FALSE),
                                                                                                    prettyCheckbox("fisher_samplelist_checkbox", label = 'Sample List', value = FALSE,inline = TRUE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                                                                                                    #bsTooltip('fisher_samplelist_checkbox','Test', placement = "top", trigger = "click"),
                                                                                                    #shinyBS::bsTooltip(id = 'fisher_samplelist_checkbox', title = 'Test', placement = 'top', trigger ='hover'),
                                                                                                    # txt file only?
                                                                                                    conditionalPanel(condition= "input.fisher_samplelist_checkbox== true",fileInput("fisher_samplelist", "Select a file including only Tumor Barcodes to include in the analysis:",multiple= FALSE, accept = 'txt', placeholder = "Choose File")),
                                                                                                    prettyCheckbox("sp_group_checkbox", "Group", value = FALSE,inline = TRUE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                                                                                                    conditionalPanel(condition= "input.sp_group_checkbox ==true",uiOutput("sp_group_choices")),
                                                                                                    #checkboxInput("fisher_samplelist_checkbox", "Sample List", value= FALSE),
                                                                                                    #checkboxInput("var2name_checkbox", "Variable 2 Name", value=FALSE),
                                                                                                    prettyCheckbox("var2name_checkbox", "Genomic Alteration II", value = FALSE,inline = TRUE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                                                                                                    conditionalPanel(condition= "input.var2name_checkbox== true",  selectizeInput("vartmp_options_select_var2","Select an additional genomic alteration to inspect:", choices= NULL, multiple= TRUE, options = list(maxItems=1, placeholder="Select a genomic alteration"), selected= NULL)),
                                                                                                    #checkboxInput("excludes_checkbox", "Exclude Variable(s)", value= FALSE),
                                                                                                    h5(strong('Exclude parameters:')),
                                                                                                    prettyCheckbox("excludes_checkbox", "Exclude Variable(s)", value = FALSE,inline = TRUE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                                                                                                    #conditionalPanel(condition= "input.excludes_checkbox== true", selectizeInput("fisher_excludes", "Select one or more variables to exclude in the fisher test:", choices =NULL, multiple =TRUE, selected= NULL)),
                                                                                                    conditionalPanel(condition= "input.excludes_checkbox== true", textInput("fisher_excludes", "Enter one or more variables to exclude from the analysis:",placeholder = 'Gender|Overall_Feature, WGD_Status|Overall_Feature')),
                                                                                                    #checkboxInput("excludes_cat_checkbox", "Exclude Source Categories", value =FALSE),
                                                                                                    prettyCheckbox("excludes_cat_checkbox", "Exclude Source Categories", value = FALSE,inline = TRUE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                                                                                                    #conditionalPanel("input.excludes_cat_checkbox ==true", selectizeInput("fisher_excludes_cat", "Select one or more source categories to exclude from the fisher test:", choices= NULL, multiple= TRUE, selected= NULL)),
                                                                                                    conditionalPanel("input.excludes_cat_checkbox ==true", textInput("fisher_excludes_cat", "Enter one or more source categories to exclude from the analysis:", placeholder = 'Overall_Feature')),
                                                                                                    h5(strong('Note: These cannot be the same variables/source categories within any genomic alterations selected, or the same as any of the Keep parameters used.')),
                                                                                                    #checkboxInput("keeps_checkbox", "Keep Variables", value =FALSE),
                                                                                                    h5(strong('Keep parameters:')),
                                                                                                    prettyCheckbox("keeps_checkbox", "Keep Variables", value = FALSE,inline = TRUE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                                                                                                    #conditionalPanel("input.keeps_checkbox ==true", selectizeInput("fisher_keeps", "Select one or more variables to keep in the fisher test:", choices= NULL, multiple= TRUE, selected= NULL)),
                                                                                                    conditionalPanel("input.keeps_checkbox ==true", textInput("fisher_keeps", "Enter one or more variables to exclude from the analysis:", placeholder = 'Kataegis|Overall_Feature')),
                                                                                                    #checkboxInput("keeps_cat_checkbox", "Keep Source Categories", value =FALSE),
                                                                                                    prettyCheckbox("keeps_cat_checkbox", "Keep Source Categories", value = FALSE,inline = TRUE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                                                                                                    #conditionalPanel("input.keeps_cat_checkbox ==true", selectizeInput("fisher_keeps_cat", "Select one or more source categories to keep in the fisher test:", choices= NULL, multiple= TRUE, selected= NULL)),
                                                                                                    conditionalPanel("input.keeps_cat_checkbox ==true", textInput("fisher_keeps_cat", "Enter one or more source categories to exclude from the analysis:", placeholder = 'Mutation_Driver')),
                                                                                                    h5(strong('Note: These cannot be the same as any of the Exclude parameters used.')),
                                                                                                    textInput("fisher_min_freq","Minimum Frequency",value= 0.03),
                                                                                                    selectInput("fisher_freq_colnames", "Frequency Column Name", choices= NULL, selected= NULL, multiple= FALSE),
                                                                                                    #radioButtons("fisher_test_type", "Select a type of test:", choices= c("fisher.test","glm"), selected= "fisher.test", inline= TRUE),
                                                                                                    #prettyRadioButtons(inputId = "fisher_test_type", label = "Select a type of test:", choices= c("fisher.test","glm"), selected = "fisher.test", inline = TRUE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                                                                                                    radioGroupButtons(
                                                                                                      inputId = "fisher_test_type",
                                                                                                      label = "Select a type of test:", 
                                                                                                      choices = c("fisher.test", "glm"),
                                                                                                      status = "primary",
                                                                                                      #selected = character(0)
                                                                                                    ),
                                                                                                   conditionalPanel(condition= "input.fisher_test_type== 'glm'", textInput("fisher_glm_input", "Enter a glm formula", value= "Var1 ~ Var2 + Gender"),strong(textOutput("fisher_glm_message"))),
                                                                                                    br(),
                                                                                                    textInput("fisher_fdr_cutoff", "FDR cutoff", value= 0.1),
                                                                                                    actionButton("calculate_fisher","Calculate", class = 'action-buttons-app')),
                                                                                                    tabsetPanel(id="fisher_results", tabPanel("Results",DT::dataTableOutput("fisher_output_table")),
                                                                                                                tabPanel("Figure 1: Volcano Plot", plotOutput("fisher_output_plots1",height=750, width = 1000), conditionalPanel(condition = 'input.calculate_fisher', downloadButton("download_fisher_volc_plot", "Download Plot", class = 'action-buttons-app'))),
                                                                                                                tabPanel("Figure 2: Bar Plot", plotOutput("fisher_output_plots3",height=500, width = 750), conditionalPanel(condition = 'input.calculate_fisher', downloadButton("download_fisher_bar_plot", "Download Plot", class = 'action-buttons-app'))))),
                                                                                                 
                                                                                                 
                                                                                                 tabPanel(value="fisher_bar", "Fisher Bar Plot", hr(style = "border-top: 1px solid #ffffff;"), box(width = 12, title = 'Instructions and Inputs', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE,
                                                                                                   selectizeInput("vartmp_options_select_bar","Genomic Alteration I", choices= NULL, multiple= FALSE, options = list(maxItems=1, placeholder="Select a genomic alteration"), selected= NULL),
                                                                                                   conditionalPanel(condition= "input.vartmp_options_select_bar",
                                                                                                        selectizeInput("vartmp_options_select_var2_bar","Genomic Alteration II", choices= NULL, multiple= FALSE, options = list(maxItems=1, placeholder="Select a genomic alteration"), selected= NULL),
                                                                                                        #checkboxInput("fisher_bar_samplelist_checkbox", "Sample List", value= FALSE),
                                                                                                        prettyCheckbox("fisher_bar_samplelist_checkbox", label = 'Sample List', value = FALSE,inline = FALSE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                                                                                                        conditionalPanel(condition= "input.fisher_bar_samplelist_checkbox== true",fileInput("fisher_bar_samplelist", "Select a file including only Tumor Barcodes to include in the analysis:",multiple= FALSE, accept = 'txt', placeholder = "Choose File")),
                                                                                                        #checkboxInput("sp_group_checkbox_bar", "Group", value=FALSE), conditionalPanel(condition= "input.sp_group_checkbox_bar ==true",uiOutput("sp_group_choices_bar")),
                                                                                                        prettyCheckbox("sp_group_checkbox_bar", label = 'Group', value = FALSE,inline = FALSE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                                                                                                        conditionalPanel(condition= "input.sp_group_checkbox_bar ==true",uiOutput("sp_group_choices_bar")),
                                                                                                        actionButton("generate_fisher_barplot", "Generate Fisher Barplot",class = 'action-buttons-app'))), 
                                                                                                        fluidRow(column(10,
                                                                                                        plotOutput("fisher_output_bar"),
                                                                                                        column(3, style = "margin-top: 75px;",conditionalPanel(condition = "input.generate_fisher_barplot", downloadButton("download_fisher_bar_plot_only","Download Plot", class = 'action-buttons-app')))))),

                                                                                                 tabPanel(value="association_test_main", "Association Testing",hr(style = "border-top: 1px solid #ffffff;"), box(width = 12, title = 'Instructions and Data Selection', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE,
                                                                                                                                                                   h4('Select a maximum of two datasets from the list below to conduct bivariable or multivariable association analyses. After selecting the data, the Data Selected box will populate with the data selected and its dimensions.'),
                                                                                                                                                                   uiOutput("data_list_association"),actionButton("load_datasets", "Load Datasets", class = 'action-buttons-app')),
                                                                                                                                                                    box(width = 12, title = 'Data Selected', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE,conditionalPanel(condition="input.load_datasets",tags$header(strong("Selected Data Dimensions")), tableOutput("assoc_datatable_dim"), tags$header(strong(textOutput("assoc_data_join_by"))), textOutput("assoc_no_tum_barcode"), uiOutput("assoc_dropdown")), DT::dataTableOutput("assoc_datatable")),
                                                                                                    tabsetPanel(id="association_testing",
                                                                                                      tabPanel(value= "bivar_analysis","Bivariable Analysis",
                                                                                                        #conditionalPanel(condition= "input.group_var_bivariable == false",
                                                                                                          uiOutput("assoc_variable_list_1"),
                                                                                                          #checkboxInput("filter_assoc_var1","Filter (>0)",value=FALSE),
                                                                                                          prettyCheckbox("filter_assoc_var1", label = "Filter (>0)", value = FALSE,inline = FALSE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                                                                                                          #checkboxInput("log2_assoc_var1", "log2", value=FALSE),
                                                                                                          prettyCheckbox("log2_assoc_var1", label = "log2", value = FALSE,inline = FALSE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),  
                                                                                                          #textInput("collapse_assoc_var1", "Collapse Level",value=""),
                                                                                                          textInput("collapse_assoc_var1", "Collapse Level",value=""),
                                                                                                          # tags$style(HTML('#check_level_var1 {margin-top: 25px}')),
                                                                                                          # splitLayout(cellWidths = c("20%","80%"),
                                                                                                          #             textInput("collapse_assoc_var1", "Collapse Level",value="")), 
                                                                                                          #             #actionButton('check_level_var1', 'Check Level', class = 'action-buttons-app')),
                                                                                                          #             #actionButton('check_level_var1_clear','Clear', class = 'action-buttons-app')), 
                                                                                                          textOutput('check_collapse_var1_output'),
                                                                                                          #conditionalPanel("input.collapse_assoc_var1 != ''",textOutput('check_collapse_var1_output')),
                                                                                                          #conditionalPanel('input.check_level_var1_clear', textOutput('check_collapse_var1_output_clear')),
                                                                                                          uiOutput("assoc_variable_list_2"),
                                                                                                          #checkboxInput("filter_assoc_var2","Filter (>0)",value=FALSE),
                                                                                                          prettyCheckbox("filter_assoc_var2", label = "Filter (>0)", value = FALSE,inline = FALSE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),
                                                                                                          #checkboxInput("log2_assoc_var2", "log2", value=FALSE),
                                                                                                          prettyCheckbox("log2_assoc_var2", label = "log2", value = FALSE,inline = FALSE, status = "primary",shape = "curve", fill = TRUE, bigger = TRUE),    
                                                                                                          #textInput("collapse_assoc_var2", "Collapse Level",value=""),
                                                                                                          textInput("collapse_assoc_var2", "Collapse Level",value=""), 
                                                                                                          # tags$style(HTML('#check_level_var2 {margin-top: 25px}')),
                                                                                                          # splitLayout(cellWidths = c("20%","80%"),
                                                                                                          #           textInput("collapse_assoc_var2", "Collapse Level",value=""), 
                                                                                                          #           actionButton('check_level_var2', 'Check Level', class = 'action-buttons-app')), 
                                                                                                          textOutput('check_collapse_var2_output'),
                                                                                                          #pickerInput("assoc_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes","fisher","skit"), multiple=FALSE),
                                                                                                          pickerInput("assoc_types_list", "Select a test type for the association test you would like to run:", choices=NULL, multiple=FALSE, options = list(style = 'picker-inputs-app')),
                                                                                                          #textInput("group_var_input_bivariable","Group Variable", value=""),
                                                                                                          tags$style(HTML('#check_group_var_bivar {margin-top: 25px}')),
                                                                                                          splitLayout(cellWidths = c("20%","80%"),
                                                                                                                    textInput("group_var_input_bivariable", "Group Variable",value=""), 
                                                                                                                    actionButton('check_group_var_bivar', 'Check Variable', class = 'action-buttons-app')),
                                                                                                          textOutput('check_group_var_bivar_output'),
                                                                                                          #checkboxInput("association_output_plot", "Save plot to computer",value = FALSE),
                                                                                                          #textInput("file_ext_assoc", "File Extension (png, svg, or jpg)", value="png"),
                                                                                                          actionButton("calculate_association","Calculate", class = 'action-buttons-app'),actionButton("reset_association", "Reset", class = 'action-buttons-app'),conditionalPanel(condition="input.group_var_input_bivariable ==''",plotOutput("bivariable_plot", width="800px")), conditionalPanel(condition="input.group_var_input_bivariable !=''", dataTableOutput("bivariable_table")),conditionalPanel('input.calculate_association', downloadButton('download_bivar_result', "Download result", class = 'action-buttons-app'))),
                                                                                                        # conditionalPanel(condition= "input.group_var_bivariable == true",uiOutput("group_var_text_box_bivariable"),uiOutput("assoc_variable_list_1_group"),
                                                                                                        #   checkboxInput("filter_assoc_var1_group","Filter (>0)",value=FALSE),
                                                                                                        #   checkboxInput("log2_assoc_var1_group", "log2", value=FALSE),
                                                                                                        #   textInput("collapse_assoc_var1_group", "Collapse Level",value=NULL),
                                                                                                        #   uiOutput("assoc_variable_list_2_group"),
                                                                                                        #   checkboxInput("filter_assoc_var2_group","Filter (>0)",value=FALSE),
                                                                                                        #   checkboxInput("log2_assoc_var2_group", "log2", value=FALSE),
                                                                                                        #   textInput("collapse_assoc_var2_group", "Collapse Level",value=NULL),
                                                                                                        #   pickerInput("assoc_types_list_group", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes","fisher","skit"), multiple=FALSE),
                                                                                                        #   actionButton("calculate_association_group","Calculate"),actionButton("reset_association_group", "Reset"),DT::dataTableOutput("bivariable_grp_result"))),
                                                                                                      
                                                                                                      tabPanel(value= "multivar_analysis","Multivariable Analysis",
                                                                                                               #textInput("regression_formula", "Regression Formula", placeholder="lm( Conpair_Concordance ~ Wave)"),
                                                                                                               tags$style(HTML('#check_multivar_formula {margin-top: 25px}')),
                                                                                                               splitLayout(cellWidths = c("20%","80%"),
                                                                                                                           textInput("regression_formula", "Regression Formula", placeholder="lm( Conpair_Concordance ~ Wave)"), 
                                                                                                                           actionButton('check_multivar_formula', 'Check Formula', class = 'action-buttons-app')),
                                                                                                               h5(strong("Supported regression model types include: lm and glm.")),
                                                                                                               textOutput('check_multivar_formula_output'), 
                                                                                                               tags$style(HTML('#check_group_var_multivar {margin-top: 25px}')),
                                                                                                               splitLayout(cellWidths = c("20%","80%"),
                                                                                                                           textInput("group_var_input", "Group Variable",value=""), 
                                                                                                                           actionButton('check_group_var_multivar', 'Check Variable', class = 'action-buttons-app')),
                                                                                                               textOutput('check_group_var_multivar_output'),
                                                                                                               #textInput(inputId= "group_var_input", label= 'Group Variable',placeholder= NULL),
                                                                                                               #checkboxInput(inputId="group_var_regression", label= "Group_Var", value= FALSE),
                                                                                                               # conditionalPanel(condition="input.group_var_regression == true", uiOutput("group_var_text_box"),
                                                                                                               #                  uiOutput("regression_model_group"),
                                                                                                               #                  uiOutput("formula_input_group"),
                                                                                                               #                  actionButton("calculate_regression_group","Calculate"),actionButton("reset_regression_group", "Reset"),DT::dataTableOutput("regression_grp_result")),
                                                                                                               # conditionalPanel(condition= "input.group_var_regression == false",uiOutput("formula_input"),
                                                                                                               #                  checkboxInput("regression_output_plot", "Save plot to computer",value = FALSE),
                                                                                                               #                  textInput("file_ext_regression", "File Extension (png, svg, or jpg)", value= "png"),
                                                                                                                                actionButton("calculate_regression","Calculate", class = 'action-buttons-app'),actionButton("reset_regression", "Reset", class = 'action-buttons-app'),conditionalPanel(condition="input.group_var_input ==''",plotOutput("regression_plot")), conditionalPanel(condition="input.group_var_input !=''", dataTableOutput("regression_table")),conditionalPanel('input.calculate_regression', downloadButton('download_multivar_result', "Download result", class = 'action-buttons-app'))))),
                                                                                                tabPanel(value= "oncoplot", "Oncoplot",hr(style = "border-top: 1px solid #ffffff;"), box(width = 12, title = 'Instructions and Data Selection', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE,
                                                                                                         radioButtons("oncoplot_input_selection", "Select an input type for generating the oncoplot.", choices = c("By genomic alteration" = "oncoplot_genalts" , "By genomic alteration text input" = "oncoplot_text_input_genalts", "By source category" = "oncoplot_categories","By source category feature frequency" = "oncoplot_cat_freq_inputs"), selected =character(0)),
                                                                                                         conditionalPanel(condition= "input.oncoplot_input_selection == 'oncoplot_genalts'", selectizeInput("oncoplot_genalts_select","Select at least one type of genomic alteration for oncoplot generation:", choices= NULL, multiple= TRUE, options = list(placeholder="Select one or more genomic alteration"), selected= NULL), checkboxInput("order_by_input1", "Order by Input", value = FALSE)),
                                                                                                         
                                                                                                         conditionalPanel(condition= "input.oncoplot_input_selection == 'oncoplot_text_input_genalts'", textAreaInput("oncoplot_user_input_genalts", "Write in a specific list of genomic alterations for oncoplot generation", placeholder="Gender|Overall_Feature \n Smoking|Overall_Feature"), checkboxInput("order_by_input2", "Order by Input", value = FALSE)),
                                                                                                         conditionalPanel(condition= "input.oncoplot_input_selection == 'oncoplot_categories'", selectizeInput("oncoplot_categories_select", "Select one or more source categories to include in the oncoplot generation:", choices= NULL, multiple= FALSE, options = list(placeholder="Select one or more source categories"), selected= NULL)),
                                                                                                         conditionalPanel(condition= "input.oncoplot_input_selection == 'oncoplot_cat_freq_inputs'", textAreaInput("oncoplot_cat_freq_select","Select altetations with one source category by frequency", placeholder="Overall_Feature,0.3")),
                                                                                                         textInput("oncoplot_min_freq","Minimum Frequency",value= 0.10),
                                                                                                         actionButton("generate_oncoplot","Generate Oncoplot", class = 'action-buttons-app')), 
                                                                                                         verticalLayout(plotOutput("data_output_plot"),
                                                                                                        conditionalPanel('input.generate_oncoplot', downloadButton("download_oncoplot", "Download Plot", class = 'action-buttons-app')))),
                                                                                                         #uiOutput("download_oncoplot"))),
                                                                                                         #conditionalPanel(condition = 'input.generate_oncoplot', uiOutput("download_oncoplot"))),
                                                                                                tabPanel(value="wordcloud", "Genomic Features Wordcloud", hr(style = "border-top: 1px solid #ffffff;"), box(width = 12, title = 'Instructions and Data Selection', status = 'primary', solidHeader = TRUE, collapsed = FALSE, collapsible = TRUE, 
                                                                                                                                                                                                            selectizeInput("wordcloud_tumor_barcode", "Select one tumor barcode from the dropdown below:", choices= NULL, multiple= FALSE), 
                                                                                                                                                                                                            actionButton("generate_wordcloud", "Generate Wordcloud", class = 'action-buttons-app')), 
                                                                                                                                                                                                            verticalLayout(plotOutput("sherlock_wordcloud_plot"),
                                                                                                                                                                                                            conditionalPanel('input.generate_wordcloud', downloadButton('download_wordcloud', 'Download Wordcloud', class = 'action-buttons-app')))))),

               tabItem(tabName = "documentation", h2("Documentation"), hr(), tabsetPanel(id= "documentation_tabs", tabPanel(title= "Data Requirement Info", value= "data_req_info", h4("Below you will find information regarding data requirements for using your data throughout the modules in the app. Data requirements vary by module."),
                                                                                                                      pickerInput('module_list', label="Modules", choices=c("Manifest Information", "NGSpurity", "Mutations","SCNA","SV","Mutational Signatures", "Genomic Landscape","Clonal Evolution","Survival Analysis", "Integrative Analysis"), options = list(style = 'picker-inputs-app')),
                                                                                                                      #actionButton("select_module_for_info",label = "Select"),
                                                                                                                      pickerInput("submodule_user_select", label="Submodules",choices=NULL, options = list(style = 'picker-inputs-app')), uiOutput("file_list"), actionButton('select_example_file', label = 'Select File', class = 'action-buttons-app'), actionButton('clear_selected_example_file', 'Reset', class = 'action-buttons-app'), textOutput('req_col_out'), uiOutput("example_file_column_names"),dataTableOutput("example_file")),
                                                                                                             tabPanel(title= "Project Info",value="project_info",h4("The following are brief descriptions of each project included in the app. In parentheses is the project code as it appears in the selection dropdown list."), includeMarkdown("Project_descriptions.Rmd")), 
                                                                                                             tabPanel(title="Module Info",value="module_info",h4("The following are brief descriptions of each module included in the app."),dataTableOutput("module_info_table")))),

               tabItem(tabName = "faq", h2("Frequently Asked Questions"), hr(), includeMarkdown("FAQ_Sherlock_Genome.Rmd")))
    ),
  )
)