
#Sherlock-Genome: A R Shiny App for Genomic Analysis and Visualization

# will need to be changed later
#setwd("/Users/kleinam/sherlock_genome/Sherlock-Genome/tmp")
#load('sherlock_landscape_v2.RData')

# specify the packages needed for the app through a character vector
# keep adding to this vector when additional packages are required to run the app successfully
packages_req <- c("shiny","shinydashboard", "markdown","shinyjs", "tibble","inspectdf","DT","dplyr","shinyWidgets","stringr","PMCMRplus","ggstatsplot","ggplot2", "ggpubr","cowplot","data.table") #ggpubr- ggqqplot()

# check for required packages and install those not installed
# lapply() function to use the packages_req vector and carry out the function written
pkg_check_fcn <-lapply(packages_req,function(x){
# if a package in the packages_req vector is not installed, install it
    if(!(x %in% installed.packages())){
      install.packages(x, dependencies = TRUE, repos=structure(c(CRAN="http://cloud.r-project.org/")))
# load the package after installing it
      library(x, character.only = TRUE)
# else statement for when the packages are already installed
    }else{
# load the package already installed into the workspace
      library(x, character.only = TRUE)
      # print(paste0(x,":" , "Package already installed"))
    }

})

# ui setup 
ui <- fluidPage(
# title at the top of the page; windowTitle = what the browser tab reads
  useShinyjs(),
  titlePanel(title = "Sherlock-Genome: A R Shiny App for Genomic Analysis and Visualization",
             windowTitle = "Sherlock-Genome"),

  # dashboard setup with left sidebar
  dashboardPage(skin = "blue",
    dashboardHeader(title = "Sherlock-Genome"),
    
# each of the following will be a tab on the left sidebar
    dashboardSidebar(sidebarMenu(
      menuItem("Data Load", tabName = "data_load"),
      menuItem("Study Overview", tabName = "study_overview"),
      menuItem("Sample QC", tabName = "sample_qc"),
      menuItem("NGSpurity", tabName = "NGSpurity"),
      menuItem("Mutations", tabName = "mutations"),
      menuItem("SCNA", tabName = "scna"),
      menuItem("SV", tabName = "sv"),
      menuItem("Mutational Signature", tabName = "mutational_signature"),
      menuItem("Genomic Landscape", tabName = "genomic_landscape"),
      menuItem("Clonal Evolution", tabName = "clonal_evolution"),
      menuItem("Survival Analysis", tabName = "survival_analysis"),
      menuItem("Integrative Analysis", tabName = "integrative_analysis")
    )),
# each of the following are what the tabs assigned page will read when the tab is clicked
    dashboardBody(
      tabItems(tabItem(tabName ="data_load", h3("Data Load"),h4("Choose one of the two methods below to upload data into the corresponding modules in the left panel:"), h5("1. Select a project code and all files available for that specific project will be populated automatically into the left panel corresponding to the appropriate module."),
                       h5("2. Select a project code, view all of the files available for the project selected, and then choose the file(s) you would like to load."),
                       selectInput("project_code","Begin by selecting a project code:", c("BRCA_HK", "mWGS","Sherlock"), selected= "Sherlock"),conditionalPanel(condition= "input.project_code",actionButton("select", "Select Project Code")),conditionalPanel(condition="input.select",uiOutput("file_list_select"),actionButton("choose_files_all","Select All Files"),actionButton("choose_files_ind","Select Individual Files")),conditionalPanel(condition="input.choose_files_all || input.choose_files_ind",actionButton(inputId="submit", label="Submit"), actionButton(inputId="reset", label="Reset")), textOutput("files_selected")),
               tabItem(tabName = "study_overview", h2("Study Overview"),uiOutput("study_overview")),
               tabItem(tabName = "sample_qc", h2("Sample QC"),uiOutput("test")),
               tabItem(tabName = "NGSpurity", h2("NGSpurity"), tabsetPanel(tabPanel("View Sample Data",h4("Select the columns you would like to view from the dropdown below:"),uiOutput("ngs_purity_header"),DT::dataTableOutput("ngs_purity_header_b")),tabPanel("Inspect Data", h4("Select filters to inspect the data based on the following options. The meaning of each is as follows:"), tags$li("cat: summary and comparison of categorical columns"), tags$li("cat_levels: summary and comparison of the most common levels in categorical columns"),
                                  tags$li("na: summary and comparison of the rate of missingness across dataframe columns"), tags$li("num: summary and comparison of numeric columns"), tags$li("types: summary and comparison of column types"),uiOutput("ngs_purity_header_inspect_tab"),uiOutput("inspect_df_test"),tags$style(type="text/css",
                                                                                                                                                                                                                                                                                                                                  ".shiny-output-error { visibility: hidden; }",
                                                                                                                                                                                                                                                                                                                                  ".shiny-output-error:before { visibility: hidden; }"),plotOutput("testprint", height="800px")),
                                  tabPanel("View Figures", h4("Using the dropdown menus below, select the figures you would like to view. The figures will be updated as choices from the dropdown menus are changed."),uiOutput("ngspurity_barcode"),uiOutput("ngspurity_battenberg"),uiOutput("ngspurity_type"),tags$head(tags$style(type="text/css","#figure img {max-width:100%; width:100%; height:auto}")),imageOutput("figure")),tabPanel("Association Testing", h4("Select two variables to be used in an association test."), uiOutput("ngspurity_variable_list_1"),uiOutput("ngspurity_variable_list_2"),actionButton(inputId="submit1",label="Submit"),verbatimTextOutput("shapiro_test"),plotOutput("scatter_corr")))),
               tabItem(tabName = "mutations", h2("Mutations")),
               tabItem(tabName = "scna", h2("SCNA")),
               tabItem(tabName = "sv", h2("SV")),
               tabItem(tabName = "mutational_signature", h2("Mutational Signatures")),
               tabItem(tabName = "genomic_landscape", h2("Genomic Landscape")),
               tabItem(tabName = "clonal_evolution", h2("Clonal Evolution")),
               tabItem(tabName = "survival_analysis", h2("Survival Analysis")),
               tabItem(tabName = "integrative_analysis", h2("Integrative Analysis"),tabsetPanel(tabPanel("Oncoplot",h3("Generate an onco")),tabPanel("Test"))))
    ),
  )
  
)

