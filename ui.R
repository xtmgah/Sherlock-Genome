
#Sherlock-Genome: A R Shiny App for Genomic Analysis and Visualization

#setwd("/Users/kleinam/sherlock_genome/Sherlock-Genome/tmp")
#load('sherlock_landscape_v2.RData')

# keep adding to this vector when additional packages are required to run the app 
packages_req <- c("shiny","shinydashboard", "markdown","shinyjs", "tibble","inspectdf","DT","dplyr","shinyWidgets","stringr","PMCMRplus","ggstatsplot","ggplot2", "ggpubr","cowplot","data.table","forcats","broom","rstudioapi","fs") #ggpubr- ggqqplot()

# check for required packages and install those not installed; load all packages required
# lapply() function to use the packages_req vector and carry out function 
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

# detect OS (Windows vs. Mac)
os_detect <- function(){
  os_name <- Sys.info()[['sysname']]
  return(os_name)
}

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
    dashboardSidebar(sidebarMenu(id = "sbmenu",
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
      menuItem("Integrative Analysis", tabName = "integrative_analysis"),
      menuItem("Documentation", tabName= "documentation"),
      menuItem("Frequently Asked Questions",tabName = "faq")
    )),
# each of the following are what the tabs assigned page will display when the tab is clicked
    dashboardBody(
      tags$script(HTML("
        var openTab = function(tabName){
          $('a', $('.sidebar')).each(function() {
            if(this.getAttribute('data-value') == tabName) {
              this.click()
            };
          });
        }
      ")),
      tabItems(tabItem(tabName ="data_load",h2("Data Load"),h4("Choose one of the two methods below to upload data into the corresponding modules in the left panel:"), h5("1. Select a project, and all of the files available for that project. This will automatically populate each file into its corresponding module in the left panel."),
                       h5("2. Select a project, specific files you wish to investigate further, and then the files you selected will be populated into their corresponding module on the left panel."), h5("Click", a("here", onclick = "openTab('documentation')", href="#"), "for documentation about each project and module available in the app."),
                       selectInput("project_code","Begin by selecting a project:", c("BRCA_HK", "mWGS","Sherlock"), selected= "Sherlock"),actionButton("select", "Select Project"), actionButton("reset_project", "Select Different Project"),conditionalPanel(condition="input.select",uiOutput("file_list_select"),actionButton("choose_files_all","Select All Files"),actionButton("choose_files_ind","Select Individual Files")),conditionalPanel(condition="input.choose_files_all || input.choose_files_ind",actionButton(inputId="submit", label="Submit"), actionButton(inputId="reset", label="Reset")), textOutput("files_selected")),
               
               tabItem(tabName = "study_overview", h2("Study Overview"),uiOutput("study_overview")),
               
               tabItem(tabName = "sample_qc", h2("Sample QC"),tabsetPanel(id="qc_tabs",tabPanel(value="qc_sample", "Sample Data", tabsetPanel(id="qc_sample_subtabs", tabPanel(value="qc_view_sample_data","View Sample Data",h4("Select the columns you would like to view from the dropdown below:"),uiOutput("qc_sample_header"),DT::dataTableOutput("qc_sample_table")),tabPanel(value="qc_inspect_sample_data","Inspect Sample Data",uiOutput("inspect_sample_qc_select_column"),uiOutput("inspect_df_qc_sample"),plotOutput("inspect_sample_plot")))),tabPanel(value="qc_subject","Subject Data",tabsetPanel(id="qc_subject_subtabs", tabPanel(value="qc_view_subject_data", "View Subject Data",h4("Select the columns you would like to view from the dropdown below:"),uiOutput("qc_subject_header"),DT::dataTableOutput("qc_subject_table")),tabPanel(value="qc_inspect_subject_data","Inspect Subject Data",uiOutput("inspect_subject_qc_select_column"),uiOutput("inspect_df_qc_subject"),plotOutput("inspect_subject_plot")))))),
               
               tabItem(tabName = "NGSpurity", h2("NGSpurity"), tabsetPanel(id="ngspurity_tabs",tabPanel(value="ngs_view_sample_data","View Sample Data",h4("Select the columns you would like to view from the dropdown below:"),uiOutput("ngs_purity_header"),DT::dataTableOutput("ngs_purity_table")),tabPanel(value="ngs_inspect_data","Inspect Data", h4("Select filters to inspect the data based on the following options. The meaning of each is as follows:"), tags$li("cat: summary and comparison of categorical columns"), tags$li("cat_levels: summary and comparison of the most common levels in categorical columns"),
                                  tags$li("na: summary and comparison of the rate of missingness across dataframe columns"), tags$li("num: summary and comparison of numeric columns"), tags$li("types: summary and comparison of column types"),uiOutput("ngs_purity_header_inspect_tab"),uiOutput("inspect_df_test"),tags$style(type="text/css",  ".shiny-output-error { visibility: hidden; }", ".shiny-output-error:before { visibility: hidden; }"),plotOutput("inspect_ngs_plot", height="800px")),

                                  if(os_detect() %in% c("Linux","Darwin")){
                                    tabPanel(value="ngs_view_figures","View Figures", h4("Using the dropdown menus below, select the figures you would like to view. The figures will be updated as choices from the dropdown menus are changed."),uiOutput("ngspurity_barcode"),uiOutput("ngspurity_battenberg"),uiOutput("ngspurity_type"),tags$head(tags$style(type="text/css","#figure img {max-width:100%; width:100%; height:auto}")),imageOutput("figure_pdf",width = "50%",height = "50%"))
                                  }else{
                                    tabPanel(value="ngs_view_figures","View Figures", h4("Using the dropdown menus below, select the figures you would like to view. The figures will be updated as choices from the dropdown menus are changed."),uiOutput("ngspurity_barcode"),uiOutput("ngspurity_battenberg"),uiOutput("ngspurity_type"),tags$head(tags$style(type="text/css","#figure img {max-width:100%; width:100%; height:auto}")),uiOutput("figure_pdf"))
                                  },

                                  tabPanel(value="ngs_association_testing","Association Testing", DT::dataTableOutput("ngs_variable_table"), tags$head(
                                  tags$style(type="text/css", "#inline label{ display: table-cell; text-align: center; vertical-align: middle; } #inline .form-group { display: table-row;}")),
                                  tabsetPanel(id="variable_assoc_analysis", 
                                  
                                  tabPanel("Multivariable Analysis","Supported regression model types include: lm and glm.",checkboxInput(inputId="group_var_regression", label= "Group_Var", value= FALSE), 
                                           
                                           conditionalPanel(condition="input.group_var_regression == true", uiOutput("group_var_text_box"), uiOutput("regression_model_group"), uiOutput("formula_input_group"),actionButton("calculate_ngs_regression_group","Calculate"),actionButton("reset_ngs_regression_group", "Reset"),DT::dataTableOutput("regression_grp_result")),
                                           conditionalPanel(condition= "input.group_var_regression == false",uiOutput("formula_input"),actionButton("calculate_ngs_regression","Calculate"),actionButton("reset_ngs_regression", "Reset"),plotOutput("ngspurity_regression_plot"))),
                                  
                                  tabPanel("Bivariable Analysis",checkboxInput(inputId="group_var_bivariable", label= "Group_Var", value= FALSE),
                                           
                                           conditionalPanel(condition= "input.group_var_bivariable == false",div(style="display: inline-block;vertical-align:top; width: 150px;",uiOutput("ngspurity_variable_list_1")),div(style="display: inline-block;vertical-align:top; width: 150px;",
                                              checkboxInput("filter_ngs_var1","Filter (>0)",value=FALSE)),div(style="display: inline-block;vertical-align:top; width: 150px;",checkboxInput("log2_ngs_var1", "log2", value=FALSE)), tags$head(tags$style("label{font-weight: 400;}")),tags$div(id = "inline",textInput("collapse_ngs_var1", "Collapse Level",value=NULL)), 
                                              uiOutput("ngspurity_variable_list_2"), checkboxInput("filter_ngs_var2","Filter (>0)",value=FALSE), checkboxInput("log2_ngs_var2", "log2", value=FALSE),textInput("collapse_ngs_var2", "Collapse Level",value=NULL),pickerInput("ngspurity_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes","fisher","skit"), multiple=FALSE),
                                              checkboxInput("ngs_output_plot", "Save plot to computer",value = FALSE),textInput("file_ext", "File Extension (png, svg, or jpg)", value="png"),actionButton("calculate_ngs_association","Calculate"),actionButton("reset_ngs_association", "Reset"),plotOutput("ngspurity_bivariable_plot", width="800px")),
                                           
                                           conditionalPanel(condition= "input.group_var_bivariable == true",uiOutput("group_var_text_box_bivariable"),div(style="display: inline-block;vertical-align:top; width: 150px;",uiOutput("ngspurity_variable_list_1_group")),div(style="display: inline-block;vertical-align:top; width: 150px;",checkboxInput("filter_ngs_var1_group","Filter (>0)",value=FALSE)),div(style="display: inline-block;vertical-align:top; width: 150px;",checkboxInput("log2_ngs_var1_group", "log2", value=FALSE)), tags$head(tags$style("label{font-weight: 400;}")),tags$div(id = "inline",textInput("collapse_ngs_var1_group", "Collapse Level",value=NULL)), 
                                                uiOutput("ngspurity_variable_list_2_group"), checkboxInput("filter_ngs_var2_group","Filter (>0)",value=FALSE), checkboxInput("log2_ngs_var2_group", "log2", value=FALSE),textInput("collapse_ngs_var2_group", "Collapse Level",value=NULL),
                                                pickerInput("ngspurity_types_list_group", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes","fisher","skit"), multiple=FALSE),actionButton("calculate_ngs_association_group","Calculate"),actionButton("reset_ngs_association_group", "Reset"),DT::dataTableOutput("bivariable_grp_result"))))))),
                             
               tabItem(tabName = "mutations", h2("Mutations")),
               
               tabItem(tabName = "scna", h2("SCNA")),
               
               tabItem(tabName = "sv", h2("SV")),
               
               tabItem(tabName = "mutational_signature", h2("Mutational Signatures")),
               
               tabItem(tabName = "genomic_landscape", h2("Genomic Landscape"), tabsetPanel(id="genomic_landscape_tabs", tabPanel(value="genome_plot_figures", "View Figures",uiOutput("genomePlot_Barcodes"), imageOutput("genomePlot_figure")))),
               
               tabItem(tabName = "clonal_evolution", h2("Clonal Evolution")),
               
               tabItem(tabName = "survival_analysis", h2("Survival Analysis")),
               
               tabItem(tabName = "integrative_analysis", h2("Integrative Analysis")),
               
               tabItem(tabName = "documentation", h2("Documentation"), tabsetPanel(id= "documentation_tabs", tabPanel(title= "Project Info",value="project_info",h4("The following are brief descriptions of each project included in the app. In parentheses is the project code as it appears in the selection dropdown list."), includeMarkdown("Project_descriptions.Rmd")), tabPanel(title="Module Info",value="module_info",h4("The following are brief descriptions of each module included in the app."),dataTableOutput("module_info_table")))),
               
               tabItem(tabName = "faq", h2("Frequently Asked Questions")))
    ),
  )
  
)

