
#Sherlock-Genome: A R Shiny App for Genomic Analysis and Visualization
#Error: Unable to retrieve package records for the following packages:

# For SKIT:
# install.packages("/Users/kleinam/Downloads/skit-0.0.2.tar.gz", repos = NULL, type="source")

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
packages_req <- c("shiny","shinydashboard", "hrbrthemes", "rlist","tidyverse", "markdown","shinyjs", "ggside","inspectdf","DT","shinyWidgets","PMCMRplus","ggsci", "ggstatsplot", "ggpubr","cowplot","data.table","broom","rstudioapi","fs","scales", "survival", "survminer", "extrafontdb","maftools") #ggpubr- ggqqplot()

# check for required packages and install those not installed; load all packages required
# lapply() function to use the packages_req vector and carry out function
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

# extrafont::loadfonts()

# detect OS (Windows vs. Mac)
os_detect <- function(){
  os_name <- Sys.info()[['sysname']]
  return(os_name)
}

inline <-  function (x) {
  tags$div(style="display:inline-block;", x)
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
      tabItems(tabItem(tabName ="data_load",h2("Data Load"),h4("Choose one of the two methods below to upload data into the corresponding modules in the left panel:"), 
                       h5("1. Select a project, and all of the files available for that project. This will automatically populate each file into its corresponding module in the left panel."),
                       h5("2. Select a project, specific files you wish to investigate further, and then the files you selected will be populated into their corresponding module on the left panel."), 
                       h5("Click", a("here", onclick = "openTab('documentation')", href="#"), "for documentation about each project and module available in the app."),
                       selectInput("project_code","Begin by selecting a project:", c("BRCA_HK", "mWGS","Sherlock"), selected= "Sherlock"),
                       actionButton("select", "Select Project"), 
                       actionButton("reset_project", "Select Different Project"),
                       conditionalPanel(condition="input.select",uiOutput("file_list_select"),actionButton(inputId="submit", label="Submit")),textOutput("file_load_message")),
               
               tabItem(tabName = "study_overview", h2("Study Overview"),uiOutput("study_overview_rmd")),
               
               tabItem(tabName = "sample_qc", h2("Sample QC"),tabsetPanel(id="qc_tabs",tabPanel(value="qc_sample", "Sample Data", tabsetPanel(id="qc_sample_subtabs", tabPanel(value="qc_view_sample_data","View Sample Data",h4("Select the columns you would like to view from the dropdown below:"),uiOutput("qc_sample_header_output"),uiOutput("filter_input_sample"), actionButton("filter_df_sample","Apply Filter(s)"), conditionalPanel(condition= "input.user_filter_input_sample == ''",DT::dataTableOutput("qc_sample_table")), conditionalPanel("input.user_filter_input_sample !=''",DT::dataTableOutput("qc_sample_table2"))),tabPanel(value="qc_inspect_sample_data","Inspect Sample Data",uiOutput("inspect_sample_qc_select_column"),uiOutput("inspect_df_qc_sample"),plotOutput("inspect_sample_plot")))),tabPanel(value="qc_subject","Subject Data", tabsetPanel(id="qc_subject_subtabs", tabPanel(value="qc_view_subject_data", "View Subject Data",h4("Select the columns you would like to view from the dropdown below:"),uiOutput("qc_subject_header"), uiOutput("filter_input_subject"),actionButton("filter_df_subject","Apply Filter(s)"), conditionalPanel(condition= "input.user_filter_input_subject == ''",DT::dataTableOutput("qc_subject_table")), conditionalPanel("input.user_filter_input_subject !=''",DT::dataTableOutput("qc_subject_table2"))),tabPanel(value="qc_inspect_subject_data","Inspect Subject Data",uiOutput("inspect_subject_qc_select_column"),uiOutput("inspect_df_qc_subject"),plotOutput("inspect_subject_plot")))))),
               
               tabItem(tabName = "NGSpurity", h2("NGSpurity"), tabsetPanel(id="ngspurity_tabs",
                        tabPanel(value="ngs_view_sample_data","View Sample Data",h4("Select the columns you would like to view from the dropdown below:"),
                          uiOutput("ngs_purity_header"),
                          uiOutput("filter_input"),
                          actionButton("filter_df","Apply Filter(s)"),
                          conditionalPanel(condition="input.user_filter_input !='' ",DT::dataTableOutput("ngs_purity_table2")),
                          conditionalPanel(condition="input.user_filter_input =='' ",DT::dataTableOutput("ngs_purity_table1"))),
                        tabPanel(value="ngs_inspect_data","Inspect Data", h4("Select filters to inspect the data based on the following options. The meaning of each is as follows:"), 
                          tags$li("cat: summary and comparison of categorical columns"), 
                          tags$li("cat_levels: summary and comparison of the most common levels in categorical columns"),
                          tags$li("na: summary and comparison of the rate of missingness across dataframe columns"), 
                          tags$li("num: summary and comparison of numeric columns"), 
                          tags$li("types: summary and comparison of column types"),
                          uiOutput("ngs_purity_header_inspect_tab"),uiOutput("inspect_df_test"),tags$style(type="text/css",  ".shiny-output-error { visibility: hidden; }", ".shiny-output-error:before { visibility: hidden; }"),
                          plotOutput("inspect_ngs_plot", height="800px")),

                                  if(os_detect() %in% c("Linux","Darwin")){
                                    tabPanel(value="ngs_view_figures","View Figures", h4("Using the dropdown menus below, select the figures you would like to view. The figures will be updated as choices from the dropdown menus are changed."),uiOutput("ngspurity_barcode"),uiOutput("ngspurity_battenberg"),uiOutput("ngspurity_type"),tags$head(tags$style(type="text/css","#figure img {max-width:100%; width:100%; height:auto}")),imageOutput("figure_pdf",width = "50%",height = "50%"))
                                  }else{
                                    tabPanel(value="ngs_view_figures","View Figures", h4("Using the dropdown menus below, select the figures you would like to view. The figures will be updated as choices from the dropdown menus are changed."),uiOutput("ngspurity_barcode"),uiOutput("ngspurity_battenberg"),uiOutput("ngspurity_type"),tags$head(tags$style(type="text/css","#figure img {max-width:100%; width:100%; height:auto}")),uiOutput("figure_pdf"))
                                  },)),
               
               tabItem(tabName = "mutations", h2("Mutations"),tabsetPanel(id="mutation_tabs", tabPanel(value="lolliplot_tab", "Lolliplot", 
                                                                                                       h4("Select the inputs below to construct a lolliplot"), 
                                                                                                       textInput("lolli_gene_name", label= "Enter a gene name", placeholder = "TP53"),
                                                                                                       checkboxInput("lolli_sp_group", "Group", value=FALSE),
                                                                                                       conditionalPanel(condition="input.lolli_sp_group == true", selectInput("lolli_group", label= "Select a group", choices= NULL, multiple=TRUE)),
                                                                                                       actionButton("get_ts_info", "Select transcript"),
                                                                                                       uiOutput("tslist_input"),
                                                                                                       #selectInput("lolli_transcript", label= "Select a transcript name", choices=NULL, multiple= FALSE),
                                                                                                       textInput("lolli_minN", label= "Enter a minN", value= 5),
                                                                                                       actionButton("lolliplot_calculate", "Generate Lolliplot"),
                                                                                                       plotOutput("lolliplot"),
                                                                                                       DT::dataTableOutput("lolliplot_table")))),
               
               tabItem(tabName = "scna", h2("SCNA")),
               
               tabItem(tabName = "sv", h2("SV")),
               
               tabItem(tabName = "mutational_signature", h2("Mutational Signatures")),
               
               tabItem(tabName = "genomic_landscape", h2("Genomic Landscape"), 
               tabsetPanel(id="genomic_landscape_tabs", 
                           tabPanel(value="genome_plot_figures", "View genomePlot", 
                              h5("The circos plots you can explore below were made using the genomePlot() function from the Signature Tools Lib R package. This is a package used for mutational signature analysis."),h5("In the plots are somatic
                  variant across the genome, displaying single nucleotide variants (SNV), small insertions and deletions (indels), copy number variants (CNV), and rearrangements."),tagList("Additional information about the 
                  Signature Tools Lib package and its functions can be found here:",a("Signature Tools Lib GitHub", href="https://github.com/Nik-Zainal-Group/signature.tools.lib"), ". The following paper also contains a lot of useful information regarding use of the package: ", a(em("A practical framework and online tool for mutational signature analyses show intertissue variation and driver dependencies"), href="https://www.nature.com/articles/s43018-020-0027-5"), "."),
                              uiOutput("genomePlot_Barcodes"),imageOutput("genomePlot_figure")))),
               
               tabItem(tabName = "clonal_evolution", h2("Clonal Evolution"), 
                       if(os_detect() %in% c("Linux","Darwin")){
                         tabPanel(value="mutationTime_view_figures","View Figures", h4("Using the dropdown menus below, select the figures you would like to view. The figures will be updated as the tumor barcode selected changes."),uiOutput("mutation_time_barcode"),tags$head(tags$style(type="text/css","#figure img {max-width:100%; width:100%; height:auto}")),imageOutput("figure_pdf_mutation_time",width = "50%",height = "50%"))
                       }else{
                         tabPanel(value="mutationTime_view_figures","View Figures", h4("Using the dropdown menus below, select the figures you would like to view. The figures will be updated as the tumor barcode selected changes"), uiOutput("mutation_time_barcode"),tags$head(tags$style(type="text/css","#figure img {max-width:100%; width:100%; height:auto}")),uiOutput("figure_pdf_mutation_time"))
                       }),
               
               tabItem(tabName = "survival_analysis", h2("Survival Analysis"), 
                       selectizeInput("vartmp_options_select_survival","Select one type of genomic alteration to inspect:", choices= NULL, multiple= FALSE, options = list(maxItems=1, placeholder="Select a genomic alteration"), selected= NULL),
                       checkboxInput("sp_group_checkbox_survival", "SP_Group", value=FALSE),conditionalPanel(condition= "input.sp_group_checkbox_survival ==true",uiOutput("sp_group_choices_survival")),
                       radioButtons("reference_survival", "Reference", choices = c("No","Yes"), selected = "No"), 
                       # checkboxInput("keyname_checkbox_survival", "Keyname", value=FALSE),
                       conditionalPanel(condition= "input.keyname_checkbox_survival ==true",selectizeInput("keyname_select_survival","Select one type of genomic alteration to inspect:", choices= NULL, multiple= FALSE, options = list(maxItems=1, placeholder="Select a genomic alteration"), selected= NULL)), 
                       # textInput("filename_survival", "Output Filename:", placeholder= "survival_plot.pdf"), 
                       actionButton("calculate_survival", "Calculate"), tableOutput("survival_value"),
                       plotOutput("survival_plot", height= 500, width= 500), uiOutput("download_survival_plot")),
               
               tabItem(tabName = "integrative_analysis", h2("Integrative Analysis"), tabsetPanel(id="integrative_analysis",
                                                                                                 tabPanel(value="fisher_test","Enrichment Analysis", h5("This tab allows you to conduct a fisher test between genomic alterations. After entering the desired inputs, the 'Results Table' tab will contain the table for the test, and the 'Result Plots' tab will include all of the different plots generated."),
                                                                                                    selectizeInput("vartmp_options_select","Select one type of genomic alteration to inspect:", choices= NULL, multiple= FALSE, options = list(maxItems=1, placeholder="Select a genomic alteration"), selected= NULL),
                                                                                                    conditionalPanel(condition= "input.vartmp_options_select",
                                                                                                    checkboxInput("sp_group_checkbox", "SP_Group", value=FALSE),
                                                                                                    conditionalPanel(condition= "input.sp_group_checkbox ==true",uiOutput("sp_group_choices")),
                                                                                                    checkboxInput("fisher_samplelist_checkbox", "Sample List", value= FALSE),
                                                                                                    # txt file only?
                                                                                                    conditionalPanel(condition= "input.fisher_samplelist_checkbox== true",fileInput("fisher_samplelist", "Select a file including only Tumor Barcodes to include in the analysis:",multiple= FALSE, accept = 'txt', placeholder = "Choose File")),
                                                                                                    checkboxInput("var2name_checkbox", "Variable 2 Name", value=FALSE),
                                                                                                    conditionalPanel(condition= "input.var2name_checkbox== true",  selectizeInput("vartmp_options_select_var2","Select one type of genomic alteration to inspect:", choices= NULL, multiple= FALSE, options = list(maxItems=1, placeholder="Select a genomic alteration"), selected= NULL)),
                                                                                                    checkboxInput("excludes_checkbox", "Exclude Variable(s)", value= FALSE),
                                                                                                    conditionalPanel(condition= "input.excludes_checkbox== true", selectizeInput("fisher_excludes", "Select one or more variables to exclude in the fisher test:", choices =NULL, multiple =TRUE, selected= NULL)),
                                                                                                    checkboxInput("excludes_cat_checkbox", "Exclude Source Categories", value =FALSE),
                                                                                                    conditionalPanel("input.excludes_cat_checkbox ==true", selectizeInput("fisher_excludes_cat", "Select one or more source categories to exclude from the fisher test:", choices= NULL, multiple= TRUE, selected= NULL)),
                                                                                                    checkboxInput("keeps_checkbox", "Keep Variables", value =FALSE),
                                                                                                    conditionalPanel("input.keeps_checkbox ==true", selectizeInput("fisher_keeps", "Select one or more variables to keep in the fisher test:", choices= NULL, multiple= TRUE, selected= NULL)),
                                                                                                    checkboxInput("keeps_cat_checkbox", "Keep Source Categories", value =FALSE),
                                                                                                    conditionalPanel("input.keeps_cat_checkbox ==true", selectizeInput("fisher_keeps_cat", "Select one or more source categories to keep in the fisher test:", choices= NULL, multiple= TRUE, selected= NULL)),
                                                                                                    textInput("fisher_min_freq","Minimum Frequency",value= 0.03),
                                                                                                    selectInput("fisher_freq_colnames", "Frequency Column Name", choices= NULL, selected= NULL, multiple= FALSE),
                                                                                                    radioButtons("fisher_test_type", "Select a type of test:", choices= c("fisher.test","glm"), selected= "fisher.test", inline= TRUE),
                                                                                                    conditionalPanel(condition= "input.fisher_test_type== 'glm'", textInput("fisher_glm_input", "Enter a glm formula", value= "Var1 ~ Var2 + Gender"),strong(textOutput("fisher_glm_message"))),
                                                                                                    br(),
                                                                                                    textInput("fisher_fdr_cutoff", "FDR cutoff", value= 0.1),
                                                                                                    actionButton("calculate_fisher","Calculate")),
                                                                                                    tabsetPanel(id="fisher_results", tabPanel("Results",DT::dataTableOutput("fisher_output_table")), 
                                                                                                                tabPanel("Figure 1: Volcano Plot", plotOutput("fisher_output_plots1",height=750, width = 1000), uiOutput("download_fisher_volcano")),
                                                                                                                tabPanel("Figure 2: Bar plot", plotOutput("fisher_output_plots3",height=500, width = 750), uiOutput("download_fisher_barplot")))),
                                                                                                 tabPanel(value="fisher_bar", "Fisher Bar Plot",
                                                                                                   selectizeInput("vartmp_options_select_bar","Select one type of genomic alteration to inspect:", choices= NULL, multiple= FALSE, options = list(maxItems=1, placeholder="Select a genomic alteration"), selected= NULL),
                                                                                                   conditionalPanel(condition= "input.vartmp_options_select_bar",
                                                                                                        selectizeInput("vartmp_options_select_var2_bar","Select one additional type of genomic alteration to inspect:", choices= NULL, multiple= FALSE, options = list(maxItems=1, placeholder="Select a genomic alteration"), selected= NULL),
                                                                                                        checkboxInput("fisher_bar_samplelist_checkbox", "Sample List", value= FALSE),
                                                                                                        conditionalPanel(condition= "input.fisher_bar_samplelist_checkbox== true",fileInput("fisher_bar_samplelist", "Select a file including only Tumor Barcodes to include in the analysis:",multiple= FALSE, accept = 'txt', placeholder = "Choose File")),
                                                                                                        checkboxInput("sp_group_checkbox_bar", "SP_Group", value=FALSE), conditionalPanel(condition= "input.sp_group_checkbox_bar ==true",uiOutput("sp_group_choices_bar")),
                                                                                                        actionButton("generate_fisher_barplot", "Generate Fisher Barplot"), plotOutput("fisher_output_bar"),uiOutput("download_fisher_barplot_separate"))),
                                                                                                 
                                                                                                 tabPanel(value="association_test_main", "Association Testing",uiOutput("data_list_association"),actionButton("load_datasets", "Load Datasets"),conditionalPanel(condition="input.load_datasets",tags$header(strong("Selected Data Dimensions")), tableOutput("assoc_datatable_dim"), tags$header(strong(textOutput("assoc_data_join_by"))), textOutput("assoc_no_tum_barcode"), uiOutput("assoc_dropdown")), DT::dataTableOutput("assoc_datatable"),
                                                                                                    tabsetPanel(id="association_testing",
                                                                                                      tabPanel(value= "bivar_analysis","Bivariable Analysis", checkboxInput(inputId="group_var_bivariable", label= "Group_Var", value= FALSE),
                                                                                                        conditionalPanel(condition= "input.group_var_bivariable == false",uiOutput("assoc_variable_list_1"),
                                                                                                          checkboxInput("filter_assoc_var1","Filter (>0)",value=FALSE),
                                                                                                          checkboxInput("log2_assoc_var1", "log2", value=FALSE),
                                                                                                          textInput("collapse_assoc_var1", "Collapse Level",value=NULL),
                                                                                                          uiOutput("assoc_variable_list_2"), 
                                                                                                          checkboxInput("filter_assoc_var2","Filter (>0)",value=FALSE), 
                                                                                                          checkboxInput("log2_assoc_var2", "log2", value=FALSE),
                                                                                                          textInput("collapse_assoc_var2", "Collapse Level",value=NULL),
                                                                                                          pickerInput("assoc_types_list", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes","fisher","skit"), multiple=FALSE),
                                                                                                          checkboxInput("association_output_plot", "Save plot to computer",value = FALSE),
                                                                                                          textInput("file_ext_assoc", "File Extension (png, svg, or jpg)", value="png"),
                                                                                                          actionButton("calculate_association","Calculate"),actionButton("reset_association", "Reset"),plotOutput("bivariable_plot", width="800px")),
                                                                                                        conditionalPanel(condition= "input.group_var_bivariable == true",uiOutput("group_var_text_box_bivariable"),uiOutput("assoc_variable_list_1_group"),
                                                                                                          checkboxInput("filter_assoc_var1_group","Filter (>0)",value=FALSE),
                                                                                                          checkboxInput("log2_assoc_var1_group", "log2", value=FALSE),
                                                                                                          textInput("collapse_assoc_var1_group", "Collapse Level",value=NULL),
                                                                                                          uiOutput("assoc_variable_list_2_group"), 
                                                                                                          checkboxInput("filter_assoc_var2_group","Filter (>0)",value=FALSE), 
                                                                                                          checkboxInput("log2_assoc_var2_group", "log2", value=FALSE),
                                                                                                          textInput("collapse_assoc_var2_group", "Collapse Level",value=NULL),
                                                                                                          pickerInput("assoc_types_list_group", "Select a test type for the association test you would like to run:", choices=c("parametric", "nonparametric", "robust","bayes","fisher","skit"), multiple=FALSE),
                                                                                                          actionButton("calculate_association_group","Calculate"),actionButton("reset_association_group", "Reset"),DT::dataTableOutput("bivariable_grp_result"))),
                                                                                                      tabPanel(value= "multivar_analysis","Multivariable Analysis","Supported regression model types include: lm and glm.",
                                                                                                               checkboxInput(inputId="group_var_regression", label= "Group_Var", value= FALSE),
                                                                                                               conditionalPanel(condition="input.group_var_regression == true", uiOutput("group_var_text_box"), 
                                                                                                                                uiOutput("regression_model_group"), 
                                                                                                                                uiOutput("formula_input_group"),
                                                                                                                                actionButton("calculate_regression_group","Calculate"),actionButton("reset_regression_group", "Reset"),DT::dataTableOutput("regression_grp_result")),
                                                                                                               conditionalPanel(condition= "input.group_var_regression == false",uiOutput("formula_input"),
                                                                                                                                checkboxInput("regression_output_plot", "Save plot to computer",value = FALSE), 
                                                                                                                                textInput("file_ext_regression", "File Extension (png, svg, or jpg)", value= "png"), 
                                                                                                                                actionButton("calculate_regression","Calculate"),actionButton("reset_regression", "Reset"),plotOutput("regression_plot"))))),
                                                                                                tabPanel(value= "oncoplot", "Oncoplot", 
                                                                                                         radioButtons("oncoplot_input_selection", "Select an input type for generating the oncoplot.", choices = c("By genomic alteration" = "oncoplot_genalts" , "By genomic alteration text input" = "oncoplot_text_input_genalts", "By source category" = "oncoplot_categories","By source category feature frequency" = "oncoplot_cat_freq_inputs"), selected =character(0)),
                                                                                                         conditionalPanel(condition= "input.oncoplot_input_selection == 'oncoplot_genalts'", selectizeInput("oncoplot_genalts_select","Select at least one type of genomic alteration tfor oncoplot generation:", choices= NULL, multiple= TRUE, options = list(placeholder="Select one or more genomic alteration"), selected= NULL), checkboxInput("order_by_input1", "Order by Input", value = FALSE)),
                                                                                                         conditionalPanel(condition= "input.oncoplot_input_selection == 'oncoplot_text_input_genalts'", textAreaInput("oncoplot_user_input_genalts", "Write in a specific list of genomic alterations for oncoplot generation", placeholder="Gender|Overall_Feature \n Smoking|Overall_Feature"), checkboxInput("order_by_input2", "Order by Input", value = FALSE)),
                                                                                                         conditionalPanel(condition= "input.oncoplot_input_selection == 'oncoplot_categories'", selectizeInput("oncoplot_categories_select", "Select one or more source categories to include in the oncoplot generation:", choices= NULL, multiple= TRUE, options = list(placeholder="Select one or more source categories"), selected= NULL)),
                                                                                                         conditionalPanel(condition= "input.oncoplot_input_selection == 'oncoplot_cat_freq_inputs'", textAreaInput("oncoplot_cat_freq_select","Select source categories by frequency", placeholder="Overall_Feature,0.3")),
                                                                                                         textInput("oncoplot_min_freq","Minimum Frequency",value= 0.10),
                                                                                                         actionButton("generate_oncoplot","Generate Oncoplot"), plotOutput("data_output_plot")))),
               
               tabItem(tabName = "documentation", h2("Documentation"), tabsetPanel(id= "documentation_tabs", tabPanel(title= "Project Info",value="project_info",h4("The following are brief descriptions of each project included in the app. In parentheses is the project code as it appears in the selection dropdown list."), includeMarkdown("Project_descriptions.Rmd")), tabPanel(title="Module Info",value="module_info",h4("The following are brief descriptions of each module included in the app."),dataTableOutput("module_info_table")))),
               
               tabItem(tabName = "faq", h2("Frequently Asked Questions")))
    ),
  )
  
)

