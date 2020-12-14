
#Sherlock-Genome: A R Shiny App for Genomic Analysis and Visualization

#updated 12/14/2020
# specify the packages needed for the app through a character vector
# keep adding to this vector when additional packages are required to run the app successfully
packages_req <- c("shiny","shinydashboard")

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
  titlePanel(title = "Sherlock-Genome: A R Shiny App for Genomic Analysis and Visualization",
             windowTitle = "Sherlock-Genome"),

  # dashboard setup with left sidebar
  dashboardPage(skin = "blue",
    dashboardHeader(title = "Sherlock-Genome"),
    
# each of the following will be a tab on the left sidebar
    dashboardSidebar(sidebarMenu(
      menuItem("Study Overview", tabName = "study_overview"),
      menuItem("Sample QC", tabName = "sample_qc"),
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
      tabItems(tabItem(tabName = "study_overview", h2("Study Overview")),
               tabItem(tabName = "sample_qc", h2("Sample QC")),
               tabItem(tabName = "mutations", h2("Mutations")),
               tabItem(tabName = "scna", h2("SCNA")),
               tabItem(tabName = "sv", h2("SV")),
               tabItem(tabName = "mutational_signature", h2("Mutational Signatures")),
               tabItem(tabName = "genomic_landscape", h2("Genomic Landscape")),
               tabItem(tabName = "clonal_evolution", h2("Clonal Evolution")),
               tabItem(tabName = "survival_analysis", h2("Survival Analysis")),
               tabItem(tabName = "integrative_analysis", h2("Integrative Analysis"),tabsetPanel(tabPanel("Oncoplot"))))
    ),
  )
  
)
