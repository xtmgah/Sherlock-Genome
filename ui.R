library(shiny)

ui <- fluidPage(
  titlePanel(title = "Sherlock-Genome: A R Shiny App for Genomic Analysis and Visualization",
             windowTitle = "Sherlock-Genome"),
  dashboardPage(skin = "blue",
    dashboardHeader(title = "Sherlock-Genome"),
    dashboardSidebar(sidebarMenu(
      menuItem("Study Overview", tabName = "study_overview"),
      menuItem("Data QC", tabName = "data_qc"),
      menuItem("Mutations", tabName = "mutations"),
      menuItem("SCNA", tabName = "sv"),
      menuItem("Mutational Signature", tabName = "mutational_signature"),
      menuItem("Genomic Landscape", tabName = "genomic_landscape"),
      menuItem("Clonal Evolution", tabName = "clonal_evolution"),
      menuItem("Survival Analysis", tabName = "survival_analysis"),
      menuItem("Integrative Analysis", tabName = "integrative_analysis")
    )),
    dashboardBody(
      tabItems(tabItem(tabName = "study_overview", h2("Study Overview")),
               tabItem(tabName = "data_qc", h2("Data Quality Control (QC)")),
               tabItem(tabName = "mutations", h2("Mutations")),
               tabItem(tabName = "sv", h2("SV")),
               tabItem(tabName = "mutational_signature", h2("Mutational Signatures")),
               tabItem(tabName = "genomic_landscape", h2("Genomic Landscape")),
               tabItem(tabName = "clonal_evolution", h2("Clonal Evolution")),
               tabItem(tabName = "survival_analysis", h2("Survival Analysis")),
               tabItem(tabName = "integrative_analysis", h2("Integrative Analysis")))
    ),
  )
  
)
