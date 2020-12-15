

server <- function(input, output, session){
  #output$selection <- eventReactive(input$select, {print(paste0("You have selected the following project: ",input$project_code))})
  #output$selection_b <- eventReactive(input$select, {print("Would you like to populate all of the modules with the corresponding files available for this project all at once?")})
  observeEvent(input$select,{
    output$selection <- renderText({paste0("You have selected the following project: ", input$project_code)})
    output$selection_b <- renderText({"Would you like to load all of the files available into their corresponding modules for this project?"})
    output$yes <- renderUI({actionButton("yes_all_files", "Yes")})
    output$no <- renderUI({actionButton("no_all_files", "No")})
  })
  
  #output$selection_yesno <- eventReactive(input$select, renderUI({actionButton("yes_all_files", "Yes")}))
  # will need to be changed depending on if the user decides to upload one-by-one or all at once for the project code selected
  output$selection_overview <- eventReactive(input$select, {includeMarkdown(paste0("/Users/kleinam/sherlock_genome/Sherlock-Genome/Genomic Data/Study_Overview_",input$project_code,".Rmd"))})
  
  # what to do if user says to upload all files at once, or to see the filenames and descriptions
  #observeEvent(input$yes_all_files, {print(paste0("Selected:",input$yes_all_files))})
  #observeEvent(input$no_all_files, {print(input$no_all_files)})
}