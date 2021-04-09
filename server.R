

server <- function(input, output, session){
  observeEvent(input$select,{
    #output$selection <- renderUI({paste0("You have selected the following project: ", input$project_code)})
    output$file_prompt <- renderUI({paste0("Below are the files available for the ", input$project_code, " project:")})
    #output$file_list <- renderUI({list.files(paste0("Genomic Data/", "Sherlock"))})
    #output$file_prompt <- renderTable({list.files(paste0("/Genomic Data/", input$project_code))})
    #output$selection_b <- renderUI({"Would you like to load all of the files available into their corresponding modules for this project?"})
    
    })

  observeEvent(input$select, {output$file_list <- renderTable({
    setwd("/Users/kleinam/sherlock_genome/Sherlock-Genome/Genomic Data")
    files <- list.files(paste0(getwd(), '/', input$project_code))
    tibble(files)
  })
    #checkboxGroupInput("file_list", label=paste0("Files Available in ", input$project_code), choices= files)})
  })
  
  
  
  # don't really need the load all files button; could just have all of the files load after clicking yes
  observeEvent(input$yes_all_files,{
  #  output$yes <- renderUI({paste0("All available files for the ", input$project_code, " project have been loaded into their corresponding module on the left.")})
    output$yes <- renderUI({actionButton("load_all_files", "Load All Files")})
  })
  
 # observeEvent(input$yes_all_files,{
 #   project_file_list <- list.files(paste0("Genomic Data/", input$project_code))s
 #   i <- 1
  #  for(each in 1:length(project_file_list)){
      #file_read_in <- readLines(con=file(paste0("Genomic Data/","Sherlock","/",project_file_list[[i]]))) %>% includeMarkdown()
  #    file_read_in <- includeMarkdown(paste0("Genomic Data/", "Sherlock/",project_file_list[[i]]))
      #print(file_read_in)
      #i <- i + 1
  #    output$yes_load <- renderUI({file_read_in})
  #  }  
  #})
  
  #observeEvent(input$no_all_files,{
  #    project_code <- input$project_code
  #    output$no <- renderUI(conditionalPanel(condition="input.no_all_files",checkboxGroupInput("file_list",label= NULL, choices=list.files(paste0("Genomic Data/", project_code)), selected= NULL),actionButton("load_selected_files", "Load selected files")))
  #})
  
  #observeEvent(input$file_list, {if("Study_Overview.Rmd" %in% input$file_list){
 #   path_overview <- file.path(paste0("Genomic Data/",input$project_code,"/Study_Overview.Rmd"))
  #} else{
    # Warning in file(con, "r") :file("") only supports open = "w+" and open = "w+b": using the former (still works)
 #   path_overview <- ""
 # }
  #  output$study_overview <- eventReactive(input$load_selected_files, {includeMarkdown(path_overview)})
  #})
  
 # observeEvent(input$file_list, {if("test.txt" %in% input$file_list){
 #   path_test <- file.path(paste0("Genomic Data/",input$project_code,"/test.txt"))
 # } else{
    # Warning in file(con, "r") :file("") only supports open = "w+" and open = "w+b": using the former (still works)
 #   path_test <- ""
 # }
   # might need to change includeMarkdown in the future when not an Rmd file
 #   output$test <- eventReactive(input$load_selected_files, {includeMarkdown(path_test)})
  #})
  
 # output$oncoplot <- renderPlot({})
  
}