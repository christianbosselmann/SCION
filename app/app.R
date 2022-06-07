#' SCN
#' Functional variant prediction for voltage-gated sodium channels

# packages
library("librarian")
librarian::shelf(tidyverse,
                 shiny,
                 shinyalert,
                 shinythemes,
                 openxlsx,
                 ontologyIndex)

# set up menu choices
vec_genes <- c("SCN1A", "SCN2A", "SCN3A", "SCN4A", "SCN5A", "SCN8A", "SCN9A", "SCN10A", "SCN11A")
vec_aa <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

ont_hpo <- get_ontology("hp.obo.txt", 
                        propagate_relationships = "is_a", 
                        extract_tags = "minimal")
list_hpo <- lapply(1:length(ont_hpo$id), function(x) paste(ont_hpo$id[[x]], ont_hpo$name[[x]], sep = " "))


# app
shinyApp(
  
  fluidPage(theme = shinytheme("flatly"),
    
    useShinyalert(), 
    
    titlePanel(
      
      tagList(
        span(h2("SCION"),
             span(actionButton("help", "", icon = icon("question")), 
                  style = "position:absolute;right:2em;top:1em"),
             h4("Predicting the functional effects of Nav variants"),
        )
      )
    ),
    
    hr(),
    fluidRow(),
    
    sidebarLayout(
      
      sidebarPanel(
        
        selectInput(inputId = "gene",
                    label = "Gene:",
                    choices = vec_genes),
        
        selectInput(inputId = "aa1",
                    label = "Reference amino acid:",
                    choices = vec_aa),
        
        numericInput(inputId = "pos",
                     label = "Protein sequence position:",
                     value = ""),
        
        selectInput(inputId = "aa2",
                    label = "Variant amino acid:",
                    choices = vec_aa),
                    
        selectizeInput(inputId = "hpo",
                    label = "Phenotypic features:",
                    choices = NULL,
                    multiple = TRUE,
                    options = list(placeholder = "Search by HPO ID or name")),
        
        actionButton(inputId = "click", label = "Predict"),
        
        helpText("For research use only. Not for use in diagnostic procedures.")
        
      ),
      
      mainPanel(
        
        h5(textOutput("flag_mkl")),
        h3(textOutput("prediction")),
        h5(textOutput("GOF")),
        h5(textOutput("LOF"))
        # dataTableOutput("table")
        # plotOutput("plot"),
        # plotOutput("plot2")
        
      )
    )
  ),
  
  function(input, output, session){
    
    # server-side selectize
    updateSelectizeInput(session, "hpo", choices = list_hpo, server = TRUE)
    
    # check if user wants to start prediction
    observeEvent(input$click, {
      
      df_in <- data.frame(gene = input$gene, 
                          aa1 = input$aa1, 
                          pos = input$pos, 
                          aa2 = input$aa2,
                          stringsAsFactors = FALSE)
      
      df_hpo <- list(input$hpo)
      
      source("master.R", local = TRUE)
      # source("viz.R", local = TRUE)
      
      output$flag_mkl <- renderText({
        if(flag_mkl == TRUE){
          "Phenotypic information provided. Predicting with multi-task multi-kernel learning."
        }else{
          "No phenotypic information provided. Predicting with multi-task learning."
        }
      })

      output$prediction <- renderText({
        paste(verb_out, " ", sep="\n")
      })
      
      output$GOF <- renderText({
        paste("Probability of gain-of-function:", round(out$GOF, 3), sep=" ")
      })
      
      output$LOF <- renderText({
        paste("Probability of loss-of-function:", round(out$LOF, 3), sep=" ")
      })
      
      # output$table <- renderDataTable(out)
      
      # output$plot <- renderPlot(p)
      
      # output$plot2 <- renderPlot(p2)
      
      # output$newline <- renderUI({
      #   HTML(paste(" ", " ", sep="<br/>"))
      # })
      
    })
    
    # check if user looks up the FAQ
    observeEvent(input$help, {
      shinyalert(title = "FAQ", 
                 text = "
                 <b> 1. What is this? </b> </br> 
                 This is prefeKt, a multi-task learning support vector machine model built to classify the functional effects of non-synonymous missense mutations in voltage-gated potassium channels, based on a dataset of >950 patch and voltage clamp experiments. For further details, please refer to the manuscript. </br></br>
                 
                 <b> 2. How do I interpret the results? </b> </br> 
                 This app displays the predicted class in plain text output. For advanced users, a class probability distribution plot (Platt scaling) and decision values as confidence measure (point distance to hyperplane) are provided. </br></br>
                 
                 <b> 3. Class prediction and probabilities don't match. Why? </b> </br> 
                 This is a known thereotical issue with the multi-class extension of Platt scaling and its implementation in LIBSVM. Please refer to the package documentation for further details. </br></br>
                 
                 <b> For suggestions and feedback, please get in touch: </b> </br> 
                 christian.bosselmann@med.uni-tuebingen / @cmbosselmann </br></br>
                 ", 
                 type = "info",
                 html = TRUE,
                 closeOnClickOutside = TRUE)
    })
    
  }
)