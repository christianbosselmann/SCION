# PREFEKT
# PREdicting the Functional Effects of Kv muTations

# packages
library("librarian")
librarian::shelf(tidyverse,
                 shiny,
                 shinyalert,
                 openxlsx)

train <- read_csv("training_data.csv")
aa_braun <- read_csv("aa_braun.csv")

# app
shinyApp(
  
  fluidPage(
    
    useShinyalert(),  # Set up shinyalert
    
    titlePanel(
      
      tagList(
        span(h2("SCN"),
             span(actionButton("help", "", icon = icon("question")), 
                  style = "position:absolute;right:2em;top:1em"),
             h4("predicting the functional effects of Nav mutations"),
        )
      )
    ),
    hr(),
    fluidRow(),
    
    sidebarLayout(
      
      sidebarPanel(
        
        selectInput(inputId = "gene",
                    label = "gene:",
                    choices = sort(unique(train$gene))),
        
        selectInput(inputId = "aa1",
                    label = "original aa:",
                    choices = sort(unique(aa_braun$aa1))),
        
        numericInput(inputId = "pos",
                     label = "position:",
                     value = ""),
        
        selectInput(inputId = "aa2",
                    label = "mutated aa:",
                    choices = sort(unique(aa_braun$aa1))),
        
        actionButton(inputId = "click", label = "Predict"),
        
        helpText("For research use only. Not for use in diagnostic procedures.")
        
      ),
      
      mainPanel(
        
        h3(textOutput("pred1")),
        h4(textOutput("pred2")),
        h4(textOutput("pred3")),
        h4(textOutput("pred4"))

      )
    )
  ),
  
  function(input, output, session){
    
    observeEvent(input$click, {
      
      df_in <- data.frame(gene = input$gene, 
                          aa1 = input$aa1, 
                          pos = input$pos, 
                          aa2 = input$aa2,
                          stringsAsFactors = FALSE)
      
      source("master.R", local = TRUE)

      output$pred1 <- renderText({
        paste(verb_out1, " ", sep="\n")
      })
      
      output$pred2 <- renderText({
        paste(verb_out2, " ", sep="\n")
      })
      
      output$pred3 <- renderText({
        paste(verb_out3, " ", sep="\n")
      })
      
      output$pred4 <- renderText({
        paste(verb_out4, " ", sep="\n")
      })

      output$newline <- renderUI({
        HTML(paste(" ", " ", sep="<br/>"))
      })
      
    })
    
    observeEvent(input$help, {
      shinyalert(title = "FAQ", 
                 text = "
                 <b> 1. What is this? </b> </br> 
                 This is an early version of a MTL-SVM to classify the functional effects of non-synonymous missense mutations in voltage-gated potassium channels, based on a dataset of 450 experiments.</br></br>
                 
                 <b> 2. How do I interpret the results? </b> </br> 
                 This app displays the predicted class in plain text output. For advanced users, a class probability distribution and decision values as confidence measure (point distance to hyperplane) are provided. </br></br>
                 
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