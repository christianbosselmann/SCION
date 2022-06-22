#' SCN
#' Functional variant prediction for voltage-gated sodium channels

# packages
library("librarian")
librarian::shelf(tidyverse,
                 shiny,
                 shinyalert,
                 shinythemes,
                 shinybusy,
                 openxlsx,
                 ontologyIndex)

# set up menu choices
vec_genes <- c("SCN1A", "SCN2A", "SCN3A", "SCN4A", "SCN5A", "SCN8A", "SCN9A", "SCN10A", "SCN11A")
vec_aa <- c("A - Ala", "R - Arg", "N - Asn", "D - Asp", "C - Cys", "E - Glu", "Q - Gln", "G - Gly", "H - His", "I - Ile", "L - Leu", "K - Lys", "M - Met", "F - Phe", "P - Pro", "S - Ser", "T - Thr", "W - Trp", "Y - Tyr", "V - Val")

ont_hpo <- get_ontology("hp.obo.txt", 
                        propagate_relationships = "is_a", 
                        extract_tags = "minimal")
list_hpo <- lapply(1:length(ont_hpo$id), function(x) paste(ont_hpo$id[[x]], ont_hpo$name[[x]], sep = " "))

# app
shinyApp(
  
  fluidPage(
    
    shinyjs::useShinyjs(),
    tags$link(rel = "stylesheet", type = "text/css", href = "custom-div.css"),
    
    theme = shinytheme("flatly"),
    
    useShinyalert(), 
    
    shinybusy::add_busy_spinner(spin = "fading-circle", color = "#2c3e50", position = "bottom-right"),
    
    tags$head(tags$link(rel="shortcut icon", href="logo.png")),
    
    titlePanel(
      windowTitle = "SCION",
      title = fluidRow(
        column(2, align="center", div(style = "height:0px;"), img(height = 116, width = 100, src = "logo.png")),
        column(9, align="center", div(style = "height:20px;"), "SCION", br(), h4(HTML("<u>S</u>odium <u>C</u>hannel functional variant predict<u>ion</u>"))),
        column(1, align="center", div(style = "height:20px;"), actionButton("help", "", icon = icon("question")))
        
      )
    ),
    
    hr(),
    fluidRow(),
    
    sidebarLayout(
      
      sidebarPanel(
        
        # inputs
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
        
        checkboxInput(inputId = "flag_exp", label = "Experimental settings", value = FALSE),
        
        column(12, # to center buttons
        actionButton(inputId = "click", label = "Predict", icon("paper-plane", lib = "font-awesome")),
        
        div(style="margin-bottom:10px"),
        
        actionButton(inputId = "reset", label = "Reset", icon("trash", lib = "font-awesome"),
                     style="color: #fff; background-color: #f39c12; border-color: #f39c12"),
       
         align = "center",
        style = "margin-bottom: 10px;",
        style = "margin-top: -5px;"
        ),
        
        # disclaimer
        helpText("For research use only. Not for use in diagnostic procedures.")
        
      ),
      
      mainPanel(
        
        # outputs
        shinyjs::hidden(
          div(id = "results",
        h5(textOutput("flag_mkl")),
        
        hr(),
        
        h3(textOutput("prediction")),
        h5(textOutput("GOF")),
        h5(textOutput("LOF"))
          )
        )
      )
    )
  ),
  
  function(input, output, session){
    
    # server-side selectize
    updateSelectizeInput(session, "hpo", choices = list_hpo, server = TRUE)
    
    # reset button, refers to div results on UI side
    observeEvent(input$reset, {
      shinyjs::reset("gene")
      shinyjs::reset("aa1")
      shinyjs::reset("pos")
      shinyjs::reset("aa2")
      shinyjs::reset("hpo")
      shinyjs::hide("results")
    })
    
    observeEvent(input$click, {
      shinyjs::show("results")
    })
    
    # check if user wants to start prediction
    observeEvent(input$click, {
      
      # prepare objects to pass to model script
      df_in <- data.frame(gene = input$gene, 
                          aa1 = input$aa1, 
                          pos = input$pos, 
                          aa2 = input$aa2,
                          stringsAsFactors = FALSE)
      
      flag_exp <- input$flag_exp
      
      df_hpo <- input$hpo
      
      # run model
      source("master.R", local = TRUE)
      
      # prepare output objects
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