#' SCN
#' Functional variant prediction for voltage-gated sodium channels

# packages
library(tidyverse)
library(tidymodels)
library(shiny)
library(shinyalert)
library(shinythemes)
library(shinybusy)
library(shinyWidgets)
library(ontologyIndex)
library(ontologySimilarity)
library(data.table)
library(yardstick)
library(caret)
library(bestNormalize)
library(kernlab)
library(e1071)
library(magic)
library(matrixcalc)
library(Matrix)
library(klic)
library(CEGO)
library(jaccard)

library(BiocManager)
options(repos = BiocManager::repositories())
library(qvalue)

# get helper fn
source("func.R")

# set seed
set.seed(42)

# get preprocessed training data set
train <- read_csv("training_data.csv") 
y <- as.factor(train$y)

# get similarity matrix
sim_matrices <- read_csv("similaritymatrix.csv")
rownames(sim_matrices) <- colnames(sim_matrices)

sim_match <- sim_matrices %>%
  rownames_to_column() %>%
  pivot_longer(cols = -c(1)) %>%
  rename(k = rowname, l = name)

# get pre-processing recipe
load("recipe.rds")

# get feature lookup tables
aa_feats <- read_tsv("aa_feats.tsv")
load("str_feats.rda")
cid_raw <- read_csv("cid.csv")

# set up menu choices
vec_genes <- c("SCN1A", "SCN2A", "SCN3A", "SCN4A", "SCN5A", "SCN8A", "SCN9A", "SCN10A", "SCN11A")
vec_aa <- c("A - Ala", "R - Arg", "N - Asn", "D - Asp", "C - Cys", "E - Glu", "Q - Gln", "G - Gly", "H - His", "I - Ile", "L - Leu", "K - Lys", "M - Met", "F - Phe", "P - Pro", "S - Ser", "T - Thr", "W - Trp", "Y - Tyr", "V - Val")

ont_hpo <- get_ontology("hp.obo.txt", 
                        propagate_relationships = "is_a", 
                        extract_tags = "minimal")
list_hpo <- lapply(1:length(ont_hpo$id), function(x) paste(ont_hpo$id[[x]], ont_hpo$name[[x]], sep = " "))

ont_omim <- read_csv("phenotype.csv")
list_omim <- ont_omim %>%
  mutate(V1 = paste(`#DatabaseID`, DiseaseName, sep = " ")) %>%
  select(V1) %>%
  distinct() %>%
  unlist(.$V1) %>%
  as.list(.$V1) %>%
  unname()

# app
shinyApp(
  
  fluidPage(
    
    # Shinyjs for reactive server-side lists
    shinyjs::useShinyjs(),
    
    # load custom css file and app logo
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom_div.css"),
      tags$link(rel = "shortcut icon", href = "logo.png")
    ),
    
    # bootstrap theme
    theme = shinytheme("flatly"),
    
    # for FAQ
    useShinyalert(), 
    
    # for loading symbol
    shinybusy::add_busy_spinner(spin = "fading-circle", color = "#2c3e50", position = "bottom-right"),
    
    # title
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
        
        tags$div(id = "gene",
                 selectizeInput(inputId = "gene",
                             label = "Channel:",
                             choices = vec_genes)),
        
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
        
        checkboxInput(inputId = "flag_exp", label = "Experimental settings", value = FALSE), # currently does not do anything
        
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
    updateSelectizeInput(session, "hpo", choices = c(list_hpo, list_omim), server = TRUE)
    
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
                 text = div(style = "overflow: auto; max-height: 50vh;", HTML(
                   "<b> 1. What is this? </b> </br> 
                 This is SCION, a multi-task multi-kernel learning support vector machine (MTMKL-SVM) built to classify the functional effects of non-synonymous missense variants in voltage-gated sodium channels, trained on data by Brunklaus et al. (DOI 10.1093/brain/awac006). </br></br>
                 
                 <b> 2. How do I interpret the results? </b> </br> 
                 The predicted class is displayed as plain text output. Class probabilities are computed via Platt scaling. </br></br>
                 
                 <b> 3. Do I have to enter phenotypic information? </b> </br> 
                 No. If no phenotypic information is provided by the user, a multi-task single-kernel learning SVM is used for prediction instead (cf. DOI 10.1101/2021.12.02.470894). </br></br>
                 
                 <b> For suggestions and feedback, please get in touch: </b> </br> 
                 christian.bosselmann@med.uni-tuebingen / @cmbosselmann </br></br>"
                 )), 
                 type = "info",
                 html = TRUE,
                 closeOnEsc = TRUE,
                 closeOnClickOutside = TRUE)
    })
  }
)