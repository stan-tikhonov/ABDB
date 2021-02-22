library(shiny)
library(tidyverse)
library(org.Mm.eg.db)
library(annotate)
library(DT)

# load the signature list
load("agingsignatures_v4.RData")
agingsignatures <- agingsignatures_v3
cols = c()
# aggregate data
for (name in names(agingsignatures)){
  colnames(agingsignatures[[name]]) = paste0(paste0(name, "_"), colnames(agingsignatures[[name]]))
  agingsignatures[[name]]$entrez = rownames(agingsignatures[[name]])
}

maintable = agingsignatures[["Brain"]][,c("Brain_logFC", "Brain_pval", "Brain_adj_pval", "entrez")]

for (name in names(agingsignatures)[-1]){
  maintable = full_join(maintable, agingsignatures[[name]][, c(paste0(name, c("_logFC", "_pval", "_adj_pval")), "entrez")])
}

maintable$presence_tissues = rowSums(!is.na(maintable[, grepl("^(Brain|Liver|Muscle)_logFC$", colnames(maintable))]))
maintable$presence_species = rowSums(!is.na(maintable[, grepl("^(Mouse|Rat|Human)_logFC$", colnames(maintable))]))
maintable$presence_compound = as.numeric(!is.na(maintable[, "All_logFC"]))
maintable$presence_total = rowSums(!is.na(maintable[, grepl("^(Mouse|Rat|Human|All|Brain|Muscle|Liver)_logFC$", colnames(maintable))]))

maintable$geommean_pval_tissues = (maintable$Brain_adj_pval * maintable$Muscle_adj_pval * maintable$Liver_adj_pval) ^ (1/3)
#maintable$geommean_pval_tissues_robust = (maintable$Brain_pval_robust * maintable$Muscle_pval_robust * maintable$Liver_pval_robust) ^ (1/3)
#$geommean_pval_tissues_LOO = (maintable$Brain_adj_pval_LOO * maintable$Muscle_adj_pval_LOO * maintable$Liver_adj_pval_LOO) ^ (1/3)

maintable$geommean_pval_species = (maintable$Human_adj_pval * maintable$Mouse_adj_pval * maintable$Rat_adj_pval) ^ (1/3)
#maintable$geommean_pval_species_robust = (maintable$Human_pval_robust * maintable$Mouse_pval_robust * maintable$Rat_pval_robust) ^ (1/3)
#maintable$geommean_pval_species_LOO = (maintable$Human_adj_pval_LOO * maintable$Mouse_adj_pval_LOO * maintable$Rat_adj_pval_LOO) ^ (1/3)

tablefordownloading = maintable

for (name in c("Brain", "Muscle", "Liver", "Mouse", "Human", "Rat", "All")){
  maintable[,paste0(name, "_logFC")] = round(maintable[,paste0(name, "_logFC")], digits = 2)
}

for (i in 1:length(rownames(maintable))){
  for (name in c("Brain", "Muscle", "Liver", "Mouse", "Human", "Rat", "All")){
    
    if (is.na(maintable[i, paste0(name, "_logFC")])){
      maintable[i, paste0(name, "_logFC")] = "NA"
    } else if (maintable[i, paste0(name, "_logFC")] == 0){
      maintable[i, paste0(name, "_logFC")] = as.character(maintable[i, paste0(name, "_logFC")])
    } else if (maintable[i, paste0(name, "_logFC")] < 0){
      maintable[i, paste0(name, "_logFC")] = paste0(as.character(maintable[i, paste0(name, "_logFC")]), ' <span style = "color:red;font-size:22px">&darr;</span>')
    } else {
      maintable[i, paste0(name, "_logFC")] = paste0(as.character(maintable[i, paste0(name, "_logFC")]), ' <span style = "color:lime;font-size:22px">&uarr;</span>')
    }
    
  }
}

presencelookup = maintable[,c("entrez", "presence_total")]

agingsignaturesapp <- agingsignatures_v3

for (name in names(agingsignaturesapp)){
  agingsignaturesapp[[name]]$entrez = rownames(agingsignaturesapp[[name]])
  agingsignaturesapp[[name]] = left_join(agingsignaturesapp[[name]], presencelookup)
  
  symbols = getSYMBOL(agingsignaturesapp[[name]]$entrez, data = "org.Mm.eg.db")
  agingsignaturesapp[[name]]$genesymbol = symbols[agingsignaturesapp[[name]]$entrez]
  
  agingsignaturesapp[[name]] = left_join(agingsignaturesapp[[name]], maintable[,c("Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC", "entrez")], by = c("entrez"))
  
  agingsignaturesapp[[name]]$logFC = round(agingsignaturesapp[[name]]$logFC, digits = 2)
  agingsignaturesapp[[name]]$pvalue = formatC(agingsignaturesapp[[name]]$pval, format = "e", digits = 2)
  agingsignaturesapp[[name]]$adjusted_pvalue = formatC(agingsignaturesapp[[name]]$adj_pval, format = "e", digits = 2)
  agingsignaturesapp[[name]]$entrez = paste0('<a href="https://www.ncbi.nlm.nih.gov/gene/?term=', agingsignaturesapp[[name]]$entrez, '">', agingsignaturesapp[[name]]$entrez, '</a>')
}


ui <- fluidPage(
  navbarPage("ABDB",
    tabPanel("Home",
            tags$h1("Aging Biomarker Database"),
            tags$p("Welcome to Aging Biomarker Database (ABDB)!", style="font-size:15px;"),
            tags$p("This is a tool where one can visualize age-related gene expression 
                   changes in different mammalian species and tissues. The tool is actually a graphical 
                    interface for 7 quantitative signatures of aging, which one can download as a .csv file at the 
                    bottom of the page. Currently, things that 
                   can be visualized only include the quantitative differential expression data (LogFC) for 
                    every gene in a form that allows easy comparison 
                    of the expression changes between the signatures.", style="font-size:15px;"),
            tags$p("The manual page contains all the information needed to understand and work with the contents of this website.", style="font-size:15px;"),
            tags$p("The differential expression data comprising all 7 signatures can be downloaded with the following link:", style="font-size:15px;"),
            downloadButton("diffexpr_download", "Download")
    ),
    navbarMenu("Signatures of tissues",
               tabPanel("Brain",
                        
                        sidebarLayout(
                          sidebarPanel(
                            sliderInput(inputId = "brain_presence",
                                        label = "In how many signatures the genes must be present?",
                                        min = 1,
                                        max = 7,
                                        value = 1
                            ),
                            selectInput(inputId = "brain_sort",
                                        label = "Sort by:",
                                        choices = c("LogFC (ascending)" = "logFC_a",
                                                    "LogFC (descending)" = "logFC_d",
                                                    "Adjusted p-value" = "adj_pval"
                                                    ),
                                        selected = "adj_pval")
                            
                          ),
                        
                          mainPanel(
                            DT::dataTableOutput("brain_maintable")
                          )
                        )
               ),
              
               tabPanel("Muscle",
                        sidebarLayout(
                          sidebarPanel(
                            sliderInput(inputId = "muscle_presence",
                                        label = "In how many signatures the genes must be present?",
                                        min = 1,
                                        max = 7,
                                        value = 1
                            ),
                            selectInput(inputId = "muscle_sort",
                                        label = "Sort by:",
                                        choices = c("LogFC (ascending)" = "logFC_a",
                                                    "LogFC (descending)" = "logFC_d",
                                                    "Adjusted p-value" = "adj_pval"
                                        ),
                                        selected = "adj_pval")
                            
                          ),
                          
                          mainPanel(
                            DT::dataTableOutput("muscle_maintable")
                          )
                        )
               ),
               tabPanel("Liver",
                        sidebarLayout(
                          sidebarPanel(
                            sliderInput(inputId = "liver_presence",
                                        label = "In how many signatures the genes must be present?",
                                        min = 1,
                                        max = 7,
                                        value = 1
                            ),
                            selectInput(inputId = "liver_sort",
                                        label = "Sort by:",
                                        choices = c("LogFC (ascending)" = "logFC_a",
                                                    "LogFC (descending)" = "logFC_d",
                                                    "Adjusted p-value" = "adj_pval"
                                        ),
                                        selected = "adj_pval")
                            
                          ),
                          
                          mainPanel(
                            DT::dataTableOutput("liver_maintable")
                          )
                        )
               )
    ),
    navbarMenu(
      "Signatures of species",
      tabPanel("Mouse",
               sidebarLayout(
                 sidebarPanel(
                   sliderInput(inputId = "mouse_presence",
                               label = "In how many signatures the genes must be present?",
                               min = 1,
                               max = 7,
                               value = 1
                   ),
                   selectInput(inputId = "mouse_sort",
                               label = "Sort by:",
                               choices = c("LogFC (ascending)" = "logFC_a",
                                           "LogFC (descending)" = "logFC_d",
                                           "Adjusted p-value" = "adj_pval"
                               ),
                               selected = "adj_pval")
                   
                 ),
                 
                 mainPanel(
                   DT::dataTableOutput("mouse_maintable")
                 )
               )
      ),
      tabPanel("Human",
               sidebarLayout(
                 sidebarPanel(
                   sliderInput(inputId = "human_presence",
                               label = "In how many signatures the genes must be present?",
                               min = 1,
                               max = 7,
                               value = 1
                   ),
                   selectInput(inputId = "human_sort",
                               label = "Sort by:",
                               choices = c("LogFC (ascending)" = "logFC_a",
                                           "LogFC (descending)" = "logFC_d",
                                           "Adjusted p-value" = "adj_pval"
                               ),
                               selected = "adj_pval")
                   
                 ),
                 
                 mainPanel(
                   DT::dataTableOutput("human_maintable")
                 )
               )
      ),
      tabPanel("Rat",
               sidebarLayout(
                 sidebarPanel(
                   sliderInput(inputId = "rat_presence",
                               label = "In how many signatures the genes must be present?",
                               min = 1,
                               max = 7,
                               value = 1
                   ),
                   selectInput(inputId = "rat_sort",
                               label = "Sort by:",
                               choices = c("LogFC (ascending)" = "logFC_a",
                                           "LogFC (descending)" = "logFC_d",
                                           "Adjusted p-value" = "adj_pval"
                               ),
                               selected = "adj_pval")
                   
                 ),
                 
                 mainPanel(
                   DT::dataTableOutput("rat_maintable")
                 )
               )
      )
    ),
    tabPanel("Compound signature",
             sidebarLayout(
               sidebarPanel(
                 sliderInput(inputId = "compound_presence",
                             label = "In how many signatures the genes must be present?",
                             min = 1,
                             max = 7,
                             value = 1
                 ),
                 selectInput(inputId = "compound_sort",
                             label = "Sort by:",
                             choices = c("LogFC (ascending)" = "logFC_a",
                                         "LogFC (descending)" = "logFC_d",
                                         "Adjusted p-value" = "adj_pval"
                             ),
                             selected = "adj_pval")
                 
               ),
               
               mainPanel(
                 DT::dataTableOutput("compound_maintable")
               )
             )
    ),
    tabPanel("Manual",
      navlistPanel(
        tabPanel(
          "Signatures of aging",
          tags$h2("Signatures of aging"),
            tags$p("This website contains a graphical interface for 7 signatures of aging. Each signature 
                 represents a general pattern of age-related gene expression changes in the corresponding datasets. The datasets chosen for 
                 the signatures of species included different tissues and sexes of the same species. The tissue signatures represent the gene expression changes 
                 that are common across species but specific for the tissues. As for the compound signature, it has been built on 
                 all of the data, including all 3 species (human, mouse, rat) and all tissues (brain, liver, muscle, lung, heart, etc.), which is approximately 100 datasets 
                 (where one dataset contains information about the age-related differential expression of genes for a specific 
                 tissue of a specific species with a specific sex).", style="font-size:15px;"),
          tags$p("One can search for a particular gene within a signature in the corresponding signature tab. For each gene, there are two identifiers: Entrez ID and GeneSymbol (both correspond to mouse genes or mouse orthologs 
                 in all of the signatures). The logFC column represents how much the expression of a given gene increases (positive values) or decreases (negative values) with age. The values in this column are not just 
                 average logFCs accross the datasets of a signature, but rather the output of a mixed-effect model, which produces an average weighted by standard errors and takes into account 
                 possible batch effects (same study, same tissue, etc.). It is important to understand that each individual logFC value in a signature is of little use because the data 
                 was normalized and the timescale was distorted as a result. These values can only be compared with each other. The table also includes columns that display statistical significance of the logFCs (adjusted pvalue 
                 corresponds to BH adjustment).", style="font-size:15px;")
        ),
        tabPanel(
          "Comparative signature interface",
          tags$h2("Comparative signature interface"),
          DT::dataTableOutput("demo_table"),
            tags$p("Comparative signature interface can be found in each of the signature tabs. Its main purpose is to show 
                   how conservative the change of each gene is across signatures. The green arrows represent that the gene is upregulated with age, and red arrows represent the downregulation of the gene. 
                   The numbers next to the arrows are the corresponding logFCs rounded to 2 decimal places (so 0 can actually be a very small up- or downregulation).", style="font-size:15px;")
        ),
        tabPanel(
          "Download data",
          tags$h2("Download data"),
          tags$p("The .csv file that can be downloaded from the home page contains one big table with all of the genes 
                 and data from all of the signatures. From each signature there are three columns: logFC, pvalue and adjusted_pvalue. These columns 
                  are the actual data, whereas the rest of the columns are all derived from it and their sole purpose is to make gene filtering easier. The presence 
                 columns represent in how many signatures the gene is present (so it's a value from 0 to 3 for presence_species and 
                 presence_tissues, a value from 1 to 7 for presence_all and a value from 0 to 1 for presence_compound). If a gene is not present in a signature, it gets NAs in the corresponding 3 columns. 
                 The table also has columns containing geometric means of adjusted pvalues accross the signatures of species and the signatures of tissues.", style="font-size:15px;")
          
        )
      )
             
    )
  )

)



server = function(input, output) {
  
  output$diffexpr_download <- downloadHandler(
    filename = "aging_diffexpr_data.csv",
    content = function(file) {
      write.table(tablefordownloading, file, row.names = FALSE, sep="\t")
    }
  )
  
  demo_table = (agingsignaturesapp[["Brain"]] %>% filter(presence_total >= 5) %>% arrange(adj_pval))[1:2, c("entrez", "genesymbol", "logFC", "pvalue", "adjusted_pvalue", "Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC")]
  
  output$demo_table = DT::renderDataTable(
    demo_table, escape = FALSE, class = "cell-border stripe", options = list(lengthMenu = c(25, 50, 100), scrollX=600, ordering=F)
  )
  
  brain_table = reactive({temp = agingsignaturesapp[["Brain"]] %>% filter(presence_total >= input$brain_presence) %>% arrange(adj_pval)
                  if (input$brain_sort == "adj_pval"){
                    temp = temp %>% arrange(adj_pval)
                  } else if (input$brain_sort == "logFC_a"){
                    temp = temp %>% arrange(logFC)
                  } else if (input$brain_sort == "logFC_d"){
                    temp = temp %>% arrange(desc(logFC))
                  }
                  return(temp[, c("entrez", "genesymbol", "logFC", "pvalue", "adjusted_pvalue", "Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC")])
  })
  
  
  output$brain_maintable = DT::renderDataTable(
    brain_table(), escape = FALSE, class = "cell-border stripe", options = list(lengthMenu = c(25, 50, 100), scrollX=600, scrollY=600, ordering=F)
  )
  
  muscle_table = reactive({temp = agingsignaturesapp[["Muscle"]] %>% filter(presence_total >= input$muscle_presence) %>% arrange(adj_pval)
  if (input$muscle_sort == "adj_pval"){
    temp = temp %>% arrange(adj_pval)
  } else if (input$muscle_sort == "logFC_a"){
    temp = temp %>% arrange(logFC)
  } else if (input$muscle_sort == "logFC_d"){
    temp = temp %>% arrange(desc(logFC))
  }
  return(temp[, c("entrez", "genesymbol", "logFC", "pvalue", "adjusted_pvalue", "Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC")])
  })
  
  
  output$muscle_maintable = DT::renderDataTable(
    muscle_table(), escape = FALSE, class = "cell-border stripe", options = list(lengthMenu = c(25, 50, 100), scrollX=600, scrollY=600, ordering=F)
  )
  
  liver_table = reactive({temp = agingsignaturesapp[["Liver"]] %>% filter(presence_total >= input$liver_presence) %>% arrange(adj_pval)
  if (input$liver_sort == "adj_pval"){
    temp = temp %>% arrange(adj_pval)
  } else if (input$liver_sort == "logFC_a"){
    temp = temp %>% arrange(logFC)
  } else if (input$liver_sort == "logFC_d"){
    temp = temp %>% arrange(desc(logFC))
  }
  return(temp[, c("entrez", "genesymbol", "logFC", "pvalue", "adjusted_pvalue", "Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC")])
  })
  
  
  output$liver_maintable = DT::renderDataTable(
    liver_table(), escape = FALSE, class = "cell-border stripe", options = list(lengthMenu = c(25, 50, 100), scrollX=600, scrollY=600, ordering=F)
  )
  
  mouse_table = reactive({temp = agingsignaturesapp[["Mouse"]] %>% filter(presence_total >= input$mouse_presence) %>% arrange(adj_pval)
  if (input$mouse_sort == "adj_pval"){
    temp = temp %>% arrange(adj_pval)
  } else if (input$mouse_sort == "logFC_a"){
    temp = temp %>% arrange(logFC)
  } else if (input$mouse_sort == "logFC_d"){
    temp = temp %>% arrange(desc(logFC))
  }
  return(temp[, c("entrez", "genesymbol", "logFC", "pvalue", "adjusted_pvalue", "Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC")])
  })
  
  
  output$mouse_maintable = DT::renderDataTable(
    mouse_table(), escape = FALSE, class = "cell-border stripe", options = list(lengthMenu = c(25, 50, 100), scrollX=600, scrollY=600, ordering=F)
  )
  
  human_table = reactive({temp = agingsignaturesapp[["Human"]] %>% filter(presence_total >= input$human_presence) %>% arrange(adj_pval)
  if (input$human_sort == "adj_pval"){
    temp = temp %>% arrange(adj_pval)
  } else if (input$human_sort == "logFC_a"){
    temp = temp %>% arrange(logFC)
  } else if (input$human_sort == "logFC_d"){
    temp = temp %>% arrange(desc(logFC))
  }
  return(temp[, c("entrez", "genesymbol", "logFC", "pvalue", "adjusted_pvalue", "Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC")])
  })
  
  
  output$human_maintable = DT::renderDataTable(
    human_table(), escape = FALSE, class = "cell-border stripe", options = list(lengthMenu = c(25, 50, 100), scrollX=600, scrollY=600, ordering=F)
  )
  
  rat_table = reactive({temp = agingsignaturesapp[["Rat"]] %>% filter(presence_total >= input$rat_presence) %>% arrange(adj_pval)
  if (input$rat_sort == "adj_pval"){
    temp = temp %>% arrange(adj_pval)
  } else if (input$rat_sort == "logFC_a"){
    temp = temp %>% arrange(logFC)
  } else if (input$rat_sort == "logFC_d"){
    temp = temp %>% arrange(desc(logFC))
  }
  return(temp[, c("entrez", "genesymbol", "logFC", "pvalue", "adjusted_pvalue", "Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC")])
  })
  
  
  output$rat_maintable = DT::renderDataTable(
    rat_table(), escape = FALSE, class = "cell-border stripe", options = list(lengthMenu = c(25, 50, 100), scrollX=600, scrollY=600, ordering=F)
  )
  
  compound_table = reactive({temp = agingsignaturesapp[["All"]] %>% filter(presence_total >= input$compound_presence) %>% arrange(adj_pval)
  if (input$compound_sort == "adj_pval"){
    temp = temp %>% arrange(adj_pval)
  } else if (input$compound_sort == "logFC_a"){
    temp = temp %>% arrange(logFC)
  } else if (input$compound_sort == "logFC_d"){
    temp = temp %>% arrange(desc(logFC))
  }
  return(temp[, c("entrez", "genesymbol", "logFC", "pvalue", "adjusted_pvalue", "Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC")])
  })
  
  
  output$compound_maintable = DT::renderDataTable(
    compound_table(), escape = FALSE, class = "cell-border stripe", options = list(lengthMenu = c(25, 50, 100), scrollX=600, scrollY=600, ordering=F)
  )
  
  
}

shinyApp(ui = ui, server = server)