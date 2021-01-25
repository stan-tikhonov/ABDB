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

maintable = agingsignatures[["Brain"]]

for (name in names(agingsignatures)[-1]){
  maintable = full_join(maintable, agingsignatures[[name]])
}

maintable$presence_tissues = rowSums(!is.na(maintable[, grepl("^(Brain|Liver|Muscle)_logFC$", colnames(maintable))]))
maintable$presence_species = rowSums(!is.na(maintable[, grepl("^(Mouse|Rat|Human)_logFC$", colnames(maintable))]))
maintable$presence_compound = as.numeric(!is.na(maintable[, "All_logFC"]))
maintable$presence_total = rowSums(!is.na(maintable[, grepl("^(Mouse|Rat|Human|All|Brain|Muscle|Liver)_logFC$", colnames(maintable))]))

maintable$geommean_pval_tissues = (maintable$Brain_adj_pval * maintable$Muscle_adj_pval * maintable$Liver_adj_pval) ^ (1/3)
maintable$geommean_pval_tissues_robust = (maintable$Brain_pval_robust * maintable$Muscle_pval_robust * maintable$Liver_pval_robust) ^ (1/3)
maintable$geommean_pval_tissues_LOO = (maintable$Brain_adj_pval_LOO * maintable$Muscle_adj_pval_LOO * maintable$Liver_adj_pval_LOO) ^ (1/3)

maintable$geommean_pval_species = (maintable$Human_adj_pval * maintable$Mouse_adj_pval * maintable$Rat_adj_pval) ^ (1/3)
maintable$geommean_pval_species_robust = (maintable$Human_pval_robust * maintable$Mouse_pval_robust * maintable$Rat_pval_robust) ^ (1/3)
maintable$geommean_pval_species_LOO = (maintable$Human_adj_pval_LOO * maintable$Mouse_adj_pval_LOO * maintable$Rat_adj_pval_LOO) ^ (1/3)

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
}


ui <- fluidPage(
  
  titlePanel("Aging Biomarkers Database"),
  
  
  navbarPage("ABDB",
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
                        sliderInput(inputId = "muscle_numsignifgenes",
                                    label = "Number of most significant genes to display",
                                    min = 1,
                                    max = 1000,
                                    value = 100
                        ),
                        sliderInput(inputId = "muscle_presence",
                                    label = "In how many signatures the genes must be present?",
                                    min = 1,
                                    max = 7,
                                    value = 1
                        ),
                        
                        tableOutput(outputId = "muscle_maintable")
               ),
               tabPanel("Liver",
                        sliderInput(inputId = "liver_numsignifgenes",
                                    label = "Number of most significant genes to display",
                                    min = 1,
                                    max = 1000,
                                    value = 100
                        ),
                        sliderInput(inputId = "liver_presence",
                                    label = "In how many signatures the genes must be present?",
                                    min = 1,
                                    max = 7,
                                    value = 1
                        ),
                        
                        tableOutput(outputId = "liver_maintable")
               )
    ),
    navbarMenu(
      "Signatures of species",
      tabPanel("Mouse",
               sliderInput(inputId = "mouse_numsignifgenes",
                           label = "Number of most significant genes to display",
                           min = 1,
                           max = 1000,
                           value = 100
               ),
               sliderInput(inputId = "mouse_presence",
                           label = "In how many signatures the genes must be present?",
                           min = 1,
                           max = 7,
                           value = 1
               ),
               
               tableOutput(outputId = "mouse_maintable")
      ),
      tabPanel("Human",
               sliderInput(inputId = "human_numsignifgenes",
                           label = "Number of most significant genes to display",
                           min = 1,
                           max = 1000,
                           value = 100
               ),
               sliderInput(inputId = "human_presence",
                           label = "In how many signatures the genes must be present?",
                           min = 1,
                           max = 7,
                           value = 1
               ),
               
               tableOutput(outputId = "human_maintable")
      ),
      tabPanel("Rat",
               sliderInput(inputId = "rat_numsignifgenes",
                           label = "Number of most significant genes to display",
                           min = 1,
                           max = 1000,
                           value = 100
               ),
               sliderInput(inputId = "rat_presence",
                           label = "In how many signatures the genes must be present?",
                           min = 1,
                           max = 7,
                           value = 1
               ),
               
               tableOutput(outputId = "rat_maintable")
      )
    ),
    tabPanel("Compound signature",
             sliderInput(inputId = "compound_numsignifgenes",
                         label = "Number of most significant genes to display",
                         min = 1,
                         max = 1000,
                         value = 100
             ),
             sliderInput(inputId = "compound_presence",
                         label = "In how many signatures the genes must be present?",
                         min = 1,
                         max = 7,
                         value = 1
             ),
             
             tableOutput(outputId = "compound_maintable")
    ),
    tabPanel("Data for parsing")
  )

)



server = function(input, output) {
  
  kok = reactive({temp = agingsignaturesapp[["Brain"]] %>% filter(presence_total >= input$brain_presence) %>% arrange(adj_pval)
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
    kok(), escape = FALSE, class = "cell-border stripe", options = list(lengthMenu = c(25, 50, 100), scrollX=600, scrollY=600, ordering=F)
  )
  
  
  output$muscle_maintable = renderTable({
    temp = (agingsignaturesapp[["Muscle"]] %>% filter(presence_total >= input$muscle_presence) %>% top_n(n = (-1 * input$muscle_numsignifgenes), wt = adj_pval) %>% arrange(adj_pval))[, c("entrez", "genesymbol", "logFC", "pvalue", "adjusted_pvalue", "Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC")]
    return(temp)
  }, sanitize.text.function = function(x) x)
  output$liver_maintable = renderTable({
    temp = (agingsignaturesapp[["Liver"]] %>% filter(presence_total >= input$liver_presence) %>% top_n(n = (-1 * input$liver_numsignifgenes), wt = adj_pval) %>% arrange(adj_pval))[, c("entrez", "genesymbol", "logFC", "pvalue", "adjusted_pvalue", "Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC")]
    return(temp)
  }, sanitize.text.function = function(x) x)
  output$human_maintable = renderTable({
    temp = (agingsignaturesapp[["Human"]] %>% filter(presence_total >= input$human_presence) %>% top_n(n = (-1 * input$human_numsignifgenes), wt = adj_pval) %>% arrange(adj_pval))[, c("entrez", "genesymbol", "logFC", "pvalue", "adjusted_pvalue", "Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC")]
    return(temp)
  }, sanitize.text.function = function(x) x)
  output$mouse_maintable = renderTable({
    temp = (agingsignaturesapp[["Mouse"]] %>% filter(presence_total >= input$mouse_presence) %>% top_n(n = (-1 * input$mouse_numsignifgenes), wt = adj_pval) %>% arrange(adj_pval))[, c("entrez", "genesymbol", "logFC", "pvalue", "adjusted_pvalue", "Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC")]
    return(temp)
  }, sanitize.text.function = function(x) x)
  output$rat_maintable = renderTable({
    temp = (agingsignaturesapp[["Rat"]] %>% filter(presence_total >= input$rat_presence) %>% top_n(n = (-1 * input$rat_numsignifgenes), wt = adj_pval) %>% arrange(adj_pval))[, c("entrez", "genesymbol", "logFC", "pvalue", "adjusted_pvalue", "Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC")]
    return(temp)
  }, sanitize.text.function = function(x) x)
  output$compound_maintable = renderTable({
    temp = (agingsignaturesapp[["All"]] %>% filter(presence_total >= input$compound_presence) %>% top_n(n = (-1 * input$compound_numsignifgenes), wt = adj_pval) %>% arrange(adj_pval))[, c("entrez", "genesymbol", "logFC", "pvalue", "adjusted_pvalue", "Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC")]
    return(temp)
  }, sanitize.text.function = function(x) x)
  
}

shinyApp(ui = ui, server = server)