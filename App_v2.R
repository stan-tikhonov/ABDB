library(shiny)
library(tidyverse)
library(org.Mm.eg.db)
library(annotate)
library(DT)
library(magrittr)
library(shinyWidgets)
library(plotly)
library(shinycssloaders)
library(htmlwidgets)

#library(shinyjs)
#library(logging)

# prepare stuff for metafor plots

load("logFCmatrixregr.RData")
load("SEmatrixregr.RData")
load("sourcedata.RData")
load("deminglist2.RData")

chosencols = list()
chosencols[["Human"]] = colnames(logFCmatrixregr)[grepl(".*Human.*", colnames(logFCmatrixregr))]
chosencols[["Rat"]] = colnames(logFCmatrixregr)[grepl(".*Rat.*", colnames(logFCmatrixregr))]
chosencols[["Mouse"]] = colnames(logFCmatrixregr)[grepl(".*Mouse.*", colnames(logFCmatrixregr))]
chosencols[["Brain"]] = colnames(logFCmatrixregr)[grepl(".*Brain.*", colnames(logFCmatrixregr)) | grepl(".*Frontalcortex.*", colnames(logFCmatrixregr)) | grepl(".*Cerebellum.*", colnames(logFCmatrixregr))]
chosencols[["Muscle"]] = colnames(logFCmatrixregr)[grepl(".*Muscle.*", colnames(logFCmatrixregr))]
chosencols[["Liver"]] = colnames(logFCmatrixregr)[grepl(".*Liver.*", colnames(logFCmatrixregr))]
chosencols[["All"]] = colnames(logFCmatrixregr)

strain_map = read.csv("all_aging_data_mod.csv", sep = ";", header = TRUE)
strain_map$GEOID = as.character(strain_map$GEOID)
strain_map$Strain = as.character(strain_map$Strain)
strain_map[which(strain_map$Strain == ""),"Strain"] = "NA"
colnames(strain_map) = c("source", "strain")
strain_map = strain_map %>% distinct(source, .keep_all = TRUE)

# define functions for metafor plots

hline <- function(y = 0, color = "red") {
  list(
    type = "line", 
    x0 = 0, 
    x1 = 1, 
    xref = "paper",
    y0 = y, 
    y1 = y, 
    line = list(color = color),
    name="Mixed effect model mean"
  )
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#prepare stuff for dataset plots
load("phenolist2.RData")
load("exprlist3.RData")
load("logFCtablelist3.RData")
load("ageunittable.RData")

names(logFCtablelist) = sub(" ", "", names(logFCtablelist))
names(sexyexprlist) = sub(" ", "", names(sexyexprlist))
names(sexyphenolist) = sub(" ", "", names(sexyphenolist))

# load the signature list
load("agingsignatures_v4_ABDB.RData")
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

symbols = getSYMBOL(maintable$entrez, data = "org.Mm.eg.db")
maintable$genesymbol = symbols[maintable$entrez]

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

# for (i in 1:length(rownames(maintable))){
#   for (name in c("Brain", "Muscle", "Liver", "Mouse", "Human", "Rat", "All")){
#     if (is.na(maintable[i, paste0(name, "_logFC")])){
#        maintable[i, paste0(name, "_logFC")] = "NA"
#     }
#   }
# }


# for (i in 1:length(rownames(maintable))){
#   for (name in c("Brain", "Muscle", "Liver", "Mouse", "Human", "Rat", "All")){
#     
#     if (is.na(maintable[i, paste0(name, "_logFC")])){
#       maintable[i, paste0(name, "_logFC")] = "NA"
#     } else if (maintable[i, paste0(name, "_logFC")] == 0){
#       maintable[i, paste0(name, "_logFC")] = as.character(maintable[i, paste0(name, "_logFC")])
#     } else if (maintable[i, paste0(name, "_logFC")] < 0){
#       maintable[i, paste0(name, "_logFC")] = paste0(as.character(maintable[i, paste0(name, "_logFC")]), ' <span style = "color:red;font-size:22px">&darr;</span>')
#     } else {
#       maintable[i, paste0(name, "_logFC")] = paste0(as.character(maintable[i, paste0(name, "_logFC")]), ' <span style = "color:lime;font-size:22px">&uarr;</span>')
#     }
#     
#   }
# }

presencelookup = maintable[,c("entrez", "presence_total")]

agingsignaturesapp <- agingsignatures_v3

signfactor = sign(maintable[,c("Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC")])
signfactor[signfactor == 0] = 1
signfactor[is.na(maintable[,c("Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC")])] = NA
maintable[,c("Brain_adj_pval_unsign", "Muscle_adj_pval_unsign", "Liver_adj_pval_unsign", "Mouse_adj_pval_unsign", "Human_adj_pval_unsign", "Rat_adj_pval_unsign", "All_adj_pval_unsign")] = maintable[,c("Brain_adj_pval", "Muscle_adj_pval", "Liver_adj_pval", "Mouse_adj_pval", "Human_adj_pval", "Rat_adj_pval", "All_adj_pval")]
maintable[,c("Brain_adj_pval", "Muscle_adj_pval", "Liver_adj_pval", "Mouse_adj_pval", "Human_adj_pval", "Rat_adj_pval", "All_adj_pval")] = maintable[,c("Brain_adj_pval", "Muscle_adj_pval", "Liver_adj_pval", "Mouse_adj_pval", "Human_adj_pval", "Rat_adj_pval", "All_adj_pval")] * signfactor

for (name in names(agingsignaturesapp)){
  agingsignaturesapp[[name]]$entrez = rownames(agingsignaturesapp[[name]])
  agingsignaturesapp[[name]] = left_join(agingsignaturesapp[[name]], presencelookup)
  
  symbols = getSYMBOL(agingsignaturesapp[[name]]$entrez, data = "org.Mm.eg.db")
  agingsignaturesapp[[name]]$genesymbol = symbols[agingsignaturesapp[[name]]$entrez]
  
  agingsignaturesapp[[name]] = left_join(agingsignaturesapp[[name]], maintable[,c("Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC", "entrez")], by = c("entrez"))
  agingsignaturesapp[[name]] = left_join(agingsignaturesapp[[name]], maintable[,c("Brain_adj_pval", "Muscle_adj_pval", "Liver_adj_pval", "Mouse_adj_pval", "Human_adj_pval", "Rat_adj_pval", "All_adj_pval", "entrez")], by = c("entrez"))
  agingsignaturesapp[[name]] = left_join(agingsignaturesapp[[name]], maintable[,c("Brain_adj_pval_unsign", "Muscle_adj_pval_unsign", "Liver_adj_pval_unsign", "Mouse_adj_pval_unsign", "Human_adj_pval_unsign", "Rat_adj_pval_unsign", "All_adj_pval_unsign", "entrez")], by = c("entrez"))
  
  agingsignaturesapp[[name]]$logFC = round(agingsignaturesapp[[name]]$logFC, digits = 2)
  agingsignaturesapp[[name]]$pvalue = agingsignaturesapp[[name]]$pval
  agingsignaturesapp[[name]]$adjusted_pvalue = agingsignaturesapp[[name]]$adj_pval
  agingsignaturesapp[[name]]$entrez = paste0('<a href="https://www.ncbi.nlm.nih.gov/gene/?term=', agingsignaturesapp[[name]]$entrez, '">', agingsignaturesapp[[name]]$entrez, '</a>')
}

agingsignaturesapp$Mouse = agingsignaturesapp$Mouse[-which(is.na(agingsignaturesapp$Mouse$genesymbol)), ]

metafor_genes = list()
for (name in names(agingsignaturesapp)){
  metafor_genes[[name]] = gsub("(<[^<>]+>)", "", agingsignaturesapp[[name]]$entrez)
  names(metafor_genes[[name]]) = agingsignaturesapp[[name]]$genesymbol
  
}


ui <- fluidPage(
  #useShinyjs(),
  tags$script("
    Shiny.addCustomMessageHandler('metafor_gene_send', function(value) {
    Shiny.setInputValue('metafor_gene', value);
    });
  "),
  navbarPage("ABDB",
    tabPanel("Home",
            tags$h1("Aging Biomarker Database"),
            tags$p("Welcome to Aging Biomarker Database (ABDB)!", style="font-size:15px;"),
            tags$p("This is a tool where one can visualize age-related gene expression 
                   changes in different mammalian species and tissues. The tool is actually a graphical 
                    interface for 7 quantitative signatures of aging, which one can download as a .csv file in the 
                    Downloads section. Currently, things that 
                   can be visualized only include the quantitative differential expression data (LogFC) for 
                    every gene in a form that allows easy comparison 
                    of the expression changes between the signatures.", style="font-size:15px;"),
            tags$p("The manual page contains all the information needed to understand and work with the contents of this website.", style="font-size:15px;")
    ),
    tabPanel("Signature comparison",
             
                        sidebarLayout(
                          sidebarPanel(
                            selectInput(inputId = "main_signature",
                                        label = "Main signature to be displayed:",
                                        choices = c("Brain" = "Brain",
                                                    "Muscle" = "Muscle",
                                                    "Liver" = "Liver",
                                                    "Human" = "Human",
                                                    "Mouse" = "Mouse",
                                                    "Rat" = "Rat",
                                                    "Compound" = "All"
                                        ),
                                        selected = "Brain"),
                            numericInput(inputId = "pval_thres",
                                         label = "Threshold for adjusted pvalue:",
                                         value = 1,
                                         max = 1,
                                         min = 0,
                                         step = 0.01),
                            sliderInput(inputId = "presence",
                                        label = "In how many signatures the genes must be present?",
                                        min = 1,
                                        max = 7,
                                        value = 1
                            ),
                            numericInput(inputId = "pval_soft_thres",
                                         label = "P-value threshold for highlighting moderately significant genes:",
                                         value = 0.05,
                                         max = 1,
                                         min = 0,
                                         step = 0.01),
                            numericInput(inputId = "pval_hard_thres",
                                         label = "P-value threshold for highlighting very significant genes:",
                                         value = 0.001,
                                         max = 1,
                                         min = 0,
                                         step = 0.01),
                            sliderInput(inputId = "num_signif",
                                        label = "In how many signatures the genes must be significant (soft threshold)?",
                                        min = 0,
                                        max = 7,
                                        value = 0
                            ),
                            pickerInput(inputId = "upregulated_in",
                                        label = "Genes shown must be upregulated in:",
                                        choices = c("Brain" = "Brain",
                                                    "Muscle" = "Muscle",
                                                    "Liver" = "Liver",
                                                    "Human" = "Human",
                                                    "Mouse" = "Mouse",
                                                    "Rat" = "Rat",
                                                    "Compound" = "All"
                                        ),
                                        options = list(
                                          "actions-box" = TRUE
                                        ),
                                        multiple=TRUE),
                            pickerInput(inputId = "downregulated_in",
                                        label = "Genes shown must be downregulated in:",
                                        choices = c("Brain" = "Brain",
                                                    "Muscle" = "Muscle",
                                                    "Liver" = "Liver",
                                                    "Human" = "Human",
                                                    "Mouse" = "Mouse",
                                                    "Rat" = "Rat",
                                                    "Compound" = "All"
                                        ),
                                        options = list(
                                          "actions-box" = TRUE
                                        ),
                                        multiple=TRUE)
                                        
                            
                          ),
                        
                          mainPanel(
                            DT::dataTableOutput("main_table")
                            #textOutput("kek")
                          )
                        )
               ),
    
    tabPanel("Gene plot",
             
             sidebarLayout(
               sidebarPanel(
                 selectInput(inputId = "metafor_signature",
                             label = "Signature to build a plot:",
                             choices = c("Brain" = "Brain",
                                         "Muscle" = "Muscle",
                                         "Liver" = "Liver",
                                         "Human" = "Human",
                                         "Mouse" = "Mouse",
                                         "Rat" = "Rat",
                                         "Compound" = "All"
                             ),
                             selected = "Brain"),
                 # multiInput(
                 #   inputId = "metafor_gene",
                 #   label = "Choose a gene:", 
                 #   choices = NULL,
                 #   options = list('limit' = 1),
                 #   choiceNames = as.character(maintable$genesymbol),
                 #   choiceValues = as.character(maintable$entrez)
                 # ),
                 
                 selectizeInput(
                   inputId = "metafor_gene",
                   label = "Choose a gene:", 
                   choices = NULL,
                   selected = NULL,
                   options = list(placeholder = "Type in gene symbol",
                                  onInitialize = I('function() { this.setValue(""); }'))
                 ),
                 
                 selectInput(inputId = "metafor_color",
                             label = "Color by:",
                             choices = c("Species" = "species",
                                         "Tissue" = "tissue",
                                         "Study" = "source",
                                         "Sex" = "sex"
                             ),
                             selected = "tissue"),
                 
                 selectInput(inputId = "metafor_shape",
                             label = "Shape by:",
                             choices = c("Species" = "species",
                                         "Sex" = "sex"
                             ),
                             selected = "species")
                 
               ),
               mainPanel(
                 withSpinner(
                  plotlyOutput(outputId = "metafor_plot")
                 ),
                 plotlyOutput(outputId = "dataset_plot", height = "100%")
               )
             )
             
          ),
              
    tabPanel("Downloads",
             tags$h2("Data for downloading"),
             tags$p("To download the table with the quantitative differential expression data, click the button below:", style="font-size:15px;"),
             downloadButton("diffexpr_download", "Download"),
             tags$p("The table is tab-separated and contains the differential expression data of all 7 signatures. From each signature there are 
                    three columns (logFC, p-value and adjusted p-value); the other columns are only present to make gene filtering easier. The complete 
                    breakdown of what each column means can be found in the Manual section.", style="font-size:15px;")),
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
            tags$p("Comparative signature interface can be found in the \"Signature comparison\" tab. Its main purpose is to show 
                   what genes are differentially expressed with age and how conservative the change of each gene is across signatures. The first columns that the user sees 
                   represent the differential expression of genes in the chosen signature (first selection menu at the left of the page). A gene's Entrez ID (clickable, is a link to NCBI for a specific gene), gene symbol (both Entrez ID and gene symbol are for mouse orthologs in all signatures), quantitative differential 
                   expression metric (logFC), it's significance before and after BH adjustment (pvalue and adjusted_pvalue) is what can be found for each gene in a selected signature.", style="font-size:15px;"),
          tags$p("One can filter out all the genes having less significance in the chosen signature than a threshold with the \"threshold for adjusted pvalue\" input. By default, all genes from the signature are shown. 
                 The presence filter is handy for filtering out the genes which are not present in a lot of signatures.", style="font-size:15px;"),
          tags$p("Now for the comparative part of the table. The green arrows represent that the gene is upregulated with age, and red arrows represent the downregulation of the gene. 
                   The numbers next to the arrows are the corresponding logFCs rounded to 2 decimal places (so 0 can actually be a very small up- or downregulation). Here, some cells are highlighted with green or red background. These cells correspond to a gene having a statistically significant up- or downregulation, respectively, in the corresponding signature. 
                 The user can specify two thresholds for adjusted p-value: the soft one for moderately significant genes (default is 0.05) and the hard one for very significant genes (default is 0.001). The genes that are more significant than the soft threshold 
                 will have lightly colored background, and those which are more significant than the hard threshold will have a brighter colored background. The genes can be filtered so that they satisfy the soft threshold in at least N signatures with a slider. One can also filter genes upregulated or downregulated in a set of signatures (selection menus at the bottom-left). 
                 These selection menus do not take into account significance, they only filter by logFC values. If several signatures are selected in them, only the intersection is displayed (genes satisfying all of the conditions).", style="font-size:15px;"),
          tags$p("Each cell in the right (\"comparative\") part of the table is clickable and will direct the user to the mixed effects model plot of the corresponding gene (table row) in the corresponding signature (table column). 
                 One can find more information on such plots in the \"Gene plots\" tab of the manual.", style="font-size:15px;")
        ),
        tabPanel(
          "Gene plots",
          tags$h2("Mixed effects model plot"),
          withSpinner(
            plotlyOutput(outputId = "metafor_demo")
          ),
          tags$p("The \"Gene plot\" tab contains plots of mixed effects models for genes within signatures. The horizontal red line on the plot is the mixed effects model mean logFC for the chosen gene in the chosen signature 
                 (the number displayed in the signature comparison table). Each point of the scatter plot represents mean logFC from a specific dataset from the signature in which the gene was present, with error bars showing the standard error for this dataset. The gene and signature 
                 for the plot can be chosen manually in the menus at the left of the table, or they will be set automatically when the user clicks on a cell in the right part of the 
                 signature comparison table. On hover, dataset information and exact logFC values are displayed for each point on the plot.", style="font-size:15px;"),
          tags$p("The \"Color by:\" and \"Shape by:\" options allow the user to tailor the plot to their needs. The default setting is to color by tissue and shape by species.", style="font-size:15px;"),
          tags$h2("Dataset plot"),
          withSpinner(
            plotlyOutput(outputId = "datasetplot_demo", height = "100%")
          ),
          tags$p("Clicking on any point on the mixed effect model plot will trigger the appearance of a dataset plot (at the bottom of the page). The dataset plot shows how the expression of the chosen gene changes across samples with different age in 
                 the dataset the clicked point corresponds to. By looking at this plot one can see where the logFC value and its standard error come from. The slope of the regression line corresponds to the logFC value, whereas the standard error (signified by the confidence interval 
                 around the line) arises from points (i.e. samples, these may be individual mice, rats or humans) being far away from the regression line. Please make note that when there are a lot of points, the confidence interval is always relatively small. On some plots 
                 the exact age of samples is not available (e.g. they are denoted as \"young\" and \"old\") and therefore the regression line cannot be fitted on the plot.", style="font-size:15px;")
          
                  ),
        tabPanel(
          "Download data",
          tags$h2("Download data"),
          tags$p("The .csv file in the Downloads section contains one big table with all of the genes 
                 and differential expression data from all of the signatures (it's tab-separated). From each signature there are three columns: logFC, pvalue and adjusted_pvalue. These columns 
                  are the actual data, whereas the rest of the columns are all derived from it and their sole purpose is to make gene filtering easier. The presence 
                 columns represent in how many signatures the gene is present (so it's a value from 0 to 3 for presence_species and 
                 presence_tissues, a value from 1 to 7 for presence_all and a value from 0 to 1 for presence_compound). If a gene is not present in a signature, it gets NAs in the corresponding 3 columns. 
                 The table also has columns containing geometric means of adjusted pvalues accross the signatures of species and the signatures of tissues.", style="font-size:15px;")
          
        )
      )
             
    )
  )

)


# basicConfig()
# 
# options(shiny.error = function() { 
#   logging::logerror(sys.calls() %>% as.character %>% paste(collapse = ", ")) })



server = function(input, output, session) {
  
  # printLogJs <- function(x, ...) {
  #   
  #   logjs(x)
  #   
  #   T
  # }
  # 
  # addHandler(printLogJs)
  # 
  # 
  output$diffexpr_download <- downloadHandler(
    filename = "aging_diffexpr_data_first100.csv",
    content = function(file) {
      write.table(tablefordownloading[1:100,], file, row.names = FALSE, sep="\t")
    }
  )
  
  
  js <- c(
    "function(row, data, displayNum, index){",
    "  var x = data[4];",
    "  $('td:eq(4)', row).html(x.toExponential(2));",
    "  var x = data[5];",
    "  $('td:eq(5)', row).html(x.toExponential(2));",
    "  for (i=6; i<13; i++){",
    "    var x = data[i];",
    "    var a = i.toString(10);",
    "    $('td:eq(' + a + ')', row).html(function(){",
    "      try{",
    "        if (x < 0){",
    "          return x + ' <span style = \"color:red;font-size:22px\">&darr;</span>';",
    "        } else if (x > 0){",
    "          return x + ' <span style = \"color:lime;font-size:22px\">&uarr;</span>';",
    "        } else if (x === null){",
    "          return 'NA';",
    "        } else {",
    "          return x;",
    "        }",
    "      } catch(error){",
    "        return x;",
    "      }",
    "//    return x.toExponential(2) + ' <span style = \"color:red;font-size:22px\">&darr;</span>';",
    "    });",
    "    $('td:eq(' + a + ')', row).css('cursor', 'pointer');",
    "    $('td:eq(' + a + ')', row).hover(function(){",
    "      $(this).css('color', 'blue');",
    "//      $(this).css('text-decoration', 'underline');",
    "      $(this).css('font-weight', 'bold');",
    "      }, function(){",
    "      $(this).css('color', 'black');",
    "//      $(this).css('text-decoration', 'none currentcolor solid');",
    "      $(this).css('font-weight', 'normal');",
    "    });",
    "  }",
    "}"
  )
  
  jss = c(
    "function(settings, json) {",
    "  var table = this.DataTable();",
    "  table.on(\"click.dt\", \"td\", function() {",
    "    var row_=table.cell(this).index().row;",
    "    var celldata = table.cell(this).data();",
    "    var geneid = table.row(this).data()[2];",
    "    const re1 = /<[^<>]*>/ig;",
    "    geneid = geneid.replaceAll(re1, \"\");",
    "    var cell = table.cell(this);",
    "    var colindex = cell.index().column;",
    "    var colname = table.column(colindex).header().innerText;",
    "    let re = new RegExp(\"Brain|Muscle|Liver|Human|Mouse|Rat|All\");",
    "    if (re.test(colname) && celldata !== null){",
    "      var signname = colname.match(re)[0];",
    "      Shiny.setInputValue(\"metafor_signature_js\", signname);",
    "      Shiny.setInputValue(\"metafor_gene_js\", geneid);",
    "      $(\".navbar-nav  li:nth-child(3) a\").click();",
    "    }",
    "  });",
    "}"
  )
  
  # "function(settings, json) {",
  # "  var table = this.DataTable();",
  
  # "    if (colname.test(\"Brain|Muscle|Liver|Human|Mouse|Rat|All\")){",
  # "      var signname = colname.match(\"Brain|Muscle|Liver|Human|Mouse|Rat|All\");",
  
  
  brks = reactive(
    if (input$pval_hard_thres < input$pval_soft_thres){
      return(c(-1*input$pval_soft_thres, -1*input$pval_hard_thres, 0, input$pval_hard_thres, input$pval_soft_thres))
    } else {
      return(c(-1*input$pval_soft_thres, 0, input$pval_soft_thres))
    }
    )
  clrs = reactive(
    if (input$pval_hard_thres < input$pval_soft_thres){
      return(c("rgb(255,255,255)", "rgb(255,240,240)", "rgb(255,200,200)", "rgb(200,255,200)", "rgb(240,255,240)", "rgb(255,255,255)"))
    } else {
      return(c("rgb(255,255,255)", "rgb(255,240,240)", "rgb(240,255,240)", "rgb(255,255,255)"))
    }
  )
  # clrs <- round(seq(100, 255, length.out = length(brks) + 1), 0) %>%
  #   {paste0("rgb(255,", ., ",", ., ")")}
  
  demo_table = (agingsignaturesapp[["Brain"]] %>% filter(presence_total >= 5) %>% arrange(adj_pval))[1:2, c("genesymbol", "entrez", "logFC", "pvalue", "adjusted_pvalue", "Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC", "Brain_adj_pval", "Muscle_adj_pval", "Liver_adj_pval", "Mouse_adj_pval", "Human_adj_pval", "Rat_adj_pval", "All_adj_pval")]
  
  output$demo_table = DT::renderDataTable(
    DT::datatable(demo_table, filter="none", extensions = "FixedColumns", escape = FALSE, class = "cell-border stripe", options = list(lengthMenu = c(25, 50, 100), scrollX=600, ordering=T, columnDefs = list(list(visible=FALSE, targets=c(13:19))), rowCallback = JS(js), fixedColumns = list(leftColumns = 2)), selection = "none") %>%
      formatStyle(c("Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC"), valueColumns = c("Brain_adj_pval", "Muscle_adj_pval", "Liver_adj_pval", "Mouse_adj_pval", "Human_adj_pval", "Rat_adj_pval", "All_adj_pval"), backgroundColor = styleInterval(brks(), clrs())) %>%
      formatStyle(c("Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC"), backgroundColor = styleEqual(c(NA), c("white"))) %>%
      formatStyle(c(0:1), 'border-right' = 'solid 1.75px') %>%
      formatStyle(c(5:5), 'border-right' = 'solid 1.75px')
  )
  
  main_table = reactive({temp = agingsignaturesapp[[input$main_signature]] %>% filter(presence_total >= input$presence) %>% filter(adj_pval < input$pval_thres) %>% arrange(adj_pval)
                  temp$signif_total = rowSums((temp[, grepl("^(Mouse|Rat|Human|All|Brain|Muscle|Liver)_adj_pval$", colnames(temp))] < input$pval_soft_thres) & (temp[, grepl("^(Mouse|Rat|Human|All|Brain|Muscle|Liver)_adj_pval$", colnames(temp))] > -1*input$pval_soft_thres), na.rm = TRUE)
                  temp = temp %>% filter(signif_total >= input$num_signif)
                  if (length(input$upregulated_in > 0)){
                    upregcols = paste0(input$upregulated_in, "_logFC")
                    for (i in 1:length(upregcols)){
                      temp = temp %>% filter(!!as.symbol(upregcols[i]) > 0)
                    }
                  }
                  if (length(input$downregulated_in > 0)){
                    downregcols = paste0(input$downregulated_in, "_logFC")
                    for (i in 1:length(downregcols)){
                      temp = temp %>% filter(!!as.symbol(downregcols[i]) < 0)
                    }
                  }
                  return(temp[, c("genesymbol", "entrez", "logFC", "pvalue", "adjusted_pvalue", "Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC", "Brain_adj_pval", "Muscle_adj_pval", "Liver_adj_pval", "Mouse_adj_pval", "Human_adj_pval", "Rat_adj_pval", "All_adj_pval")])
  })
  
  
  #output$kek = reactive(length(colnames(main_table())))
  output$main_table = DT::renderDataTable(
    DT::datatable(main_table(), filter="none", extensions = "FixedColumns", escape = FALSE, class = "cell-border stripe", options = list(lengthMenu = c(25, 50, 100), scrollX=600, scrollY=600, ordering=T, columnDefs = list(list(visible=FALSE, targets=c((length(colnames(main_table()))-6):(length(colnames(main_table())))))), initComplete = JS(jss), rowCallback = JS(js), fixedColumns = list(leftColumns = 2)), selection = "none") %>% 
      formatStyle(c("Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC"), valueColumns = c("Brain_adj_pval", "Muscle_adj_pval", "Liver_adj_pval", "Mouse_adj_pval", "Human_adj_pval", "Rat_adj_pval", "All_adj_pval"), backgroundColor = styleInterval(brks(), clrs())) %>%
      formatStyle(c("Brain_logFC", "Muscle_logFC", "Liver_logFC", "Mouse_logFC", "Human_logFC", "Rat_logFC", "All_logFC"), backgroundColor = styleEqual(c(NA), c("white"))) %>%
      formatStyle(c(0:1), 'border-right' = 'solid 1.75px') %>%
      formatStyle(c(5:5), 'border-right' = 'solid 1.75px')
  )
  
  observeEvent(input$metafor_signature_js, 
               {updateSelectInput(session, inputId = "metafor_signature", selected = input$metafor_signature_js)
                })
  
  observeEvent(input$metafor_signature, 
               {temp = input$metafor_gene
                updateSelectizeInput(session, inputId = "metafor_gene", selected = temp, choices = metafor_genes[[input$metafor_signature]])})
  
  # observeEvent(input$metafor_gene_js, 
  #              {updateMultiInput(session, inputId = "metafor_gene", selected = input$metafor_gene_js)})
  observeEvent(input$metafor_gene_js, {
    #session$sendCustomMessage("metafor_gene_send", input$metafor_gene_js)
    updateSelectizeInput(session, inputId = "metafor_gene", selected = input$metafor_gene_js)
  })
  
  output$metafor_demo = renderPlotly({
    name = "All"
    
    # filter datasets for the individual signature:
    logFCmatrixchosen = logFCmatrixregr[, chosencols[[name]]]
    SEmatrixchosen = SEmatrixregr[, chosencols[[name]]]
    
    
    # get the deming coefs:
    minimums = c()
    for (i in 1:10){
      #deminglist[[name]][[i]] = deming_minimizer(logFCmatrixchosen)
      minimums = c(minimums, deminglist[[name]][[i]]$minimum)
    }
    kres = deminglist[[name]][[which.min(minimums)]]$coefs
    
    # normalize by deming coefficients:
    
    for (i in 1:length(colnames(logFCmatrixchosen))){
      SEmatrixchosen[,i] = SEmatrixchosen[,i] / kres[i]
      logFCmatrixchosen[,i] = logFCmatrixchosen[,i] / kres[i]
    }
    
    # discard bad boys:
    logFCmatrixchosen$NACount = rowSums(is.na(logFCmatrixchosen))
    if (name != "Liver"){
      #ggplot(logFCmatrixchosen, aes(x = NACount)) + geom_density() + geom_vline(xintercept = floor(length(colnames(logFCmatrixchosen))/2))
      goodboys = subset(rownames(logFCmatrixchosen), logFCmatrixchosen$NACount < floor(length(colnames(logFCmatrixchosen))/2))
    } else {
      #ggplot(logFCmatrixchosen, aes(x = NACount)) + geom_density() + geom_vline(xintercept = 4)
      goodboys = subset(rownames(logFCmatrixchosen), logFCmatrixchosen$NACount < 4)
    }
    logFCmatrixchosen = subset(logFCmatrixchosen, rownames(logFCmatrixchosen) %in% goodboys)
    SEmatrixchosen = subset(SEmatrixchosen, rownames(SEmatrixchosen) %in% goodboys)
    logFCmatrixchosen$NACount = NULL
    
    helpertable = as.data.frame(t(logFCmatrixchosen["72930",]), stringsAsFactors = FALSE)
    rownames(helpertable) = colnames(logFCmatrixchosen)
    colnames(helpertable) = c("logFC")
    helpertable$SE = as.numeric(as.character(t(SEmatrixchosen["72930",])))
    helpertable$source = as.character(sourcedata[rownames(helpertable),"dataset"])
    helpertable$dataset = rownames(helpertable)
    helpertable$species = sub("^([^_]+)_(.*)", "\\1", rownames(helpertable))
    helpertable$tissue = sub("^([^_]+)_([^_]+)_([^_]+)_(.*)", "\\3", rownames(helpertable))
    helpertable$tissue = sub("Frontalcortex", "Brain", helpertable$tissue)
    helpertable$tissue = sub("Cerebellum", "Brain", helpertable$tissue)
    helpertable$tissue = sub("LimbMuscle", "Muscle", helpertable$tissue)
    helpertable$sex = sub("^([^_]+)_([^_]+)_([^_]+)_([^_]+)(.*)", "\\4", rownames(helpertable))
    helpertable$sex = sub("NoSex", "Unknown", helpertable$sex)
    helpertable = left_join(helpertable, strain_map, by="source")
    if (sum(str_count(helpertable$dataset, "_")) != 0){
      helpertable[which(str_count(helpertable$dataset, "_") == 4),"strain"] = sub("^([^_]+)_([^_]+)_([^_]+)_([^_]+)_(.*)", "\\5", helpertable[which(str_count(helpertable$dataset, "_") == 4),"dataset"])
    }
    helpertable$colorby = helpertable[,"tissue"]
    helpertable$shapeby = helpertable[,"species"]
    helpertable = na.omit(helpertable)
    #helpertable = data.frame(helpertable, stringsAsFactors = FALSE)
    helpertable = helpertable %>% arrange(dataset)
    #helpertablelist[[input$metafor_gene]] = helpertable
    border = max(abs(helpertable$logFC)) + max(helpertable$SE)
    
    # fig <- plot_ly(data = helpertable, x = ~dataset, y = ~logFC, type = 'scatter', mode = 'markers',
    #                color= ~tissue, symbol= ~species, colors = gg_color_hue(length(unique(helpertable$tissue))),
    #                error_y = ~list(array = SE),
    #                marker = list(size=10))
    # 
    # 
    helpertable$dummy = rep("", length(rownames(helpertable)))
    
    rainbow_hues = gg_color_hue(length(unique(helpertable$colorby)))
    
    fig <- plot_ly(data = helpertable, colors = c("black", rainbow_hues)) %>%
      add_markers(data = helpertable, x = ~dataset, y = ~logFC, type = 'scatter',
                  color= ~colorby, legendgroup = "Tissue", hoverinfo = "none") %>%
      add_markers(data = helpertable, x = ~dataset, y = ~logFC, type = 'scatter',
                  symbol= ~shapeby, color = ~dummy, legendgroup = "Species", hoverinfo = "none")
    for (i in 1:length(unique(helpertable$colorby))){
      one_color = unique(helpertable$colorby)[i]
      for (j in 1:length(unique(helpertable[which(helpertable$colorby == one_color),]$shapeby))){
        one_shape = unique(helpertable[which(helpertable$colorby == one_color),]$shapeby)[j]
        fig = fig %>% add_markers(data = helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),], x = ~dataset, y = ~logFC, type = 'scatter',
                                  error_y = ~list(array = SE),
                                  text = as.character(formatC(helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$SE, digits = 2, format="f")),
                                  color = ~colorby,
                                  colors = helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$colorby,
                                  symbol = ~shapeby,
                                  symbols = helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$shapeby,
                                  marker = list(size=10),
                                  showlegend=FALSE,
                                  hovertemplate = paste(
                                    '<br><b>LogFC</b>: %{y:.2f}&plusmn;', formatC(helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$SE, digits = 2, format = "f"), '<extra></extra>',
                                    '<br><b>Species</b>: ', helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$species,
                                    '<br><b>Study</b>: ', helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$source,
                                    '<br><b>Tissue</b>: ', helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$tissue,
                                    '<br><b>Sex</b>: ', helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$sex,
                                    '<br><b>Strain</b>: ', helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$strain)
        )
      }
    }
    
    fig %>% layout(title = paste0("Mixed effects model for gene ", maintable[which(maintable$entrez == "72930"),'genesymbol'], " in ", "Compound", " signature"),  shapes = list(hline(agingsignatures_v3[[name]]["72930", "logFC"])), xaxis= list(showticklabels = FALSE, categoryorder = "array", categoryarray = helpertable[,"colorby"]))
  })  

  output$metafor_plot <- renderPlotly({
    name = input$metafor_signature
      
    # filter datasets for the individual signature:
    logFCmatrixchosen = logFCmatrixregr[, chosencols[[name]]]
    SEmatrixchosen = SEmatrixregr[, chosencols[[name]]]
    
    
    # get the deming coefs:
    minimums = c()
    for (i in 1:10){
      #deminglist[[name]][[i]] = deming_minimizer(logFCmatrixchosen)
      minimums = c(minimums, deminglist[[name]][[i]]$minimum)
    }
    kres = deminglist[[name]][[which.min(minimums)]]$coefs
    
    # normalize by deming coefficients:
    
    for (i in 1:length(colnames(logFCmatrixchosen))){
      SEmatrixchosen[,i] = SEmatrixchosen[,i] / kres[i]
      logFCmatrixchosen[,i] = logFCmatrixchosen[,i] / kres[i]
    }
    
    # discard bad boys:
    logFCmatrixchosen$NACount = rowSums(is.na(logFCmatrixchosen))
    if (name != "Liver"){
      #ggplot(logFCmatrixchosen, aes(x = NACount)) + geom_density() + geom_vline(xintercept = floor(length(colnames(logFCmatrixchosen))/2))
      goodboys = subset(rownames(logFCmatrixchosen), logFCmatrixchosen$NACount < floor(length(colnames(logFCmatrixchosen))/2))
    } else {
      #ggplot(logFCmatrixchosen, aes(x = NACount)) + geom_density() + geom_vline(xintercept = 4)
      goodboys = subset(rownames(logFCmatrixchosen), logFCmatrixchosen$NACount < 4)
    }
    logFCmatrixchosen = subset(logFCmatrixchosen, rownames(logFCmatrixchosen) %in% goodboys)
    SEmatrixchosen = subset(SEmatrixchosen, rownames(SEmatrixchosen) %in% goodboys)
    logFCmatrixchosen$NACount = NULL
    
    if (input$metafor_gene == ""){
      fig <- plot_ly() %>%
        layout(xaxis = list("visible" = FALSE), yaxis = list("visible" = FALSE), annotations = list("text" = "Please select a gene",
                                                                                                    "xref" = "paper",
                                                                                                    "yref" = "paper",
                                                                                                    "showarrow" = FALSE,
                                                                                                    "font" = list("size" = 28)
        )) %>% config(displayModeBar = F)
    } else if (!(input$metafor_gene %in% rownames(logFCmatrixchosen)) | is.na(agingsignatures_v3[[name]][input$metafor_gene, "logFC"])){
      
      fig <- plot_ly() %>%
        layout(xaxis = list("visible" = FALSE), yaxis = list("visible" = FALSE), annotations = list("text" = "Gene is not present in this signature",
                                                                                                    "xref" = "paper",
                                                                                                    "yref" = "paper",
                                                                                                    "showarrow" = FALSE,
                                                                                                    "font" = list("size" = 28)
        )) %>% config(displayModeBar = F)
      
    } else {
      
      helpertable = as.data.frame(t(logFCmatrixchosen[input$metafor_gene,]), stringsAsFactors = FALSE)
      rownames(helpertable) = colnames(logFCmatrixchosen)
      colnames(helpertable) = c("logFC")
      helpertable$SE = as.numeric(as.character(t(SEmatrixchosen[input$metafor_gene,])))
      helpertable$source = as.character(sourcedata[rownames(helpertable),"dataset"])
      helpertable$dataset = rownames(helpertable)
      helpertable$species = sub("^([^_]+)_(.*)", "\\1", rownames(helpertable))
      helpertable$tissue = sub("^([^_]+)_([^_]+)_([^_]+)_(.*)", "\\3", rownames(helpertable))
      helpertable$tissue = sub("Frontalcortex", "Brain", helpertable$tissue)
      helpertable$tissue = sub("Cerebellum", "Brain", helpertable$tissue)
      helpertable$tissue = sub("LimbMuscle", "Muscle", helpertable$tissue)
      helpertable$sex = sub("^([^_]+)_([^_]+)_([^_]+)_([^_]+)(.*)", "\\4", rownames(helpertable))
      helpertable$sex = sub("NoSex", "Unknown", helpertable$sex)
      helpertable = left_join(helpertable, strain_map, by="source")
      if (sum(str_count(helpertable$dataset, "_")) != 0){
        helpertable[which(str_count(helpertable$dataset, "_") == 4),"strain"] = sub("^([^_]+)_([^_]+)_([^_]+)_([^_]+)_(.*)", "\\5", helpertable[which(str_count(helpertable$dataset, "_") == 4),"dataset"])
      }
      helpertable$colorby = helpertable[,input$metafor_color]
      helpertable$shapeby = helpertable[,input$metafor_shape]
      helpertable = na.omit(helpertable)
      #helpertable = data.frame(helpertable, stringsAsFactors = FALSE)
      helpertable = helpertable %>% arrange(dataset)
      #helpertablelist[[input$metafor_gene]] = helpertable
      border = max(abs(helpertable$logFC)) + max(helpertable$SE)
      
      # fig <- plot_ly(data = helpertable, x = ~dataset, y = ~logFC, type = 'scatter', mode = 'markers',
      #                color= ~tissue, symbol= ~species, colors = gg_color_hue(length(unique(helpertable$tissue))),
      #                error_y = ~list(array = SE),
      #                marker = list(size=10))
      # 
      # 
      helpertable$dummy = rep("", length(rownames(helpertable)))
      
      rainbow_hues = gg_color_hue(length(unique(helpertable$colorby)))
      
      fig <- plot_ly(data = helpertable, colors = c("black", rainbow_hues)) %>%
        add_markers(data = helpertable, x = ~dataset, y = ~logFC, type = 'scatter',
                    color= ~colorby, legendgroup = "Tissue", hoverinfo = "none") %>%
        add_markers(data = helpertable, x = ~dataset, y = ~logFC, type = 'scatter',
                    symbol= ~shapeby, color = ~dummy, legendgroup = "Species", hoverinfo = "none") %>%
        onRender("
                 function(el, x) {
                   el.on('plotly_click', function(d) {
                     var point = d.points[0];
                     //console.log('Click: ', d);
                     var mypoint = point.x;
                     Shiny.setInputValue(\"dataset_id_js\", mypoint);
                   });
                 }") %>% config(displayModeBar = F)
        

      for (i in 1:length(unique(helpertable$colorby))){
        one_color = unique(helpertable$colorby)[i]
        for (j in 1:length(unique(helpertable[which(helpertable$colorby == one_color),]$shapeby))){
          one_shape = unique(helpertable[which(helpertable$colorby == one_color),]$shapeby)[j]
          fig = fig %>% add_markers(data = helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),], x = ~dataset, y = ~logFC, type = 'scatter',
                                    error_y = ~list(array = SE),
                                    text = as.character(formatC(helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$SE, digits = 2, format="f")),
                                    color = ~colorby,
                                    colors = helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$colorby,
                                    symbol = ~shapeby,
                                    symbols = helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$shapeby,
                                    marker = list(size=10),
                                    showlegend=FALSE,
                                    hovertemplate = paste(
                                      '<br><b>LogFC</b>: %{y:.2f}&plusmn;', formatC(helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$SE, digits = 2, format = "f"), '<extra></extra>',
                                      '<br><b>Species</b>: ', helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$species,
                                      '<br><b>Study</b>: ', helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$source,
                                      '<br><b>Tissue</b>: ', helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$tissue,
                                      '<br><b>Sex</b>: ', helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$sex,
                                      '<br><b>Strain</b>: ', helpertable[which(helpertable$colorby == one_color & helpertable$shapeby == one_shape),]$strain)
          )
        }
      }
      
      if (name == "All"){
        fig %>% layout(title = paste0("Mixed effects model for gene ", maintable[which(maintable$entrez == input$metafor_gene),'genesymbol'], " in ", "Compound", " signature"),  shapes = list(hline(agingsignatures_v3[[name]][input$metafor_gene, "logFC"])), xaxis= list(showticklabels = FALSE, categoryorder = "array", categoryarray = helpertable[,"colorby"]))
      } else {
        fig %>% layout(title = paste0("Mixed effects model for gene ", maintable[which(maintable$entrez == input$metafor_gene),'genesymbol'], " in ", name, " signature"),  shapes = list(hline(agingsignatures_v3[[name]][input$metafor_gene, "logFC"])), xaxis= list(showticklabels = FALSE, categoryorder = "array", categoryarray = helpertable[,"colorby"]))
      }
      
    }
  })
  
  output$dataset_plot = renderPlotly({
    #event.data <- event_data(event = "plotly_click", source = "metafor_plot")
    #View(event.data)
    if (length(input$dataset_id_js) != 0){

      # EXPR PLOT
      genename = input$metafor_gene
      dataset_id = input$dataset_id_js


      if (str_count(dataset_id, "_") == 3){
        strain = strain_map[which(strain_map$source == sub("([^_]+)_([^_]+)_.*", "\\2", dataset_id)), "strain"]
        sex = sub("([^_]+)_([^_]+)_([^_]+)_([^_]+)", "\\4", dataset_id)
      } else if (str_count(dataset_id, "_") == 4){
        sex = sub("([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)", "\\4", dataset_id)
        strain = sub("([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)", "\\5", dataset_id)
      }

      topgenes = logFCtablelist[[dataset_id]]
      if ((is.numeric(sexyphenolist[[dataset_id]]$Age)) | (sum(grepl("^[A-Za-z]+$", sexyphenolist[[dataset_id]]$Age)) == 0)){
        copps = as.data.frame(cbind(colnames(sexyexprlist[[dataset_id]]), as.numeric(sexyexprlist[[dataset_id]][genename,]), as.integer(sexyphenolist[[dataset_id]]$Age)))
        copps$V3 = as.integer(as.character(copps$V3))
      } else {
        copps = as.data.frame(cbind(colnames(sexyexprlist[[dataset_id]]), as.numeric(sexyexprlist[[dataset_id]][genename,]), casefold(sexyphenolist[[dataset_id]]$Age, upper=F)))
        copps$V3 = relevel(factor(copps$V3), ref="young")
      }
      copps$V2 = as.double(as.character(copps$V2))
      # copps$V2 = as.numeric(sub(",", ".", copps$V2))
      copps = copps %>% arrange(V2)
      
      if (sum(is.na(copps$V2)) == 0){
        kek = ggplot(copps, aes(x = copps$V3, y = copps$V2)) + geom_point(color = "blue", aes(text=NULL)) + labs(x = paste0("age in ", ageunittable[sub("([^_]+)_([^_]+)_.*", "\\2", dataset_id), "unit"]), y = "normalized expression (log scale)",title = paste0("Gene symbol: ", maintable[which(maintable$entrez == genename),'genesymbol'], "\n",
                                                                                                                                                                                                                                                                   "Species: ", sub("([^_]+)_([^_]+)_([^_]+)_.*", "\\1", dataset_id), "\n",
                                                                                                                                                                                                                                                                   "Study: ", sub("([^_]+)_([^_]+)_([^_]+)_.*", "\\2", dataset_id), "\n",
                                                                                                                                                                                                                                                                   "Tissue: ", sub("([^_]+)_([^_]+)_([^_]+)_.*", "\\3", dataset_id), "\n",
                                                                                                                                                                                                                                                                   "Sex: ", sex, "\n",
                                                                                                                                                                                                                                                                   "Strain: ", strain, "\n")) +
          theme_minimal() + geom_smooth(method='lm', formula= y~x) + theme(plot.title = element_text(size=9), plot.margin = margin(t=100, l=30, b=5))
        
        fig <- ggplotly(kek) %>% config(displayModeBar = F)
        style(fig, hoverinfo = "none")
      } else {
        fig <- plot_ly() %>%
          layout(xaxis = list("visible" = FALSE), yaxis = list("visible" = FALSE), annotations = list("text" = "Plot not available",
                                                                                                      "xref" = "paper",
                                                                                                      "yref" = "paper",
                                                                                                      "showarrow" = FALSE,
                                                                                                      "font" = list("size" = 28)
          )) %>% config(displayModeBar = F)
      }
      
      
      # fig %>% layout(title=list("y" = 0.9, "font" = list(
      #
      #   "size" = 12
      #
      #   )), margin = list(
      #     "t" = 125
      #   )
      # )

    } else {
      fig <- plot_ly() %>%
        layout(xaxis = list("visible" = FALSE), yaxis = list("visible" = FALSE), annotations = list("text" = "Please click on a dataset",
                                                                                                    "xref" = "paper",
                                                                                                    "yref" = "paper",
                                                                                                    "showarrow" = FALSE,
                                                                                                    "font" = list("size" = 28)
        )) %>% config(displayModeBar = F)
      #ggplot() + theme_void()
    }


  })
  
  output$datasetplot_demo = renderPlotly({
    #event.data <- event_data(event = "plotly_click", source = "metafor_plot")
    #View(event.data)
      
    # EXPR PLOT
    genename = "12268"
    dataset_id = "Mouse_GSE132040_Brain_Male"
    
    
    if (str_count(dataset_id, "_") == 3){
      strain = strain_map[which(strain_map$source == sub("([^_]+)_([^_]+)_.*", "\\2", dataset_id)), "strain"]
      sex = sub("([^_]+)_([^_]+)_([^_]+)_([^_]+)", "\\4", dataset_id)
    } else if (str_count(dataset_id, "_") == 4){
      sex = sub("([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)", "\\4", dataset_id)
      strain = sub("([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)", "\\5", dataset_id)
    }
    
    topgenes = logFCtablelist[[dataset_id]]
    if ((is.numeric(sexyphenolist[[dataset_id]]$Age)) | (sum(grepl("^[A-Za-z]+$", sexyphenolist[[dataset_id]]$Age)) == 0)){
      copps = as.data.frame(cbind(colnames(sexyexprlist[[dataset_id]]), as.numeric(sexyexprlist[[dataset_id]][genename,]), as.integer(sexyphenolist[[dataset_id]]$Age)))
      copps$V3 = as.integer(as.character(copps$V3))
    } else {
      copps = as.data.frame(cbind(colnames(sexyexprlist[[dataset_id]]), as.numeric(sexyexprlist[[dataset_id]][genename,]), casefold(sexyphenolist[[dataset_id]]$Age, upper=F)))
      copps$V3 = relevel(factor(copps$V3), ref="young")
    }
    copps$V2 = as.double(as.character(copps$V2))
    # copps$V2 = as.numeric(sub(",", ".", copps$V2))
    copps = copps %>% arrange(V2)
    
    kek = ggplot(copps, aes(x = copps$V3, y = copps$V2)) + geom_point(color = "blue", aes(text=NULL)) + labs(x = paste0("age in ", ageunittable[sub("([^_]+)_([^_]+)_.*", "\\2", dataset_id), "unit"]), y = "normalized expression (log scale)",title = paste0("Gene symbol: ", maintable[which(maintable$entrez == genename),'genesymbol'], "\n",
                                                                                                                                                                                                                                                               "Species: ", sub("([^_]+)_([^_]+)_([^_]+)_.*", "\\1", dataset_id), "\n",
                                                                                                                                                                                                                                                               "Study: ", sub("([^_]+)_([^_]+)_([^_]+)_.*", "\\2", dataset_id), "\n",
                                                                                                                                                                                                                                                               "Tissue: ", sub("([^_]+)_([^_]+)_([^_]+)_.*", "\\3", dataset_id), "\n",
                                                                                                                                                                                                                                                               "Sex: ", sex, "\n",
                                                                                                                                                                                                                                                               "Strain: ", strain, "\n")) +
      theme_minimal() + geom_smooth(method='lm', formula= y~x) + theme(plot.title = element_text(size=9), plot.margin = margin(t=100, l=30, b=5))
    
    fig <- ggplotly(kek) %>% config(displayModeBar = F)
    style(fig, hoverinfo = "none")
    
  })
  
  
  # 
  # output$dataset_plot = renderPlot({
  #   #event.data <- event_data(event = "plotly_click", source = "metafor_plot")
  #   #View(event.data)
  #   
  #   # EXPR PLOT
  #   genename = input$metafor_gene
  #   dataset_id = input$dataset_id_js
  #   topgenes = logFCtablelist[[dataset_id]]
  #   copps = as.data.frame(cbind(colnames(sexyexprlist[[dataset_id]]), as.numeric(sexyexprlist[[dataset_id]][genename,]), as.integer(sexyphenolist[[dataset_id]]$Age)))
  #   copps$V3 = as.integer(as.character(copps$V3))
  #   copps$V2 = as.double(as.character(copps$V2))
  #   # copps$V2 = as.numeric(sub(",", ".", copps$V2))
  #   copps = copps %>% arrange(V2)
  #   ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age", x = "age", y = "Expression",title = paste("entrez ID", genename, sep = " "))
  #   
  # })
  
  #options(shiny.error = NULL)
  
}

shinyApp(ui = ui, server = server)