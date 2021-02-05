library(tidyverse)

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

ggplot(maintable, aes(x = presence_species)) + geom_histogram()

ggplot(maintable, aes(x = presence_tissues)) + geom_histogram()

ggplot(maintable, aes(x = presence_total)) + geom_histogram()

sum(na.omit(maintable$geommean_pval_tissues)< 0.05)

sum(na.omit(maintable$geommean_pval_species)< 0.05)



# make logFCmatrix and SEmatrix and plot mixed model examples

library(biomaRt)
library(tidyverse)
library(deming)
library(reshape2)
library(metafor)

logFClist = readRDS(file = "logFClist.rds")

# create entrez mappings between mouse and rat, and mouse and human

# default:
ensembl <- useMart("ensembl")
# in case www.ensembl.org is under maintenance:
ensembl = useEnsembl("ensembl", host = "uswest.ensembl.org")

datasets <- listDatasets(ensembl)
rat_dataset = useDataset("rnorvegicus_gene_ensembl", mart=ensembl)
human_dataset = useDataset("hsapiens_gene_ensembl", mart=ensembl)
mouse_dataset = useDataset("mmusculus_gene_ensembl", mart=ensembl)

# default:
rat = useMart("ensembl","rnorvegicus_gene_ensembl")
human = useMart("ensembl", "hsapiens_gene_ensembl")
mouse = useMart("ensembl","mmusculus_gene_ensembl")

# in case www.ensembl.org is under maintenance:
rat = useEnsembl("ensembl","rnorvegicus_gene_ensembl", host = "uswest.ensembl.org")
human = useEnsembl("ensembl", "hsapiens_gene_ensembl", host = "uswest.ensembl.org")
mouse = useEnsembl("ensembl","mmusculus_gene_ensembl", host = "uswest.ensembl.org")


Rat_to_mouse_orthologs <- getLDS(attributes=c("entrezgene_id"),
                                 mart=rat_dataset,attributesL=c("entrezgene_id"), martL=mouse_dataset)
colnames(Rat_to_mouse_orthologs) <- c("Rat_Entrez","Mouse_Entrez")
Mouse_to_rat_orthologs <- getLDS(attributes=c("entrezgene_id"), filters="entrezgene_id",
                                 values=as.character(unique(Rat_to_mouse_orthologs$Mouse_Entrez)),
                                 mart=mouse_dataset,attributesL=c("entrezgene_id"), martL=rat_dataset)
colnames(Mouse_to_rat_orthologs) <- c("Mouse_Entrez","Rat_Entrez")
Rat_to_mouse_orthologs1 = Rat_to_mouse_orthologs %>% group_by(Rat_Entrez) %>% filter(n() == 1)
Rat_to_mouse_orthologs1 = na.omit(Rat_to_mouse_orthologs1)
Mouse_to_rat_orthologs1 = subset(Mouse_to_rat_orthologs, Mouse_Entrez %in% Rat_to_mouse_orthologs1$Mouse_Entrez)
Mouse_to_rat_orthologs1 = Mouse_to_rat_orthologs1 %>% group_by(Mouse_Entrez) %>% filter(n() == 1)
mouse_rat_entrez_map = na.omit(Mouse_to_rat_orthologs1)
mouse_rat_entrez_map = as.data.frame(mouse_rat_entrez_map)
mouse_rat_entrez_map = mouse_rat_entrez_map %>% mutate_all(as.character)

Human_to_mouse_orthologs <- getLDS(attributes=c("entrezgene_id"),
                                   mart=human_dataset,attributesL=c("entrezgene_id"), martL=mouse_dataset)
colnames(Human_to_mouse_orthologs) <- c("Human_Entrez","Mouse_Entrez")
Mouse_to_human_orthologs <- getLDS(attributes=c("entrezgene_id"), filters="entrezgene_id",
                                   values=as.character(unique(Human_to_mouse_orthologs$Mouse_Entrez)),
                                   mart=mouse_dataset,attributesL=c("entrezgene_id"), martL=human_dataset)
colnames(Mouse_to_human_orthologs) <- c("Mouse_Entrez","Human_Entrez")
Human_to_mouse_orthologs1 = Human_to_mouse_orthologs %>% group_by(Human_Entrez) %>% filter(n() == 1)
Human_to_mouse_orthologs1 = na.omit(Human_to_mouse_orthologs1)
Mouse_to_human_orthologs1 = subset(Mouse_to_human_orthologs, Mouse_Entrez %in% Human_to_mouse_orthologs1$Mouse_Entrez)
Mouse_to_human_orthologs1 = Mouse_to_human_orthologs1 %>% group_by(Mouse_Entrez) %>% filter(n() == 1)
mouse_human_entrez_map = na.omit(Mouse_to_human_orthologs1)
mouse_human_entrez_map = as.data.frame(mouse_human_entrez_map)
mouse_human_entrez_map = mouse_human_entrez_map %>% mutate_all(as.character)


#remove bad datasets:
logFClist1 = logFClist
logFClist1$Mouse$GSE53959 = NULL
logFClist1$Human$GSE40645 = NULL
logFClist1$Human$GSE5086 = NULL
#logFClist1$Mouse$`E-MTAB-3374` = NULL
#logFClist1$Human$GSE9103 = NULL

# convert to orthologs and make totalrownames:
totalrownames = union(rownames(logFClist$Mouse$GSE6591$Lung$Male$DBA2J), rownames(logFClist$Mouse$GSE6591$Lung$Male$C57BL6J))
logFClist1$Mouse$GSE6591 = NULL
for (species in names(logFClist1)){
  for (dataset in names(logFClist1[[species]])){
    for (tissue in names(logFClist1[[species]][[dataset]])){
      for (sex in names(logFClist1[[species]][[dataset]][[tissue]])){
        if (species == "Rat"){
          exd <- logFClist1[[species]][[dataset]][[tissue]][[sex]]
          rownames(exd) = as.character(rownames(exd))
          exd$Rat_Entrez = rownames(exd)
          exd <- left_join(exd, mouse_rat_entrez_map)
          exd = na.omit(exd, cols=Mouse_Entrez)
          logFClist1[[species]][[dataset]][[tissue]][[sex]] <- exd %>% remove_rownames() %>% column_to_rownames(var = "Mouse_Entrez")
        }
        if (species == "Human"){
          exd <- logFClist1[[species]][[dataset]][[tissue]][[sex]]
          rownames(exd) = as.character(rownames(exd))
          exd$Human_Entrez = rownames(exd)
          exd <- left_join(exd, mouse_human_entrez_map)
          exd = na.omit(exd, cols=Mouse_Entrez)
          logFClist1[[species]][[dataset]][[tissue]][[sex]] <- exd %>% remove_rownames() %>% column_to_rownames(var = "Mouse_Entrez")
        }
        totalrownames = union(totalrownames, rownames(logFClist1[[species]][[dataset]][[tissue]][[sex]]))
        #        logFCmatrix = merge(logFCmatrix, logFClist1[[species]][[dataset]][[tissue]][[sex]]["logFC"], by=0, all=TRUE)
        #        logFCmatrix = logFCmatrix %>% column_to_rownames(var = "Row.names")
        #        colnames(logFCmatrix)[j + 2] = paste0(species, "_", dataset, "_", tissue, "_", sex)
        #        j = j + 1
      }
    }
  }
}

# create logFCunlisted:
logFCunlisted = list()
logFCunlisted[["Mouse_GSE6591_Lung_Male_DBA2J"]] = logFClist$Mouse$GSE6591$Lung$Male$DBA2J
logFCunlisted[["Mouse_GSE6591_Lung_Male_C57BL6J"]] = logFClist$Mouse$GSE6591$Lung$Male$C57BL6J
for (species in names(logFClist1)){
  for (dataset in names(logFClist1[[species]])){
    for (tissue in names(logFClist1[[species]][[dataset]])){
      for (sex in names(logFClist1[[species]][[dataset]][[tissue]])){
        logFCunlisted[[paste(species, dataset, tissue, sex, sep = "_")]] = logFClist1[[species]][[dataset]][[tissue]][[sex]]
      }
    }
  }
}

# create logFCmatrix:
logFCmatrix = matrix(nrow = length(totalrownames), ncol = length(logFCunlisted))
rownames(logFCmatrix) = totalrownames
colnames(logFCmatrix) = names(logFCunlisted)

SEmatrixregr = matrix(nrow = length(totalrownames), ncol = length(logFCunlisted))
rownames(SEmatrixregr) = totalrownames
colnames(SEmatrixregr) = names(logFCunlisted)

pvalmatrix = matrix(nrow = length(totalrownames), ncol = length(logFCunlisted))
rownames(pvalmatrix) = totalrownames
colnames(pvalmatrix) = names(logFCunlisted)

for (name in names(logFCunlisted)){
  pvalmatrix[rownames(logFCunlisted[[name]]), name] = logFCunlisted[[name]]$adj.P.Val
  logFCmatrix[rownames(logFCunlisted[[name]]), name] = logFCunlisted[[name]]$logFC
  SEmatrixregr[rownames(logFCunlisted[[name]]), name] = logFCunlisted[[name]]$SE
}

pvalmatrix = as.data.frame(pvalmatrix)
logFCmatrix = as.data.frame(logFCmatrix)
SEmatrixregr = as.data.frame(SEmatrixregr)
logFCmatrixregr = logFCmatrix

# fix spaces in logFCmatrixregr colnames:
colnames(logFCmatrixregr) = sub(" ", "", colnames(logFCmatrixregr))
colnames(SEmatrixregr) = sub(" ", "", colnames(SEmatrixregr))

names(logFCunlisted) = sub(" ", "", names(logFCunlisted))

# normalize by sd:

for (i in 1:length(colnames(logFCmatrixregr))){
  SEmatrixregr[,i] = SEmatrixregr[,i] / sd(na.omit(logFCmatrixregr[,i]))
  logFCmatrixregr[,i] = logFCmatrixregr[,i] / sd(na.omit(logFCmatrixregr[,i]))
}

# get source table:
sourcedata = as.data.frame(colnames(logFCmatrixregr))
rownames(sourcedata) = colnames(logFCmatrixregr)
colnames(sourcedata) = "kekkekkek"
sourcedata = sourcedata %>% separate(kekkekkek, c(NA, "dataset", NA), sep = "_")

helpertablelist = list()
# plot examples:
geneids = rownames(agingsignatures_v3[["All"]])
for (i in 1:length(geneids)){
  helpertable = as.data.frame(t(logFCmatrixregr[geneids[i],]))
  rownames(helpertable) = colnames(logFCmatrixregr)
  colnames(helpertable) = c("logFC")
  helpertable$SE = t(SEmatrixregr[geneids[i],])
  helpertable$source = as.factor(sourcedata[rownames(helpertable),"dataset"])
  helpertable$dataset = rownames(helpertable)
  helpertable$species = sub("^([^_]+)_(.*)", "\\1", rownames(helpertable))
  helpertable$tissue = sub("^([^_]+)_([^_]+)_([^_]+)_(.*)", "\\3", rownames(helpertable))
  helpertable$tissue = sub("Frontalcortex", "Brain", helpertable$tissue)
  helpertable$tissue = sub("Cerebellum", "Brain", helpertable$tissue)
  helpertable = na.omit(helpertable)
  helpertable = helpertable %>% arrange(dataset)
  helpertablelist[[geneids[i]]] = helpertable
  border = max(abs(helpertable$logFC)) + max(helpertable$SE)
  ggheatmap = ggplot(helpertable, aes(x = dataset, y = logFC, color = tissue, shape = species, label = source)) + geom_pointrange(aes(ymin = logFC - SE, ymax = logFC + SE)) + geom_hline(yintercept = agingsignatures_v3[[name]][geneids[i], "logFC"], colour = "red") +
    geom_hline(yintercept = 0) + ylim(-border, border) + theme_minimal()# + geom_text_repel(aes(label=source))
  #print(ggheatmap)
  pdf(paste0("./newplots/signatureplots1/", name, "/mixedmodelexample_shapes", geneids[i], ".pdf"))
  print(ggheatmap)
  dev.off()
  print(paste0("I'm ", as.character(i / length(geneids) * 100), "% done."))
}

