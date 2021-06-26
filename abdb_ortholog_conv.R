load("exprlist2.RData")
load("phenolist2.RData")
load("logFCtablelist2.RData")

library(biomaRt)
library(tidyverse)
library(deming)
library(reshape2)
library(metafor)

# default:
ensembl <- useMart("ensembl")
# in case www.ensembl.org is under maintenance:
ensembl = useEnsembl("ensembl", host = "uswest.ensembl.org")

#use archived version:
ensembl = useEnsembl("ensembl", host="jan2020.archive.ensembl.org")

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
# archived version:
rat = useEnsembl("ensembl","rnorvegicus_gene_ensembl", host = "jan2020.archive.ensembl.org")
human = useEnsembl("ensembl", "hsapiens_gene_ensembl", host = "jan2020.archive.ensembl.org")
mouse = useEnsembl("ensembl","mmusculus_gene_ensembl", host = "jan2020.archive.ensembl.org")


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


# convert to orthologs and make totalrownames:

sexyexprlistid = sexyexprlist
logFCtablelistid = logFCtablelist

#totalrownames = union(rownames(logFClist$Mouse$GSE6591$Lung$Male$DBA2J), rownames(logFClist$Mouse$GSE6591$Lung$Male$C57BL6J))
#logFClist1$Mouse$GSE6591 = NULL
for (name in names(sexyexprlist)){
  species = sub("([^_]+)_.*", "\\1", name)
  if (species == "Rat"){
    exd <- sexyexprlist[[name]]
    rownames(exd) = as.character(rownames(exd))
    exd$Rat_Entrez = rownames(exd)
    exd <- left_join(exd, mouse_rat_entrez_map)
    exd = na.omit(exd, cols=Mouse_Entrez)
    exd$Rat_Entrez = NULL
    sexyexprlistid[[name]] <- exd %>% remove_rownames() %>% column_to_rownames(var = "Mouse_Entrez")
  }
  if (species == "Human"){
    exd <- sexyexprlist[[name]]
    rownames(exd) = as.character(rownames(exd))
    exd$Human_Entrez = rownames(exd)
    exd <- left_join(exd, mouse_human_entrez_map)
    exd = na.omit(exd, cols=Mouse_Entrez)
    exd$Human_Entrez = NULL
    sexyexprlistid[[name]] <- exd %>% remove_rownames() %>% column_to_rownames(var = "Mouse_Entrez")
  }
  #totalrownames = union(totalrownames, rownames(logFClist1[[species]][[dataset]][[tissue]][[sex]]))
  #        logFCmatrix = merge(logFCmatrix, logFClist1[[species]][[dataset]][[tissue]][[sex]]["logFC"], by=0, all=TRUE)
  #        logFCmatrix = logFCmatrix %>% column_to_rownames(var = "Row.names")
  #        colnames(logFCmatrix)[j + 2] = paste0(species, "_", dataset, "_", tissue, "_", sex)
  #        j = j + 1
}
  

for (name in names(logFCtablelist)){
  species = sub("([^_]+)_.*", "\\1", name)
  if (species == "Rat"){
    exd <- logFCtablelist[[name]]
    rownames(exd) = as.character(rownames(exd))
    exd$Rat_Entrez = rownames(exd)
    exd <- left_join(exd, mouse_rat_entrez_map)
    exd = na.omit(exd, cols=Mouse_Entrez)
    exd$Rat_Entrez = NULL
    logFCtablelistid[[name]] <- exd %>% remove_rownames() %>% column_to_rownames(var = "Mouse_Entrez")
  }
  if (species == "Human"){
    exd <- logFCtablelist[[name]]
    rownames(exd) = as.character(rownames(exd))
    exd$Human_Entrez = rownames(exd)
    exd <- left_join(exd, mouse_human_entrez_map)
    exd = na.omit(exd, cols=Mouse_Entrez)
    exd$Human_Entrez = NULL
    logFCtablelistid[[name]] <- exd %>% remove_rownames() %>% column_to_rownames(var = "Mouse_Entrez")
  }
  #totalrownames = union(totalrownames, rownames(logFClist1[[species]][[dataset]][[tissue]][[sex]]))
  #        logFCmatrix = merge(logFCmatrix, logFClist1[[species]][[dataset]][[tissue]][[sex]]["logFC"], by=0, all=TRUE)
  #        logFCmatrix = logFCmatrix %>% column_to_rownames(var = "Row.names")
  #        colnames(logFCmatrix)[j + 2] = paste0(species, "_", dataset, "_", tissue, "_", sex)
  #        j = j + 1
}

logFCtablelist = logFCtablelistid
sexyexprlist = sexyexprlistid

save(logFCtablelist, file = "logFCtablelist3.RData")
save(sexyexprlist, file = "exprlist3.RData")

ageunittable <-as.data.frame(matrix(nrow=0,ncol=2))

ageunittable = rbind(ageunittable, c("GSE9103", "years"))
ageunittable = rbind(ageunittable, c("E-MTAB-3374", "months"))
ageunittable = rbind(ageunittable, c("GSE123981", "months"))
ageunittable = rbind(ageunittable, c("GSE3150", "months"))
ageunittable = rbind(ageunittable, c("GSE6591", "months"))
ageunittable = rbind(ageunittable, c("GSE74463", "weeks"))
ageunittable = rbind(ageunittable, c("PRJNA281127", "months"))
ageunittable = rbind(ageunittable, c("GSE53960", "weeks"))
ageunittable = rbind(ageunittable, c("PRJNA516151", "months"))
ageunittable = rbind(ageunittable, c("GSE66715", "months"))
ageunittable = rbind(ageunittable, c("GSE11291", "months"))
ageunittable = rbind(ageunittable, c("GSE34378", "weeks"))
ageunittable = rbind(ageunittable, c("GSE27625", "months"))
ageunittable = rbind(ageunittable, c("GSE12480", "months"))
ageunittable = rbind(ageunittable, c("GSE36192", "years"))
ageunittable = rbind(ageunittable, c("GSE1572", "years"))
ageunittable = rbind(ageunittable, c("GSE28422", "years"))
ageunittable = rbind(ageunittable, c("GSE25941", "years"))
ageunittable = rbind(ageunittable, c("GSE53890", "years"))
ageunittable = rbind(ageunittable, c("GSE38718", "years"))
ageunittable = rbind(ageunittable, c("GSE674", "years"))
ageunittable = rbind(ageunittable, c("GSE17612", "years"))
ageunittable = rbind(ageunittable, c("GSE21935", "years"))
ageunittable = rbind(ageunittable, c("GSE362", "years"))
ageunittable = rbind(ageunittable, c("GSE132040", "months"))

colnames(ageunittable) = c("dataset", "unit")
ageunittable = ageunittable %>% column_to_rownames("dataset")

save(ageunittable, file = "ageunittable.RData")
