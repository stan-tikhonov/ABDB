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

ggplot(maintable, aes(x = presence_species)) + geom_histogram()

ggplot(maintable, aes(x = presence_tissues)) + geom_histogram()

ggplot(maintable, aes(x = presence_total)) + geom_histogram()

sum(na.omit(maintable$geommean_pval_tissues)< 0.05)

sum(na.omit(maintable$geommean_pval_species)< 0.05)
