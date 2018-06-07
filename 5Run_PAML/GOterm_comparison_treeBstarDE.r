# ch. 2 GO term comparison from codeml PSGs
# does not include genes that could not be annotated (just have zebra finch Ensembl IDs for those)

library(GOplot)

setwd('C:/Users/mcwlim/Dropbox/Marisacompfiles/Transcriptome files/Jan2018_codemlreanalysis')

# ---- Create GO dfs ------
tree_go_terms <- read.csv('Pantherdb_GO_alltrees.csv', header=TRUE)
#head(tree_go_terms)

tree_B <- tree_go_terms[tree_go_terms$Tree == 'B' & tree_go_terms$Use.this.result == 'Yes',]
tree_star <- tree_go_terms[tree_go_terms$Tree == 'star' & tree_go_terms$Use.this.result == 'Yes',]
tree_D <- tree_go_terms[tree_go_terms$Tree == 'D' & tree_go_terms$Use.this.result == 'Yes',]
tree_E <- tree_go_terms[tree_go_terms$Tree == 'E' & tree_go_terms$Use.this.result == 'Yes',]

# Biological process
tree_B_BP <- tree_B[, c(1,3,10)]
tree_B_BP1 <- strsplit(as.character(tree_B_BP$GO.database.BP.complete), ";")
tree_B_BP2 <- data.frame('Tree'=rep(tree_B_BP$Tree, sapply(tree_B_BP1, length)), 
                         'Mapped.IDs'=rep(tree_B_BP$Mapped.IDs, sapply(tree_B_BP1, length)),
                         'GO_term'=unlist(tree_B_BP1))
tree_B_BP2$GO_cat <- 'BP'

tree_star_BP <- tree_star[, c(1,3,10)]
tree_star_BP1  <- strsplit(as.character(tree_star_BP$GO.database.BP.complete), ';')
tree_star_BP2 <- data.frame('Tree'=rep(tree_star_BP$Tree, sapply(tree_star_BP1, length)),
                            'Mapped.IDs'=rep(tree_star_BP$Mapped.IDs, sapply(tree_star_BP1, length)),
                            'GO_term'=unlist(tree_star_BP1))
tree_star_BP2$GO_cat <- 'BP'

tree_D_BP <- tree_D[, c(1,3,10)]
tree_D_BP1  <- strsplit(as.character(tree_D_BP$GO.database.BP.complete), ';')
tree_D_BP2 <- data.frame('Tree'=rep(tree_D_BP$Tree, sapply(tree_D_BP1, length)),
                            'Mapped.IDs'=rep(tree_D_BP$Mapped.IDs, sapply(tree_D_BP1, length)),
                            'GO_term'=unlist(tree_D_BP1))
tree_D_BP2$GO_cat <- 'BP'

tree_E_BP <- tree_E[, c(1,3,10)]
tree_E_BP1  <- strsplit(as.character(tree_E_BP$GO.database.BP.complete), ';')
tree_E_BP2 <- data.frame('Tree'=rep(tree_E_BP$Tree, sapply(tree_E_BP1, length)),
                            'Mapped.IDs'=rep(tree_E_BP$Mapped.IDs, sapply(tree_E_BP1, length)),
                            'GO_term'=unlist(tree_E_BP1))
tree_E_BP2$GO_cat <- 'BP'

# Molecular function
tree_B_MF <- tree_B[, c(1,3, 9)]
tree_B_MF1  <- strsplit(as.character(tree_B_MF$GO.database.MF.complete), ';')
tree_B_MF2 <- data.frame('Tree'=rep(tree_B_MF$Tree, sapply(tree_B_MF1, length)),
                            'Mapped.IDs'=rep(tree_B_MF$Mapped.IDs, sapply(tree_B_MF1, length)),
                            'GO_term'=unlist(tree_B_MF1))
tree_B_MF2$GO_cat <- 'MF'

tree_star_MF <- tree_star[, c(1,3, 9)]
tree_star_MF1  <- strsplit(as.character(tree_star_MF$GO.database.MF.complete), ';')
tree_star_MF2 <- data.frame('Tree'=rep(tree_star_MF$Tree, sapply(tree_star_MF1, length)),
                         'Mapped.IDs'=rep(tree_star_MF$Mapped.IDs, sapply(tree_star_MF1, length)),
                         'GO_term'=unlist(tree_star_MF1))
tree_star_MF2$GO_cat <- 'MF'

tree_D_MF <- tree_D[, c(1,3, 9)]
tree_D_MF1  <- strsplit(as.character(tree_D_MF$GO.database.MF.complete), ';')
tree_D_MF2 <- data.frame('Tree'=rep(tree_D_MF$Tree, sapply(tree_D_MF1, length)),
                         'Mapped.IDs'=rep(tree_D_MF$Mapped.IDs, sapply(tree_D_MF1, length)),
                         'GO_term'=unlist(tree_D_MF1))
tree_D_MF2$GO_cat <- 'MF'

tree_E_MF <- tree_E[, c(1,3, 9)]
tree_E_MF1  <- strsplit(as.character(tree_E_MF$GO.database.MF.complete), ';')
tree_E_MF2 <- data.frame('Tree'=rep(tree_E_MF$Tree, sapply(tree_E_MF1, length)),
                         'Mapped.IDs'=rep(tree_E_MF$Mapped.IDs, sapply(tree_E_MF1, length)),
                         'GO_term'=unlist(tree_E_MF1))
tree_E_MF2$GO_cat <- 'MF'

# Cellular component
tree_B_CC <- tree_B[, c(1,3,11)]
tree_B_CC1  <- strsplit(as.character(tree_B_CC$GO.database.CC.complete), ';')
tree_B_CC2 <- data.frame('Tree'=rep(tree_B_CC$Tree, sapply(tree_B_CC1, length)),
                         'Mapped.IDs'=rep(tree_B_CC$Mapped.IDs, sapply(tree_B_CC1, length)),
                         'GO_term'=unlist(tree_B_CC1))
tree_B_CC2$GO_cat <- 'CC'

tree_star_CC <- tree_star[, c(1,3,11)]
tree_star_CC1  <- strsplit(as.character(tree_star_CC$GO.database.CC.complete), ';')
tree_star_CC2 <- data.frame('Tree'=rep(tree_star_CC$Tree, sapply(tree_star_CC1, length)),
                         'Mapped.IDs'=rep(tree_star_CC$Mapped.IDs, sapply(tree_star_CC1, length)),
                         'GO_term'=unlist(tree_star_CC1))
tree_star_CC2$GO_cat <- 'CC'

tree_D_CC <- tree_D[, c(1,3,11)]
tree_D_CC1  <- strsplit(as.character(tree_D_CC$GO.database.CC.complete), ';')
tree_D_CC2 <- data.frame('Tree'=rep(tree_D_CC$Tree, sapply(tree_D_CC1, length)),
                         'Mapped.IDs'=rep(tree_D_CC$Mapped.IDs, sapply(tree_D_CC1, length)),
                         'GO_term'=unlist(tree_D_CC1))
tree_D_CC2$GO_cat <- 'CC'

tree_E_CC <- tree_E[, c(1,3,11)]
tree_E_CC1  <- strsplit(as.character(tree_E_CC$GO.database.CC.complete), ';')
tree_E_CC2 <- data.frame('Tree'=rep(tree_E_CC$Tree, sapply(tree_E_CC1, length)),
                         'Mapped.IDs'=rep(tree_E_CC$Mapped.IDs, sapply(tree_E_CC1, length)),
                         'GO_term'=unlist(tree_E_CC1))
tree_E_CC2$GO_cat <- 'CC'

# ---- Analysis 1: compare GO terms within Tree_star PSGs -----
star_bp_overlaps <- as.vector(tree_star_BP2[c(which(duplicated(tree_star_BP2$GO_term))),]$GO_term)
tree_star_BP2[tree_star_BP2$GO_term %in% star_bp_overlaps,]
length(star_bp_overlaps)
length(unique(tree_star_BP2[tree_star_BP2$GO_term %in% star_bp_overlaps,]$Mapped.IDs))

star_mf_overlaps <- as.vector(tree_star_MF2[c(which(duplicated(tree_star_MF2$GO_term))),]$GO_term)
tree_star_MF2[tree_star_MF2$GO_term %in% star_mf_overlaps,]
length(star_mf_overlaps)
length(unique(tree_star_MF2[tree_star_MF2$GO_term %in% star_mf_overlaps,]$Mapped.IDs))

star_cc_overlaps <- as.vector(tree_star_CC2[c(which(duplicated(tree_star_CC2$GO_term))),]$GO_term)
tree_star_CC2[tree_star_CC2$GO_term %in% star_cc_overlaps,]
length(star_cc_overlaps)
length(unique(tree_star_CC2[tree_star_CC2$GO_term %in% star_cc_overlaps,]$Mapped.IDs))

# ---- Analysis 1 vs. 2: compare GO terms Tree_star vs. Tree_E PSGs -----
A2_bp_overlap <- intersect(tree_star_BP2$GO_term, tree_E_BP2$GO_term)
tree_star_BP2[tree_star_BP2$GO_term %in% A2_bp_overlap,]
tree_E_BP2[tree_E_BP2$GO_term %in% A2_bp_overlap,]

A2_mf_overlap <- intersect(tree_star_MF2$GO_term, tree_E_MF2$GO_term)
tree_star_MF2[tree_star_MF2$GO_term %in% A2_mf_overlap,]
tree_E_MF2[tree_E_MF2$GO_term %in% A2_mf_overlap,]

A2_cc_overlap <- intersect(tree_star_CC2$GO_term, tree_E_CC2$GO_term)
tree_star_CC2[tree_star_CC2$GO_term %in% A2_cc_overlap,]
tree_E_CC2[tree_E_CC2$GO_term %in% A2_cc_overlap,]

# ---- Analysis 1 vs. 1a: compare GO terms Tree_star vs. Tree_B PSGs ------
A1a_bp_overlap <- intersect(tree_star_BP2$GO_term, tree_B_BP2$GO_term)
tree_star_BP2[tree_star_BP2$GO_term %in% A1a_bp_overlap,]
tree_B_BP2[tree_B_BP2$GO_term %in% A1a_bp_overlap,]

A1a_mf_overlap <- intersect(tree_star_MF2$GO_term, tree_B_MF2$GO_term)
tree_star_MF2[tree_star_MF2$GO_term %in% A1a_mf_overlap,]
tree_B_MF2[tree_B_MF2$GO_term %in% A1a_mf_overlap,]

A1a_cc_overlap <- intersect(tree_star_CC2$GO_term, tree_B_CC2$GO_term)
tree_star_CC2[tree_star_CC2$GO_term %in% A1a_cc_overlap,]
tree_B_CC2[tree_B_CC2$GO_term %in% A1a_cc_overlap,]

# ---- Analysis 2 vs. 2a: compare GO terms Tree_E vs. Tree_D PSGs ------
A2a_bp_overlap <- intersect(tree_D_BP2$GO_term, tree_E_BP2$GO_term)
tree_D_BP2[tree_D_BP2$GO_term %in% A2a_bp_overlap,]
tree_E_BP2[tree_E_BP2$GO_term %in% A2a_bp_overlap,]

A2a_mf_overlap <- intersect(tree_D_MF2$GO_term, tree_E_MF2$GO_term)
tree_D_MF2[tree_D_MF2$GO_term %in% A2a_mf_overlap,]
tree_E_MF2[tree_E_MF2$GO_term %in% A2a_mf_overlap,]

A2ab_cc_overlap <- intersect(tree_D_CC2$GO_term, tree_E_CC2$GO_term)
tree_D_CC2[tree_D_CC2$GO_term %in% A2ab_cc_overlap,]
tree_E_CC2[tree_E_CC2$GO_term %in% A2ab_cc_overlap,]

# ---- Save tables ----
BP_df <- rbind(tree_B_BP2, tree_star_BP2, tree_D_BP2, tree_E_BP2)
MF_df <- rbind(tree_B_MF2, tree_star_MF2, tree_D_MF2, tree_E_MF2)
CC_df <- rbind(tree_B_CC2, tree_star_CC2, tree_D_CC2, tree_E_CC2)

write.table(BP_df, 'GO_BP_treeBstarDE.txt', row.names=FALSE, quote=FALSE, sep='\t')
write.table(MF_df, 'GO_MF_treeBstarDE.txt', row.names=FALSE, quote=FALSE, sep='\t')
write.table(CC_df, 'GO_CC_treeBstarDE.txt', row.names=FALSE, quote=FALSE, sep='\t')

# ---- Plots -----
A2_df <- rbind(tree_star_BP2, tree_E_BP2, tree_star_MF2, tree_E_MF2, tree_star_CC2, tree_E_CC2)
A1a_df <- rbind(tree_star_BP2, tree_B_BP2, tree_star_MF2, tree_B_MF2, tree_star_CC2, tree_B_CC2)
A2a_df <- rbind(tree_E_BP2, tree_D_BP2, tree_E_MF2, tree_D_MF2, tree_E_CC2, tree_D_CC2)


GOHeat(table(A2a_df[A2_df$GO_cat=='MF',]$Mapped.IDs, 
             A2a_df[A2_df$GO_cat=='MF',]$GO_term), nlfc=0)


  
