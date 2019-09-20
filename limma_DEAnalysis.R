rm(list=ls())
graphics.off()


##########################################################
############# Differential Expression Analysis  ##########
##########################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR", version = "3.8")
BiocManager::install("EnhancedVolcano", version = "3.8")
BiocManager::install("biomaRt")
BiocManager::install("KEGGREST", version = "3.8")
require(edgeR)
require(EnhancedVolcano)
require(biomaRt)

X <- read.table("htseq_output.txt", header = F, row.names =1, col.names =c("genes","w0_R1","w0_R2","w0_R3","w1_R1","w1_R2","w1_R3"))
X <- as.matrix(X)
# Cutting last 4 rows of heseq data(X)
X1 <- X[-c((nrow(X)-4):nrow(X)),]

design <- model.matrix(~ 0+factor(c(1,1,1,2,2,2))) #design matrix
colnames(design) <- c("w0","w1") 
contrast.matrix <- makeContrasts(w1-w0, levels = design) #creating contrast matrix

# creating a DGElist object(dge)
dge <- DGEList(counts=X1)
#removing rows with consistant zeros
keep <- filterByExpr(dge, design, min.count = 1)
dge1 <- dge[keep,,keep.lib.sizes=FALSE]
normalised_TMM <- calcNormFactors(dge1)

# Differential Expression plot : voom
v <- voom(normalised_TMM, design, plot = TRUE)

#lima pipeline for differential expression
X_fit <- lmFit(v, design)
X_fit <- contrasts.fit(X_fit, contrast.matrix)
X_fit2 <- eBayes(X_fit)

c_bind <- cbind(X_fit2$coefficients, X_fit2$p.value)
colnames(c_bind) <- c("logFC","p.value")



#############################################################
################# PCA & HCA Analysis ########################
#############################################################

# PCA using factoextra
install.packages("factoextra")
library(factoextra)

#computing PCA
pca <- prcomp(t(log2(1+dge1$counts)), scale = F, center = T, retx = T)


# graph of similar individual PCSs  grouped together
fviz_pca_ind(pca, col.ind = "cos2",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE,title = "PC Analysis of htseq data")

clust  <- hclust(dist(pca$x[,]))
plot(clust, xlab = "PCA Clusters of Replicates", main = "HCA Dendrogram")


#################################################################
############ top 100 differentially expressed genes &    ########
############ volcano plot (top 10 DE genes Higilighted)  ########
#################################################################

#top 100 DEG table
top_100_DEG <- topTable(X_fit2, coef = 1, adjust.method = "BH", p.value = 0.01, lfc = 0.5, sort.by = "logFC", number = 100)
# top 10 DE genes for volcano plot
top10table <- topTable(X_fit2, coef=1, sort.by = "logFC") # picks top 10 by default
# Enhanced volcano plot with top 10 differetially expressed genes highlighted
EnhancedVolcano(c_bind, rownames(c_bind), "logFC", "p.value", selectLab = rownames(top10table), DrawConnectors = T ,title = " Differentially Expressed genes (Top 10 Highlighted)")

################################################################
############# Counting OVER and UNDER expressed genes ##########
################################################################

# counting Overexpresed genes
c("Total number of Overexpressed Genes:", sum(c_bind[,1]>0))
Over_expressed_genes <- sum(c_bind[,1]>0)
# counting underxpressed genes
c("Total number of Underexpressed Genes:", sum(c_bind[,1]<0))
under_expressed_genes <- sum(c_bind[,1]<0)


#############################################################################
################### Annotated list of top 100 DE Genes ######################
#############################################################################

ensembl <- useEnsembl(biomart = "plants_mart", host = "plants.ensembl.org")
#listDatasets(ensembl)

potato_dataset_ensembl <- useMart(biomart = "plants_mart",dataset = "stuberosum_eg_gene", host = "plants.ensembl.org")

#Remove the # from the next line of code to execute it to see the list of attributes
#View(listAttributes(potato_dataset_ensembl, page = "feature_page"))

#getting attributes
biomart_annotation <- getBM(attributes = c("ensembl_gene_id","start_position","end_position","entrezgene_id", "strand","name_1006","description"), mart = potato_dataset_ensembl, filters = "ensembl_gene_id", values = rownames(top_100_DEG)[1:100])
# Biomart Annotation 2 for KEGG Pathways
biomart_annotation2 <- getBM(attributes = "entrezgene_id", mart = potato_dataset_ensembl, filters = "ensembl_gene_id", values = rownames(top_100_DEG)[1:100])
#top100 table with annotations
annot_table <- merge(x=top_100_DEG, y=biomart_annotation, by.x = "row.names", by.y = "ensembl_gene_id", all.x = T, all.y = T)
annot_table <- annot_table[order(abs(annot_table$logFC), decreasing = TRUE),]

# function for collapsing repeating GO Terms
collapse.rows = function(df, sep = ", ") {
  df2 = aggregate(df[, 2:length(df)], list(df[, 1]), function(x) paste0(unique(x), 
                                                                        collapse = sep))
  df2[df2 == "NA"] = NA
  colnames(df2) = colnames(df)
  return(df2)
}


top100.final <- collapse.rows(annot_table, '; ')
top100.final <- top100.final[order(abs(as.numeric(top100.final$logFC)), decreasing = TRUE),]

# writing to CSV file
write.csv(top100.final, file = "potato_annot.csv", row.names = T)


#########################################################
############ Integration of Gene Expression ############¢
############ with KEGG Metabolic Pathways   #############
#########################################################

library(KEGGREST)

omit_NA <- na.omit(biomart_annotation2)
omit_NA <- omit_NA[,1]
kegg_ID <- paste("sot:", omit_NA, sep = "")
kegg_ID <- keggLink("pathway",kegg_ID)

#names <- (names(kegg_ID[1:14]))

# seperating kegg_ID into to due to server limiting it to 10
keggget_1 <- keggGet(kegg_ID[1:10])
keggget_2 <- keggGet(kegg_ID[11:14])
# merging keggget_1 & keggget_2 together
keggget_complete <- c(keggget_1,keggget_2)

#Installing reqired packages "magicfor"
install.packages("magicfor")
library(magicfor)
magic_for(print,silent=T)
# for loop
for(i in 1:14){
  print(keggget_complete[[i]]$NAME[1])
}

kegg_Paths <- magic_result_as_dataframe(iter = FALSE) 

path <- as.matrix(kegg_ID)
Kegg_pathways_results <- cbind(path,kegg_Paths)
Kegg_pathways_results <- data.frame(names = row.names(Kegg_pathways_results),Kegg_pathways_results)
colnames(Kegg_pathways_results) <- c("Gene ID","Path ID", "Pathway")
rownames(Kegg_pathways_results) <- c(1:14)

# writing to .csv file
write.csv(Kegg_pathways_results, "Kegg_pathways_results.csv")

