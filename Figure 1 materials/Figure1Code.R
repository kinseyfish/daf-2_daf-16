library(reshape2)
library(ggplot2)
library(edgeR)

daf2_daf16_counts<-read.csv("/Users/kinseyfisher/Documents/daf2_daf16_RNAseq/daf2_daf16_counts.csv",header = T)

WS273_geneNames <- read.csv("/Users/kinseyfisher/Documents/daf2_daf16_RNAseq/WS273_geneNames.csv",header = T) #this has the extra information beyond the wormbase IDs (WB_id)


#Mapping by biotype
daf2_daf16_counts_merge <- merge (daf2_daf16_counts, WS273_geneNames, by.x = "gene_id", by.y = "WB_id") 
daf2_daf16_counts_byBiotype <- daf2_daf16_counts_merge[,-1]
head(daf2_daf16_counts_byBiotype)
daf2_daf16_counts_byBiotype<- daf2_daf16_counts_byBiotype[,-21:-23]
daf2_daf16_counts_byBiotype_melt <- melt(daf2_daf16_counts_byBiotype)

daf2_daf16_byBiotype<- aggregate (daf2_daf16_counts_byBiotype_melt$value,list(variable=daf2_daf16_counts_byBiotype_melt$variable,type=daf2_daf16_counts_byBiotype_melt$type),sum)

Mapping_by_Biotype<-(ggplot(daf2_daf16_byBiotype,aes(x=variable,y=x,fill=type))+
                       geom_bar(stat="identity")+labs(y="Total counts",x="Condition and replicate")+
                       theme_classic(base_size = 15)+ggtitle("Number of reads per biotype")+
                       theme(axis.text.x = element_text(angle = 90, size = 7)) +
                       theme(legend.key.size = unit(.2, "cm")) +
                       theme(aspect.ratio = 1)) 


Mapping_by_Biotype

#Subsetting and prepping the data for DE analysis
daf2_daf16_proteinCoding<-subset(daf2_daf16_counts_merge,type=="protein_coding_gene") 
daf2_daf16_counts2<-daf2_daf16_proteinCoding[,1:21]
rownames(daf2_daf16_counts2)<-daf2_daf16_counts2$gene_id #make the rownames the gene_id column
daf2_daf16_counts2 <- daf2_daf16_counts2[,-1]

counts_groups <- c("N2fed", "N2fed", "N2fed", "N2fed","N2starved", "N2starved", "N2starved", "N2starved", "Daf-2", "Daf-2", "Daf-2", "Daf-2", "Daf16;Daf2", "Daf16;Daf2", "Daf16;Daf2", "Daf16;Daf2", "Daf-16", "Daf-16", "Daf-16", "Daf-16")


#Filter genes that are expressed in >4 samples
d2<-DGEList(counts=daf2_daf16_counts2,group=factor(counts_groups))
dim(d2) #check the dimensions of the d2 without any filtering
keep_filter<-rowSums(cpm(d2)>1)>=4 #decide how you want to filter the data
d2<-d2[keep_filter,] #filter d2 to only include genes that passed filtering
dim(d2) #check the new dimensions
cpm_d2<-cpm(d2,normalized.lib.sizes = TRUE) #make a counts per million object containing normalized CPM
cpm_d2<-as.data.frame(cpm_d2)
cpm_d2_melt<-melt(cpm_d2)


#PCA analysis
conditions<-counts_groups
cpm_d2_df<-data.frame(cpm_d2)
cpm_d2_df$mean<-rowMeans(cpm_d2_df)
cpm_d2_df2<-cpm_d2_df[,1:20]/cpm_d2_df$mean #mean normalize
cpm_d2_df2<-log2(cpm_d2_df2+1) #log2 transform 
pca = prcomp(t(cpm_d2_df2)) #principal component analysis (PCA) on the log2 mean normalized CPM values
summary(pca)
pca_genes<-pca$x
pca_genes_dataframe<-as.data.frame(pca_genes)
pca_genes_dataframe<-data.frame(conditions,pca_genes_dataframe)

replicates<-c("rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3", "rep4","rep1","rep2","rep3", "rep4","rep1","rep2","rep3", "rep4")


PCA1<-(ggplot(pca_genes_dataframe,aes(x=PC1,y=PC2,colour=conditions))+
         geom_point(size=5)+
         ggtitle("PCA of daf2/daf16, 95% CI")+
         labs(x="PC1 (66.5% of variance)",y="PC2 (15.56% of variance)")+
         stat_ellipse(level = 0.95)+theme_classic(base_size = 15)+
         theme(aspect.ratio = 1))
PCA1


#Differential expression analysis
d2<-calcNormFactors(d2) #calculate normalization factors
d2<-estimateCommonDisp(d2) #calculate dispersion 
d2<-estimateTagwiseDisp(d2)

#N2 starved/N2 fed comparison:
de.tag<-exactTest(d2,d2$tagwise.dispersion,pair=c("N2fed","N2starved"))
de.tag_sort<-topTags(de.tag,n=nrow(de.tag$table))$table
de.tag_top<-rownames(de.tag_sort)[de.tag_sort$FDR<=0.05]
de.tag_merge<-merge(WS273_geneNames,de.tag_sort,by.x="WB_id",by.y = 0)
plotSmear(de.tag,de.tags = de.tag_top,main="N2starved/N2fed")
write.table(de.tag_merge,"N2starvedvsN2fed.txt",quote=F,sep="\t")


#daf-16;daf-2/daf-16 comparison:
de.tag<-exactTest(d2,d2$tagwise.dispersion,pair=c("Daf-16","Daf16;Daf2")) #make sure to have these in the right order. The control, or denominator, should come first.
de.tag_sort<-topTags(de.tag,n=nrow(de.tag$table))$table
de.tag_top<-rownames(de.tag_sort)[de.tag_sort$FDR<=0.05]
de.tag_merge<-merge(WS273_geneNames,de.tag_sort,by.x="WB_id",by.y = 0)
plotSmear(de.tag,de.tags = de.tag_top,main="daf-2;daf-16/daf-16")
write.table(de.tag_merge,"daf-16;daf-2vsdaf-16.txt",quote=F,sep="\t") #this writes our .txt file
#Can swap out for different comparisons






#GLMQL
counts <- read.csv("daf2_daf16_counts_starved.csv", header = T, row.names ="gene_id" )
design <- read.csv("glm_design.csv", header = T, row.names = "Sample")

WS273_geneNames <- read.csv("WS273_geneNames.csv",header = T)
counts$gene_id <- row.names(counts)
counts_merge <- merge (counts, WS273_geneNames, by.x = "gene_id", by.y = "WB_id") 
counts_proteinCoding<-subset(counts_merge,type=="protein_coding_gene") 
counts2<-counts_proteinCoding[,1:17]

row.names(counts2) <- counts2$gene_id
counts2 <- counts2[,-1]

dge <- DGEList(counts = counts2, group = design$Genotype)
keep <- filterByExpr(dge, group = design$Genotype)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

designmatrix <- model.matrix(~Genotype, data=design)
colnames(designmatrix) <- levels(design$Genotype)

dge <- estimateDisp(dge, designmatrix)
fit <- glmQLFit(dge, designmatrix)
qlf <- glmQLFTest(fit)
top_genes <- topTags(qlf)
print(top_genes$table)


all_genes <- as.data.frame(qlf$table)
all_genes$FDR <- p.adjust(all_genes$PValue, method = "BH")
write.csv(all_genes, "GLMQL_results.csv")


#Heatmap
library(pheatmap)
library(RColorBrewer)

data <- read.csv("GLM_logFC.csv", row.names = 1, check.names = FALSE)
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(100)

# Convert to matrix
mat <- as.matrix(data)

# Define custom color limits
max_abs <- max(abs(mat), na.rm = TRUE)

# Create more breaks around 0 (nonlinear scale)
# Example: denser near zero using a square root transform
n_breaks <- 101  # must be odd so 0 is included
symmetric_seq <- seq(-sqrt(max_abs), sqrt(max_abs), length.out = n_breaks)
breaks <- sign(symmetric_seq) * symmetric_seq^2  # squash outer values

# Create color palette with same number of colors as breaks - 1
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(length(breaks) - 1)

# Plot heatmap
pheatmap(
  mat = mat,
  color = my_palette,
  breaks = breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  scale = "none",
  fontsize_col = 12
)



#VENN DIAGRAMS
library(biomaRt)

# Step 2: Connect to the Ensembl Mart for WormBase genes
ensembl <- useMart("ensembl", dataset = "celegans_gene_ensembl")

kaplanandtepper <- read.csv("kaplanandtepper.csv", header = T)
# Example gene names (replace this with your actual dataframe column)

kaplanup <- kaplanandtepper$Kaplan.up
kaplandown <- kaplanandtepper$Kaplan.down
classI <- kaplanandtepper$Tepper.class.I
classII <- kaplanandtepper$Tepper.class.II

# Step 3: Map gene names to WormBase IDs
mapped_genes <- getBM(
  attributes = c("external_gene_name", "wormbase_gene"),
  filters = "external_gene_name",
  values = classII,
  mart = ensembl
)

kaplanup_wbIDlist <- as.list(mapped_genes$wormbase_gene)
kaplandown_wbIDlist <- as.list(mapped_genes$wormbase_gene)
classI <- as.list(kaplanandtepper$Tepper.class.I)
classII <- as.list(kaplanandtepper$Tepper.class.II)
double_daf2_down_wb <- as.list(kaplanandtepper$double_daf2_down)
double_daf2_up_wb <- as.list(kaplanandtepper$double_daf2_up)

my_lists <- list(
  double_daf2_down = double_daf2_down_wb,
  double_daf2_up = double_daf2_up_wb,
  classI = classI,
  classII = classII
  #kaplan_down = kaplandown_wbIDlist,
  #kaplan_up = kaplanup_wbIDlist
)


# Plot Venn diagram
ggVennDiagram(my_lists)

#Hypergeometric p-values
phyper(718, 2371, 11951, 1664,lower.tail = FALSE,log.p = FALSE) #down in double/daf-2, Class I
phyper(115, 2282, 12040, 1664,lower.tail = FALSE,log.p = FALSE) #up in double/daf-2, Class I 
phyper(130, 2371, 11951, 1733,lower.tail = FALSE,log.p = FALSE) #down in double/daf-2, Class II
phyper(479, 2282, 12040, 1733,lower.tail = FALSE,log.p = FALSE) #up in double/daf-2, Class II

phyper(396, 2371, 11951, 886,lower.tail = FALSE,log.p = FALSE) #down in double/daf-2, up in Kaplan
phyper(51, 2282, 12040, 886,lower.tail = FALSE,log.p = FALSE) #up in double/daf-2, up in Kaplan
phyper(10, 2371, 11951, 457,lower.tail = FALSE,log.p = FALSE) #down in double/daf-2, down in Kaplan
phyper(334, 2282, 12040, 457,lower.tail = FALSE,log.p = FALSE) #up in double/daf-2, down in Kaplan




#GO-term plot
all <- read.csv("Allcomparisons_logFC_FDR.csv", header = T) #dataframe that has logFC and FDR for each comparison for all genes
#GOTERMS
background <- all$symbol
double_daf2_GO <- subset(all, FDR_double_daf2 < 0.05)
double_daf2_GO_up <- subset(double_daf2_GO, logFC_double_daf2 > 0)
double_daf2_GO_down <- subset(double_daf2_GO, logFC_double_daf2 < 0)

background_list <- as.list(background)
double_daf2_GO_up_list <- as.list(double_daf2_up$symbol)
double_daf2_GO_down_list <- as.list(double_daf2_down$symbol)


library(clusterProfiler)
library(org.Ce.eg.db)
library(enrichplot)
library(GOSemSim)

# Map SYMBOL to ENTREZID
gene_up_entrez <- bitr(double_daf2_GO_up_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Ce.eg.db)$ENTREZID
gene_down_entrez <- bitr(double_daf2_GO_down_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Ce.eg.db)$ENTREZID
background_entrez <- bitr(background_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Ce.eg.db)$ENTREZID


ego_up <- enrichGO(
  gene = gene_up_entrez,
  universe = background_entrez,
  OrgDb = org.Ce.eg.db,
  ont = "BP",
  #pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

ego_down <- enrichGO(
  gene = gene_down_entrez,
  universe = background_entrez,
  OrgDb = org.Ce.eg.db,
  ont = "BP",
  #pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)
combined_go <- rbind(
  ego_up@result %>% dplyr::mutate(direction = "Up"),
  ego_down@result %>% dplyr::mutate(direction = "Down")
)


library(clusterProfiler)

# Simplify UP terms
ego_up_simple <- simplify(ego_up, cutoff = 0.4, by = "p.adjust", select_fun = min)

# Simplify DOWN terms
ego_down_simple <- simplify(ego_down, cutoff = 0.4, by = "p.adjust", select_fun = min)

# Add direction column
ego_up_simple@result$direction <- "Up"
ego_down_simple@result$direction <- "Down"

# Combine results
combined_go <- rbind(ego_up_simple@result, ego_down_simple@result)

# Proceed with semantic similarity matrix, MDS, and plotting (as before)
# Optional: filter to keep terms with adjusted p-value < 0.05 to keep the plot cleaner
combined_go <- combined_go %>% filter(p.adjust < 0.05)

# 4. Prepare GO IDs vector
go_ids <- unique(combined_go$ID)

# 5. Load semantic data for C. elegans biological process (BP)
semdata <- godata('org.Ce.eg.db', ont = "BP")

# 6. Compute semantic similarity matrix using Wang method
sim_matrix <- matrix(NA, nrow = length(go_ids), ncol = length(go_ids),
                     dimnames = list(go_ids, go_ids))

for(i in seq_along(go_ids)) {
  for(j in seq_along(go_ids)) {
    sim_matrix[i,j] <- termSim(go_ids[i], go_ids[j], semData = semdata, method = "Wang")
  }
}

# 7. Perform multidimensional scaling (MDS) to reduce dimensions to 2D
mds_coords <- cmdscale(as.dist(1 - sim_matrix), k = 2)

# 8. Merge MDS coordinates back with combined GO results for plotting
plot_data <- combined_go %>%
  filter(ID %in% rownames(mds_coords)) %>%
  mutate(
    Dim1 = mds_coords[ID, 1],
    Dim2 = mds_coords[ID, 2],
    log10p = -log10(p.adjust)
  )

# 9. Plot with ggplot2

x_range <- range(plot_data$Dim1)
y_range <- range(plot_data$Dim2)

padding_x <- diff(x_range) * 0.1  # 10% padding on x
padding_y <- diff(y_range) * 0.1  # 10% padding on y

ggplot(plot_data, aes(x = Dim1, y = Dim2, color = direction, size = log10p)) +
  geom_point(alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = Description), size = 3, max.overlaps = 100, color = "black") +
  scale_color_manual(values = c("Up" = "yellow", "Down" = "steelblue")) +
  scale_size_continuous(range = c(4, 10))+  # instead of (12, 30)
  coord_fixed(xlim = c(x_range[1] - padding_x, x_range[2] + padding_x),
              ylim = c(y_range[1] - padding_y, y_range[2] + padding_y)) +
  labs(title = "Simplified GO Terms Semantic Similarity Plot",
       x = "Semantic Dimension 1",
       y = "Semantic Dimension 2",
       size = "-log10(adj. p-value)") +
  theme_minimal()

