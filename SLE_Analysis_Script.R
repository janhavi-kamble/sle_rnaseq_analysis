library(readxl)
install.packages("readxl")
library(readxl)
raw_data <- read_xlsx("GSE122459_normalized_genes_all_samples.xlsx")
setwd("C:/Users/Admin/OneDrive/Desktop/SLE_Project")
raw_data <- read_xlsx("GSE122459_normalized_genes_all_samples.xlsx")
install.packages("Rcpp")
"C:/Users/Admin/OneDrive/Desktop/SLE_Project/GSE122459_normalized_genes_all_samples.xlsx"
raw_data <- readxl::read_xlsx("C:/Users/Admin/OneDrive/Desktop/SLE_Project/GSE122459_normalized_genes_all_samples.xlsx")
"C:/Users/Admin/OneDrive/NGS/SLE_RNAseq_Project/GSE122459_normalized_genes_all_samples.xlsx.gz"
if (!require("R.utils")) install.packages("R.utils")
"C:/Users/Admin/OneDrive/NGS/SLE_RNAseq_Project/GSE122459_normalized_genes_all_samples.xlsx"
R.utils::gunzip(zipped_path, destname = unzipped_path, remove = FALSE, overwrite = TRUE)
"C:/Users/Admin/OneDrive/NGS/SLE_RNAseq_Project/GSE122459_normalized_genes_all_samples.xlsx.gz"
if (!require("R.utils")) install.packages("R.utils")
R.utils::gunzip(zipped_path, destname = unzipped_path, remove = FALSE, overwrite = TRUE)
library(readxl)
raw_data <- read_excel(unzipped_path)
# 1. Define file paths
zipped_path <- "C:/Users/Admin/OneDrive/NGS/SLE_RNAseq_Project/GSE122459_normalized_genes_all_samples.xlsx.gz"
unzipped_path <- "C:/Users/Admin/OneDrive/NGS/SLE_RNAseq_Project/GSE122459_normalized_genes_all_samples.xlsx"
# 2. Install and load necessary tools
if (!require("R.utils")) install.packages("R.utils")
if (!require("readxl")) install.packages("readxl")
library(R.utils)
library(readxl)
# 3. Decompress and Read
if (file.exists(zipped_path)) {
gunzip(zipped_path, destname = unzipped_path, remove = FALSE, overwrite = TRUE)
raw_data <- read_excel(unzipped_path)
# 4. Final Check - This will print the dimensions of your data
print(paste("Rows:", nrow(raw_data)))
print(paste("Columns:", ncol(raw_data)))
head(raw_data)
} else {
print("The .gz file was not found at that path. Please check the folder name.")
}
colnames(raw_data)[1] <- "Ensembl_ID"
head(raw_data)
# 1. Remove rows with all zeros (genes that aren't expressed at all)
raw_data <- raw_data[rowSums(raw_data[, 3:28]) > 0, ]
# 2. Split into two groups for comparison
sle_samples <- raw_data[, grep("SLE", colnames(raw_data))]
healthy_samples <- raw_data[, grep("Healthy", colnames(raw_data))]
# 3. Calculate the mean expression for each group
raw_data$mean_SLE <- rowMeans(sle_samples)
raw_data$mean_Healthy <- rowMeans(healthy_samples)
# 4. View the result
head(raw_data[, c("Ensembl_ID", "name", "mean_SLE", "mean_Healthy")])
# 1. Calculate Log2 Fold Change
# We add 1 to avoid dividing by zero or taking log of zero (log-shifting)
raw_data$log2FC <- log2( (raw_data$mean_SLE + 1) / (raw_data$mean_Healthy + 1) )
# 2. Filter for Highly Up-regulated genes in SLE (Fold change > 1)
up_regulated <- raw_data[raw_data$log2FC > 1, ]
# 3. Filter for Highly Down-regulated genes in SLE (Fold change < -1)
down_regulated <- raw_data[raw_data$log2FC < -1, ]
# 4. See how many genes are significantly different
print(paste("Up-regulated genes:", nrow(up_regulated)))
print(paste("Down-regulated genes:", nrow(down_regulated)))
# 5. Look at the top up-regulated genes
head(up_regulated[order(-up_regulated$log2FC), c("name", "log2FC")])
# 1. Install and load ggplot2 for plotting
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)
# 2. Create the Volcano Plot
ggplot(raw_data, aes(x = log2FC, y = mean_SLE)) +
geom_point(aes(color = abs(log2FC) > 1), alpha = 0.5) +
scale_color_manual(values = c("black", "red")) +
theme_minimal() +
labs(title = "Gene Expression: SLE vs Healthy",
x = "Log2 Fold Change",
y = "Mean Expression (SLE)",
color = "Significant (|L2FC| > 1)") +
geom_text(data = subset(raw_data, log2FC > 5),
aes(label = name), vjust = 1, hjust = 1, check_overlap = TRUE)
library(ggplot2)
# Create the plot with threshold lines and labels
ggplot(raw_data, aes(x = log2FC, y = log10(mean_SLE + 1))) +
geom_point(aes(color = abs(log2FC) > 1), alpha = 0.4) +
geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
scale_color_manual(values = c("grey", "red")) +
theme_classic() +
labs(title = "SLE Gene Expression Analysis",
subtitle = "Red points indicate |Log2FC| > 1",
x = "Log2 Fold Change (SLE / Healthy)",
y = "Log10 Mean Expression (SLE)",
color = "Fold Change > 2") +
geom_text(data = subset(raw_data, log2FC > 4),
aes(label = name),
nudge_y = 0.1,
check_overlap = TRUE,
color = "darkred")
if (!require("pheatmap")) install.packages("pheatmap")
library(pheatmap)
# 1. Force install from a reliable source
install.packages("pheatmap", repos = "http://cran.us.r-project.org")
# 2. Try loading it again
library(pheatmap)
# 3. If the library loaded, run the heatmap code
top_genes <- up_regulated[order(-up_regulated$log2FC), ][1:30, ]
plot_data <- as.matrix(top_genes[, 3:28])
rownames(plot_data) <- top_genes$name
pheatmap(plot_data,
scale = "row",
clustering_distance_rows = "correlation",
main = "Expression Heatmap: Top 30 SLE Genes",
color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
# 1. Ensure ggplot2 is loaded
library(ggplot2)
# 2. Add a 'Significance' column for coloring the dots
raw_data$diffexpressed <- "No"
raw_data$diffexpressed[raw_data$log2FC > 1] <- "Up-regulated"
raw_data$diffexpressed[raw_data$log2FC < -1] <- "Down-regulated"
# 3. Create the Volcano Plot
ggplot(data = raw_data, aes(x = log2FC, y = log10(mean_SLE + 1), color = diffexpressed)) +
geom_point(alpha = 0.4, size = 1.5) +
theme_minimal() +
# Add vertical and horizontal lines for thresholds
geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") +
# Custom colors: Red for Up, Blue for Down, Grey for No Change
scale_color_manual(values = c("Down-regulated" = "blue", "No" = "grey", "Up-regulated" = "red")) +
# Label the top genes (adjusting 'log2FC > 4' to catch your specific high-flyers)
geom_text(data = subset(raw_data, log2FC > 4),
aes(label = name),
size = 3,
vjust = -1,
check_overlap = TRUE,
show.legend = FALSE) +
# Labels and Titles
labs(title = "Volcano Plot: SLE vs Healthy Controls",
subtitle = "Highlighting genes with |Log2FC| > 1",
x = "Log2 Fold Change (Effect Size)",
y = "Log10 Mean Expression (Magnitude)",
color = "Expression Status")
ggsave("SLE_Volcano_Plot.png", width = 8, height = 6, dpi = 300)
# Save the current plot in your project folder
ggsave("SLE_Volcano_Plot.png", width = 8, height = 6, dpi = 300)
# 1. Open a new PNG file
png("SLE_Heatmap.png", width = 1200, height = 1500, res = 150)
# 2. Run the heatmap command (this "draws" into the PNG)
pheatmap(plot_data,
scale = "row",
annotation_col = annotation_col,
annotation_colors = ann_colors,
clustering_distance_rows = "correlation",
main = "Annotated Heatmap: SLE vs Healthy",
color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
fontsize_row = 8)
# 1. Redefine the annotation data (Functional requirement for the plot)
annotation_col <- data.frame(
Group = factor(c(rep("SLE", 20), rep("Healthy", 6)))
)
rownames(annotation_col) <- colnames(plot_data)
# 2. Redefine the colors
ann_colors = list(
Group = c(SLE = "firebrick", Healthy = "navy")
)
# 3. Open the PNG file
png("SLE_Heatmap.png", width = 1200, height = 1500, res = 150)
# 4. Run the heatmap command
pheatmap(plot_data,
scale = "row",
annotation_col = annotation_col,
annotation_colors = ann_colors,
clustering_distance_rows = "correlation",
main = "Annotated Heatmap: SLE vs Healthy",
color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
fontsize_row = 8)
# 5. Close and save
dev.off()
# Save everything to a file named 'SLE_Analysis_Session.RData'
save.image(file = "SLE_Analysis_Session.RData")
# Save the full table with all your new calculations
write.csv(raw_data, "SLE_Full_Results_Analysis.csv", row.names = FALSE)
# 1. Run this multiple times until you see "Error: device 1 (the null device)"
# This ensures all '0 byte' stuck files are released.
dev.off()
# 2. Re-run the save command
png("SLE_Heatmap.png", width = 1200, height = 1500, res = 150)
pheatmap(plot_data,
scale = "row",
annotation_col = annotation_col,
annotation_colors = ann_colors,
clustering_distance_rows = "correlation",
main = "Annotated Heatmap: SLE vs Healthy",
color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
fontsize_row = 8)
# 3. This is the most important line - do not skip it!
dev.off()
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
library(clusterProfiler)
library(org.Hs.eg.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
# Extract gene symbols from your top_genes object
up_genes <- row.names(top_genes)
# Run the GO Enrichment (Biological Process category)
ego <- enrichGO(gene          = up_genes,
OrgDb         = org.Hs.eg.db,
keyType       = 'SYMBOL',
ont           = "BP",
pAdjustMethod = "BH",
pvalueCutoff  = 0.05,
readable      = TRUE)
# Remove NA values and empty strings
up_genes_clean <- up_genes[!is.na(up_genes) & up_genes != ""]
# Show the first 10 genes to double-check their format
head(up_genes_clean)
# This checks how many of your genes match the database
valid_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
sum(up_genes_clean %in% valid_symbols)
ego <- enrichGO(gene          = up_genes_clean,
OrgDb         = org.Hs.eg.db,
keyType       = 'SYMBOL',
ont           = "BP",
pAdjustMethod = "BH",
pvalueCutoff  = 0.05,
readable      = TRUE)
head(up_genes)
table(is.na(up_genes))
print(up_genes)
colnames(top_genes)
# Pull the actual gene names from the 'name' column
up_genes <- top_genes$name
# Quick check: this should now show names like IFI27, not numbers
head(up_genes)
ego <- enrichGO(gene          = up_genes,
OrgDb         = org.Hs.eg.db,
keyType       = 'SYMBOL',
ont           = "BP",
pAdjustMethod = "BH",
pvalueCutoff  = 0.05,
readable      = TRUE)
# View the biological insights
head(ego)
library(enrichplot)
# Create the Dotplot
final_dotplot <- dotplot(ego, showCategory=10) +
ggtitle("Biological Pathway Enrichment: SLE vs Healthy") +
theme(plot.title = element_text(hjust = 0.5))
# Save the final piece of your portfolio
ggsave("SLE_Enrichment_Final.png", plot = final_dotplot, width = 10, height = 7, dpi = 300)
# Display the plot
print(final_dotplot)
