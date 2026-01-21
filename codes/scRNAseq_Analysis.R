################################################
# scRNA-seq
################################################

# ---- Setup ----
.libPaths("~/Rlibs")

# ================================= LIBRARIES ==================================

library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(openxlsx)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(patchwork)
library(gprofiler2)
library(scales)

library(future)
plan("multicore", workers = 8, future.seed = TRUE)
options(future.globals.maxSize = 50 * 1024^3)

rm(list = ls())

# ================================ PARAMETERS ==================================

input_file = "SeuratFinal.rds"
output_folder = "Results/"
image_folder = "result_images/"

# Create directories
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(image_folder, showWarnings = FALSE, recursive = TRUE)

# Define color scheme
CTL_COLOR <- "#4DAF4A"      # Green for Control
PINK1_COLOR <- "#E91E8C"    # Pink for PINK1(PD-mutant)
DAY_COLORS <- viridis(5, option = "plasma")  # For differentiation days

# Color palettes for different visualizations
CONDITION_COLORS <- c("CTL" = CTL_COLOR, "PINK1" = PINK1_COLOR)
DIRECTION_COLORS <- c("Up in CTL" = CTL_COLOR, "Up in PINK1" = PINK1_COLOR)

# ================================ LOAD DATA ===================================

cat("Loading data...\n")
S = readRDS(input_file)

# Metadata setup
S$Day = factor(S$Day, levels = c(0,18,25,37,57))
S$Condition = str_replace_all(S$Condition, "ND","PINK1")
S$Condition = str_replace_all(S$Condition,"WT","CTL")
S$ConditionDay = paste(S$Condition,S$Day,sep = " ")
S$ConditionDay = factor(S$ConditionDay, levels = c("CTL 0","CTL 18","CTL 25","CTL 37","CTL 57",
                                                   "PINK1 0","PINK1 18","PINK1 25","PINK1 37","PINK1 57"))
S$Condition = factor(S$Condition, levels = c("CTL","PINK1"))

# Normalize
S = NormalizeData(S, assay = "RNA", verbose = FALSE)
S = NormalizeData(S, assay = "RNA_imputed", verbose = FALSE)

# PLOT 1: SPLIT UMAP - CTL vs PINK1 SIDE-BY-SIDE

cat("\n=== Creating Plot 1: Split UMAP ===\n")

DefaultAssay(S) = "integrated"
S = ScaleData(S, verbose = FALSE)
S = RunPCA(S, npcs = 50, verbose = FALSE)
set.seed(18)
S = RunUMAP(S, dims = 1:50, assay = "integrated", n.neighbors = 50, seed.use = 18, verbose = FALSE)

# Create plot
p1 <- DimPlot(S, reduction = "umap", group.by = "Day", split.by = "Condition", ncol = 2) +
  scale_color_manual(values = DAY_COLORS) +
  labs(title = "UMAP: Differentiation Trajectory - CTL vs PINK1",
       x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

ggsave(filename = paste0(image_folder, "plot1_split_umap.png"), 
       plot = p1, width = 12, height = 6, dpi = 300, bg = "white")

# PLOT 2: CELL COMPOSITION BAR PLOT

cat("\n=== Creating Plot 2: Cell Composition ===\n")

cell_counts <- S[[]] %>%
  group_by(Condition, Day) %>%
  summarise(n_cells = n(), .groups = 'drop')

p2 <- ggplot(cell_counts, aes(x = Day, y = n_cells, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = n_cells), position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = CONDITION_COLORS) +
  theme_minimal() +
  labs(title = "Cell Numbers Across Differentiation",
       subtitle = "Quantitative comparison of cell representation",
       x = "Differentiation Day", y = "Number of Cells",
       fill = "Genotype") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "top",
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

ggsave(filename = paste0(image_folder, "plot2_cell_composition.png"), 
       plot = p2, width = 10, height = 6, dpi = 300, bg = "white")

# RUN DEG ANALYSIS (for plots 3-4)

cat("\n=== Running DEG Analysis ===\n")

DefaultAssay(S) = "RNA"
S = SetIdent(object = S, value = S$ConditionDay)

for(pink in c("PINK1")){
  for(day in c(0,18,25,37,57)) {
    cat("Analyzing Day", day, "...\n")
    temp = FindMarkers(S, 
                       assay = "RNA", 
                       ident.1 = paste("CTL",day,sep=" "), 
                       ident.2 = paste(pink,day,sep=" "),
                       logfc.threshold = 0.25,
                       min.pct = 0.05,
                       densify = TRUE,
                       test.use = "wilcox",
                       layer = "data",
                       verbose = FALSE)
    temp$Day = day
    temp$Gene = rownames(temp)
    temp$ident.1 = "CTL"
    temp$ident.2 = pink
    rownames(temp) = paste(pink,day,c(1:dim(temp)[1]),sep="_")
    if(day==0 & pink == "PINK1"){
      DEG = temp
    }else{
      DEG = rbind(DEG,temp)
    }
  }
}

DEG = DEG %>% 
  select(Gene,Day,pct.1,pct.2,avg_log2FC,p_val,p_val_adj,ident.1,ident.2)  %>%
  arrange(Day,Gene,ident.2)

write.xlsx(DEG, paste0(output_folder,"DEG_complete.xlsx"), colNames = TRUE, rowNames = FALSE)

# PLOT 3: DEG TEMPORAL DYNAMICS

cat("\n=== Creating Plot 3: DEG Temporal Dynamics ===\n")

deg_summary <- DEG %>%
  mutate(Direction = case_when(
    p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "Up in CTL",
    p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Up in PINK1",
    TRUE ~ "Not Significant"
  )) %>%
  filter(Direction != "Not Significant") %>%
  group_by(Day, Direction) %>%
  summarise(Count = n(), .groups = 'drop')

p3 <- ggplot(deg_summary, aes(x = factor(Day), y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = DIRECTION_COLORS) +
  theme_minimal() +
  labs(title = "Number of Differentially Expressed Genes Over Time",
       subtitle = "Progressive divergence between CTL and PINK1 (padj < 0.05, |FC| > 0.5)",
       x = "Differentiation Day", y = "Number of DEGs",
       fill = "Direction") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "top",
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

ggsave(filename = paste0(image_folder, "plot3_deg_temporal_dynamics.png"), 
       plot = p3, width = 10, height = 6, dpi = 300, bg = "white")


# PLOT 4: GENE EXPRESSION TRAJECTORIES

cat("\n=== Creating Plot 4: Gene Expression Trajectories ===\n")

key_genes <- c("TH", "DDC", "NR4A2", "POU5F1", "SOX2", "LMX1A")
key_genes_present <- key_genes[key_genes %in% rownames(S)]

DefaultAssay(S) = "RNA_imputed"

gene_trajectories <- lapply(key_genes_present, function(gene) {
  exp_data <- GetAssayData(S, assay = "RNA_imputed", layer = "data")[gene, ]
  
  trajectory <- S[[]] %>%
    mutate(Expression = exp_data) %>%
    group_by(Condition, Day) %>%
    summarise(Mean_Exp = mean(Expression),
              SE = sd(Expression) / sqrt(n()),
              .groups = 'drop') %>%
    mutate(Gene = gene)
  
  return(trajectory)
})

gene_trajectories_df <- bind_rows(gene_trajectories)

p4 <- ggplot(gene_trajectories_df, 
             aes(x = as.numeric(as.character(Day)), y = Mean_Exp, 
                 color = Condition, group = Condition)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean_Exp - SE, ymax = Mean_Exp + SE), 
                width = 2, alpha = 0.5, linewidth = 0.8) +
  facet_wrap(~Gene, scales = "free_y", ncol = 3) +
  scale_color_manual(values = CONDITION_COLORS) +
  theme_minimal() +
  labs(title = "Gene Expression Trajectories During Differentiation",
       subtitle = "Mean ± SE across cells",
       x = "Differentiation Day", y = "Average Expression",
       color = "Genotype") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "top",
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 11))

ggsave(filename = paste0(image_folder, "plot4_gene_trajectories.png"), 
       plot = p4, width = 12, height = 8, dpi = 300, bg = "white")

# PLOT 5: GO ENRICHMENT - MITOCHONDRIAL PATHWAYS

cat("\n=== Creating Plot 5: GO Enrichment Analysis ===\n")

go_results_list <- list()

for(day in c(0,18,25,37,57)) {
  
  deg_day <- DEG %>% filter(Day == day, p_val_adj < 0.05)
  
  if(nrow(deg_day) < 10) next
  
  # Genes up in CTL
  genes_up_ctl <- deg_day %>% filter(avg_log2FC > 0.5) %>% pull(Gene)
  
  if(length(genes_up_ctl) >= 5) {
    gost_ctl <- gost(genes_up_ctl, 
                     organism = "hsapiens",
                     ordered_query = FALSE,
                     significant = TRUE,
                     user_threshold = 0.05,
                     correction_method = "fdr",
                     sources = c("GO:BP", "GO:CC"))
    
    if(!is.null(gost_ctl$result)) {
      gost_ctl$result$Day <- day
      gost_ctl$result$Direction <- "CTL"
      go_results_list[[paste0("Day", day, "_CTL")]] <- gost_ctl$result
    }
  }
  
  # Genes up in PINK1
  genes_up_pink1 <- deg_day %>% filter(avg_log2FC < -0.5) %>% pull(Gene)
  
  if(length(genes_up_pink1) >= 5) {
    gost_pink1 <- gost(genes_up_pink1, 
                       organism = "hsapiens",
                       ordered_query = FALSE,
                       significant = TRUE,
                       user_threshold = 0.05,
                       correction_method = "fdr",
                       sources = c("GO:BP", "GO:CC"))
    
    if(!is.null(gost_pink1$result)) {
      gost_pink1$result$Day <- day
      gost_pink1$result$Direction <- "PINK1"
      go_results_list[[paste0("Day", day, "_PINK1")]] <- gost_pink1$result
    }
  }
}

if(length(go_results_list) > 0) {
  GO_combined <- bind_rows(go_results_list)
  write.xlsx(GO_combined, paste0(output_folder, "GO_enrichment_daywise.xlsx"))
  
  # Plot mitochondrial pathways
  mito_keywords <- c("mitochondri", "oxidative", "respiratory", "electron transport", "ATP")
  
  go_mito <- GO_combined %>%
    filter(grepl(paste(mito_keywords, collapse = "|"), term_name, ignore.case = TRUE)) %>%
    group_by(Day, Direction) %>%
    arrange(p_value) %>%
    slice_head(n = 3) %>%
    ungroup()
  
  if(nrow(go_mito) > 0) {
    p5 <- ggplot(go_mito, 
                 aes(x = factor(Day), y = term_name, 
                     size = intersection_size, color = -log10(p_value))) +
      geom_point(alpha = 0.8) +
      scale_color_viridis(option = "plasma") +
      scale_size_continuous(range = c(3, 10)) +
      facet_wrap(~Direction, scales = "free_y", ncol = 1) +
      theme_minimal() +
      labs(title = "Mitochondrial Pathway Enrichment Across Differentiation",
           subtitle = "Top mitochondrial-related GO terms per day and condition",
           x = "Differentiation Day", y = "GO Term",
           size = "Gene Count", color = "-log10(p-value)") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 10),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 10),
            strip.text = element_text(face = "bold", size = 12),
            legend.position = "right")
    
    ggsave(filename = paste0(image_folder, "plot5_go_mitochondrial.png"), 
           plot = p5, width = 12, height = 10, dpi = 300, bg = "white")
  }
}

# PLOT 6: CLUSTERING UMAP

cat("\n=== Creating Plot 6: Clustering UMAP ===\n")

DefaultAssay(S) = "integrated"

# Find neighbors and clusters using the PCA
S = FindNeighbors(S, reduction = "pca", dims = 1:20, verbose = FALSE)
S = FindClusters(S, resolution = 0.4, verbose = FALSE)

# Define colors for clusters
cluster_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(S$seurat_clusters)))

p6 <- DimPlot(S, reduction = "umap", group.by = "seurat_clusters", 
              label = TRUE, label.size = 5, repel = TRUE) +
  scale_color_manual(values = cluster_colors) +
  labs(title = "UMAP: Cell Clusters (Resolution 0.4)",
       subtitle = paste0(length(unique(S$seurat_clusters)), " clusters identified"),
       x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

ggsave(filename = paste0(image_folder, "plot6_clustering_umap.png"), 
       plot = p6, width = 10, height = 8, dpi = 300, bg = "white")

# PLOT 7: CLUSTER MARKERS HEATMAP

cat("\n=== Creating Plot 7: Cluster Markers Heatmap ===\n")

# Find markers for each cluster
DefaultAssay(S) = "RNA"
Idents(S) = "seurat_clusters"

cluster_markers <- FindAllMarkers(S, 
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25,
                                  verbose = FALSE,
                                  random.seed = 42)

write.xlsx(cluster_markers, paste0(output_folder, "Cluster_markers.xlsx"))

# Get top 10 markers per cluster
top10_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  ungroup()

if(nrow(top10_markers) > 0) {
  S = ScaleData(S, features = unique(top10_markers$gene), verbose = FALSE)
  
  p7 <- DoHeatmap(S, features = unique(top10_markers$gene),
                  group.by = "seurat_clusters",
                  assay = "RNA", slot = "scale.data",
                  size = 3) +
    scale_fill_viridis(option = "viridis") +
    labs(title = "Top 10 Marker Genes per Cluster") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text.y = element_text(size = 6))
  
  ggsave(filename = paste0(image_folder, "plot7_cluster_heatmap.png"), 
         plot = p7, width = 14, height = 10, dpi = 300, bg = "white")
}

# PLOT 8: CLUSTER MARKERS DOT PLOT

cat("\n=== Creating Plot 8: Cluster Markers Dot Plot ===\n")

# Get top 5 markers per cluster
top5_markers <- cluster_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 5) %>%
  ungroup()

if(nrow(top5_markers) > 0) {
  top_genes <- unique(top5_markers$gene)
  
  p8 <- DotPlot(S, features = top_genes, group.by = "seurat_clusters",
                assay = "RNA", dot.scale = 8) +
    coord_flip() +
    scale_color_viridis(option = "plasma") +
    theme_minimal() +
    labs(title = "Top 5 Marker Genes per Cluster",
         x = "Gene", y = "Cluster") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 10))
  
  ggsave(filename = paste0(image_folder, "plot8_cluster_dotplot.png"), 
         plot = p8, width = 12, height = 10, dpi = 300, bg = "white")
}

# SUMMARY

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("All plots saved to:", image_folder, "\n")
cat("Excel files saved to:", output_folder, "\n")
cat("\nGenerated plots:\n")
cat("  1. plot1_split_umap.png - Split UMAP showing CTL vs PINK1\n")
cat("  2. plot2_cell_composition.png - Cell numbers across differentiation\n")
cat("  3. plot3_deg_temporal_dynamics.png - DEG counts over time\n")
cat("  4. plot4_gene_trajectories.png - Key gene expression trajectories\n")
cat("  5. plot5_go_mitochondrial.png - Mitochondrial pathway enrichment\n")
cat("  6. plot6_clustering_umap.png - Cell clusters visualization\n")
cat("  7. plot7_cluster_heatmap.png - Cluster marker heatmap\n")
cat("  8. plot8_cluster_dotplot.png - Cluster marker dot plot\n")

# PLOT 9: VOLCANO PLOTS FOR EACH DAY

cat("\n=== Creating Plot 9: Combined Volcano Plots ===\n")

volcano_plots <- list()

for(day in c(0,18,25,37,57)) {
  
  deg_subset <- DEG %>% filter(Day == day)
  
  # Add significance categories
  deg_subset$Significance <- "Not Significant"
  deg_subset$Significance[deg_subset$p_val_adj < 0.05 & deg_subset$avg_log2FC > 0.5] <- "Up in CTL"
  deg_subset$Significance[deg_subset$p_val_adj < 0.05 & deg_subset$avg_log2FC < -0.5] <- "Up in PINK1"
  
  deg_subset$Significance <- factor(deg_subset$Significance, 
                                     levels = c("Up in CTL", "Up in PINK1", "Not Significant"))
  
  # Top genes to label (most significant)
  top_genes <- deg_subset %>%
    filter(p_val_adj < 0.001, abs(avg_log2FC) > 1) %>%
    arrange(p_val_adj) %>%
    slice_head(n = 8) %>%
    pull(Gene)
  
  # Count significant genes
  n_up_ctl <- sum(deg_subset$Significance == "Up in CTL")
  n_up_pink1 <- sum(deg_subset$Significance == "Up in PINK1")
  
  p <- ggplot(deg_subset, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray40", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", linewidth = 0.3) +
    scale_color_manual(values = c("Up in CTL" = CTL_COLOR, 
                                   "Up in PINK1" = PINK1_COLOR,
                                   "Not Significant" = "gray70")) +
    geom_text_repel(data = deg_subset %>% filter(Gene %in% top_genes),
                    aes(label = Gene), size = 2, max.overlaps = 15,
                    box.padding = 0.3, point.padding = 0.2,
                    segment.color = "gray30", segment.size = 0.2,
                    min.segment.length = 0) +
    theme_minimal() +
    labs(title = paste0("Day ", day),
         subtitle = paste0("↑CTL: ", n_up_ctl, " | ↑PINK1: ", n_up_pink1),
         x = "log2 Fold Change (CTL/PINK1)",
         y = "-log10(Adj. P-value)",
         color = "") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
          legend.position = "none",
          plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 9)) +
    xlim(c(-max(abs(deg_subset$avg_log2FC), na.rm = TRUE) - 0.5, 
           max(abs(deg_subset$avg_log2FC), na.rm = TRUE) + 0.5))
  
  volcano_plots[[paste0("Day_", day)]] <- p
}

# Combine all volcano plots
combined_volcano <- wrap_plots(volcano_plots, ncol = 3) +
  plot_annotation(title = "Volcano Plots: Differential Gene Expression Across Differentiation",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)))

ggsave(filename = paste0(image_folder, "plot9_volcano_combined.png"), 
       plot = combined_volcano, width = 15, height = 10, dpi = 300, bg = "white")

# PLOT 10: TOP DEGs HEATMAP (Day-wise comparison)

cat("\n=== Creating Plot 10: Top DEGs Heatmap ===\n")

# Get top DEGs from each day
top_degs_per_day <- DEG %>%
  filter(p_val_adj < 0.01) %>%
  group_by(Day) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  pull(Gene) %>%
  unique()

if(length(top_degs_per_day) > 5) {
  
  DefaultAssay(S) = "RNA_imputed"
  S = ScaleData(S, assay = "RNA_imputed", features = top_degs_per_day, verbose = FALSE)
  
  # Extract scaled data
  scaled_data <- GetAssayData(S, assay = "RNA_imputed", layer = "scale.data")
  
  # Average expression per condition-day
  avg_exp <- sapply(levels(S$ConditionDay), function(cond) {
    cells <- WhichCells(S, expression = ConditionDay == cond)
    if(length(cells) > 0) {
      rowMeans(scaled_data[top_degs_per_day, cells, drop = FALSE])
    } else {
      rep(0, length(top_degs_per_day))
    }
  })
  rownames(avg_exp) <- top_degs_per_day
  
  # Create annotation
  annotation_col <- data.frame(
    Condition = c(rep("CTL", 5), rep("PINK1", 5)),
    Day = rep(c("0","18","25","37","57"), 2)
  )
  rownames(annotation_col) <- colnames(avg_exp)
  
  ann_colors <- list(
    Condition = c(CTL = CTL_COLOR, PINK1 = PINK1_COLOR),
    Day = setNames(DAY_COLORS, c("0","18","25","37","57"))
  )
  
  # Save pheatmap using grDevices (works without X11)
  # Open device
  grDevices::png(filename = paste0(image_folder, "plot10_top_degs_heatmap.png"), 
                 width = 3600, height = 4200, res = 300, type = "cairo-png")
  
  pheatmap(avg_exp,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = seq(-2, 2, length.out = 101),
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 7,
           fontsize_col = 10,
           main = "Top DEGs Across Differentiation (Top 10 per Day)",
           border_color = NA)
  
  dev.off()
}

# PLOT 11: COMPREHENSIVE GO ENRICHMENT - BIOLOGICAL PROCESS

cat("\n=== Creating Plot 11: GO Biological Process Enrichment ===\n")

go_bp_results_list <- list()

for(day in c(0,18,25,37,57)) {
  
  deg_day <- DEG %>% filter(Day == day, p_val_adj < 0.05)
  
  if(nrow(deg_day) < 10) next
  
  # Genes up in CTL
  genes_up_ctl <- deg_day %>% filter(avg_log2FC > 0.5) %>% pull(Gene)
  
  if(length(genes_up_ctl) >= 10) {
    gost_ctl <- gost(genes_up_ctl, 
                     organism = "hsapiens",
                     ordered_query = FALSE,
                     significant = TRUE,
                     user_threshold = 0.05,
                     correction_method = "fdr",
                     sources = c("GO:BP"))
    
    if(!is.null(gost_ctl$result)) {
      gost_ctl$result$Day <- day
      gost_ctl$result$Direction <- "CTL"
      go_bp_results_list[[paste0("Day", day, "_CTL")]] <- gost_ctl$result
    }
  }
  
  # Genes up in PINK1
  genes_up_pink1 <- deg_day %>% filter(avg_log2FC < -0.5) %>% pull(Gene)
  
  if(length(genes_up_pink1) >= 10) {
    gost_pink1 <- gost(genes_up_pink1, 
                       organism = "hsapiens",
                       ordered_query = FALSE,
                       significant = TRUE,
                       user_threshold = 0.05,
                       correction_method = "fdr",
                       sources = c("GO:BP"))
    
    if(!is.null(gost_pink1$result)) {
      gost_pink1$result$Day <- day
      gost_pink1$result$Direction <- "PINK1"
      go_bp_results_list[[paste0("Day", day, "_PINK1")]] <- gost_pink1$result
    }
  }
}

if(length(go_bp_results_list) > 0) {
  GO_BP_combined <- bind_rows(go_bp_results_list)
  write.xlsx(GO_BP_combined, paste0(output_folder, "GO_BP_enrichment_daywise.xlsx"))
  
  # Get top 5 terms per day/direction
  go_bp_top <- GO_BP_combined %>%
    group_by(Day, Direction) %>%
    arrange(p_value) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    mutate(term_name = str_trunc(term_name, width = 60))
  
  if(nrow(go_bp_top) > 0) {
    p11 <- ggplot(go_bp_top, 
                  aes(x = factor(Day), y = reorder(term_name, -p_value), 
                      size = intersection_size, color = -log10(p_value))) +
      geom_point(alpha = 0.7) +
      scale_color_viridis(option = "plasma") +
      scale_size_continuous(range = c(3, 10)) +
      facet_wrap(~Direction, scales = "free_y", ncol = 1) +
      theme_minimal() +
      labs(title = "GO Biological Process Enrichment Across Differentiation",
           subtitle = "Top 5 GO:BP terms per day and condition",
           x = "Differentiation Day", y = "GO Biological Process Term",
           size = "Gene Count", color = "-log10(p-value)") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 11),
            axis.text.y = element_text(size = 7),
            axis.text.x = element_text(size = 10),
            strip.text = element_text(face = "bold", size = 12),
            legend.position = "right")
    
    ggsave(filename = paste0(image_folder, "plot11_go_bp_enrichment.png"), 
           plot = p11, width = 14, height = 12, dpi = 300, bg = "white")
  }
}

# PLOT 12: CELLULAR COMPONENT GO ENRICHMENT

cat("\n=== Creating Plot 12: GO Cellular Component Enrichment ===\n")

go_cc_results_list <- list()

for(day in c(0,18,25,37,57)) {
  
  deg_day <- DEG %>% filter(Day == day, p_val_adj < 0.05)
  
  if(nrow(deg_day) < 10) next
  
  # Genes up in CTL
  genes_up_ctl <- deg_day %>% filter(avg_log2FC > 0.5) %>% pull(Gene)
  
  if(length(genes_up_ctl) >= 10) {
    gost_ctl <- gost(genes_up_ctl, 
                     organism = "hsapiens",
                     ordered_query = FALSE,
                     significant = TRUE,
                     user_threshold = 0.05,
                     correction_method = "fdr",
                     sources = c("GO:CC"))
    
    if(!is.null(gost_ctl$result)) {
      gost_ctl$result$Day <- day
      gost_ctl$result$Direction <- "CTL"
      go_cc_results_list[[paste0("Day", day, "_CTL")]] <- gost_ctl$result
    }
  }
  
  # Genes up in PINK1
  genes_up_pink1 <- deg_day %>% filter(avg_log2FC < -0.5) %>% pull(Gene)
  
  if(length(genes_up_pink1) >= 10) {
    gost_pink1 <- gost(genes_up_pink1, 
                       organism = "hsapiens",
                       ordered_query = FALSE,
                       significant = TRUE,
                       user_threshold = 0.05,
                       correction_method = "fdr",
                       sources = c("GO:CC"))
    
    if(!is.null(gost_pink1$result)) {
      gost_pink1$result$Day <- day
      gost_pink1$result$Direction <- "PINK1"
      go_cc_results_list[[paste0("Day", day, "_PINK1")]] <- gost_pink1$result
    }
  }
}

if(length(go_cc_results_list) > 0) {
  GO_CC_combined <- bind_rows(go_cc_results_list)
  write.xlsx(GO_CC_combined, paste0(output_folder, "GO_CC_enrichment_daywise.xlsx"))
  
  # Get top 5 terms per day/direction
  go_cc_top <- GO_CC_combined %>%
    group_by(Day, Direction) %>%
    arrange(p_value) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    mutate(term_name = str_trunc(term_name, width = 60))
  
  if(nrow(go_cc_top) > 0) {
    p12 <- ggplot(go_cc_top, 
                  aes(x = factor(Day), y = reorder(term_name, -p_value), 
                      size = intersection_size, color = -log10(p_value))) +
      geom_point(alpha = 0.7) +
      scale_color_viridis(option = "plasma") +
      scale_size_continuous(range = c(3, 10)) +
      facet_wrap(~Direction, scales = "free_y", ncol = 1) +
      theme_minimal() +
      labs(title = "GO Cellular Component Enrichment Across Differentiation",
           subtitle = "Top 5 GO:CC terms per day and condition",
           x = "Differentiation Day", y = "GO Cellular Component Term",
           size = "Gene Count", color = "-log10(p-value)") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 11),
            axis.text.y = element_text(size = 7),
            axis.text.x = element_text(size = 10),
            strip.text = element_text(face = "bold", size = 12),
            legend.position = "right")
    
    ggsave(filename = paste0(image_folder, "plot12_go_cc_enrichment.png"), 
           plot = p12, width = 14, height = 12, dpi = 300, bg = "white")
  }
}

# PLOT 13: DEG OVERLAP ACROSS DAYS

cat("\n=== Creating Plot 13: DEG Overlap Heatmap ===\n")

# Get DEGs for each day
deg_lists <- list()
for(day in c(0,18,25,37,57)) {
  deg_day <- DEG %>% 
    filter(Day == day, p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
    pull(Gene)
  deg_lists[[paste0("Day_", day)]] <- deg_day
}

# Create presence/absence matrix
all_genes <- unique(unlist(deg_lists))
overlap_matrix <- matrix(0, nrow = length(all_genes), ncol = length(deg_lists))
rownames(overlap_matrix) <- all_genes
colnames(overlap_matrix) <- names(deg_lists)

for(i in seq_along(deg_lists)) {
  overlap_matrix[deg_lists[[i]], i] <- 1
}

# Calculate overlap counts
overlap_counts <- matrix(0, nrow = length(deg_lists), ncol = length(deg_lists))
rownames(overlap_counts) <- names(deg_lists)
colnames(overlap_counts) <- names(deg_lists)

for(i in 1:length(deg_lists)) {
  for(j in 1:length(deg_lists)) {
    overlap_counts[i,j] <- length(intersect(deg_lists[[i]], deg_lists[[j]]))
  }
}

# Convert to long format for ggplot
overlap_df <- as.data.frame(overlap_counts) %>%
  rownames_to_column("Day1") %>%
  pivot_longer(-Day1, names_to = "Day2", values_to = "Count") %>%
  mutate(Day1 = factor(Day1, levels = names(deg_lists)),
         Day2 = factor(Day2, levels = names(deg_lists)))

p13 <- ggplot(overlap_df, aes(x = Day2, y = Day1, fill = Count)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Count), color = "black", size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "white", mid = "#FFA500", high = "#8B0000", 
                       midpoint = max(overlap_df$Count)/2) +
  theme_minimal() +
  labs(title = "DEG Overlap Across Differentiation Days",
       subtitle = "Number of shared differentially expressed genes between time points",
       x = "Day", y = "Day",
       fill = "Shared\nDEGs") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "right",
        panel.grid = element_blank())

ggsave(filename = paste0(image_folder, "plot13_deg_overlap_heatmap.png"), 
       plot = p13, width = 10, height = 8, dpi = 300, bg = "white")

# PLOT 14: TOP UP/DOWN REGULATED GENES (BAR PLOT)

cat("\n=== Creating Plot 14: Top Up/Down Regulated Genes ===\n")

# Get top 15 up and down genes across all days
top_up_ctl <- DEG %>%
  filter(p_val_adj < 0.01, avg_log2FC > 0) %>%
  group_by(Gene) %>%
  summarise(mean_fc = mean(avg_log2FC),
            max_fc = max(avg_log2FC),
            min_pval = min(p_val_adj),
            .groups = 'drop') %>%
  arrange(desc(mean_fc)) %>%
  slice_head(n = 15) %>%
  mutate(Direction = "Up in CTL")

top_up_pink1 <- DEG %>%
  filter(p_val_adj < 0.01, avg_log2FC < 0) %>%
  group_by(Gene) %>%
  summarise(mean_fc = mean(avg_log2FC),
            max_fc = max(abs(avg_log2FC)),
            min_pval = min(p_val_adj),
            .groups = 'drop') %>%
  arrange(mean_fc) %>%
  slice_head(n = 15) %>%
  mutate(Direction = "Up in PINK1")

top_genes_combined <- bind_rows(top_up_ctl, top_up_pink1) %>%
  mutate(Gene = factor(Gene, levels = Gene[order(mean_fc)]))

p14 <- ggplot(top_genes_combined, aes(x = mean_fc, y = Gene, fill = Direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = DIRECTION_COLORS) +
  theme_minimal() +
  labs(title = "Top 15 Consistently Up/Down Regulated Genes",
       subtitle = "Genes showing strongest average fold change across all time points (padj < 0.01)",
       x = "Mean log2 Fold Change",
       y = "Gene",
       fill = "Direction") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 11),
        legend.position = "top")

ggsave(filename = paste0(image_folder, "plot14_top_updown_genes.png"), 
       plot = p14, width = 11, height = 9, dpi = 300, bg = "white")

# Export the gene lists
write.xlsx(list("Top_Up_CTL" = top_up_ctl, "Top_Up_PINK1" = top_up_pink1),
           paste0(output_folder, "Top_consistently_regulated_genes.xlsx"))

# PLOT 15: EXPRESSION PATTERNS OF PARKINSON'S DISEASE GENES

cat("\n=== Creating Plot 15: PD-Related Gene Expression ===\n")

# Key Parkinson's disease related genes
pd_genes <- c("SNCA", "LRRK2", "PARK7", "PRKN", "PINK1", "DJ1", 
              "UCHL1", "GBA", "VPS35", "ATP13A2", "FBXO7", "PLA2G6")

# Check which are present
pd_genes_present <- pd_genes[pd_genes %in% rownames(S)]

if(length(pd_genes_present) > 2) {
  
  DefaultAssay(S) = "RNA_imputed"
  
  # Calculate average expression per condition/day
  pd_trajectories <- lapply(pd_genes_present, function(gene) {
    exp_data <- GetAssayData(S, assay = "RNA_imputed", layer = "data")[gene, ]
    
    trajectory <- S[[]] %>%
      mutate(Expression = exp_data) %>%
      group_by(Condition, Day) %>%
      summarise(Mean_Exp = mean(Expression),
                SE = sd(Expression) / sqrt(n()),
                .groups = 'drop') %>%
      mutate(Gene = gene)
    
    return(trajectory)
  })
  
  pd_trajectories_df <- bind_rows(pd_trajectories)
  
  p15 <- ggplot(pd_trajectories_df, 
                aes(x = as.numeric(as.character(Day)), y = Mean_Exp, 
                    color = Condition, group = Condition)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = Mean_Exp - SE, ymax = Mean_Exp + SE), 
                  width = 2, alpha = 0.4, linewidth = 0.6) +
    facet_wrap(~Gene, scales = "free_y", ncol = 4) +
    scale_color_manual(values = CONDITION_COLORS) +
    theme_minimal() +
    labs(title = "Parkinson's Disease-Related Gene Expression Trajectories",
         subtitle = "Expression dynamics of key PD-associated genes during differentiation",
         x = "Differentiation Day", y = "Average Expression (log-normalized)",
         color = "Genotype") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          legend.position = "top",
          strip.text = element_text(face = "bold", size = 10),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10))
  
  ggsave(filename = paste0(image_folder, "plot15_pd_gene_trajectories.png"), 
         plot = p15, width = 14, height = 10, dpi = 300, bg = "white")
}

# EXPORT GENE LISTS FOR STRING ANALYSIS

cat("\n=== Exporting Gene Lists for STRING Analysis ===\n")

# Export top DEGs from Day 25 (peak divergence) for STRING
day25_degs <- DEG %>%
  filter(Day == 25, p_val_adj < 0.01, abs(avg_log2FC) > 1) %>%
  arrange(p_val_adj)

# Up in CTL
day25_up_ctl <- day25_degs %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 100) %>%
  select(Gene, avg_log2FC, p_val_adj)

write.table(day25_up_ctl$Gene, 
            file = paste0(output_folder, "STRING_Day25_Up_CTL_genes.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Up in PINK1
day25_up_pink1 <- day25_degs %>%
  filter(avg_log2FC < -1) %>%
  slice_head(n = 100) %>%
  select(Gene, avg_log2FC, p_val_adj)

write.table(day25_up_pink1$Gene, 
            file = paste0(output_folder, "STRING_Day25_Up_PINK1_genes.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Export all significant DEGs for Day 25
write.xlsx(day25_degs, paste0(output_folder, "Day25_DEGs_for_STRING.xlsx"))

cat("\n=== Gene lists exported for STRING analysis ===\n")    # For string db
cat("  1. STRING_Day25_Up_CTL_genes.txt (Top 100 genes up in CTL)\n")
cat("  2. STRING_Day25_Up_PINK1_genes.txt (Top 100 genes up in PINK1)\n")

# SUMMARY

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Additional plots saved to:", image_folder, "\n")
cat("  9. plot9_volcano_combined.png - Volcano plots for all 5 days\n")
cat(" 10. plot10_top_degs_heatmap.png - Heatmap of top DEGs across days\n")
cat(" 11. plot11_go_bp_enrichment.png - GO Biological Process enrichment\n")
cat(" 12. plot12_go_cc_enrichment.png - GO Cellular Component enrichment\n")
cat(" 13. plot13_deg_overlap_heatmap.png - DEG overlap between days\n")
cat(" 14. plot14_top_updown_genes.png - Top consistently regulated genes\n")
cat(" 15. plot15_pd_gene_trajectories.png - PD-related gene expression\n")
cat("\nGene lists for STRING analysis:\n")
cat("  - STRING_Day25_Up_CTL_genes.txt\n")
cat("  - STRING_Day25_Up_PINK1_genes.txt\n")
cat("  - Day25_DEGs_for_STRING.xlsx\n")