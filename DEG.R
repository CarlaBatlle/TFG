## DESeq 2

# load libraries 
library(DESeq2)
library(dplyr)
library(org.Hs.eg.db)
library(annotate)
library(clusterProfiler)
library(ggpubr)
library(pheatmap)
library(tibble)
library(tidyr)
library(factoextra)
library(FactoMineR)
library(caret)
library(umap)
library(enrichplot)
library(ggplot2)
library(forcats)
library(ggupset)

# set working directory
setwd("~/Documents/practiques/data")


# FUNCIO PCA METADATA
pca_donant_analysis <- function(file_path = "donor.csv") {
  library(ggplot2)
  
  metaData <- read.csv(file_path, header = TRUE, sep = ",")
  
  vars_pca <- metaData[, c("record_id", "donorage", "donorsex", "donationtype", 
                           "donorheight", "donorweight", "bodymassindex", "diagnosis_computed")]
  
  # Eliminar NA
  vars_pca <- na.omit(vars_pca)
  
  # Eliminar nens (menors de 10 anys)
  vars_pca <- vars_pca[vars_pca$donorage >= 10, ]
  
  # Convertir sexe a binari
  vars_pca$donorsex <- ifelse(vars_pca$donorsex == "Male", 1, 0)
  
  # Variables dummy per donationtype
  vars_pca <- cbind(vars_pca[, !names(vars_pca) %in% "donationtype"],
                    model.matrix(~ donationtype - 1, data = vars_pca))
  
  labels <- vars_pca$record_id
  colors <- vars_pca$diagnosis_computed
  vars_numeric <- vars_pca[, !(names(vars_pca) %in% c("record_id", "diagnosis_computed"))]
  
  pca_result <- prcomp(vars_numeric, scale. = TRUE)
  pca_df <- as.data.frame(pca_result$x)
  pca_df$record_id <- labels
  pca_df$diagnosis_computed <- colors
  
  # Visualització PCA
  p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = diagnosis_computed, label = record_id)) +
    geom_point(size = 3) +
    geom_text(vjust = -1, size = 3.5) +
    scale_color_manual(values = c("None" = "#1E90FF", "Type2" = "#FF69B4", "Pre.T2D" = "orange")) +
    theme_minimal() +
    labs(title = "PCA de variables de donant", x = "PC1", y = "PC2") +
    theme(legend.position = "top", legend.title = element_blank())
  
  print(p1)
  
  # Barplots de contribucions
  loadings_list <- list()
  for (pc in c("PC1", "PC2")) {
    loadings_sorted <- sort(abs(pca_result$rotation[, pc]), decreasing = TRUE)
    df_loadings <- data.frame(Variable = names(loadings_sorted), Loading = loadings_sorted)
    loadings_list[[pc]] <- df_loadings
    
    p <- ggplot(df_loadings, aes(x = reorder(Variable, Loading), y = Loading)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(title = paste("Contribució de cada variable a", pc),
           x = "Variable", y = "Pes (absolut)") +
      theme_minimal()
    
    print(p)
  }
  
  # Retornar resultats
  invisible(list(
    pca_result = pca_result,
    pca_data = pca_df,
    loadings = loadings_list
  ))
}

full_expression_boxplots <- function(count_file, tpm_file) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  
  # Llegeix els fitxers d’entrada
  countData <- read.csv(count_file)
  countData_2 <- read.csv(count_file, row.names = 1)
  metaData <- read.csv('donor.csv', header = TRUE, sep = ",")
  metaData <- as.data.frame(metaData)
  
  # Filtrar només els pacients presents a countData
  filtered_df <- subset(metaData, record_id %in% colnames(countData))
  tpmData <- read.csv(tpm_file, row.names = 1)
  
  # Transformació log10(TPM + 1)
  tpmData_log <- as.data.frame(lapply(tpmData, function(col) {
    suppressWarnings(log10(as.numeric(as.character(col)) + 1))
  }))
  rownames(tpmData_log) <- rownames(tpmData)
  
  # Merge TPM + diagnosis
  tpm_data_m <- tpmData_log %>% rownames_to_column(var = "record_id")
  metaData$record_id <- as.character(metaData$record_id)
  tpmData_merged <- merge(tpm_data_m, metaData[, c("record_id", "diagnosis_computed", "donorage")], by = "record_id", all.x = TRUE)
  tpmData_merged$diagnosis <- as.factor(tpmData_merged$diagnosis_computed)
  
  # Filtrar per diagnòstics d’interès
  tpmData_merged <- tpmData_merged %>% filter(diagnosis_computed %in% c("None", "Pre.T2D", "Type2"))
  table(tpmData_merged$diagnosis_computed)
  
  # 🔹 Afegim el boxplot de la mitjana d'expressió de tots els gens
  genes_of_interest <- c("INS", "GCG", "MKI67", "ARX")
  tpmData_merged$mean_expression <- rowMeans(tpmData_merged[, genes_of_interest], na.rm = TRUE)
  
  print(
    ggplot(tpmData_merged, aes(x = factor(record_id, levels = unique(record_id)), y = mean_expression, fill = diagnosis_computed)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(aes(color = diagnosis_computed), width = 0.3, size = 2) +
      scale_fill_manual(values = c("None" = "blue", "Pre.T2D" = "orange", "Type2" = "pink")) +
      scale_color_manual(values = c("None" = "blue", "Pre.T2D" = "orange", "Type2" = "pink")) +
      labs(title = "Boxplot de la mitjana d'expressió dels gens (log₁₀ TPM + 1)",
           x = "Pacients", y = "Mitjana log₁₀(TPM + 1)") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
      scale_x_discrete(guide = guide_axis(n.dodge = 2))
  )
  
  # 🔸 FILTRATGES EXTRA 🔸
  
  # 1. Eliminar pacients menors de 10 anys
  tpmData_merged <- tpmData_merged %>% filter(donorage >= 10)
  
  # 2. Eliminar mostres amb >80% gens amb expressió 0 per condició
  zero_prop_per_sample <- colSums(tpmData_log == 0, na.rm = TRUE) / nrow(tpmData_log)
  sample_condition <- metaData$diagnosis_computed
  names(sample_condition) <- metaData$record_id
  
  samples_to_remove <- c()
  for (cond in unique(sample_condition)) {
    cond_samples <- names(sample_condition)[sample_condition == cond]
    cond_zeros <- zero_prop_per_sample[cond_samples]
    to_remove <- names(cond_zeros[cond_zeros > 0.8])
    samples_to_remove <- c(samples_to_remove, to_remove)
  }
  
  # 3. Eliminar pacients específics manualment
  manual_exclude <- c("R032", "R036", "R228", "R233", "R234", "R235", "R246", "R240", "R064")
  samples_to_remove <- unique(c(samples_to_remove, manual_exclude))
  
  # Aplicar el filtratge global
  tpmData_merged <- tpmData_merged[!tpmData_merged$record_id %in% samples_to_remove, ]
  
  # 🔹 Funció interna per a boxplot per gen
  plot_gene_expression <- function(gene) {
    ggplot(tpmData_merged, aes(x = diagnosis_computed, y = .data[[gene]], fill = diagnosis)) +
      geom_boxplot() +
      labs(title = paste("Expressió de", gene, "(log₁₀ TPM + 1)"),
           x = "Diagnòstic",
           y = "log₁₀(TPM + 1)") +
      theme_minimal() +
      scale_fill_manual(values = c("None" = "blue", "Pre.T2D" = "orange", "Type2" = "red"))
  }
  
  # 🔹 Mostrar boxplots globals per gen
  lapply(genes_of_interest, function(gene) print(plot_gene_expression(gene)))
  
  # 🔹 Boxplots individuals per pacient
  tpmData_genes <- tpmData_merged[, c("record_id", genes_of_interest, "diagnosis_computed")]
  
  for (gene in genes_of_interest) {
    print(
      ggplot(tpmData_genes, aes(x = factor(record_id, levels = unique(record_id)), y = .data[[gene]], fill = diagnosis_computed)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.6) +
        geom_jitter(aes(color = diagnosis_computed), width = 0.3, size = 2) +
        scale_fill_manual(values = c("None" = "blue", "Pre.T2D" = "orange", "Type2" = "pink")) +
        scale_color_manual(values = c("None" = "blue", "Pre.T2D" = "orange", "Type2" = "pink")) +
        labs(title = paste("Boxplot de l'expressió gènica de", gene, "(log₁₀ TPM + 1)"),
             x = "Pacients", y = "log₁₀(TPM + 1)") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
        scale_x_discrete(guide = guide_axis(n.dodge = 2))
    )
  }
}





distrib_dades <- function(count_file, meta_file = "donor.csv", pacients_a_eliminar = c("R032", "R036", "R228", "R233", "R234", "R235", "R246", "R240", "R064")) {
  # Llegim dades de comptatge i filtrem pacients
  countData <- read.csv(count_file, row.names = 1)
  countData <- countData[, !colnames(countData) %in% pacients_a_eliminar]
  
  # Llegim metadades i filtrem per les mostres presents
  metaData <- read.csv(meta_file, header = TRUE, sep = ",")
  metaData_filtered <- subset(metaData, record_id %in% colnames(countData))
  rownames(metaData_filtered) <- metaData_filtered$record_id
  metaData_filtered$record_id <- NULL
  
  # 🔴 Filtratge per eliminar pacients nens (donorage < 10)
  metaData_filtered <- subset(metaData_filtered, donorage >= 18)
  countData <- countData[, rownames(metaData_filtered)]
  
  # Creem la categoria de BMI
  metaData_filtered$BMI_category <- cut(metaData_filtered$bodymassindex,
                                        breaks = c(-Inf, 18.5, 24.9, 29.9, Inf),
                                        labels = c("Prim", "Normal", "Sobrepes", "Obesitat"),
                                        right = TRUE)
  
  # Assegurem que diagnosis_computed és un factor
  metaData_filtered$diagnosis_computed <- factor(metaData_filtered$diagnosis_computed)
  
  # Transformació logarítmica de les dades de comptatge
  countData_log <- log2(countData + 1)
  
  # Histograma de les dades log-transformades abans del filtratge
  countData_log_melted <- as.data.frame(countData_log)
  countData_log_melted$genes <- rownames(countData_log)
  countData_log_melted <- reshape2::melt(countData_log_melted, id.vars = "genes")
  
  p_hist_before <- ggplot(countData_log_melted, aes(x = value)) +
    geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black", alpha = 0.7) +
    theme_minimal() +
    labs(title = "Distribució de les dades log-transforms (Abans del filtratge)", x = "Log2(Count + 1)", y = "Frequència")
  
  # Mostrar histograma abans del filtratge
  print(p_hist_before)
  
  # Calcular la mitjana i el coeficient de variació (CV) per cada gen
  gene_means <- rowMeans(countData_log)  # Mitjana per gen
  gene_sds <- apply(countData_log, 1, sd)  # Desviació estàndard per gen
  gene_cv <- gene_sds / gene_means  # Coeficient de variació per gen
  
  # Crear el gràfic mitjana vs. CV
  cv_df <- data.frame(gene = rownames(countData_log), mean = gene_means, cv = gene_cv)
  p_cv <- ggplot(cv_df, aes(x = mean, y = cv)) +
    geom_point(color = "darkgreen", alpha = 0.7) +
    theme_minimal() +
    labs(title = "Mitjana vs. Coeficient de Variació per Gen", x = "Mitjana de l'expressió (Log2)", y = "Coeficient de Variació (CV)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Mostrar gràfic mitjana vs. CV
  print(p_cv)
  
  # Filtratge pel 80% d'expressió només en mostres SAS (None) o Type2
  sas_type2_samples <- rownames(metaData_filtered)[metaData_filtered$diagnosis_computed %in% c("None", "Type2")]
  countData_log_sas_type2 <- countData_log[, sas_type2_samples]
  
  # Calcular la proporció de valors >0 per gen en aquest subgrup
  expr_prop <- apply(countData_log_sas_type2, 1, function(x) mean(x > 0))
  
  # Mantenim només els gens que s'expressen en almenys el 20% de les mostres SAS/Type2
  countData_log_filtered_zeros <- countData_log[expr_prop >= 0.2, ]
  
  # Eliminar gens amb variància zero
  variances <- apply(countData_log_filtered_zeros, 1, var)
  countData_log_filtered <- countData_log_filtered_zeros[variances > 0, ]
  
  # Generem el histograma després del filtratge
  countData_log_melted_filtered <- as.data.frame(countData_log_filtered)
  countData_log_melted_filtered$genes <- rownames(countData_log_filtered)
  countData_log_melted_filtered <- reshape2::melt(countData_log_melted_filtered, id.vars = "genes")
  
  p_hist_after <- ggplot(countData_log_melted_filtered, aes(x = value)) +
    geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black", alpha = 0.7) +
    theme_minimal() +
    labs(title = "Distribució de les dades log-transforms (Després del filtratge)", x = "Log2(Count + 1)", y = "Frequència")
  
  # Mostrar el histograma després del filtratge
  print(p_hist_after)
  
  # Realitzem la PCA sobre les mostres filtrades
  pca_res <- prcomp(t(countData_log_filtered), scale. = TRUE)  # Transposem per aplicar PCA a les mostres
  pca_df <- as.data.frame(pca_res$x)  # Obtenim els resultats de la PCA
  pca_df$sample <- rownames(pca_df)
  pca_df$diagnosis_computed <- metaData_filtered$diagnosis_computed
  pca_df$donationtype <- metaData_filtered$donationtype
  pca_df$BMI_category <- metaData_filtered$BMI_category  # Afegim BMI_category a PCA
  pca_df$donorsex <- metaData_filtered$donorsex  # Afegim donorsex a PCA
  
  # Gràfiques PCA per diverses categories
  p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = diagnosis_computed)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "PCA segons diagnosis_computed")
  
  p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = donationtype)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "PCA segons donationtype")
  
  p3 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = BMI_category)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "PCA segons BMI_category")
  
  p4 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = donorsex)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "PCA segons donorsex")
  
  # Mostrar gràfiques PCA
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  
  # Realitzar UMAP sobre les mostres filtrades
  library(umap)
  umap_res <- umap(t(countData_log_filtered))
  umap_df <- as.data.frame(umap_res$layout)  # Obtenim els resultats de UMAP
  umap_df$sample <- rownames(umap_df)
  umap_df$diagnosis_computed <- metaData_filtered$diagnosis_computed
  umap_df$donationtype <- metaData_filtered$donationtype
  umap_df$BMI_category <- metaData_filtered$BMI_category  # Afegim BMI_category a UMAP
  umap_df$donorsex <- metaData_filtered$donorsex  # Afegim donorsex a UMAP
  
  # Gràfiques UMAP
  u1 <- ggplot(umap_df, aes(x = V1, y = V2, color = diagnosis_computed)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "UMAP segons diagnosis_computed")
  
  u2 <- ggplot(umap_df, aes(x = V1, y = V2, color = donationtype)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "UMAP segons donationtype")
  
  u3 <- ggplot(umap_df, aes(x = V1, y = V2, color = BMI_category)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "UMAP segons BMI_category")
  
  u4 <- ggplot(umap_df, aes(x = V1, y = V2, color = donorsex)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "UMAP segons donorsex")
  
  # Mostrar gràfiques UMAP
  print(u1)
  print(u2)
  print(u3)
  print(u4)
  
  # Retornar les dades i resultats
  return(list(
    countData = countData,
    countData_log = countData_log,
    countData_log_filtered = countData_log_filtered,
    metaData_filtered = metaData_filtered,
    pca = pca_res,
    umap = umap_res
  ))
}


run_deg_analysis <- function(count_file, meta_file, pacients_a_eliminar = c("R032", "R036", "R228", "R233", "R234", "R235", "R246", "R240", "R064"), log2FC_threshold = 1) {
  
  # ─── 1. Llegir i netejar les dades ────────────────────────────────────────────────
  
  # Llegir matriu de comptatge i eliminar pacients específics
  countData <- read.csv(count_file, row.names = 1)
  countData <- countData[, !colnames(countData) %in% pacients_a_eliminar]
  
  # Llegir metadades i filtrar només mostres presents a countData
  metaData <- read.csv(meta_file, header = TRUE, sep = ",")
  metaData_filtered <- subset(metaData, record_id %in% colnames(countData))
  rownames(metaData_filtered) <- metaData_filtered$record_id
  metaData_filtered$record_id <- NULL
  
  # 🔁 Filtrar per edat: només adults (donorage >= 10 anys)
  metaData_filtered <- metaData_filtered[metaData_filtered$donorage >= 10, ]
  countData <- countData[, rownames(metaData_filtered)]
  
  # 🔁 Eliminar pacients amb diagnosis_computed == PreType2
  metaData_filtered <- metaData_filtered[metaData_filtered$diagnosis_computed != "PreType2", ]
  countData <- countData[, rownames(metaData_filtered)]
  
  # 🔁 Eliminar mostres amb discrepància entre diagnosis i diagnosis_computed
  matching_diag <- metaData_filtered$diagnosis == metaData_filtered$diagnosis_computed
  metaData_filtered <- metaData_filtered[matching_diag, ]
  countData <- countData[, rownames(metaData_filtered)]
  
  # ─── 2. Selecció per DEG (None vs Type2) ──────────────────────────────────────────
  
  meta_deg <- metaData_filtered[metaData_filtered$diagnosis_computed %in% c("None", "Type2"), ]
  count_deg <- countData[, rownames(meta_deg)]
  meta_deg$diagnosis_computed <- factor(meta_deg$diagnosis_computed, levels = c("None", "Type2"))
  
  # 🔁 Filtratge per gens amb >80% de zeros PER CONDICIÓ
  samples_none <- rownames(meta_deg[meta_deg$diagnosis_computed == "None", ])
  samples_type2 <- rownames(meta_deg[meta_deg$diagnosis_computed == "Type2", ])
  
  zero_none <- apply(count_deg[, samples_none], 1, function(x) sum(x == 0) / length(x))
  zero_type2 <- apply(count_deg[, samples_type2], 1, function(x) sum(x == 0) / length(x))
  
  genes_to_keep <- names(which(zero_none <= 0.8 & zero_type2 <= 0.8))
  count_deg_filtered <- count_deg[genes_to_keep, ]
  
  cat("Nombre de gens eliminats per >80% de zeros (per condició):", nrow(count_deg) - length(genes_to_keep), "\n")
  
  # ─── 3. DESeq2 ────────────────────────────────────────────────────────────────────
  
  dds <- DESeqDataSetFromMatrix(countData = round(count_deg_filtered),
                                colData = meta_deg,
                                design = ~ diagnosis_computed)
  dds <- DESeq(dds)
  res <- results(dds)
  res <- res[order(res$padj), ]
  res2<- res[order(res$pvalue), ]
  res <- res[!is.na(res$log2FoldChange) & !is.na(res$pvalue), ]
  
  # Filtrar per padj i log2FC
  res_padj <- subset(res, padj < 0.05)
  res_pvalue<-subset(res, pvalue < 0.05)
  res_filtered <- subset(res_padj, abs(log2FoldChange) > log2FC_threshold)
  res_filtered_pval <- subset(res_pvalue, abs(log2FoldChange) > log2FC_threshold)
  
  # ─── 4. Transformació log2 +1 i dades per visualitzar ────────────────────────────
  
  countData_log <- log2(countData + 1)
  countData_logDEG <- countData_log[rownames(countData_log) %in% rownames(res_filtered), ]
  
  # 🔁 Obtenir els dos gens amb més canvi d’expressió
  top_2_genes <- head(res_filtered[order(abs(res_filtered$log2FoldChange), decreasing = TRUE), ], 2)
  top_2_genes_names <- rownames(top_2_genes)
  
  countData_2_genes <- countData_log[top_2_genes_names, colnames(countData_logDEG), drop = FALSE]
  
  # ─── 5. Output ───────────────────────────────────────────────────────────────────
  
  cat("Comparació DEG entre:", levels(meta_deg$diagnosis_computed)[1], "i", levels(meta_deg$diagnosis_computed)[2], "\n")
  cat("Nombre de gens DEG amb |log2FC| >", log2FC_threshold, "i padj < 0.05:", nrow(res_filtered), "\n")
  
  return(list(
    res = res,
    res_deg = res_filtered,
    res_deg_pval = res_filtered_pval,
    countData_logDEG = countData_logDEG,
    countData_2_genes = countData_2_genes,
    metaData_filtered = metaData_filtered,
    top_2_genes = top_2_genes,
    top_2_genes_names = top_2_genes_names
  ))
}




visualize_deg_results <- function(results_list, log2FC_threshold = 1, tag = "analisi") {
  res <- results_list$res
  res_deg <- results_list$res_deg
  countData_logDEG <- results_list$countData_logDEG
  countData_2_genes <- results_list$countData_2_genes
  metaData_filtered <- results_list$metaData_filtered
  top_2_genes <- results_list$top_2_genes
  top_2_genes_names <- rownames(top_2_genes)
  
  res_volcano <- res[!is.na(res$padj), ]  # Filtra només els gens amb padj disponible
  
  res_volcano$diffexpressed <- "NO"
  res_volcano$diffexpressed[res_volcano$log2FoldChange > log2FC_threshold & res_volcano$padj < 0.05] <- "UP"
  res_volcano$diffexpressed[res_volcano$log2FoldChange < -log2FC_threshold & res_volcano$padj < 0.05] <- "DOWN"
  
  volcano_plot <- ggplot(res_volcano, aes(x = log2FoldChange, y = -log10(pvalue), color = diffexpressed)) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("UP" = "red3", "DOWN" = "blue4", "NO" = "grey80")) +
    labs(title = paste("Volcano Plot -", tag), x = "Log2 Fold Change", y = "-Log10 P-value", color = "Expressió") +
    theme_minimal(base_size = 12)
  print(volcano_plot)
  
  # Filtrar "PreType2" per eliminar-la
  metaData_filtered <- metaData_filtered[metaData_filtered$diagnosis_computed != "PreType2", ]
  
  # Heatmap (sense PreType2)
  condition_vector <- metaData_filtered$diagnosis_computed[match(colnames(countData_logDEG), rownames(metaData_filtered))]
  annotation_col <- data.frame(Condition = factor(condition_vector, levels = c("None", "Type2")))
  rownames(annotation_col) <- colnames(countData_logDEG)
  ann_colors <- list(Condition = c("None" = "lightblue", "Type2" = "pink"))
  
  pheatmap(countData_logDEG,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           cluster_rows = TRUE, cluster_cols = TRUE, scale = "row",
           fontsize_row = 8, fontsize_col = 10,
           main = paste("Heatmap DEG -", tag),
           annotation_col = annotation_col, annotation_colors = ann_colors,
           labels_row = NULL, treeheight_row = 0, show_colnames = FALSE, show_rownames = FALSE)
  
  # PCA (sense PreType2)
  pca_result <- prcomp(t(countData_logDEG), scale. = TRUE)
  pca_data <- as.data.frame(pca_result$x)
  pca_data$pacient <- rownames(pca_data)
  pca_data$condicio <- metaData_filtered$diagnosis_computed[match(pca_data$pacient, rownames(metaData_filtered))]
  
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condicio)) +
    geom_point(size = 4) +
    scale_color_manual(values = c("None" = "lightblue", "Type2" = "pink")) +
    labs(title = paste("PCA DEG -", tag), x = "PC1", y = "PC2") +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = "top")
  print(pca_plot)
  
  # UMAP (sense PreType2)
  umap_result <- umap(t(countData_logDEG))
  umap_data <- as.data.frame(umap_result$layout)
  umap_data$pacient <- rownames(umap_data)
  umap_data$condicio <- metaData_filtered$diagnosis_computed[match(umap_data$pacient, rownames(metaData_filtered))]
  
  umap_plot <- ggplot(umap_data, aes(x = V1, y = V2, color = condicio)) +
    geom_point(size = 4) +
    scale_color_manual(values = c("None" = "lightblue", "Type2" = "pink")) +
    labs(title = paste("UMAP DEG -", tag), x = "UMAP1", y = "UMAP2") +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = "top")
  print(umap_plot)
  
  # Heatmap top 2 gens (sense PreType2)
  annotation_col_top2 <- annotation_col[match(colnames(countData_2_genes), rownames(annotation_col)), , drop = FALSE]
  pheatmap(countData_2_genes,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           cluster_rows = TRUE, cluster_cols = TRUE, scale = "row",
           fontsize_row = 8, fontsize_col = 10,
           main = paste("Heatmap top2 gens -", tag),
           annotation_col = annotation_col_top2, annotation_colors = ann_colors,
           labels_row = NULL, treeheight_row = 0, show_colnames = FALSE, show_rownames = FALSE)
  
  # PCA top 2 gens (sense PreType2)
  pca2_result <- prcomp(t(countData_2_genes), scale. = TRUE)
  pca2_data <- as.data.frame(pca2_result$x)
  pca2_data$pacient <- colnames(countData_2_genes)
  pca2_data$condicio <- metaData_filtered$diagnosis_computed[match(pca2_data$pacient, rownames(metaData_filtered))]
  
  pca2_plot <- ggplot(pca2_data, aes(x = PC1, y = PC2, color = condicio)) +
    geom_point(size = 4) +
    scale_color_manual(values = c("None" = "lightblue", "Type2" = "pink")) +
    labs(title = paste("PCA top2 gens -", tag), x = "PC1", y = "PC2") +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = "top")
  print(pca2_plot)
  
  # Scatter Plot top 2 gens (sense PreType2)
  gene1 <- top_2_genes_names[1]
  gene2 <- top_2_genes_names[2]
  scatter_df <- data.frame(
    pacient = colnames(countData_2_genes),
    expr_gene1 = as.numeric(countData_2_genes[gene1, ]),
    expr_gene2 = as.numeric(countData_2_genes[gene2, ]),
    condicio = metaData_filtered$diagnosis_computed[match(colnames(countData_2_genes), rownames(metaData_filtered))]
  )
  
  scatter_plot <- ggplot(scatter_df, aes(x = expr_gene1, y = expr_gene2, color = condicio)) +
    geom_point(size = 4) +
    scale_color_manual(values = c("None" = "lightblue", "Type2" = "pink")) +
    labs(
      title = paste("Scatter Plot top2 gens DEG -", tag),
      x = paste("Expressió log2 de", gene1),
      y = paste("Expressió log2 de", gene2)
    ) +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = "top")
  print(scatter_plot)
}


get_deg_results <- function(deg_results, pvalue_cutoff = 0.05, log2FC_cutoff = 1) {
  # Convertir DESeqResults a data.frame
  deg_results_df <- as.data.frame(deg_results)
  
  # Afegir identificador de gen com columna
  deg_results_df$gene <- rownames(deg_results_df)
  
  # Filtrar els gens DEG per padj i log2FoldChange
  deg_table <- deg_results_df %>%
    filter(pvalue < pvalue_cutoff, abs(log2FoldChange) > log2FC_cutoff) %>%
    select(gene, log2FoldChange, padj, pvalue) # Seleccionar la columna 'gene'
  
  # Convertir ENSEMBL a SYMBOL i ENTREZID mantenint l'ordre original
  conversion <- clusterProfiler::bitr(
    deg_results_df$gene,  # Usar la columna gene
    fromType = "ENSEMBL",
    toType = c("SYMBOL", "ENTREZID"),
    OrgDb = org.Hs.eg.db
  )
  
  # Eliminar duplicats en ENSEMBL (evitar múltiples SYMBOLs/ENTREZ per un gen)
  conversion <- conversion[!duplicated(conversion$ENSEMBL), ]
  
  # Filtrar res_filtered perquè només contingui gens amb conversió vàlida
  res_filtered <- deg_results_df[deg_results_df$gene %in% conversion$ENSEMBL, ]
  
  # Afegir SYMBOLs com rownames a res_filtered
  rownames(res_filtered) <- conversion$SYMBOL[match(res_filtered$gene, conversion$ENSEMBL)]
  
  # Afegir SYMBOL i ENTREZID a deg_table
  deg_table$SYMBOL <- conversion$SYMBOL[match(deg_table$gene, conversion$ENSEMBL)]
  deg_table$ENTREZID <- conversion$ENTREZID[match(deg_table$gene, conversion$ENSEMBL)]
  
  # Retornar les dues taules: res_filtered (original amb rownames=SYMBOL) i deg_table (amb SYMBOL i ENTREZID)
  return(list(res_filtered = res_filtered, deg_table = deg_table))
}


countData <- read.csv("unproc_pbrna_Alpha.csv", row.names = 1)


pca_donant_analysis("donor.csv")

full_expression_boxplots("unproc_rnaseq.csv", "TPM_merged.csv")

distrib_bulk <- distrib_dades("unproc_rnaseq.csv")
distrib_beta <- distrib_dades("unproc_pbrna_Beta.csv")
distrib_alpha <- distrib_dades("unproc_pbrna_Alpha.csv")

deg_bulk <- run_deg_analysis(
  count_file = "unproc_rnaseq.csv",
  meta_file = "donor.csv",
  log2FC_threshold = 1
)

deg_beta <- run_deg_analysis(
  count_file = "unproc_pbrna_Beta.csv",
  meta_file = "donor.csv",
  log2FC_threshold = 1
)

deg_alpha <- run_deg_analysis(
  count_file = "unproc_pbrna_Alpha.csv",
  meta_file = "donor.csv",
  log2FC_threshold = 1
)

deg_alpha$top_2_genes_names

visualize_deg_results(deg_bulk, log2FC_threshold = 1, tag = "BULK")
visualize_deg_results(deg_beta, log2FC_threshold = 1, tag = "BETA")
visualize_deg_results(deg_alpha, log2FC_threshold = 1, tag = "ALPHA")


result_beta<- get_deg_results(deg_beta$res_deg_pval)
res_filtered_beta <- result_beta$res_filtered
deg_table_beta <- result_beta$deg_table
gene_list_beta <-deg_table_beta$SYMBOL
clean_list_beta <- na.omit(gene_list_beta)
cat(clean_list_beta, sep = "\n")

result_beta<- get_deg_results(deg_beta$res_deg)
res_filtered_beta <- result_beta$res_filtered
deg_table_beta <- result_beta$deg_table
gene_list_beta <-deg_table_beta$SYMBOL
clean_list_beta <- na.omit(gene_list_beta)
cat(clean_list_beta, sep = "\n")


result_alpha <- get_deg_results(deg_alpha$res_deg)
res_filtered_alpha <- result_alpha$res_filtered
deg_table_alpha <- result_alpha$deg_table
gene_list_alpha <- deg_table_alpha$SYMBOL
clean_list_alpha <- na.omit(gene_list_alpha)
cat(clean_list_alpha, sep = "\n")



































distrib_dades <- function(count_file, meta_file = "donor.csv", pacients_a_eliminar = c("R032", "R036", "R228", "R233", "R234", "R235", "R246", "R240", "R064"), n_top_genes = 2000) {
  # Llegim dades de comptatge i filtrem pacients
  countData <- read.csv(count_file, row.names = 1)
  countData <- countData[, !colnames(countData) %in% pacients_a_eliminar]
  
  # Llegim metadades i filtrem per les mostres presents
  metaData <- read.csv(meta_file, header = TRUE, sep = ",")
  metaData_filtered <- subset(metaData, record_id %in% colnames(countData))
  rownames(metaData_filtered) <- metaData_filtered$record_id
  metaData_filtered$record_id <- NULL
  
  # Creem la categoria de BMI
  metaData_filtered$BMI_category <- cut(metaData_filtered$bodymassindex,
                                        breaks = c(-Inf, 18.5, 24.9, 29.9, Inf),
                                        labels = c("Prim", "Normal", "Sobrepes", "Obesitat"),
                                        right = TRUE)
  
  # Assegurem que diagnosis_computed és un factor
  metaData_filtered$diagnosis_computed <- factor(metaData_filtered$diagnosis_computed)
  
  # Transformació logarítmica de les dades de comptatge
  countData_log <- log2(countData + 1)
  
  # Calcular la mitjana i el coeficient de variació (CV) per cada gen
  gene_means <- rowMeans(countData_log)  # Mitjana per gen
  gene_sds <- apply(countData_log, 1, sd)  # Desviació estàndard per gen
  gene_cv <- gene_sds / gene_means  # Coeficient de variació per gen
  
  # Crear el dataframe amb el coeficient de variació
  cv_df <- data.frame(gene = rownames(countData_log), mean = gene_means, cv = gene_cv)
  
  # Filtratge per coeficient de variació (top 2000 gens amb més variabilitat)
  cv_df_sorted <- cv_df[order(cv_df$cv, decreasing = TRUE), ]  # Ordenem per CV (de major a menor)
  top_genes_cv <- head(cv_df_sorted$gene, n_top_genes)  # Seleccionem els 2000 gens amb major CV
  
  # Filtratge de les dades per mantenir només els top gens seleccionats
  countData_log_filtered_top_genes <- countData_log[top_genes_cv, ]
  
  # Generar el gràfic del Coeficient de Variació per als 2000 gens seleccionats
  cv_selected <- cv_df[cv_df$gene %in% top_genes_cv, ]
  
  p_cv <- ggplot(cv_selected, aes(x = mean, y = cv)) +
    geom_point(color = "darkgreen", alpha = 0.7) +
    theme_minimal() +
    labs(title = "Mitjana vs. Coeficient de Variació (Top 2000 gens)", x = "Mitjana de l'expressió (Log2)", y = "Coeficient de Variació (CV)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Mostrar gràfic CV
  print(p_cv)
  
  # Generem el histograma amb les dades després del filtratge
  countData_log_melted_filtered <- as.data.frame(countData_log_filtered_top_genes)
  countData_log_melted_filtered$genes <- rownames(countData_log_filtered_top_genes)
  countData_log_melted_filtered <- reshape2::melt(countData_log_melted_filtered, id.vars = "genes")
  
  p_hist_after <- ggplot(countData_log_melted_filtered, aes(x = value)) +
    geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black", alpha = 0.7) +
    theme_minimal() +
    labs(title = "Distribució de les dades log-transforms (Després del filtratge)", x = "Log2(Count + 1)", y = "Frequència")
  
  # Mostrar histograma després del filtratge
  print(p_hist_after)
  
  # Retornar les dades i resultats només amb els top gens seleccionats
  return(list(
    countData = countData,
    countData_log = countData_log,
    countData_log_filtered_top_genes = countData_log_filtered_top_genes,
    metaData_filtered = metaData_filtered,
    p_cv = p_cv,
    p_hist_after = p_hist_after
  ))
}





























# 1. Definir els gens significatius
signif_genes <- rownames(genes_deg_a[genes_deg_a$padj < 0.05, ])

# 2. Mapar SYMBOL → ENTREZID
mapped_genes <- bitr(signif_genes, fromType = "SYMBOL", 
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# 3. Enriquiment GO per processos biològics (BP)
ego <- enrichGO(gene          = mapped_genes$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",  # prova també "MF" o "CC" si vols altres termes GO
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.1,  # més permissiu
                qvalueCutoff  = 0.2,
                readable      = TRUE)

# Comprovar els resultats
head(ego@result)
nrow(ego@result)
summary(ego@result$p.adjust)

# Barplot
barplot(ego_significant, showCategory = 10) + 
  theme(axis.text.y = element_text(size = 8))

# Dotplot
dotplot(ego, showCategory = 10) + 
  theme(axis.text.y = element_text(size = 8))

# 5. Càlcul del Rich Factor
ego3 <- mutate(ego@result, 
               richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

# Gràfic del Rich Factor
ggplot(ego3, aes(richFactor, fct_reorder(Description, richFactor))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(color = p.adjust, size = Count)) +
  scale_color_gradientn(colours = c("#f7ca64", "#46bac2", "#7e62a3"), 
                        trans = "log10",
                        guide = guide_colorbar(reverse = TRUE, order = 1)) +
  scale_size_continuous(range = c(2, 10)) +
  xlab("Rich Factor") +
  labs(y = NULL) +
  ggtitle("Biological Processes")

# 6. Gene Set Enrichment Analysis (GSEA)

geneList <- genes_deg_a$log2FoldChange
names(geneList) <- rownames(genes_deg_a)
geneList <- geneList[order(geneList, decreasing = TRUE)]

ggo <- gseGO(geneList, keyType = "SYMBOL",
             OrgDb = org.Hs.eg.db, pvalueCutoff = 1)
ggo2 <- simplify(ggo, cutoff = 0.7, by = "p.adjust", select_fun = min)

# Similitud de termes
ggo2 <- pairwise_termsim(ggo2)

# 7. Gene-Concept Network

cnetplot(ego, foldChange = geneList, showCategory = 2, 
         cex_label_gene = 0.2, cex_label_category = 0.5)

# 8. Heatplot de GSEA
heatplot(ggo2, foldChange = geneList, showCategory = 10) +
  theme(axis.text.x = element_text(angle = 90, size = 5))

# 9. UpSet plot per veure la superposició de termes GO
upsetplot(ego, n = 5)




# REACTOME
library(ReactomePA)

# Realitza l'anàlisi funcional amb ReactomePA
enrich_result <- enrichPathway(gene = genes_interessants, organism = "human", pvalueCutoff = 0.05)

# Visualitza els resultats amb un gràfic de barres
barplot(enrich_result, showCategory = 10)  # Mostra les 10 vies més significatives

# També pots visualitzar els resultats en una taula
summary(enrich_result)




Change <- clusterProfiler:: bitr(rownames(res_filtered),
                                 fromType = "ENSEMBLE",
                                 toType = "SYMBOL",
                                 OrgDb = org.Hs.eg.db,
                                 drop = TRUE)

rownames(res_filtered) <- Change$SYMBOL

























# import data
pacients_a_eliminar <- c("R032", "R036", "R228", "R233", "R234", "R235", "R246", "R240", "R064")

countData <- read.csv("unproc_pbrna_Alpha.csv", row.names = 1)
countData <- read.csv("unproc_rnaseq.csv", row.names = 1)
countData <- countData[, !colnames(countData) %in% pacients_a_eliminar]


metaData <- read.csv('donor.csv', header = TRUE, sep = ",")

metaData_filtered <- subset(metaData, record_id %in% colnames(countData))
rownames(metaData_filtered) <- metaData_filtered$record_id
metaData_filtered$record_id <- NULL

metaData_filtered <- metaData_filtered[
  metaData_filtered$diagnosis == metaData_filtered$diagnosis_computed &
    metaData_filtered$diagnosis %in% c("None", "Type2"),
]
samples_to_keep <- rownames(metaData_filtered)[metaData_filtered$diagnosis_computed != "Pre.T2D"]
countData <- countData[, samples_to_keep, drop = FALSE]
countData_t <- t(countData)
countData_log <- log2(countData + 1)
countData_ri <- cbind(ID = rownames(countData), countData)
rownames(countData_ri) <- NULL


metaData_filtered <- metaData_filtered[samples_to_keep, , drop = FALSE]


tpmData <- read.csv("TPM_results.csv", row.names = 1)
tpmData <- tpmData[!rownames(tpmData) %in% pacients_a_eliminar, ]
tpmData <- tpmData[, samples_to_keep, drop = FALSE]
tpmData_log <- log2(tpmData + 1)

#------------------------------------------------------------------------------------
# DESeq 2 
metaData_filtered$diagnosis_computed <- as.factor(metaData_filtered$diagnosis_computed)

dds <- DESeqDataSetFromMatrix(countData=countData_ri, 
                              colData=metaData_filtered, 
                              design= ~ diagnosis_computed, tidy = TRUE)
        
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]


## FILTRATGE Padj, Log2FC, expresion
res_padj<- subset(res, padj<0.05)
res_filtered <- subset(res_padj, log2FoldChange < -1 | log2FoldChange > 1)
res_filtered_na <- subset(res_padj, log2FoldChange < -1 | log2FoldChange > 1)
summary(res_filtered)

# Annotate genes
Change <- clusterProfiler:: bitr(rownames(res_filtered),
                                 fromType = "ENTERZID",
                                 toType = "SYMBOL",
                                 OrgDb = org.Hs.eg.db,
                                 drop = TRUE)

rownames(res_filtered) <- Change$SYMBOL

genes_deg_a <- rownames(res_filtered)
genes_deg_na <- rownames(res_filtered_na)

countData_logDEG <- countData_log[genes_deg_na, , drop = FALSE]
tpmData_logDEG <- tpmData_log[genes_deg_a, , drop = FALSE]
tpmData_logDEG<- na.omit(tpmData_logDEG)
#----------------------------------------------------------------------------------------
## VOLCANOPLOT

res_padj$diffexpressed <- "NO"
res_padj$diffexpressed[res_padj$log2FoldChange > 3 & res_padj$pvalue < 0.05] <- "UP"
res_padj$diffexpressed[res_padj$log2FoldChange < -3 & res_padj$pvalue < 0.05] <- "DOWN"

res_padj <- res_padj[!is.na(res_padj$log2FoldChange) & !is.na(res_padj$pvalue), ]

# Gràfic de volcà amb anotació
ggplot(data = res_padj, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  geom_text(data = res_padj, 
            aes(label = rownames(res_padj)), 
            vjust = 1.5, hjust = 0.5, size = 3, check_overlap = TRUE) +
  theme_minimal() +
  ggtitle("Volcano Plot")  

#----------------------------------------------------------------------------------------
## HEATMAP DEG

condition_vector <- metaData$diagnosis_computed[match(colnames(countData_logDEG), metaData$record_id)]
annotation_col <- data.frame(Condition = factor(condition_vector, levels = c("None", "Type2")))
rownames(annotation_col) <- colnames(countData_logDEG)
ann_colors <- list(Condition = c("None" = "lightblue", "Type2" = "pink"))

#COUNTS
pheatmap(countData_logDEG, 
         color = colorRampPalette(c("blue", "white", "red"))(50),  
         cluster_rows = TRUE,  # Manté el clustering de les files (gens)
         cluster_cols = TRUE,  # Manté el clustering de les columnes (pacients)
         scale = "row",    
         fontsize_row = 8,  
         fontsize_col = 10,  
         main = "Beta Heatmap de l'expressió gènica (counts) L2FC 3",
         annotation_col = annotation_col,  
         annotation_colors = ann_colors,
         labels_row = NULL,  # Elimina els labels de les files
         treeheight_row = 0,  # Elimina la línia de clustering de les files
         show_colnames = FALSE,
         show_rownames = FALSE)


## TPM
#rownames(annotation_col) <- colnames(tpmData_logDEG)
#pheatmap(tpmData_logDEG, 
         color = colorRampPalette(c("blue", "white", "red"))(50),  
         cluster_rows = TRUE,  
         cluster_cols = TRUE,  
         scale = "row",    
         fontsize_row = 8,  
         fontsize_col = 10,  
         main = "Heatmap de l'expressió gènica per pacient tpm",
         annotation_col = annotation_col,  
         annotation_colors = ann_colors,
         labels_row = NULL,  
         treeheight_row = 0,  
         show_colnames = FALSE,
         show_rownames = FALSE)  

#---------------------------------------------------------------------------------------
# PCA 

## COUNTS
pca_result <- prcomp(t(countData_logDEG), scale. = TRUE)
pca_data <- as.data.frame(pca_result$x)
pca_data$pacient <- colnames(countData_log)  


pca_data$condicio <- metaData$diagnosis_computed[match(pca_data$pacient, metaData$record_id)]
pca_data$color <- ifelse(pca_data$condicio == "Type2", "#FF69B4", "#1E90FF")  


ggplot(pca_data, aes(x = PC1, y = PC2, color = condicio)) +
  geom_point(size = 4) + 
  scale_color_manual(values = c("Type2" = "#FF69B4", "None" = "#1E90FF")) +  # Assignar colors a les condicions
  labs(title = "PCA dels Pacients counts BETA L2FC 3", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "top")  # Afegir llegenda al top


## TOP 2 GENS
top_2_genes <- head(res_filtered_na[order(abs(res_filtered_na$log2FoldChange), decreasing = TRUE), ], 2)
top_2_genes_names <- rownames(top_2_genes)
countData_2_genes <- countData_logDEG[top_2_genes_names, , drop = FALSE]
countData_2_genes

# HEATMAP 2
condition_vector <- metaData$diagnosis_computed[match(colnames(countData_2_genes), metaData$record_id)]
annotation_col <- data.frame(Condition = factor(condition_vector, levels = c("None", "Type2")))
rownames(annotation_col) <- colnames(countData_2_genes)
ann_colors <- list(Condition = c("None" = "lightblue", "Type2" = "pink"))

pheatmap(countData_2_genes, 
         color = colorRampPalette(c("blue", "white", "red"))(50),  
         cluster_rows = TRUE,  # Manté el clustering de les files (gens)
         cluster_cols = TRUE,  # Manté el clustering de les columnes (pacients)
         scale = "row",    
         fontsize_row = 8,  
         fontsize_col = 10,  
         main = "Heatmap de l'expressió gènica per pacient top2",
         annotation_col = annotation_col,  
         annotation_colors = ann_colors,
         labels_row = NULL,  # Elimina els labels de les files
         treeheight_row = 0,  # Elimina la línia de clustering de les files
         show_colnames = FALSE,
         show_rownames = FALSE)

#PCA 2
pca_result <- prcomp(t(countData_2_genes), scale. = TRUE)
pca_data <- as.data.frame(pca_result$x)
pca_data$pacient <- colnames(countData_2_genes)  # Afegir noms de pacients


pca_data$condicio <- metaData$diagnosis_computed[match(pca_data$pacient, metaData$record_id)]
pca_data$color <- ifelse(pca_data$condicio == "Type2", "#FF69B4", "#1E90FF")  


ggplot(pca_data, aes(x = PC1, y = PC2, color = condicio)) +
  geom_point(size = 4) + 
  scale_color_manual(values = c("Type2" = "#FF69B4", "None" = "#1E90FF")) +  # Assignar colors a les condicions
  labs(title = "PCA dels Pacients counts top2", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "top")  # Afegir llegenda al top


#----------------------------------------------------------
# PCA pacients tpm

pca_result <- prcomp(t(tpmData_log), scale. = TRUE)
pca_data <- as.data.frame(pca_result$x)

pca_data$pacient <- colnames(countData_log)  
pca_data$condicio <- metaData_2$diagnosis_computed[match(pca_data$pacient, metaData_2$record_id)]
pca_data$color <- ifelse(pca_data$condicio == "Type2", "#FF69B4", "#1E90FF")  


ggplot(pca_data, aes(x = PC1, y = PC2, color = condicio)) +
  geom_point(size = 4) + 
  scale_color_manual(values = c("Type2" = "#FF69B4", "None" = "#1E90FF")) +  # Assignar colors a les condicions
  labs(title = "PCA dels Pacients tpm", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "top")  # Afegir llegenda al top

# PCA GENES

pca_result_genes <- prcomp(tpmData_log, scale. = TRUE)
pca_data_genes <- as.data.frame(pca_result_genes$x)
pca_data_genes$gen <- rownames(tpmData_log)  

ggplot(pca_data_genes, aes(x = PC1, y = PC2, label = gen)) +
  geom_point(aes(color = PC1), size = 4) +  # Assignar colors en funció del valor de PC1
  scale_color_gradient(low = "blue", high = "red") +  # Assignar escala de colors (de blau a vermell)
  labs(title = "PCA dels Gens tpm", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "top")  # Afegir llegenda al top

#-----------------------------------------------------------------------------------------------
## Scatterplot dos gens més deg


top_10_gene_names <- rownames(top_10_genes)  
top_10_gene_names[top_2_gene_names %in% rownames(tpmData_log)]
top_10_gene_names
#no estan --> "CCDC197" "FAM230I"
setdiff(rownames(res_filtered), rownames(tpmData_log))

tpmData_top_2_geneslong <- data.frame(
  record_id = rep(colnames(tpmData_top_2_genes), 2),
  expression = c(as.numeric(tpmData_top_2_genes[1, ]), as.numeric(tpmData_top_2_genes[2, ])),  # Assegurar-nos que sigui numèric
  gene = rep(top_2_gene_names, each = ncol(tpmData_top_2_genes)),
  diagnosis = rep(condition_vector_tpm, 2)
)

ggplot(tpmData_top_2_geneslong, aes(x = gene, y = expression, color = diagnosis)) +  
  geom_point(aes(shape = diagnosis), size = 3) +  # Punts per als pacients
  scale_color_manual(values = c("None" = "blue", "Type2" = "pink")) +  # Colors per diagnòstic
  labs(title = paste("Scatterplot dels dos gens més diferencialment expressats: ", 
                     top_2_gene_names[1], "i ", top_2_gene_names[2]),
       x = "Gens", y = "Expressió gènica (log2TPM)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Millorar la visualització de les etiquetes

#  PCA dos gens mes DEG
pca_result_2 <- prcomp(t(tpmData_top_2_genes), scale. = TRUE)

pca_data_2 <- as.data.frame(pca_result_2$x)

pca_data_2$pacient <- colnames(tpmData_top_2_genes)  # Afegir noms de pacients
pca_data_2$condicio <- metaData_2$diagnosis[match(pca_data_2$pacient, metaData_2$record_id)]
pca_data_2$color <- ifelse(pca_data_2$condicio == "Type2", "#FF69B4", "#1E90FF")  

# Crear el gràfic de PCA
ggplot(pca_data_2, aes(x = PC1, y = PC2, color = condicio)) +
  geom_point(size = 4) + 
  scale_color_manual(values = c("Type2" = "#FF69B4", "None" = "#1E90FF")) +  # Assignar colors a les condicions
  labs(title = "PCA dels dos gens més diferencials", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "top")  # Afegir llegenda al top

#----------------------------------------------------------------------------------------          
# Mirar si coincideixen amb els de HI
hi_results <- read.csv("dea_results_2025-03-10T11-43-06.csv")
hi_results

intersect(rownames(res_filtered), hi_results$Feature)

hi_results_adj<- subset(hi_results, Adjusted.p_value<0.05)

intersect<-intersect(rownames(res_filtered), hi_results_adj$Feature)
summary(intersect)

#107

#------------------------------------------------------------------------------------------

## FUNCTIONAL ANALYSIS GO


signif_genes <- rownames(res_filtered)
head(signif_genes)

ego <- enrichGO(gene = signif_genes,
                keyType = "SYMBOL",
                OrgDb = org.Hs.eg.db,
                ont = "BP", # BP: Biological Process
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = F
)
# use simplify() to remove redundant terms
ego2 <- simplify(ego, 
                 cutoff = 0.7,  #llindar del minim redundant
                 by = "p.adjust", 
                 select_fun = min)
# GRÀFICS
barplot(ego2, showCategory = 10) +
  theme(axis.text.y = element_text(size = 8))

dotplot(ego2, showCategory = 10) +
  theme(axis.text.y = element_text(size = 8))

#RICH FACTOR
ego3 <- mutate(ego2, 
               richFactor = Count/as.numeric(
                 sub("/\\d+", "", BgRatio)
               )
)

ggplot(ego3, showCategory = 10,
       aes(richFactor, fct_reorder(Description, richFactor))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(color = p.adjust, size = Count)) +
  scale_color_gradientn(colours = c("#f7ca64", "#46bac2", "#7e62a3"), 
                        trans = "log10",
                        guide = guide_colorbar(reverse = T, order = 1)) +
  scale_size_continuous(range = c(2,10)) +
  #theme_classic(text = element_text(size = 14)) +
  xlab("Rich Factor") +
  labs(y = NULL) +
  ggtitle("Biological Processes")


#Gene Set Enrichment Analysis
geneList <- res_filtered$log2FoldChange
names(geneList) <- rownames(res_filtered)
head(geneList)
geneList <- geneList[order(geneList, decreasing = T)]
ggo <- gseGO(geneList, keyType = "SYMBOL",
             OrgDb = org.Hs.eg.db, pvalueCutoff = 1)
ggo2 <- simplify(ggo, cutoff = 0.7, by = "p.adjust", select_fun = min)

ggo2 <- pairwise_termsim(ggo2)

#Gene-Concept Network
cnetplot(ego2, foldChange = geneList, showCategory = 2,
         cex_label_gene = 0.2,
         cex_label_category = 0.5)

#Heatplot
heatplot(ggo2, foldChange = geneList, showCategory = 10) +
  theme(axis.text.x = element_text(angle = 90, size = 5))

#Upseat plot
upsetplot(ego, n = 5)


# REACTOME
library(ReactomePA)

# Realitza l'anàlisi funcional amb ReactomePA
enrich_result <- enrichPathway(gene = genes_interessants, organism = "human", pvalueCutoff = 0.05)

# Visualitza els resultats amb un gràfic de barres
barplot(enrich_result, showCategory = 10)  # Mostra les 10 vies més significatives

# També pots visualitzar els resultats en una taula
summary(enrich_result)

