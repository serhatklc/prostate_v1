############################################################
# Full Bioinformatics Workflow (per manuscript)
# Datasets: GSE3325, GSE6919, GSE55945, GSE26910 (discovery), GSE46602 (validation)
# Steps:
# 1) Download GEO data & extract ExpressionSets
# 2) Probe -> Gene Symbol mapping (platform-aware) and collapse to genes
# 3) Merge datasets on common genes
# 4) Log2 transform (if needed) + Quantile normalization
# 5) Batch correction with ComBat
# 6) Missing value imputation (kNN)
# 7) Quality Control (boxplot/density/PCA before & after)
# 8) Differential Expression (limma, Tumor vs Normal) with BH-FDR
# 9) Optional housekeeping filter (e.g., ACTB/GAPDH) sensitivity check
# 10) Enrichment analyses (GO, KEGG) via clusterProfiler
# 11) Save machine-learning ready matrices: combined expression + labels
# 12) Process validation set (GSE46602) with identical steps and save
############################################################

# ----------------------------
# Package setup
# ----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("GEOquery","limma","sva","impute","WGCNA","clusterProfiler",
          "org.Hs.eg.db","ggplot2","RColorBrewer")
for(p in pkgs){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE)))){
    BiocManager::install(p, ask = FALSE, update = FALSE)
    library(p, character.only = TRUE)
  }
}

# ----------------------------
# Helpers
# ----------------------------
safelog2 <- function(m){
  # If already roughly log-scale (median < 100), skip; else log2(x+1)
  if(median(m, na.rm = TRUE) > 100){
    return(log2(m + 1))
  } else {
    return(m)
  }
}

pick_symbol_col <- function(fd){
  # Try common column names for Gene Symbol in GEO platform featureData
  cand <- c("Gene Symbol","Gene symbol","GENE_SYMBOL","Symbol","SYMBOL","gene_assignment","GeneSymbol","Entrez Gene Symbol")
  has <- cand[cand %in% colnames(fd)]
  if(length(has) == 0) return(NULL)
  return(has[1])
}

extract_symbol_from_assignment <- function(x){
  # For fields like "GENE1 // description // ..."
  if(is.na(x) || x=="") return(NA_character_)
  sp <- unlist(strsplit(as.character(x), " // |///|;|\\|", perl=TRUE))
  if(length(sp)==0) return(NA_character_)
  return(trimws(sp[1]))
}

collapse_to_genes <- function(eset){
  ex <- exprs(eset)
  fd <- fData(eset)
  symcol <- pick_symbol_col(fd)
  if(is.null(symcol)){
    # fallback: try to parse from gene_assignment-like field
    if("gene_assignment" %in% colnames(fd)){
      syms <- vapply(fd$gene_assignment, extract_symbol_from_assignment, character(1))
    } else {
      stop("No gene symbol column found in platform feature data.")
    }
  } else {
    syms <- fd[[symcol]]
  }
  # Clean symbols
  syms <- toupper(trimws(as.character(syms)))
  syms[syms=="" | syms=="NA"] <- NA
  # Remove rows without symbols
  keep <- !is.na(syms)
  ex <- ex[keep, , drop=FALSE]
  syms <- syms[keep]
  # For probes mapping to multiple symbols like "GENE1 /// GENE2" keep first
  syms <- vapply(syms, function(s){
    sp <- unlist(strsplit(s, "///|;|,|\\s+/\\s+|\\|", perl=TRUE))
    if(length(sp)==0) return(NA_character_)
    trimws(sp[1])
  }, character(1))
  keep <- !is.na(syms) & syms!="" & syms!="NA"
  ex <- ex[keep, , drop=FALSE]; syms <- syms[keep]
  # Collapse duplicate symbols by highest IQR (probe with best variability)
  # Compute IQR per probe
  iqrs <- apply(ex, 1, IQR, na.rm=TRUE)
  ord <- order(syms, -iqrs)  # sort by symbol, then descending IQR
  ex <- ex[ord, , drop=FALSE]
  syms <- syms[ord]
  # Keep first occurrence of each symbol
  dup <- duplicated(syms)
  ex_unique <- ex[!dup, , drop=FALSE]
  rownames(ex_unique) <- syms[!dup]
  return(ex_unique)
}

detect_labels <- function(ph){
  # Try to infer Tumor/Normal from common phenotype fields
  txt <- paste(apply(ph, 2, function(col) paste(col, collapse = " ")), collapse = " ")
  # This is just for pattern search scope; real decision uses each sample's rows
  keys <- tolower(apply(ph, 2, as.character))
  keys[is.na(keys)] <- ""
  # per-sample decision
  lab <- apply(keys, 1, function(row){
    s <- paste(row, collapse=" ")
    s <- tolower(s)
    if(grepl("tumou?r|prostate cancer|adenocarcinoma|cancerous|primary tumor|malignant", s)) return("Tumor")
    if(grepl("normal|benign|adjacent normal|non-tumou?r|noncancer", s)) return("Normal")
    return(NA_character_)
  })
  return(lab)
}

qc_plots <- function(mat, prefix){
  # Boxplot
  png(sprintf("%s_boxplot.png", prefix), width=1200, height=700)
  par(mar=c(8,5,2,1))
  boxplot(mat, las=2, outline=FALSE, col=brewer.pal(8,"Set3"), main=paste0(prefix, " - Boxplot"))
  dev.off()
  # Density
  png(sprintf("%s_density.png", prefix), width=1000, height=700)
  matplot(density(mat[,1])$x, t(apply(mat, 2, function(x) density(x, na.rm=TRUE)$y)),
          type="l", lty=1, main=paste0(prefix," - Density"), xlab="Expression", ylab="Density")
  dev.off()
  # PCA
  png(sprintf("%s_pca.png", prefix), width=900, height=800)
  pca <- prcomp(t(mat), scale.=TRUE)
  plot(pca$x[,1], pca$x[,2], pch=19, col="steelblue", xlab=sprintf("PC1 (%.1f%%)", 100*summary(pca)$importance[2,1]),
       ylab=sprintf("PC2 (%.1f%%)", 100*summary(pca)$importance[2,2]), main=paste0(prefix," - PCA"))
  dev.off()
}

# ----------------------------
# 1) Download discovery datasets and preprocess each
# ----------------------------
discovery_ids <- c("GSE3325","GSE6919","GSE55945","GSE26910")
expr_list <- list()
labels_list <- list()
batch_vec <- c()

for(id in discovery_ids){
  message("Downloading ", id)
  gset <- getGEO(id, GSEMatrix = TRUE)
  eset <- gset[[1]]
  ex <- exprs(eset)
  # Probe -> Gene mapping & collapse
  ex_gene <- collapse_to_genes(eset)
  # Ensure numeric matrix
  storage.mode(ex_gene) <- "double"
  # Save raw QC per dataset
  qc_plots(ex_gene, paste0(id,"_RAW_gene"))
  # Keep
  expr_list[[id]] <- ex_gene
  # Labels
  ph <- pData(eset)
  labs <- detect_labels(ph)
  # Fallback to characteristics_ch1.1 etc if needed
  if(all(is.na(labs)) && "characteristics_ch1" %in% colnames(ph)){
    labs <- ifelse(grepl("tumou?r|cancer", tolower(ph$characteristics_ch1)), "Tumor",
                   ifelse(grepl("normal|benign", tolower(ph$characteristics_ch1)), "Normal", NA))
  }
  labels_list[[id]] <- labs
  batch_vec <- c(batch_vec, rep(id, ncol(ex_gene)))
}

# ----------------------------
# 2) Intersect common genes & merge
# ----------------------------
common_genes <- Reduce(intersect, lapply(expr_list, rownames))
expr_merged <- do.call(cbind, lapply(expr_list, function(m) m[common_genes, , drop=FALSE]))

# safety: replace infs/nas prior to log/normalization
expr_merged[!is.finite(expr_merged)] <- NA

# ----------------------------
# 3) Log2 (if needed) + Quantile normalization
# ----------------------------
expr_log <- safelog2(expr_merged)
qc_plots(expr_log, "DISCOVERY_preNorm")
expr_norm <- normalizeBetweenArrays(expr_log, method="quantile")
qc_plots(expr_norm, "DISCOVERY_postNorm")

# ----------------------------
# 4) Batch correction (ComBat)
# ----------------------------
batch <- factor(batch_vec)
expr_combat <- ComBat(dat = as.matrix(expr_norm), batch = batch, par.prior = TRUE, prior.plots = FALSE)
qc_plots(expr_combat, "DISCOVERY_postComBat")

# ----------------------------
# 5) Missing value imputation (kNN)
# ----------------------------
expr_imputed <- impute.knn(expr_combat)$data
qc_plots(expr_imputed, "DISCOVERY_postImpute")

# ----------------------------
# 6) Build label vector
# ----------------------------
labels <- unlist(labels_list, use.names = FALSE)
if(length(labels) != ncol(expr_imputed)){
  warning("Label length mismatch. Unlabeled samples will be set to NA and removed for DEG/ML.")
}
labels <- as.factor(labels)
keep_samples <- !is.na(labels) & labels %in% c("Tumor","Normal")
expr_final <- expr_imputed[, keep_samples, drop=FALSE]
labels_final <- droplevels(labels[keep_samples])

# Save ML-ready matrices
write.csv(t(expr_final), "combined_expression_matrix.csv") # samples x genes
write.csv(data.frame(class = ifelse(labels_final=="Tumor",1,0)), "labels.csv", row.names = FALSE)

# ----------------------------
# 7) DEG analysis with limma (Tumor vs Normal) including batch as covariate
# ----------------------------
batch_final <- batch[keep_samples]
design <- model.matrix(~ 0 + labels_final + batch_final)
colnames(design) <- gsub("labels_final","",colnames(design))
fit <- lmFit(expr_final, design)
contrast.matrix <- makeContrasts(Tumor - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
write.csv(deg, "DEG_full_results.csv")

# Significant set & volcano
deg_sig <- subset(deg, abs(logFC) > 1 & adj.P.Val <= 0.05)
write.csv(deg_sig, "DEG_significant.csv")

png("Volcano_DEG.png", width=1000, height=800)
plot(deg$logFC, -log10(deg$adj.P.Val), pch=20, col=rgb(0.7,0,0,0.5),
     xlab="log2 Fold Change", ylab="-log10(FDR)", main="Volcano Plot (Tumor vs Normal)")
abline(v=c(-1,1), lty=2, col="gray40"); abline(h=-log10(0.05), lty=2, col="gray40")
dev.off()

# ----------------------------
# 8) OPTIONAL: Housekeeping sensitivity (remove ACTB/GAPDH & recompute)
# ----------------------------
hk <- c("ACTB","GAPDH","RPLP0","B2M","HPRT1")
expr_nohk <- expr_final[!(rownames(expr_final) %in% hk), , drop=FALSE]
fit_hk <- lmFit(expr_nohk, design)
fit_hk <- eBayes(contrasts.fit(fit_hk, contrast.matrix))
deg_hk <- topTable(fit_hk, coef = 1, number = Inf, adjust.method = "BH")
write.csv(deg_hk, "DEG_without_housekeeping.csv")

# ----------------------------
# 9) Enrichment (GO, KEGG)
# ----------------------------
gene_symbols <- rownames(deg_sig)
entrez <- bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego <- enrichGO(gene=entrez$ENTREZID, OrgDb=org.Hs.eg.db, ont="ALL", pAdjustMethod="BH",
                qvalueCutoff=0.05, readable=TRUE)
ekegg <- enrichKEGG(gene=entrez$ENTREZID, organism='hsa', pAdjustMethod="BH", qvalueCutoff=0.05)
write.csv(as.data.frame(ego), "GO_enrichment.csv")
write.csv(as.data.frame(ekegg), "KEGG_enrichment.csv")

# ----------------------------
# 10) Process Validation: GSE46602 (same steps)
# ----------------------------
message("Processing validation dataset GSE46602")
gset_val <- getGEO("GSE46602", GSEMatrix = TRUE)
eset_val <- gset_val[[1]]
ex_val_gene <- collapse_to_genes(eset_val)
storage.mode(ex_val_gene) <- "double"
qc_plots(ex_val_gene, "GSE46602_RAW_gene")

# Align genes to discovery common set
common_to_train <- intersect(rownames(expr_final), rownames(ex_val_gene))
ex_val_gene <- ex_val_gene[common_to_train, , drop=FALSE]
ex_val <- safelog2(ex_val_gene)
ex_val <- normalizeBetweenArrays(ex_val, method="quantile")
# Note: No ComBat across discovery+validation (keep validation independent!)
ex_val <- impute.knn(ex_val)$data
qc_plots(ex_val, "GSE46602_postNormImpute")

# Validation labels
ph_val <- pData(eset_val)
labs_val <- detect_labels(ph_val)
if(all(is.na(labs_val)) && "characteristics_ch1" %in% colnames(ph_val)){
  labs_val <- ifelse(grepl("tumou?r|cancer", tolower(ph_val$characteristics_ch1)), "Tumor",
                     ifelse(grepl("normal|benign", tolower(ph_val$characteristics_ch1)), "Normal", NA))
}
labs_val <- as.factor(labs_val)
keep_val <- !is.na(labs_val) & labs_val %in% c("Tumor","Normal")
ex_val <- ex_val[, keep_val, drop=FALSE]
labs_val <- droplevels(labs_val[keep_val])

# Save ML-ready validation matrices (samples x genes)
write.csv(t(ex_val), "validation_expression_matrix.csv")
write.csv(data.frame(class = ifelse(labs_val=="Tumor",1,0)), "validation_labels.csv", row.names = FALSE)

message("Discovery and Validation preprocessing complete.")
