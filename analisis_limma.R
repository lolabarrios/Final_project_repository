# =============================================================================
# Análisis de expresión diferencial con limma (R/Bioconductor)
# Validación cruzada del pipeline Python para el TFG de Lola Barrios Naranjo
# Dataset: GSE28735 (Pancreatic ductal adenocarcinoma, Affymetrix Hu Gene 1.0 ST)
#
# Plataforma recomendada para ejecución:
#   - Posit Cloud (https://posit.cloud) — gratuito, sin instalación local
#   - RStudio Desktop con R 4.0+
# =============================================================================

# -----------------------------------------------------------------------------
# 1. INSTALACIÓN DE PAQUETES (solo si no están)
# -----------------------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

required_packages <- c("limma", "GEOquery", "hugene10sttranscriptcluster.db")
for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
        library(pkg, character.only = TRUE)
    }
}
library(limma)
library(GEOquery)
library(hugene10sttranscriptcluster.db)

# -----------------------------------------------------------------------------
# 2. DESCARGA DEL DATASET DESDE GEO
# -----------------------------------------------------------------------------
gse_id <- "GSE28735"
cat("Descargando", gse_id, "desde GEO...\n")
gse <- getGEO(gse_id, GSEMatrix = TRUE, destdir = ".")
eset <- gse[[1]]
expr_matrix <- exprs(eset)
pheno       <- pData(eset)
cat("Matriz de expresión:", dim(expr_matrix)[1], "sondas x", dim(expr_matrix)[2], "muestras\n")

# -----------------------------------------------------------------------------
# 3. EXTRACCIÓN DE METADATOS (paciente y condición)
# -----------------------------------------------------------------------------
titles <- as.character(pheno$title)

# Condición: Tumor vs Normal
condition <- ifelse(grepl("nontumor", titles, ignore.case = TRUE), "Normal", "Tumor")
condition <- factor(condition, levels = c("Normal", "Tumor"))

# Paciente: extraemos el ID numérico al final del título
patient_id <- sub(".*?(\\d+)$", "\\1", titles)
patient    <- factor(patient_id)

cat("Distribución por condición:\n"); print(table(condition))
cat("Pacientes únicos:", nlevels(patient), "\n")

# -----------------------------------------------------------------------------
# 4. MATRIZ DE DISEÑO: BLOQUE POR PACIENTE + EFECTO DE CONDICIÓN
# -----------------------------------------------------------------------------
# Este es el diseño correcto para datos pareados
design <- model.matrix(~ patient + condition)
cat("Matriz de diseño:", dim(design)[1], "muestras x", dim(design)[2], "coeficientes\n")

# -----------------------------------------------------------------------------
# 5. AJUSTE DEL MODELO LINEAL CON SHRINKAGE BAYESIANO (eBayes)
# -----------------------------------------------------------------------------
fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)

# Tabla completa de DEGs ordenados por p-valor ajustado
deg_table <- topTable(fit,
                       coef          = "conditionTumor",
                       number        = Inf,
                       adjust.method = "BH",
                       sort.by       = "P")

cat("Sondas testadas:", nrow(deg_table), "\n")

# -----------------------------------------------------------------------------
# 6. ANOTACIÓN DE SONDAS A SÍMBOLOS GÉNICOS HGNC
# -----------------------------------------------------------------------------
deg_table$ProbeID <- rownames(deg_table)
gene_symbols <- mapIds(hugene10sttranscriptcluster.db,
                        keys      = deg_table$ProbeID,
                        column    = "SYMBOL",
                        keytype   = "PROBEID",
                        multiVals = "first")
deg_table$GeneSymbol <- gene_symbols

# Filtramos sondas sin anotación
deg_filtered <- deg_table[!is.na(deg_table$GeneSymbol), ]

# Colapsamos por gen, conservando la sonda más significativa por gen
deg_sorted <- deg_filtered[order(deg_filtered$adj.P.Val), ]
deg_unique <- deg_sorted[!duplicated(deg_sorted$GeneSymbol), ]

cat("Genes únicos tras anotación y colapso:", nrow(deg_unique), "\n")

# -----------------------------------------------------------------------------
# 7. EXPORTACIÓN DE RESULTADOS
# -----------------------------------------------------------------------------
output_cols <- c("GeneSymbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
write.csv(deg_unique[, output_cols],
           file      = "DEGs_limma_R.csv",
           row.names = FALSE)
cat("Resultados guardados en: DEGs_limma_R.csv\n")

# -----------------------------------------------------------------------------
# 8. SANITY CHECK: GENES CANÓNICOS DE PDAC EN AMBOS UMBRALES
# -----------------------------------------------------------------------------
genes_canonicos <- c("KRAS","TP53","CDKN2A","SMAD4","MYC","BRCA2",
                      "MUC1","S100A4","S100P","MMP7","LOXL2",
                      "KRT19","REG1A","CEACAM5","CEACAM6",
                      "DPEP1","TPX2","PRR11","DCBLD2","ANLN",
                      "CDO1","C7","ALDH1A1","NR3C2")

# DEGs según umbrales
strict     <- subset(deg_unique, adj.P.Val < 0.05 & abs(logFC) > 1.0)
permissive <- subset(deg_unique, adj.P.Val < 0.05 & abs(logFC) > 0.585)

cat("\n=== RESUMEN DE DEGs ===\n")
cat("  Estricto  (|logFC|>1.0, FDR<0.05):", nrow(strict),     "DEGs\n")
cat("  Permisivo (|logFC|>0.585, FDR<0.05):", nrow(permissive), "DEGs\n")

cat("\n=== COBERTURA DE GENES CANÓNICOS PDAC ===\n")
cat("  Estricto  :",
    paste(intersect(genes_canonicos, strict$GeneSymbol), collapse = ", "), "\n")
cat("  Permisivo :",
    paste(intersect(genes_canonicos, permissive$GeneSymbol), collapse = ", "), "\n")

cat("\n=== TABLA DETALLADA DE GENES CANÓNICOS ===\n")
canonicos_data <- subset(deg_unique, GeneSymbol %in% genes_canonicos)
print(canonicos_data[, c("GeneSymbol", "logFC", "P.Value", "adj.P.Val")])

# -----------------------------------------------------------------------------
# 9. COMPARACIÓN CON RESULTADOS DE PYTHON (opcional)
# -----------------------------------------------------------------------------
# Si has exportado los DEGs del cuaderno Python a un CSV, descomenta:
#
# python_degs <- read.csv("DEGs_python_permissive.csv")  # del cuaderno v5
# overlap_strict     <- intersect(strict$GeneSymbol,     python_degs$GeneSymbol)
# overlap_permissive <- intersect(permissive$GeneSymbol, python_degs$GeneSymbol)
# cat("\n=== SOLAPAMIENTO PYTHON ∩ R limma ===\n")
# cat("  Estricto  :", length(overlap_strict),     "/", length(strict$GeneSymbol),     "\n")
# cat("  Permisivo :", length(overlap_permissive), "/", length(permissive$GeneSymbol), "\n")

# -----------------------------------------------------------------------------
# 10. INFORMACIÓN DE SESIÓN (reproducibilidad)
# -----------------------------------------------------------------------------
cat("\n=== sessionInfo() ===\n")
sessionInfo()
