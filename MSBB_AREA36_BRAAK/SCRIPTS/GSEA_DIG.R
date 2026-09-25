#V2 - GSEA
# GSEA PARA CENTRALIDADE DIFERENCIAL 

suppressPackageStartupMessages({
  library(vroom)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

dir_out <- "C:/Users/Tannus/Desktop/msbb_lioness/"

# 1. CARREGAMENTO E ANÁLISE LIMMA 
lioness_degree <- readRDS(file.path(dir_out, "lioness_gene_degree_matrix.rds"))
meta <- vroom(file.path(dir_out, "metadata_braak.csv"), delim = ",")

# Filtragem de amostras
meta_sub <- meta %>% filter(BraakGroup %in% c("Low", "High"))
amostras_comuns <- intersect(colnames(lioness_degree), meta_sub$specimenID)

matriz_centralidade <- as.matrix(lioness_degree[, amostras_comuns])
meta <- meta_sub[match(colnames(matriz_centralidade), meta_sub$specimenID), ]

# Transformação log2
matriz_transformada <- log2(matriz_centralidade + 1)

variancias <- apply(matriz_transformada, 1, var)
matriz_filtrada <- matriz_transformada[variancias > 0, ]

# Limma
meta$BraakGroup <- factor(meta$BraakGroup, levels = c("Low", "High")) # Low = ref, High = caso
design <- model.matrix(~ BraakGroup, data = meta)

fit <- lmFit(matriz_filtrada, design)
fit <- eBayes(fit)

resultados <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")

resultados$ensembl_gene_id <- gsub("\\..*", "", rownames(resultados))

# 2. MAPEAMENTO DE IDS  ----------

mapeamento <- bitr(
  resultados$ensembl_gene_id,
  fromType = "ENSEMBL",
  toType   = "SYMBOL",
  OrgDb    = org.Hs.eg.db
)

resultados_com_symbol <- merge(
  resultados, 
  mapeamento, 
  by.x = "ensembl_gene_id", 
  by.y = "ENSEMBL", 
  all.x = TRUE
)

# 3. PREPARAÇÃO DA LISTA RANQUEADA PARA O GSEA ----
gsea_prep <- resultados_com_symbol %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  filter(!is.na(t)) %>%
  arrange(desc(abs(t))) %>%  # Trata duplicatas 
  distinct(SYMBOL, .keep_all = TRUE)

gene_list_centrality <- gsea_prep$t
names(gene_list_centrality) <- gsea_prep$SYMBOL
gene_list_centrality <- sort(gene_list_centrality, decreasing = TRUE)


# 4. EXECUÇÃO DO GSEA 

gsea_centrality_res <- gseGO(
  geneList      = gene_list_centrality,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "CC",         # BP, CC 
  pvalueCutoff  = 1,            
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 500,
  verbose       = FALSE
)

# Exportar resultados em tabela
if (!is.null(gsea_centrality_res) && nrow(gsea_centrality_res@result) > 0) {
  
  file_csv <- file.path(dir_out, "DIG_GSEA_Centrality_GO_CC_results.csv")
  write.csv(as.data.frame(gsea_centrality_res@result), file = file_csv, row.names = FALSE)
  
  
  # Filtrar apenas vias estatisticamente significativas para os gráficos
  vias_sig <- gsea_centrality_res@result %>% filter(p.adjust < 0.05)
  message(paste("vias significativas (p.adjust < 0.05):", nrow(vias_sig)))
  }
# 5. GERAR E SALVAR GRÁFICOS DO GSEA
gsea_sim <- tryCatch({
  simplify(gsea_centrality_res, cutoff = 0.6, by = "p.adjust", select_fun = min)
}, error = function(e) {
  return(gsea_centrality_res)
})

#  ordenação  pelo menor p.adjust
gsea_sim@result <- gsea_sim@result[order(gsea_sim@result$p.adjust), ]

# Plota as top  vias globais
p_dot_centrality <- dotplot(gsea_sim, showCategory = 8) +
  ggtitle("DIG_GSEA Centrality: Gene Ontology - CC (High vs Low Braak)") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title  = element_text(face = "bold", hjust = 0.5)
  )

ggsave(
  filename = file.path(dir_out, "1Dig_GSEA_centrality_dotplot_CC.png"),
  plot     = p_dot_centrality,
  width    = 6.5, 
  height   = 8.5, 
  dpi      = 300
)


#----alternativa de grafico
# 1. dados do DIG 
df_dig <- gsea_centrality_res@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(
    Direction  = ifelse(NES > 0, "Increased centrality", "Decreased centrality"),
    log10_padj = -log10(p.adjust)
  )

# Top 10 vias com maior ganho de centralidade (NES > 0)
top_up_dig <- df_dig %>%
  filter(Direction == "Increased centrality") %>%
  arrange(p.adjust) %>%
  head(10)

# Top 10 vias com maior perda de centralidade (NES < 0)
top_down_dig <- df_dig %>%
  filter(Direction == "Decreased centrality") %>%
  arrange(p.adjust) %>%
  head(10)

# ajustar eixos direcionais
plot_data_dig <- bind_rows(top_up_dig, top_down_dig) %>%
  mutate(
    # Barras Decreased vão para a esquerda (-) e Increased para a direita (+)
    plot_val = ifelse(Direction == "Increased centrality", log10_padj, -log10_padj),
    Description_wrap = str_wrap(Description, width = 45)
  )

# Increased no topo e Decreased embaixo
plot_data_dig <- plot_data_dig %>%
  arrange(plot_val) %>%
  mutate(Description_wrap = factor(Description_wrap, levels = Description_wrap))

# 2. Definir limites simétricos para o eixo X
max_limit_dig <- ceiling(max(abs(plot_data_dig$plot_val))) * 1.4

# 3. Construir o gráfico
p_diverging_dig <- ggplot(plot_data_dig, aes(x = plot_val, y = Description_wrap, fill = Direction)) +
  # Barras
  geom_col(width = 0.75) +
  
  # Textos das vias posicionados no lado oposto de cada barra
  geom_text(
    aes(
      x = ifelse(Direction == "Increased centrality", -0.4, 0.4),
      label = Description_wrap,
      hjust = ifelse(Direction == "Increased centrality", 1, 0)
    ),
    size = 3.8,
    color = "black"
  ) +
  
  scale_fill_manual(
    values = c(
      "Increased centrality" = "#8B1A1A", 
      "Decreased centrality" = "#275396"
    )
  ) +
  
  # Eixo X com valores absolutos 
  scale_x_continuous(
    labels = abs,
    limits = c(-max_limit_dig, max_limit_dig),
    breaks = seq(-30, 30, by = 5)
  ) +
  
  # Linha vertical central divisória
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  
  labs(
    x = "-log10(P.adjust)",
    y = NULL,
    fill = NULL
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.line          = element_blank(),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.ticks.x       = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x        = element_text(size = 11, color = "black", face = "bold"),
    axis.title.x       = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    
    # Legenda no canto superior esquerdo
    legend.position    = c(0.18, 0.90),
    legend.background  = element_blank(),
    legend.text        = element_text(size = 10, face = "bold"),
    legend.key.size    = unit(0.6, "cm")
  )

# 4. Salvar
ggsave(
  filename = file.path(dir_out, "dig_diverging_bar_plot.png"),
  plot     = p_diverging_dig,
  width    = 9.0,
  height   = 7.5,
  dpi      = 300,
  units    = "in"
)
