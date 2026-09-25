# PIPELINE: DESeq2 + GSEA

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(vroom)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(msigdbr)
  library(enrichplot)
})

# 1. CAMINHOS E DADOS

path_metadata <- "C:/Users/Tannus/Desktop/msbb_lioness/metadata_braak.csv"
path_counts   <- "C:/Users/Tannus/Desktop/IC LUISA/synapse/msbb_filtered_counts_(greater_than_1cpm).tsv"
path_14k_txt  <- "C:/Users/Tannus/Desktop/msbb_lioness/lista_genes_selecionados.txt" 

dir_out <- "C:/Users/Tannus/Desktop/msbb_lioness/"

# Metadata (Apenas grupos High e Low)
metadata <- vroom(path_metadata, delim = ",") %>%
  filter(BraakGroup %in% c("High", "Low")) %>%
  as.data.frame()

rownames(metadata) <- metadata$specimenID

# Lista de ~14k Genes do LIONESS
genes_14k <- readLines(path_14k_txt) %>%
  trimws()
genes_14k <- genes_14k[genes_14k != ""]
genes_14k <- gsub("\\..*", "", genes_14k) # Remove versão de transcript (.X)

# Matriz de counts 
gene_count <- read.csv(path_counts, sep = "\t", check.names = FALSE) %>%
  column_to_rownames(var = "feature")

rownames(gene_count) <- gsub("\\..*", "", rownames(gene_count))

# 2. ALINHAMENTO E FILTRAGEM GENES DO LIONESS

amostras_comuns <- intersect(colnames(gene_count), rownames(metadata))
gene_count <- gene_count[, amostras_comuns, drop = FALSE]
metadata   <- metadata[amostras_comuns, , drop = FALSE]

gene_count_14k <- gene_count[rownames(gene_count) %in% genes_14k, ]

# 3. DESeq2

metadata$BraakGroup <- factor(metadata$BraakGroup)
metadata$BraakGroup <- relevel(metadata$BraakGroup, ref = "Low") # "Low" é a referência (Controle)

dds <- DESeqDataSetFromMatrix(
  countData = gene_count_14k,
  colData   = metadata,
  design    = ~ BraakGroup
)

dds <- DESeq(dds)

# Extrair resultados (Contraste: High vs Low)
res_dge <- results(dds, contrast = c("BraakGroup", "High", "Low"))
res_dge_df <- as.data.frame(res_dge)

# 4.LISTA RANQUEADA PARA O GSEA

res_dge_df$ENSEMBL <- rownames(res_dge_df)

# Mapear Ensembl IDs para Gene Symbols
res_dge_df$symbol <- mapIds(
  x         = org.Hs.eg.db,
  keys      = res_dge_df$ENSEMBL,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

# Adicionar a estatística Wald ('stat') como métrica de ranking
res_dge_df <- res_dge_df %>%
  mutate(
    padj = ifelse(is.na(padj), 1, padj),
    rank_metric = ifelse(!is.na(stat), stat, log2FoldChange),
    signif = case_when(
      padj <= 0.05 & log2FoldChange > 0 ~ "Upregulated",
      padj <= 0.05 & log2FoldChange < 0 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

# Salvar a tabela completa do DEG 
write.csv(res_dge_df, file = file.path(dir_out, "DEG_results_MSBB_14k_complete.csv"), row.names = FALSE)

# Criar vetor numérico ordenado para o GSEA
gsea_data <- res_dge_df %>%
  filter(!is.na(symbol) & symbol != "") %>%
  filter(!is.na(rank_metric)) %>%
  arrange(desc(abs(rank_metric))) %>%
  distinct(symbol, .keep_all = TRUE) # Trata Gene Symbols duplicados mantendo o maior valor absoluto

gene_list <- gsea_data$rank_metric
names(gene_list) <- gsea_data$symbol
gene_list <- sort(gene_list, decreasing = TRUE)

# 5. EXECUÇÃO DO GSEA 

gsea_res <- gseGO(
  geneList      = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",     
  ont           = "CC",         
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

# 6. Gerar grafico e salvar 

#simplificar termos
# gsea_res <- simplify(gsea_res, cutoff = 0.7, by = "p.adjust", select_fun = min)

# Exportar tabela de resultados do GSEA
if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
  write.csv(
    as.data.frame(gsea_res@result), 
    file = file.path(dir_out, "DEG_GSEA_GO_CC_results.csv"), 
    row.names = FALSE
  )
} else {
  cat("Nenhuma via atingiu o corte pvalueCutoff = 0.05.\n")
}

if (nrow(gsea_res@result) > 0) {
  
  # 1.  ordenação  pelo menor p.adjust 
  gsea_res@result <- gsea_res@result[order(gsea_res@result$p.adjust), ]
  
  # 2. Criar o Dotplot vertical global (top 10 sem separar por sinal)
  p_dot <- dotplot(gsea_res, showCategory = 10) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 28)) + # quebra termos longos
    ggtitle("DEG GSEA: Gene Ontology - CC (High vs Low Braak)") +
    theme_bw(base_size = 14) + 
    theme(
      axis.text.y   = element_text(size = 11, color = "black"),
      axis.text.x   = element_text(size = 11, color = "black"),
      axis.title    = element_text(size = 12, face = "bold"),
      plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.title  = element_text(size = 11, face = "bold"),
      legend.text   = element_text(size = 10),
      legend.box    = "vertical",
      legend.margin = margin(l = 2, r = 2)
    )
  
  # 3. Salvar em formato vertical 
  ggsave(
    filename = file.path(dir_out, "deg_GSEA_dotplot - CC.png"),  
    plot     = p_dot, 
    width    = 6.5,   
    height   = 8.5,   
    units    = "in",
    dpi      = 300
  )
  
  message("--> Dotplot de DEG salvo com sucesso no padrão vertical!")
} else {
  message("--> Nenhuma via significativa encontrada para DEG.")
}

#----------alternativa grafico 
# 1. Preparar e filtrar os top 10 Up e top 10 Down por p.adjust
df_gsea <- gsea_res@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(
    Direction = ifelse(NES > 0, "Up-regulated pathways", "Down-regulated pathways"),
    log10_padj = -log10(p.adjust)
  )

# Top 10 de cada sentido
top_up <- df_gsea %>%
  filter(Direction == "Up-regulated pathways") %>%
  arrange(p.adjust) %>%
  head(10)

top_down <- df_gsea %>%
  filter(Direction == "Down-regulated pathways") %>%
  arrange(p.adjust) %>%
  head(10)

# Unificar e configurar valores direcionais para o gráfico
plot_data <- bind_rows(top_up, top_down) %>%
  mutate(
    # Barras Down vão para a esquerda (-) e Up para a direita (+)
    plot_val = ifelse(Direction == "Up-regulated pathways", log10_padj, -log10_padj),
    # Quebra nomes muito longos se necessário
    Description_wrap = str_wrap(Description, width = 45)
  )

# Ordenação para posicionar os Up no topo e Down embaixo
plot_data <- plot_data %>%
  arrange(plot_val) %>%
  mutate(Description_wrap = factor(Description_wrap, levels = Description_wrap))

# 2. Definir limites simétricos para o eixo X
max_limit <- ceiling(max(abs(plot_data$plot_val))) * 1.4

# 3. Construir o gráfico
p_diverging <- ggplot(plot_data, aes(x = plot_val, y = Description_wrap, fill = Direction)) +
  # Barras
  geom_col(width = 0.75) +
  
  # Textos das vias posicionados no lado oposto de cada barra
  geom_text(
    aes(
      x = ifelse(Direction == "Up-regulated pathways", -0.4, 0.4),
      label = Description_wrap,
      hjust = ifelse(Direction == "Up-regulated pathways", 1, 0)
    ),
    size = 3.8,
    color = "black"
  ) +
  
  # Cores personalizadas (Vermelho escuro e Azul escuro idênticos ao exemplo)
  scale_fill_manual(
    values = c(
      "Up-regulated pathways"   = "#8B1A1A", 
      "Down-regulated pathways" = "#275396"
    )
  ) +
  
  # Eixo X com valores absolutos (positivos de ambos os lados)
  scale_x_continuous(
    labels = abs,
    limits = c(-max_limit, max_limit),
    breaks = seq(-30, 30, by = 5)
  ) +
  
  # Linha vertical central divisória
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  
  labs(
    x = "-log10(P.adjust)",
    y = NULL,
    fill = NULL
  ) +
  
  # Tema com borda preta em volta e fundo limpo
  theme_classic(base_size = 13) +
  theme(
    # Borda externa preta
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.line          = element_blank(),
    
    # Oculta texto e ticks do eixo Y tradicional (pois usamos geom_text)
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

# 4. Salvar com proporções adequadas
ggsave(
  filename = file.path(dir_out, "deg_diverging_bar_plot.png"),
  plot     = p_diverging,
  width    = 9.0,
  height   = 7.5,
  dpi      = 300,
  units    = "in"
)

message("--> Gráfico divergente salvo com sucesso!")