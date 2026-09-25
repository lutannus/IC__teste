library(lionessR)
library(DESeq2)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(vroom)
library(dplyr)
library(tibble)
library(parallel)

# PARAMETRIZAÇÃO
CONFIG <- list(

  #Entrada
  caminho_metadata     = "/home/ltlima/dig_ad/rnaseq_harmonization_msbb_combined_metadata.csv",
  caminho_cpm          = "/home/ltlima/dig_ad/msbb_filtered_counts_(greater_than_1cpm).tsv",
  caminho_string_alias = "/home/ltlima/dig_ad/9606.protein.aliases.v12.0.txt",
  caminho_string_links = "/home/ltlima/dig_ad/9606.protein.links.detailed.v12.0.txt",
 
  # Saida
  diretorio_resultados = "/home/ltlima/dig_ad/results",

  # Filtros
  area        = "36",
  n_amostras_grupo     = 80,
  semente_random       = 73,
  n_genes_variaveis    = 3
)

cat(">>> Verificando e criando saída de resultados...\n")
if(!dir.exists(CONFIG$diretorio_resultados)) {
  dir.create(CONFIG$diretorio_resultados, recursive = TRUE)
  cat("    Saída criada em:", CONFIG$diretorio_resultados, "\n")
} else {
  cat("    Saída já existente.\n")
}

# CARREGAMENTO DE ARQUIVOS
cat("\n--- [PASSO 1] Carregando arquivos  ---\n")
df <- vroom(CONFIG$caminho_metadata, show_col_types = FALSE)
cpm_matrix <- vroom(CONFIG$caminho_cpm, show_col_types = FALSE)

cat(" [MONITORAMENTO] Metadados brutos:", nrow(df), "linhas.\n")
cat(" [MONITORAMENTO] Matriz CPM bruta:", nrow(cpm_matrix), "genes x", ncol(cpm_matrix) - 1, "amostras.\n")

# Limpeza inicial dos nomes das colunas e unificação de IDs Ensembl
colnames(cpm_matrix) <- trimws(colnames(cpm_matrix))
cpm_matrix <- cpm_matrix %>%
  mutate(feature = sub("\\..*", "", feature)) %>%
  distinct(feature, .keep_all = TRUE) %>%
  column_to_rownames("feature")

# FILTRAGEM DE AMOSTRAS VÁLIDAS E CATEGORIZAÇÃO BRAAK
cat("\n--- [PASSO 2] Alinhando Metadados e Categorizando Braak ---\n")
coluna_id_metadado <- "specimenID"

# Assegurar que a coluna 'braak' seja numérica
df <- df %>% mutate(braak = as.numeric(Braak))

df_filtrado <- df %>% 
  filter(is.na(exclude) | exclude != TRUE) %>%
  filter(grepl(CONFIG$area, tolower(BrodmannArea))) %>%
  filter(.data[[coluna_id_metadado]] %in% colnames(cpm_matrix)) %>%
  filter(!is.na(braak)) %>%
  # MODIFICAÇÃO BRAAK: Definição dos grupos Low, Mid e High
  mutate(BraakGroup = case_when(
    braak >= 0 & braak <= 2 ~ "Low",
    braak >= 3 & braak <= 4 ~ "Mid",
    braak >= 5 & braak <= 6 ~ "High",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(BraakGroup))

cat(" [MONITORAMENTO] Amostras válidas com estágio Braak definido:", nrow(df_filtrado), "\n")
print(table(df_filtrado$BraakGroup))

# Sorteio equilibrado por grupo de Braak
set.seed(CONFIG$semente_random)

df_low  <- df_filtrado %>% filter(BraakGroup == "Low")  %>% sample_n(min(CONFIG$n_amostras_grupo, nrow(.)))
df_mid  <- df_filtrado %>% filter(BraakGroup == "Mid")  %>% sample_n(min(CONFIG$n_amostras_grupo, nrow(.)))
df_high <- df_filtrado %>% filter(BraakGroup == "High") %>% sample_n(min(CONFIG$n_amostras_grupo, nrow(.)))

metadata_subgrupo <- bind_rows(df_low, df_mid, df_high)
amostras_sorteadas <- metadata_subgrupo[[coluna_id_metadado]]

cpm_subgrupo <- cpm_matrix[, amostras_sorteadas, drop = FALSE]

# FILTRO PROTEIN-CODING 
cat("\n--- [PASSO 3] Filtrando Genes Protein-Coding via org.Hs.eg.db ---\n")
keys_db <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
gene_info <- AnnotationDbi::select(org.Hs.eg.db, 
                                   keys = keys_db, 
                                   columns = c("ENSEMBL", "GENETYPE"), 
                                   keytype = "ENSEMBL")

r_protein_genes <- gene_info %>%
  filter(GENETYPE == "protein-coding") %>%
  pull(ENSEMBL) %>%
  na.omit() %>%
  unique()

genes_na_matriz_validos <- intersect(rownames(cpm_subgrupo), r_protein_genes)
cpm_codificantes <- cpm_subgrupo[genes_na_matriz_validos, , drop = FALSE]

rm(keys_db, gene_info, r_protein_genes, genes_na_matriz_validos); gc()

# CÁLCULO DE VARIABILIDADE
cat("\n--- [PASSO 4] Seleção de Genes mais variáveis E Protein-Coding ---\n")
matriz_log_provisoria <- log2(as.matrix(cpm_codificantes) + 1)
vars_log <- apply(matriz_log_provisoria, 1, var)

top_genes <- names(sort(vars_log, decreasing = TRUE)[1:min(CONFIG$n_genes_variaveis, length(vars_log))])
log2_cpm_matrix <- matriz_log_provisoria[top_genes, , drop = FALSE]

writeLines(top_genes, file.path(CONFIG$diretorio_resultados, "lista_genes_selecionados.txt"))
write.csv(metadata_subgrupo[, c(coluna_id_metadado, "braak", "BraakGroup")],
          file.path(CONFIG$diretorio_resultados, "metadata_braak.csv"), row.names = FALSE)

rm(matriz_log_provisoria, cpm_matrix, cpm_subgrupo, cpm_codificantes); gc()

# LIONESS
cat("\n--- [PASSO 5] Inicializando Estruturas do LIONESS ---\n")
rowData <- DataFrame(row.names = rownames(log2_cpm_matrix), gene = rownames(log2_cpm_matrix))
colData <- DataFrame(row.names = metadata_subgrupo[[coluna_id_metadado]],
                     sample = metadata_subgrupo[[coluna_id_metadado]],
                     braak_group = metadata_subgrupo$BraakGroup)

se_filtered <- SummarizedExperiment(assays = list(log_cpm = log2_cpm_matrix),
                                   colData = colData, rowData = rowData)

lioness_cor_optimized <- function(se) {
  X <- assay(se)
  n <- ncol(X); p <- nrow(X)
  S <- rowSums(X); SS <- rowSums(X^2); SP <- X %*% t(X)

  num_full <- SP - outer(S, S) / n
  den_full <- sqrt((SS - S^2 / n) %o% (SS - S^2 / n))
  C_full <- num_full / den_full

  ut_idx <- which(upper.tri(C_full), arr.ind = TRUE)
  genes <- rownames(X)

  reg_names <- genes[ut_idx[, 1]]
  tar_names <- genes[ut_idx[, 2]]
  v_full <- C_full[ut_idx]

  out <- matrix(NA_real_, nrow = length(reg_names), ncol = n)
  colnames(out) <- colnames(X)

  for (q in seq_len(n)) {
    if(q %% 10 == 0) cat("    Processando amostra (LOO):", q, "/", n, "\n")
    xq <- X[, q]
    S_q <- S - xq; SS_q <- SS - xq^2; SP_q <- SP - xq %*% t(xq)
    num_q <- SP_q - outer(S_q, S_q) / (n - 1)
    den_q <- sqrt((SS_q - S_q^2 / (n - 1)) %o% (SS_q - S_q^2 / (n - 1)))
    C_q <- num_q / den_q
    out[, q] <- n * (v_full - C_q[ut_idx]) + C_q[ut_idx]
  }
  return(list(W = out, reg = reg_names, tar = tar_names))
}

cat(" [STATUS] Iniciando LIONESS...\n")
lioness_out <- lioness_cor_optimized(se_filtered)
W <- lioness_out$W
reg_names <- lioness_out$reg
tar_names <- lioness_out$tar
rm(lioness_out); gc()

# STRING INTERACTION FILTER
cat("\n--- [PASSO 6] Cruzando Arestas com Interações do STRING ---\n")
# 1. Carregar e tratar o dicionário de Aliases do STRING
cat(" [STATUS] Carregando arquivo de aliases do STRING...\n")
alias_df <- vroom(CONFIG$caminho_string_alias, show_col_types = FALSE)

# Identifica dinamicamente as colunas do alias_df
col_protein <- colnames(alias_df)[1] # Geralmente #string_protein_id
col_alias   <- "alias"

# Cria o mapeamento ENSP -> ENSG limpando o prefixo 9606.
ensp_to_ensg <- alias_df %>%
  mutate(protein_clean = gsub("9606\\.", "", .data[[col_protein]])) %>%
  filter(grepl("^ENSG", .data[[col_alias]])) %>%
  dplyr::select(protein_clean, .data[[col_alias]]) %>%
  distinct(protein_clean, .keep_all = TRUE) %>%
  deframe()

cat(" [MONITORAMENTO] Dicionário ENSP -> ENSG construído com", length(ensp_to_ensg), "mapeamentos.\n")

# 2. Carregar e filtrar os links experimentais do STRING
links <- vroom(CONFIG$caminho_string_links, delim = " ", show_col_types = FALSE) %>%
  filter(experimental > 0) %>%
  mutate(protein1 = gsub("9606\\.", "", protein1),
         protein2 = gsub("9606\\.", "", protein2)) %>%
  dplyr::select(protein1, protein2)

# 3. Traduzir os IDs de proteína para IDs de gene
links$gene1 <- ensp_to_ensg[links$protein1]
links$gene2 <- ensp_to_ensg[links$protein2]

# Filtra apenas interações onde AMBOS os genes foram traduzidos com sucesso
links <- links %>% filter(!is.na(gene1) & !is.na(gene2))

cat(" [MONITORAMENTO] Links do STRING traduzidos com sucesso para IDs Ensembl:", nrow(links), "\n")

# 4. Cruzamento com a rede LIONESS
s_keys <- unique(paste(pmin(links$gene1, links$gene2), pmax(links$gene1, links$gene2), sep = "_"))
rm(links, alias_df, ensp_to_ensg); gc()

l_keys <- paste(pmin(reg_names, tar_names), pmax(reg_names, tar_names), sep = "_")
valid_mask <- l_keys %in% s_keys

cat(" [MONITORAMENTO] Arestas validadas pelo STRING encontradas na rede LIONESS:", sum(valid_mask), "\n")

if (sum(valid_mask) == 0) {
  stop("Erro Crítico: Nenhuma aresta bateu com a rede STRING. Verifique se os IDs da matriz do R estão no formato ENSG000...")
}

W <- W[valid_mask, , drop = FALSE]
reg_all <- reg_names[valid_mask]
tar_all <- tar_names[valid_mask]
rownames(W) <- paste(reg_all, tar_all, sep = "_")
rm(l_keys, s_keys, valid_mask); gc()

# DEGREE MATRIX
cat("\n--- [PASSO 7] Construindo Matriz de Grau ---\n")
all_genes_in_network <- unique(c(reg_all, tar_all))

deg_list <- lapply(seq_len(ncol(W)), function(i) {
  col_w <- abs(W[, i])
  thr <- mean(col_w, na.rm = TRUE) + (2 * sd(col_w, na.rm = TRUE))
  keep <- !is.na(col_w) & col_w > thr
  active_genes <- c(reg_all[keep], tar_all[keep])
  gene_counts <- table(active_genes)

  sample_deg <- setNames(integer(length(all_genes_in_network)), all_genes_in_network)
  sample_deg[names(gene_counts)] <- as.integer(gene_counts)
  return(sample_deg)
})

deg_mat <- as.data.frame(do.call(cbind, deg_list))
colnames(deg_mat) <- colnames(W)
rownames(deg_mat) <- all_genes_in_network
saveRDS(deg_mat, file.path(CONFIG$diretorio_resultados, "lioness_gene_degree_matrix.rds"))
rm(deg_list, deg_mat); gc()

# Aplicando Thresholding Individual
for (i in 1:ncol(W)) {
  w_col <- abs(W[, i])
  thr <- mean(w_col, na.rm=TRUE) + 2 * sd(w_col, na.rm=TRUE)
  W[w_col <= thr, i] <- NA
}

# ANÁLISE DIFERENCIAL PAIRWISE (LOW x MID, MID x HIGH, LOW x HIGH)
cat("\n--- [PASSO 8] Análise de Conectividade Diferencial por Estágios de Braak ---\n")

# Função auxiliar para realizar as comparações pairwise
executar_teste_pairwise <- function(grupo_A, grupo_B, nome_comparacao, matriz_W, metadados) {
  cat(paste0(" [MONITORAMENTO] Executando comparação: ", nome_comparacao, "\n"))
  
  idx_A <- metadados$BraakGroup == grupo_A
  idx_B <- metadados$BraakGroup == grupo_B
  
  pvals <- apply(matriz_W, 1, function(x) {
    vec_A <- x[idx_A]; vec_A <- vec_A[!is.na(vec_A)]
    vec_B <- x[idx_B]; vec_B <- vec_B[!is.na(vec_B)]
    
    if (length(vec_A) < 2 || length(vec_B) < 2) return(NA_real_)
    var_A <- var(vec_A); var_B <- var(vec_B)
    if (is.na(var_A) || is.na(var_B) || var_A == 0 || var_B == 0) return(NA_real_)
    
    tryCatch({ return(t.test(vec_A, vec_B)$p.value) }, error = function(e) { return(NA_real_) })
  })
  
  padj <- p.adjust(pvals, method = "BH")
  
  res <- data.frame(
    edge = rownames(matriz_W),
    gene1 = reg_all,
    gene2 = tar_all,
    pval = pvals,
    padj = padj,
    stringsAsFactors = FALSE
  )
  
  sig <- res %>% filter(!is.na(padj) & padj < 0.05) %>% arrange(padj)
  
  cat(paste0("    -> Total de arestas significativas (FDR < 0.05) em ", nome_comparacao, ": ", nrow(sig), "\n"))
  
  # Salvando resultados individuais por comparação
  write.csv(sig, file.path(CONFIG$diretorio_resultados, paste0("sig_edges_", nome_comparacao, ".csv")), row.names = FALSE)
  saveRDS(res, file.path(CONFIG$diretorio_resultados, paste0("full_results_", nome_comparacao, ".rds")))
  
  return(res)
}

# Executando as 3 comparações solicitadas
res_low_mid  <- executar_teste_pairwise("Low", "Mid", "Low_vs_Mid", W, metadata_subgrupo)
res_mid_high <- executar_teste_pairwise("Mid", "High", "Mid_vs_High", W, metadata_subgrupo)
res_low_high <- executar_teste_pairwise("Low", "High", "Low_vs_High", W, metadata_subgrupo)

cat("\n--- [PASSO 9] Salvando Matriz LIONESS de Pesos Global ---\n")
saveRDS(W, file.path(CONFIG$diretorio_resultados, "lioness_weight_matrix_full.rds"))
cat(">>> [FIM] Pipeline Completo Executado com Sucesso! <<<\n")