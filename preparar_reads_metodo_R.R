#!/usr/bin/env Rscript

# ============================================================
# PANHERPES DASHBOARD
# EXPORTAÇÃO READ-LEVEL BASEADA NO MÉTODO ORIGINAL EM R
#
# Método preservado:
#   - FASTQ após Chopper / FASTQ canônico deduplicado
#   - filtro 180-400 pb
#   - no máximo 300 reads por execução
#   - set.seed(123)
#   - normalização forward/reverse-complement
#   - distância de Levenshtein com adist()
#   - divergência = distância / maior comprimento do par * 100
#   - cmdscale(..., k=2, eig=TRUE, add=TRUE)
#
# Uso:
# Rscript preparar_reads_metodo_R.R \
#   "/Users/navio01/Downloads/PanHerpes_dashboard_data" \
#   "./data"
# ============================================================

suppressPackageStartupMessages({
  library(ShortRead)
  library(Biostrings)
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    paste0(
      "Uso:\n",
      "Rscript preparar_reads_metodo_R.R ",
      "\"/caminho/PanHerpes_dashboard_data\" \"./data\"\n"
    )
  )
}

dados_dir <- normalizePath(
  args[1],
  winslash = "/",
  mustWork = TRUE
)

saida_dir <- if (length(args) >= 2) {
  args[2]
} else {
  "./data"
}

dir.create(
  saida_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

heatmap_dir <- file.path(
  saida_dir,
  "heatmaps"
)

dir.create(
  heatmap_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# CONFIGURAÇÃO ORIGINAL
# ============================================================

min_compr <- 180
max_compr <- 400
n_reads_diversidade <- 300

set.seed(123)


# ============================================================
# FUNÇÕES
# ============================================================

calc_moda <- function(x) {

  if (length(x) == 0) {
    return(NA_real_)
  }

  freq <- table(x)

  as.numeric(
    names(freq)[which.max(freq)]
  )
}


estatisticas_completas <- function(comprimentos) {

  if (length(comprimentos) == 0) {

    return(
      data.frame(
        total = 0,
        media = NA_real_,
        mediana = NA_real_,
        moda = NA_real_,
        dp = NA_real_,
        minimo = NA_real_,
        maximo = NA_real_,
        q25 = NA_real_,
        q75 = NA_real_
      )
    )
  }

  data.frame(
    total = length(comprimentos),
    media = mean(comprimentos),
    mediana = median(comprimentos),
    moda = calc_moda(comprimentos),
    dp = if (length(comprimentos) > 1) sd(comprimentos) else 0,
    minimo = min(comprimentos),
    maximo = max(comprimentos),
    q25 = as.numeric(quantile(comprimentos, 0.25)),
    q75 = as.numeric(quantile(comprimentos, 0.75))
  )
}


orientar_reads <- function(seqs) {

  seqs_char <- as.character(seqs)

  tamanhos <- nchar(seqs_char)

  comprimento_mediano <- median(tamanhos)

  indice_referencia <- which.min(
    abs(
      tamanhos -
        comprimento_mediano
    )
  )

  referencia <- seqs_char[
    indice_referencia
  ]

  seqs_rc <- as.character(
    reverseComplement(
      DNAStringSet(
        seqs_char
      )
    )
  )

  dist_forward <- as.numeric(
    adist(
      seqs_char,
      referencia
    )
  )

  dist_reverse <- as.numeric(
    adist(
      seqs_rc,
      referencia
    )
  )

  usar_reverse <- (
    dist_reverse <
      dist_forward
  )

  seqs_orientadas <- seqs_char

  seqs_orientadas[
    usar_reverse
  ] <- seqs_rc[
    usar_reverse
  ]

  tabela_orientacao <- data.frame(
    indice = seq_along(seqs_char),
    comprimento = tamanhos,
    distancia_forward = dist_forward,
    distancia_reverse = dist_reverse,
    orientacao_escolhida = ifelse(
      usar_reverse,
      "reverse_complement",
      "forward"
    ),
    stringsAsFactors = FALSE
  )

  list(
    sequencias = DNAStringSet(seqs_orientadas),
    tabela = tabela_orientacao,
    indice_referencia = indice_referencia,
    numero_reverse = sum(usar_reverse),
    percentual_reverse = mean(usar_reverse) * 100
  )
}


write_tsv_gz <- function(df, path) {

  con <- gzfile(
    path,
    open = "wt"
  )

  on.exit(
    close(con),
    add = TRUE
  )

  write.table(
    df,
    con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}


plot_data_with_execution <- function(plot_object, execution) {

  x <- ggplot_build(
    plot_object
  )$data[[1]]

  x$execution <- execution

  x
}


# ============================================================
# MANIFESTO DOS FASTQs CANÔNICOS/DEDUPLICADOS
# ============================================================

manifest_candidates <- c(
  file.path(
    dados_dir,
    "dataset_valido",
    "competitive_inputs",
    "competitive_input_manifest.tsv"
  ),
  file.path(
    dados_dir,
    "dataset_valido",
    "canonical_fastq_manifest.tsv"
  )
)

manifest_file <- manifest_candidates[
  file.exists(
    manifest_candidates
  )
][1]

if (
  length(manifest_file) == 0 ||
  is.na(manifest_file)
) {
  stop(
    "Manifesto dos FASTQs canônicos não encontrado."
  )
}

manifest <- read.delim(
  manifest_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

find_col <- function(candidatos) {

  x <- candidatos[
    candidatos %in% names(manifest)
  ]

  if (length(x) == 0) {
    return(NA_character_)
  }

  x[1]
}

exec_col <- find_col(
  c(
    "execution",
    "sample",
    "sample_id"
  )
)

fastq_col <- find_col(
  c(
    "fastq_for_competitive_mapping",
    "canonical_fastq",
    "fastq_path",
    "trimmed_fastq"
  )
)

status_col <- find_col(
  c(
    "input_status",
    "status"
  )
)

if (
  is.na(exec_col) ||
  is.na(fastq_col)
) {
  stop(
    paste0(
      "Colunas incompatíveis no manifesto. Encontradas: ",
      paste(
        names(manifest),
        collapse = ", "
      )
    )
  )
}

manifest <- manifest[
  order(
    manifest[[exec_col]]
  ),
  ,
  drop = FALSE
]


# ============================================================
# OBJETOS DE SAÍDA
# ============================================================

resumo_qc <- data.frame()
resumo_div <- data.frame()
orientacao_geral <- data.frame()
pcoa_geral <- data.frame()
hist_comprimento_geral <- data.frame()
densidade_comprimento_geral <- data.frame()
hist_filtrado_geral <- data.frame()
hist_divergencia_geral <- data.frame()
erros <- data.frame()


cat("\n")
cat("============================================================\n")
cat("PANHERPES - EXPORTAÇÃO READ-LEVEL PARA O DASHBOARD\n")
cat("============================================================\n")
cat("Manifesto:", manifest_file, "\n")
cat("Execuções:", nrow(manifest), "\n")
cat("Filtro:", min_compr, "-", max_compr, "pb\n")
cat("Reads diversidade:", n_reads_diversidade, "\n")
cat("Semente: 123\n")
cat("============================================================\n")


# ============================================================
# LOOP
# ============================================================

for (i in seq_len(nrow(manifest))) {

  execution <- manifest[[exec_col]][i]

  if (
    !is.na(status_col) &&
    nzchar(manifest[[status_col]][i]) &&
    toupper(manifest[[status_col]][i]) != "OK"
  ) {

    cat(
      "[SKIP]",
      execution,
      "status =",
      manifest[[status_col]][i],
      "\n"
    )

    next
  }

  fastq_file <- manifest[[fastq_col]][i]

  cat(
    "\n[",
    i,
    "/",
    nrow(manifest),
    "] ",
    execution,
    "\n",
    sep = ""
  )

  tryCatch(
    {

      fastq_data <- readFastq(
        fastq_file
      )

      sequencias <- sread(
        fastq_data
      )

      comprimentos <- width(
        sequencias
      )

      stats_orig <- estatisticas_completas(
        comprimentos
      )


      # --------------------------------------------------------
      # HISTOGRAMA ORIGINAL - mesmo geom_histogram(bins = 100)
      # --------------------------------------------------------

      df_original <- data.frame(
        comprimento = comprimentos
      )

      p_hist <- ggplot(
        df_original,
        aes(
          x = comprimento
        )
      ) +
        geom_histogram(
          bins = 100
        )

      hist_tmp <- plot_data_with_execution(
        p_hist,
        execution
      )

      hist_comprimento_geral <- bind_rows(
        hist_comprimento_geral,
        hist_tmp
      )


      # --------------------------------------------------------
      # DENSIDADE ORIGINAL
      # --------------------------------------------------------

      if (
        length(comprimentos) > 1 &&
        length(unique(comprimentos)) > 1
      ) {

        p_dens <- ggplot(
          df_original,
          aes(
            x = comprimento
          )
        ) +
          geom_density()

        dens_tmp <- plot_data_with_execution(
          p_dens,
          execution
        )

        densidade_comprimento_geral <- bind_rows(
          densidade_comprimento_geral,
          dens_tmp
        )
      }


      # --------------------------------------------------------
      # FILTRO 180-400
      # --------------------------------------------------------

      filtro <- (
        comprimentos >= min_compr &
          comprimentos <= max_compr
      )

      fastq_filtrado <- fastq_data[
        filtro
      ]

      comprimentos_filtrados <- comprimentos[
        filtro
      ]

      total_original <- length(
        fastq_data
      )

      total_filtrado <- length(
        fastq_filtrado
      )

      percentual <- if (
        total_original > 0
      ) {
        (
          total_filtrado /
            total_original
        ) * 100
      } else {
        NA_real_
      }

      stats_filt <- estatisticas_completas(
        comprimentos_filtrados
      )


      # --------------------------------------------------------
      # HISTOGRAMA FILTRADO - mesmo binwidth = 1
      # --------------------------------------------------------

      if (
        total_filtrado > 0
      ) {

        df_filtrado <- data.frame(
          comprimento =
            comprimentos_filtrados
        )

        p_filtrado <- ggplot(
          df_filtrado,
          aes(
            x = comprimento
          )
        ) +
          geom_histogram(
            binwidth = 1
          )

        filt_tmp <- plot_data_with_execution(
          p_filtrado,
          execution
        )

        hist_filtrado_geral <- bind_rows(
          hist_filtrado_geral,
          filt_tmp
        )
      }


      # --------------------------------------------------------
      # QC SUMMARY
      # --------------------------------------------------------

      linha_qc <- data.frame(
        execution = execution,
        reads_totais = total_original,
        media_original = stats_orig$media,
        mediana_original = stats_orig$mediana,
        moda_original = stats_orig$moda,
        dp_original = stats_orig$dp,
        minimo_original = stats_orig$minimo,
        maximo_original = stats_orig$maximo,
        q25_original = stats_orig$q25,
        q75_original = stats_orig$q75,
        reads_filtradas = total_filtrado,
        reads_removidas =
          total_original -
          total_filtrado,
        percentual_mantido = percentual,
        media_filtrada = stats_filt$media,
        mediana_filtrada = stats_filt$mediana,
        moda_filtrada = stats_filt$moda,
        dp_filtrado = stats_filt$dp,
        minimo_filtrado = stats_filt$minimo,
        maximo_filtrado = stats_filt$maximo,
        q25_filtrado = stats_filt$q25,
        q75_filtrado = stats_filt$q75,
        min_compr = min_compr,
        max_compr = max_compr,
        stringsAsFactors = FALSE
      )

      resumo_qc <- bind_rows(
        resumo_qc,
        linha_qc
      )


      # --------------------------------------------------------
      # DIVERSIDADE
      # --------------------------------------------------------

      if (
        total_filtrado >= 3
      ) {

        n_div <- min(
          n_reads_diversidade,
          total_filtrado
        )

        if (
          total_filtrado >
          n_reads_diversidade
        ) {

          indices_div <- sample(
            seq_len(
              total_filtrado
            ),
            size = n_div,
            replace = FALSE
          )

        } else {

          indices_div <- seq_len(
            total_filtrado
          )
        }


        seqs_div_original <- sread(
          fastq_filtrado[
            indices_div
          ]
        )

        ids_fastq <- as.character(
          id(
            fastq_filtrado[
              indices_div
            ]
          )
        )

        nomes_reads <- paste0(
          execution,
          "_read",
          seq_len(
            n_div
          )
        )

        names(
          seqs_div_original
        ) <- nomes_reads


        resultado_orientacao <- orientar_reads(
          seqs_div_original
        )

        seqs_div <- resultado_orientacao$sequencias

        names(
          seqs_div
        ) <- nomes_reads

        n_reverse <- resultado_orientacao$numero_reverse

        pct_reverse <- resultado_orientacao$percentual_reverse


        tabela_orientacao <- resultado_orientacao$tabela

        tabela_orientacao$execution <- execution

        tabela_orientacao$read <- nomes_reads

        tabela_orientacao$id_fastq <- ids_fastq

        tabela_orientacao$sequence_oriented <- as.character(
          seqs_div
        )

        orientacao_geral <- bind_rows(
          orientacao_geral,
          tabela_orientacao
        )


        n_seqs_unicas <- length(
          unique(
            as.character(
              seqs_div
            )
          )
        )

        pct_seqs_unicas <- (
          n_seqs_unicas /
            n_div
        ) * 100


        # ------------------------------------------------------
        # LEVENSHTEIN + DIVERGÊNCIA NORMALIZADA
        # ------------------------------------------------------

        seqs_char <- as.character(
          seqs_div
        )

        matriz_dist <- adist(
          seqs_char,
          seqs_char
        )

        rownames(
          matriz_dist
        ) <- nomes_reads

        colnames(
          matriz_dist
        ) <- nomes_reads

        tamanhos_div <- width(
          seqs_div
        )

        denominador <- outer(
          tamanhos_div,
          tamanhos_div,
          FUN = pmax
        )

        matriz_divergencia <- (
          matriz_dist /
            denominador
        ) * 100

        diag(
          matriz_divergencia
        ) <- 0


        # ------------------------------------------------------
        # SALVAR MATRIZ PARA HEATMAP INTERATIVO
        # ------------------------------------------------------

        heatmap_file <- file.path(
          heatmap_dir,
          paste0(
            execution,
            "_divergence_matrix.tsv.gz"
          )
        )

        con_heat <- gzfile(
          heatmap_file,
          open = "wt"
        )

        write.table(
          matriz_divergencia,
          con_heat,
          sep = "\t",
          quote = FALSE,
          row.names = TRUE,
          col.names = NA
        )

        close(
          con_heat
        )


        divergencias <- matriz_divergencia[
          upper.tri(
            matriz_divergencia
          )
        ]


        div_media <- mean(
          divergencias,
          na.rm = TRUE
        )

        div_mediana <- median(
          divergencias,
          na.rm = TRUE
        )

        div_dp <- sd(
          divergencias,
          na.rm = TRUE
        )

        div_q25 <- as.numeric(
          quantile(
            divergencias,
            0.25,
            na.rm = TRUE
          )
        )

        div_q75 <- as.numeric(
          quantile(
            divergencias,
            0.75,
            na.rm = TRUE
          )
        )

        pct_menor_1 <- mean(
          divergencias <= 1,
          na.rm = TRUE
        ) * 100

        pct_menor_3 <- mean(
          divergencias <= 3,
          na.rm = TRUE
        ) * 100

        pct_menor_5 <- mean(
          divergencias <= 5,
          na.rm = TRUE
        ) * 100


        # ------------------------------------------------------
        # HISTOGRAMA DIVERGÊNCIA - mesmos 80 bins
        # ------------------------------------------------------

        df_div <- data.frame(
          divergencia =
            divergencias
        )

        p_div <- ggplot(
          df_div,
          aes(
            x = divergencia
          )
        ) +
          geom_histogram(
            bins = 80
          )

        div_hist_tmp <- plot_data_with_execution(
          p_div,
          execution
        )

        hist_divergencia_geral <- bind_rows(
          hist_divergencia_geral,
          div_hist_tmp
        )


        # ------------------------------------------------------
        # PCoA - MESMO cmdscale(add = TRUE)
        # ------------------------------------------------------

        var1 <- NA_real_
        var2 <- NA_real_
        cailliez_additive_constant <- NA_real_

        if (
          n_seqs_unicas > 1 &&
          max(
            matriz_divergencia,
            na.rm = TRUE
          ) > 0
        ) {

          pcoa <- tryCatch(
            {
              cmdscale(
                as.dist(
                  matriz_divergencia
                ),
                k = 2,
                eig = TRUE,
                add = TRUE
              )
            },
            error = function(e) {
              NULL
            }
          )

          if (
            !is.null(
              pcoa
            )
          ) {

            eig_pos <- pcoa$eig[
              pcoa$eig > 0
            ]

            if (
              length(eig_pos) >= 2
            ) {

              var1 <- (
                eig_pos[1] /
                  sum(
                    eig_pos
                  )
              ) * 100

              var2 <- (
                eig_pos[2] /
                  sum(
                    eig_pos
                  )
              ) * 100
            }

            if (
              !is.null(
                pcoa$ac
              )
            ) {
              cailliez_additive_constant <- pcoa$ac
            }


            df_pcoa <- data.frame(
              execution = execution,
              read = nomes_reads,
              id_fastq = ids_fastq,
              PCoA1 = pcoa$points[, 1],
              PCoA2 = pcoa$points[, 2],
              orientacao_escolhida =
                tabela_orientacao$orientacao_escolhida,
              comprimento =
                tabela_orientacao$comprimento,
              sequence_oriented =
                as.character(
                  seqs_div
                ),
              stringsAsFactors = FALSE
            )

            pcoa_geral <- bind_rows(
              pcoa_geral,
              df_pcoa
            )
          }
        }


        linha_div <- data.frame(
          execution = execution,
          reads_filtradas = total_filtrado,
          reads_amostradas = n_div,
          reads_reverse_complementadas = n_reverse,
          percentual_reverse_complementadas = pct_reverse,
          sequencias_unicas = n_seqs_unicas,
          percentual_sequencias_unicas = pct_seqs_unicas,
          numero_comparacoes = length(divergencias),
          divergencia_media = div_media,
          divergencia_mediana = div_mediana,
          divergencia_dp = div_dp,
          divergencia_q25 = div_q25,
          divergencia_q75 = div_q75,
          pares_menor_igual_1pct = pct_menor_1,
          pares_menor_igual_3pct = pct_menor_3,
          pares_menor_igual_5pct = pct_menor_5,
          PCoA1_var_pct = var1,
          PCoA2_var_pct = var2,
          cmdscale_additive_constant =
            cailliez_additive_constant,
          stringsAsFactors = FALSE
        )

        resumo_div <- bind_rows(
          resumo_div,
          linha_div
        )
      }


      rm(
        fastq_data,
        sequencias
      )

      gc()

    },
    error = function(e) {

      cat(
        "  ERRO:",
        conditionMessage(e),
        "\n"
      )

      erros <<- bind_rows(
        erros,
        data.frame(
          execution = execution,
          fastq = fastq_file,
          error = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      )
    }
  )
}


# ============================================================
# EXPORTAR
# ============================================================

write.table(
  resumo_qc,
  file.path(
    saida_dir,
    "read_qc_R.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

write.table(
  resumo_div,
  file.path(
    saida_dir,
    "diversity_summary_R.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

write.table(
  erros,
  file.path(
    saida_dir,
    "read_method_R_errors.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

write_tsv_gz(
  orientacao_geral,
  file.path(
    saida_dir,
    "read_orientation_R.tsv.gz"
  )
)

write_tsv_gz(
  pcoa_geral,
  file.path(
    saida_dir,
    "pcoa_points_R.tsv.gz"
  )
)

write_tsv_gz(
  hist_comprimento_geral,
  file.path(
    saida_dir,
    "plot_hist_length_R.tsv.gz"
  )
)

write_tsv_gz(
  densidade_comprimento_geral,
  file.path(
    saida_dir,
    "plot_density_length_R.tsv.gz"
  )
)

write_tsv_gz(
  hist_filtrado_geral,
  file.path(
    saida_dir,
    "plot_hist_filtered_R.tsv.gz"
  )
)

write_tsv_gz(
  hist_divergencia_geral,
  file.path(
    saida_dir,
    "plot_hist_divergence_R.tsv.gz"
  )
)


cat("\n")
cat("============================================================\n")
cat("FINAL\n")
cat("============================================================\n")
cat("Execuções QC:", nrow(resumo_qc), "\n")
cat("Execuções diversidade:", nrow(resumo_div), "\n")
cat("Pontos PCoA:", nrow(pcoa_geral), "\n")
cat("Erros:", nrow(erros), "\n")
cat("Saída:", normalizePath(saida_dir, mustWork = FALSE), "\n")
cat("============================================================\n")
