# PanHerpes Dashboard — método de diversidade do script R

Esta versão usa, para a camada read-level, o mesmo método do script R fornecido:

1. FASTQ após Chopper / FASTQ canônico deduplicado;
2. filtro 180–400 pb;
3. amostragem máxima de 300 reads;
4. `set.seed(123)`;
5. orientação relativa forward/reverse-complement;
6. distância de Levenshtein (`adist`);
7. divergência normalizada = distância / maior comprimento do par × 100;
8. PCoA com `cmdscale(..., eig=TRUE, add=TRUE)`.

## Plots disponíveis

- histograma de comprimento;
- densidade de comprimento;
- histograma de reads filtradas;
- histograma de divergência;
- PCoA;
- heatmap da matriz de divergência;
- gráfico de orientação forward/reverse-complement;
- agrupamentos exploratórios opcionais sobre a PCoA;
- read representante por agrupamento.

## Preparação

### 1. Tabelas agregadas

```bash
cd ~/Downloads/panherpes-dashboard-metodo-R

python preparar_dados.py \
  --dados "/Users/navio01/Downloads/PanHerpes_dashboard_data"
```

### 2. Camada read-level usando o método R

Use um ambiente R que já tenha:

```r
ShortRead
Biostrings
ggplot2
dplyr
```

Depois:

```bash
Rscript preparar_reads_metodo_R.R \
  "/Users/navio01/Downloads/PanHerpes_dashboard_data" \
  "./data"
```

Arquivos esperados:

```text
data/
├── read_qc_R.tsv
├── diversity_summary_R.tsv
├── read_orientation_R.tsv.gz
├── pcoa_points_R.tsv.gz
├── plot_hist_length_R.tsv.gz
├── plot_density_length_R.tsv.gz
├── plot_hist_filtered_R.tsv.gz
├── plot_hist_divergence_R.tsv.gz
├── read_method_R_errors.tsv
└── heatmaps/
    └── LibXbarcodeYY_divergence_matrix.tsv.gz
```

### 3. Teste local

```bash
streamlit run app.py
```

### 4. Deploy

No Streamlit Community Cloud:

- main file: `app.py`
- Python: `3.13`

O Streamlit Cloud não precisa executar R. O R é usado somente localmente para
pré-calcular as tabelas e matrizes que serão armazenadas no GitHub.
