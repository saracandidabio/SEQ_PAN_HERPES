from pathlib import Path
from collections import deque

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st


st.set_page_config(
    page_title="PanHerpes Dashboard",
    page_icon="🧬",
    layout="wide",
)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
HEATMAPS = DATA / "heatmaps"

TARGET_LABELS = {
    "montagem_LIB_3_HSV6": "HSV6",
    "montagem_LIB_3_HSV7": "HSV7",
    "montagens_Alouatta_macconnelli_cytomegalovirus": "Alouatta macconnelli CMV",
    "montagens_Alouatta_palliata_cytomegalovirus": "Alouatta palliata CMV",
    "montagens_CALHV3": "CalHV3",
    "montagens_Cebus_albifrons_lymphocryptovirus_1": "Cebus albifrons LCV1",
    "montagens_consenso_H1": "H1 — consenso interno",
}

FILES = {
    "exec": "execution_evidence.tsv",
    "target": "target_evidence.tsv",
    "aa": "amac_apal.tsv",
    "qc": "qc_summary.tsv",
    "tax_audit": "taxonomy_source_audit.tsv",
    "ambiguity": "execution_ambiguity.tsv",
    "pair_global": "pair_ambiguity_global.tsv",
    "dedup": "dedup_summary.tsv",
    "refs": "references.tsv",
    "comp_exec": "competitive_execution_summary.tsv",
    "comp_target": "competitive_target_summary.tsv",
    "comp_pair": "competitive_pairwise_overlap.tsv",
    "read_qc_R": "read_qc_R.tsv",
    "div_R": "diversity_summary_R.tsv",
    "orientation_R": "read_orientation_R.tsv.gz",
    "pcoa_R": "pcoa_points_R.tsv.gz",
    "hist_length_R": "plot_hist_length_R.tsv.gz",
    "density_length_R": "plot_density_length_R.tsv.gz",
    "hist_filtered_R": "plot_hist_filtered_R.tsv.gz",
    "hist_div_R": "plot_hist_divergence_R.tsv.gz",
}


def _read_readlevel_tsv(path):
    """
    Leitor robusto para os TSV read-level.

    Alguns IDs originais do FASTQ podem conter caracteres que interferem
    com a interpretação convencional de TSV pelo pandas. Como sabemos a
    posição da coluna id_fastq, reconstruímos essa coluna sem descartar reads.
    """

    opener = gzip.open if str(path).endswith(".gz") else open

    rows = []

    with opener(
        path,
        "rt",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as fh:

        header_line = fh.readline()

        if not header_line:
            return pd.DataFrame()

        header = header_line.rstrip("\r\n").split("\t")

        if "id_fastq" not in header:
            raise ValueError(
                f"Arquivo read-level sem coluna id_fastq: {path.name}"
            )

        expected_columns = len(header)

        id_idx = header.index(
            "id_fastq"
        )

        # Quantas colunas existem depois de id_fastq.
        tail_count = (
            expected_columns -
            id_idx -
            1
        )

        for line_number, raw_line in enumerate(
            fh,
            start=2,
        ):

            line = raw_line.rstrip(
                "\r\n"
            )

            parts = line.split(
                "\t"
            )

            # ----------------------------------------------------
            # Caso normal
            # ----------------------------------------------------

            if len(parts) == expected_columns:

                pass

            # ----------------------------------------------------
            # Caso com TAB dentro do id_fastq
            #
            # Reúne novamente os fragmentos pertencentes ao ID
            # sem alterar as demais colunas.
            # ----------------------------------------------------

            elif len(parts) > expected_columns:

                if tail_count > 0:

                    id_parts = parts[
                        id_idx:
                        len(parts) - tail_count
                    ]

                    tail = parts[
                        len(parts) - tail_count:
                    ]

                else:

                    id_parts = parts[
                        id_idx:
                    ]

                    tail = []

                reconstructed_id = " ".join(
                    id_parts
                )

                parts = (
                    parts[:id_idx]
                    +
                    [reconstructed_id]
                    +
                    tail
                )

            # ----------------------------------------------------
            # Menos colunas que o esperado indica realmente
            # uma linha incompleta.
            # ----------------------------------------------------

            else:

                raise ValueError(
                    f"{path.name}: linha {line_number} possui "
                    f"{len(parts)} colunas; eram esperadas "
                    f"{expected_columns}."
                )

            if len(parts) != expected_columns:

                raise ValueError(
                    f"{path.name}: não foi possível reconstruir "
                    f"a linha {line_number}. "
                    f"Obtidas {len(parts)} colunas; "
                    f"esperadas {expected_columns}."
                )

            # Limpa caracteres problemáticos somente do ID.
            parts[id_idx] = (
                parts[id_idx]
                .replace('"', "'")
                .replace("\r", " ")
                .replace("\n", " ")
                .strip()
            )

            rows.append(
                parts
            )

    return pd.DataFrame(
        rows,
        columns=header,
    )


def read_tsv(filename):

    p = DATA / filename

    if not p.exists():
        return pd.DataFrame()

    # Estes dois arquivos carregam IDs originais de FASTQ
    # e usam o leitor robusto.
    if filename in {
        "pcoa_points_R.tsv.gz",
        "read_orientation_R.tsv.gz",
    }:

        return _read_readlevel_tsv(
            p
        )

    # Demais tabelas continuam usando o leitor convencional.
    return pd.read_csv(
        p,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )


@st.cache_data(show_spinner=False)
def load_small():
    keys = [
        "exec", "target", "aa", "qc", "tax_audit", "ambiguity",
        "pair_global", "dedup", "refs", "comp_exec", "comp_target",
        "comp_pair", "read_qc_R", "div_R",
    ]
    return {k: read_tsv(FILES[k]) for k in keys}


@st.cache_data(show_spinner=False)
def load_big(key):
    return read_tsv(FILES[key])


@st.cache_data(show_spinner=False)
def load_heatmap(execution):
    p = HEATMAPS / f"{execution}_divergence_matrix.tsv.gz"
    if not p.exists():
        return pd.DataFrame()
    return pd.read_csv(p, sep="\t", index_col=0)


def numeric(df, cols):
    df = df.copy()
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def fmt_int(x):
    try:
        return f"{int(float(x)):,}".replace(",", ".")
    except Exception:
        return "—"


def fmt_pct(x, digits=1):
    try:
        return f"{float(x):.{digits}f}%"
    except Exception:
        return "—"


def nice_target(x):
    return TARGET_LABELS.get(x, x)


def dbscan_xy(xy, eps=0.35, min_samples=5):
    n = len(xy)
    if n == 0:
        return np.array([], dtype=int)

    mu = np.nanmean(xy, axis=0)
    sd = np.nanstd(xy, axis=0)
    sd[sd == 0] = 1.0
    z = (xy - mu) / sd

    diff = z[:, None, :] - z[None, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=2))
    neighbors = [np.where(dist[i] <= eps)[0] for i in range(n)]

    UNVISITED = -99
    NOISE = -1
    labels = np.full(n, UNVISITED, dtype=int)
    cluster = 0

    for i in range(n):
        if labels[i] != UNVISITED:
            continue
        nbs = neighbors[i]
        if len(nbs) < min_samples:
            labels[i] = NOISE
            continue

        labels[i] = cluster
        queue = deque(int(x) for x in nbs if x != i)

        while queue:
            j = queue.popleft()

            if labels[j] == NOISE:
                labels[j] = cluster

            if labels[j] != UNVISITED:
                continue

            labels[j] = cluster
            jn = neighbors[j]

            if len(jn) >= min_samples:
                for k in jn:
                    k = int(k)
                    if labels[k] in (UNVISITED, NOISE):
                        queue.append(k)

        cluster += 1

    return labels


def cluster_representatives(df):
    rows = []
    for cid, g in df[df["cluster_id"] >= 0].groupby("cluster_id"):
        cx = g["PCoA1"].mean()
        cy = g["PCoA2"].mean()
        d2 = (g["PCoA1"] - cx) ** 2 + (g["PCoA2"] - cy) ** 2
        r = g.loc[d2.idxmin()]
        rows.append(
            {
                "cluster": f"Cluster {int(cid)+1}",
                "n_reads": len(g),
                "read": r["read"],
                "id_fastq": r.get("id_fastq", ""),
                "comprimento": r.get("comprimento", ""),
                "orientacao": r.get("orientacao_escolhida", ""),
                "sequence_oriented": r.get("sequence_oriented", ""),
            }
        )
    return pd.DataFrame(rows)


dfs = load_small()

exec_df = numeric(
    dfs["exec"],
    [
        "input_reads_competitive", "mapped_any_target", "unmapped",
        "reads_one_reported_target", "reads_multiple_reported_targets",
        "reads_unique_best_AS", "reads_tied_best_AS",
        "dedup_records_removed", "dedup_fraction_removed",
    ],
)

target_df = numeric(
    dfs["target"],
    [
        "reads_with_reported_alignment", "reads_unique_best_AS",
        "reads_tied_best_AS_including_target", "median_AS",
        "median_identity", "median_query_coverage",
    ],
)

aa_df = numeric(
    dfs["aa"],
    [
        "amac_only_reported", "apal_only_reported",
        "both_amac_best_AS", "both_apal_best_AS", "both_tied_best_AS",
        "shared_reads", "consensus_joint_acgt_positions",
        "consensus_joint_identical", "consensus_joint_different",
        "diagnostic_sites_called_in_both_consensuses",
        "diagnostic_sites_consensus_same",
        "diagnostic_sites_consensus_different",
    ],
)

read_qc = numeric(
    dfs["read_qc_R"],
    [
        "reads_totais", "media_original", "mediana_original", "moda_original",
        "dp_original", "minimo_original", "maximo_original", "q25_original",
        "q75_original", "reads_filtradas", "reads_removidas",
        "percentual_mantido", "media_filtrada", "mediana_filtrada",
        "moda_filtrada", "dp_filtrado", "minimo_filtrado",
        "maximo_filtrado", "q25_filtrado", "q75_filtrado",
    ],
)

div_summary = numeric(
    dfs["div_R"],
    [
        "reads_filtradas", "reads_amostradas",
        "reads_reverse_complementadas",
        "percentual_reverse_complementadas", "sequencias_unicas",
        "percentual_sequencias_unicas", "numero_comparacoes",
        "divergencia_media", "divergencia_mediana", "divergencia_dp",
        "divergencia_q25", "divergencia_q75",
        "pares_menor_igual_1pct", "pares_menor_igual_3pct",
        "pares_menor_igual_5pct", "PCoA1_var_pct", "PCoA2_var_pct",
        "cmdscale_additive_constant",
    ],
)

if exec_df.empty:
    st.error("Dados agregados ausentes. Rode `python preparar_dados.py`.")
    st.stop()


st.sidebar.title("🧬 PanHerpes")
page = st.sidebar.radio(
    "Página",
    [
        "Visão geral",
        "Reads — comprimento",
        "Reads — divergência",
        "Reads — PCoA & agrupamentos",
        "Reads — heatmap",
        "Reads — orientação",
        "Evidência por execução",
        "Explorar execução",
        "Alvos e referências",
        "Ambiguidade Amac × Apal",
        "QC e auditoria",
        "Downloads",
    ],
)

st.sidebar.caption(
    "Diversidade read-level: método original em R com orientação relativa, "
    "Levenshtein normalizado e cmdscale(add=TRUE)."
)


if page == "Visão geral":
    st.title("PanHerpes — visão geral")

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Execuções", fmt_int(exec_df["execution"].nunique()))
    c2.metric("Alvos válidos", fmt_int(target_df["target_folder"].nunique()))
    c3.metric("Reads de entrada", fmt_int(exec_df["input_reads_competitive"].sum()))
    c4.metric("Reads mapeadas ≥1 alvo", fmt_int(exec_df["mapped_any_target"].sum()))

    left, right = st.columns(2)

    with left:
        st.subheader("Evidência original")
        if "original_evidence_category" in exec_df:
            st.bar_chart(
                exec_df["original_evidence_category"]
                .replace("", "NÃO INFORMADO")
                .value_counts()
            )

    with right:
        st.subheader("Táxons Kraken")
        if "original_top_organism" in exec_df:
            st.bar_chart(
                exec_df["original_top_organism"]
                .replace("", "NÃO INFORMADO")
                .value_counts()
                .head(15)
            )

    if not read_qc.empty:
        st.subheader("Reads após filtro 180–400 pb")
        by_exec = read_qc[["execution", "reads_totais", "reads_filtradas"]].copy()
        by_exec = by_exec.set_index("execution")
        st.bar_chart(by_exec)


elif page == "Reads — comprimento":
    st.title("Reads — distribuição de comprimento")
    st.caption(
        "Equivale aos gráficos 01, 02 e 03 do script R: histograma original, "
        "densidade e histograma após filtro 180–400 pb."
    )

    if read_qc.empty:
        st.warning("Execute `Rscript preparar_reads_metodo_R.R ...`.")
        st.stop()

    execution = st.selectbox(
        "Execução",
        sorted(read_qc["execution"].unique()),
    )

    r = read_qc[read_qc["execution"] == execution].iloc[0]

    a, b, c, d = st.columns(4)
    a.metric("Reads totais", fmt_int(r["reads_totais"]))
    b.metric("Reads 180–400 pb", fmt_int(r["reads_filtradas"]))
    c.metric("Mantidas", fmt_pct(r["percentual_mantido"]))
    d.metric("Mediana filtrada", f"{r['mediana_filtrada']:.0f} pb")

    stats = pd.DataFrame(
        {
            "Métrica": [
                "Média", "Mediana", "Moda", "DP",
                "Mínimo", "Máximo", "Q25", "Q75",
            ],
            "Original": [
                r["media_original"], r["mediana_original"], r["moda_original"],
                r["dp_original"], r["minimo_original"], r["maximo_original"],
                r["q25_original"], r["q75_original"],
            ],
            "Filtrado 180–400 pb": [
                r["media_filtrada"], r["mediana_filtrada"], r["moda_filtrada"],
                r["dp_filtrado"], r["minimo_filtrado"], r["maximo_filtrado"],
                r["q25_filtrado"], r["q75_filtrado"],
            ],
        }
    )
    st.dataframe(stats, width="stretch", hide_index=True)

    hist = numeric(
        load_big("hist_length_R"),
        ["x", "xmin", "xmax", "count"],
    )
    h = hist[hist["execution"] == execution].copy()

    if not h.empty:
        fig = go.Figure(
            go.Bar(
                x=h["x"],
                y=h["count"],
                width=(h["xmax"] - h["xmin"]),
                name="Original",
            )
        )
        fig.update_layout(
            title=f"{execution} — Distribuição dos tamanhos",
            xaxis_title="Tamanho (pb)",
            yaxis_title="Frequência",
        )
        st.plotly_chart(fig, width="stretch")

    dens = numeric(
        load_big("density_length_R"),
        ["x", "density"],
    )
    d = dens[dens["execution"] == execution].copy()

    if not d.empty:
        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                x=d["x"],
                y=d["density"],
                mode="lines",
                name="Densidade",
                fill="tozeroy",
            )
        )
        fig.add_vline(
            x=float(r["media_original"]),
            line_dash="dash",
            annotation_text="Média",
        )
        fig.add_vline(
            x=float(r["mediana_original"]),
            line_dash="dot",
            annotation_text="Mediana",
        )
        fig.update_layout(
            title=f"{execution} — Densidade dos tamanhos",
            xaxis_title="Tamanho (pb)",
            yaxis_title="Densidade",
        )
        st.plotly_chart(fig, width="stretch")

    filt = numeric(
        load_big("hist_filtered_R"),
        ["x", "xmin", "xmax", "count"],
    )
    f = filt[filt["execution"] == execution].copy()

    if not f.empty:
        fig = go.Figure(
            go.Bar(
                x=f["x"],
                y=f["count"],
                width=(f["xmax"] - f["xmin"]),
                name="180–400 pb",
            )
        )
        fig.update_layout(
            title=f"{execution} — Reads filtradas 180–400 pb",
            xaxis_title="Tamanho (pb)",
            yaxis_title="Número de reads",
        )
        st.plotly_chart(fig, width="stretch")


elif page == "Reads — divergência":
    st.title("Reads — divergência entre sequências")
    st.caption(
        "Distância de edição de Levenshtein normalizada pelo maior comprimento "
        "do par. Esta porcentagem não é identidade BLAST."
    )

    if div_summary.empty:
        st.warning("Execute `Rscript preparar_reads_metodo_R.R ...`.")
        st.stop()

    execution = st.selectbox(
        "Execução",
        sorted(div_summary["execution"].unique()),
    )

    r = div_summary[div_summary["execution"] == execution].iloc[0]

    a, b, c, d = st.columns(4)
    a.metric("Reads analisadas", fmt_int(r["reads_amostradas"]))
    b.metric("Comparações", fmt_int(r["numero_comparacoes"]))
    c.metric("Divergência mediana", fmt_pct(r["divergencia_mediana"], 2))
    d.metric("Sequências únicas", fmt_int(r["sequencias_unicas"]))

    e, f, g = st.columns(3)
    e.metric("Pares ≤1%", fmt_pct(r["pares_menor_igual_1pct"], 1))
    f.metric("Pares ≤3%", fmt_pct(r["pares_menor_igual_3pct"], 1))
    g.metric("Pares ≤5%", fmt_pct(r["pares_menor_igual_5pct"], 1))

    hist = numeric(
        load_big("hist_div_R"),
        ["x", "xmin", "xmax", "count"],
    )
    h = hist[hist["execution"] == execution].copy()

    if not h.empty:
        fig = go.Figure(
            go.Bar(
                x=h["x"],
                y=h["count"],
                width=(h["xmax"] - h["xmin"]),
                name="Pares de reads",
            )
        )
        fig.add_vline(
            x=float(r["divergencia_mediana"]),
            line_dash="dash",
            annotation_text="Mediana",
        )
        fig.update_layout(
            title=f"{execution} — Divergência entre reads",
            xaxis_title="Distância de edição normalizada (%)",
            yaxis_title="Número de pares de reads",
        )
        st.plotly_chart(fig, width="stretch")

    summary_table = pd.DataFrame(
        {
            "Métrica": [
                "Divergência média", "Mediana", "DP",
                "Q25", "Q75", "% sequências exatamente únicas",
            ],
            "Valor": [
                r["divergencia_media"], r["divergencia_mediana"],
                r["divergencia_dp"], r["divergencia_q25"],
                r["divergencia_q75"], r["percentual_sequencias_unicas"],
            ],
        }
    )
    st.dataframe(summary_table, width="stretch", hide_index=True)


elif page == "Reads — PCoA & agrupamentos":
    st.title("Reads — PCoA")
    st.caption(
        "PCoA calculada no R com `cmdscale(as.dist(matriz_divergencia), "
        "k=2, eig=TRUE, add=TRUE)`, a partir da divergência de Levenshtein normalizada."
    )

    p = numeric(
        load_big("pcoa_R"),
        ["PCoA1", "PCoA2", "comprimento"],
    )

    if p.empty or div_summary.empty:
        st.warning("Execute `Rscript preparar_reads_metodo_R.R ...`.")
        st.stop()

    executions = sorted(p["execution"].unique())
    default = executions.index("Lib2barcode76") if "Lib2barcode76" in executions else 0

    execution = st.selectbox(
        "Execução",
        executions,
        index=default,
    )

    sub = p[p["execution"] == execution].copy()
    vr = div_summary[div_summary["execution"] == execution]

    if not vr.empty:
        vr = vr.iloc[0]
        v1 = vr["PCoA1_var_pct"]
        v2 = vr["PCoA2_var_pct"]

        a, b, c, d = st.columns(4)
        a.metric("Reads PCoA", fmt_int(vr["reads_amostradas"]))
        b.metric("Reverse-complementadas", fmt_int(vr["reads_reverse_complementadas"]))
        c.metric("PCoA1", fmt_pct(v1, 1))
        d.metric("PCoA2", fmt_pct(v2, 1))
    else:
        v1 = v2 = np.nan

    show_clusters = st.toggle(
        "Mostrar agrupamentos exploratórios sobre a PCoA",
        value=True,
    )

    if show_clusters:
        c1, c2 = st.columns(2)

        eps = c1.slider(
            "Tolerância entre pontos (eps)",
            0.05,
            2.00,
            0.35,
            0.05,
        )

        min_samples = c2.slider(
            "Mínimo de reads por grupo",
            2,
            30,
            5,
            1,
        )

        labels = dbscan_xy(
            sub[["PCoA1", "PCoA2"]].to_numpy(float),
            eps=eps,
            min_samples=min_samples,
        )

        sub["cluster_id"] = labels
        sub["cluster"] = [
            "Ruído / não agrupada" if x < 0 else f"Cluster {x+1}"
            for x in labels
        ]

        color = "cluster"

        c1, c2 = st.columns(2)
        c1.metric("Clusters", len(set(labels) - {-1}))
        c2.metric("Reads não agrupadas", int(np.sum(labels == -1)))

    else:
        sub["cluster_id"] = -1
        sub["cluster"] = "Reads"
        color = None

    xlab = (
        f"PCoA1 ({v1:.1f}%)"
        if np.isfinite(v1)
        else "PCoA1"
    )
    ylab = (
        f"PCoA2 ({v2:.1f}%)"
        if np.isfinite(v2)
        else "PCoA2"
    )

    fig = px.scatter(
        sub,
        x="PCoA1",
        y="PCoA2",
        color=color,
        hover_data=[
            "read",
            "id_fastq",
            "comprimento",
            "orientacao_escolhida",
        ],
        labels={
            "PCoA1": xlab,
            "PCoA2": ylab,
        },
        title=f"{execution} — Similaridade entre reads",
    )

    fig.update_traces(
        marker={
            "size": 8,
            "opacity": 0.7,
        }
    )

    st.plotly_chart(
        fig,
        width="stretch",
    )

    if show_clusters:
        st.info(
            "Os clusters são uma camada exploratória adicionada ao dashboard. "
            "Eles não fazem parte do script R original e não equivalem, por si só, "
            "a haplótipos ou táxons."
        )

        reps = cluster_representatives(
            sub
        )

        st.subheader("Read representante de cada agrupamento")

        if reps.empty:
            st.write("Nenhum cluster definido com estes parâmetros.")
        else:
            display_reps = reps.drop(
                columns=["sequence_oriented"],
                errors="ignore",
            )
            st.dataframe(
                display_reps,
                width="stretch",
                hide_index=True,
            )

            fasta = []

            for _, r in reps.iterrows():
                fasta.append(
                    f">{execution}|{r['cluster']}|{r['id_fastq']}\n"
                    f"{r['sequence_oriented']}"
                )

            st.download_button(
                "Baixar representantes em FASTA",
                "\n".join(fasta) + "\n",
                file_name=f"{execution}_representantes_clusters.fasta",
                mime="text/plain",
            )


elif page == "Reads — heatmap":
    st.title("Reads — heatmap das distâncias")
    st.caption(
        "Matriz simétrica da divergência de Levenshtein normalizada entre as "
        "reads amostradas na análise de diversidade."
    )

    if div_summary.empty:
        st.warning("Execute `Rscript preparar_reads_metodo_R.R ...`.")
        st.stop()

    execution = st.selectbox(
        "Execução",
        sorted(div_summary["execution"].unique()),
    )

    matrix = load_heatmap(
        execution
    )

    if matrix.empty:
        st.warning(
            f"Matriz de distância não encontrada para {execution}."
        )
        st.stop()

    matrix = matrix.apply(
        pd.to_numeric,
        errors="coerce",
    )

    fig = px.imshow(
        matrix.to_numpy(),
        aspect="auto",
        labels={
            "x": "Reads",
            "y": "Reads",
            "color": "Divergência (%)",
        },
        title=f"{execution} — Distância entre reads",
    )

    fig.update_xaxes(
        showticklabels=False
    )
    fig.update_yaxes(
        showticklabels=False
    )

    st.plotly_chart(
        fig,
        width="stretch",
    )

    st.caption(
        f"Matriz: {matrix.shape[0]} × {matrix.shape[1]} reads."
    )


elif page == "Reads — orientação":
    st.title("Reads — orientação forward/reverse-complement")

    orientation = numeric(
        load_big("orientation_R"),
        [
            "indice", "comprimento",
            "distancia_forward", "distancia_reverse",
        ],
    )

    if orientation.empty or div_summary.empty:
        st.warning("Execute `Rscript preparar_reads_metodo_R.R ...`.")
        st.stop()

    execution = st.selectbox(
        "Execução",
        sorted(orientation["execution"].unique()),
    )

    sub = orientation[
        orientation["execution"] == execution
    ].copy()

    vr = div_summary[
        div_summary["execution"] == execution
    ]

    if not vr.empty:
        vr = vr.iloc[0]

        a, b, c = st.columns(3)
        a.metric(
            "Reads analisadas",
            fmt_int(vr["reads_amostradas"]),
        )
        b.metric(
            "Reverse-complementadas",
            fmt_int(vr["reads_reverse_complementadas"]),
        )
        c.metric(
            "% reverse-complementadas",
            fmt_pct(vr["percentual_reverse_complementadas"]),
        )

    counts = (
        sub["orientacao_escolhida"]
        .value_counts()
        .rename_axis("orientacao")
        .reset_index(name="reads")
    )

    st.bar_chart(
        counts.set_index("orientacao")
    )

    show_seq = st.checkbox(
        "Mostrar sequência orientada"
    )

    cols = [
        "read",
        "id_fastq",
        "comprimento",
        "distancia_forward",
        "distancia_reverse",
        "orientacao_escolhida",
    ]

    if show_seq:
        cols.append(
            "sequence_oriented"
        )

    st.dataframe(
        sub[cols],
        width="stretch",
        hide_index=True,
    )


elif page == "Evidência por execução":
    st.title("Evidência por execução")

    filt = exec_df.copy()

    f1, f2, f3 = st.columns(3)

    libs = sorted(
        x
        for x in filt.get(
            "meta_library",
            pd.Series(dtype=str),
        ).unique()
        if x
    )

    evid = sorted(
        x
        for x in filt.get(
            "original_evidence_category",
            pd.Series(dtype=str),
        ).unique()
        if x
    )

    taxa = sorted(
        x
        for x in filt.get(
            "original_top_organism",
            pd.Series(dtype=str),
        ).unique()
        if x
    )

    a = f1.multiselect("Biblioteca", libs)
    b = f2.multiselect("Evidência", evid)
    c = f3.multiselect("Táxon Kraken", taxa)

    if a:
        filt = filt[
            filt["meta_library"].isin(a)
        ]

    if b:
        filt = filt[
            filt["original_evidence_category"].isin(b)
        ]

    if c:
        filt = filt[
            filt["original_top_organism"].isin(c)
        ]

    st.dataframe(
        filt,
        width="stretch",
        hide_index=True,
    )


elif page == "Explorar execução":
    st.title("Explorar execução")

    execution = st.selectbox(
        "Execução",
        sorted(
            exec_df["execution"].unique()
        ),
    )

    row = exec_df[
        exec_df["execution"] == execution
    ].iloc[0]

    a, b, c = st.columns(3)
    a.metric("Reads entrada", fmt_int(row["input_reads_competitive"]))
    b.metric("Mapeadas", fmt_int(row["mapped_any_target"]))
    c.metric("Multi-alvo", fmt_int(row["reads_multiple_reported_targets"]))

    td = target_df[
        target_df["execution"] == execution
    ].copy()

    if not td.empty:
        td["alvo"] = td["target_folder"].map(nice_target)

        st.dataframe(
            td,
            width="stretch",
            hide_index=True,
        )

        st.bar_chart(
            td[
                [
                    "alvo",
                    "reads_with_reported_alignment",
                    "reads_unique_best_AS",
                ]
            ].set_index(
                "alvo"
            )
        )


elif page == "Alvos e referências":
    st.title("Alvos e referências")

    refs = dfs["refs"].copy()

    if not refs.empty:
        refs["alvo"] = refs["target_folder"].map(nice_target)

        st.dataframe(
            refs,
            width="stretch",
            hide_index=True,
        )


elif page == "Ambiguidade Amac × Apal":
    st.title("Ambiguidade AmacCMV × ApalCMV")

    if not aa_df.empty:
        st.dataframe(
            aa_df,
            width="stretch",
            hide_index=True,
        )

        st.bar_chart(
            aa_df[
                [
                    "execution",
                    "amac_only_reported",
                    "apal_only_reported",
                    "both_amac_best_AS",
                    "both_apal_best_AS",
                    "both_tied_best_AS",
                ]
            ].set_index(
                "execution"
            )
        )


elif page == "QC e auditoria":
    st.title("QC e auditoria")

    if not dfs["dedup"].empty:
        st.subheader("Deduplicação")
        st.dataframe(
            dfs["dedup"],
            width="stretch",
            hide_index=True,
        )

    if not dfs["tax_audit"].empty:
        st.subheader("Fonte taxonômica")
        st.dataframe(
            dfs["tax_audit"],
            width="stretch",
            hide_index=True,
        )


elif page == "Downloads":
    st.title("Downloads")

    for key, filename in FILES.items():
        p = DATA / filename

        if p.exists():
            st.download_button(
                f"Baixar — {filename}",
                p.read_bytes(),
                file_name=filename,
                mime="application/octet-stream",
                key=key,
            )
