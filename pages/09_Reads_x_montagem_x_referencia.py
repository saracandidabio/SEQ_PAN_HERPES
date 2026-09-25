from pathlib import Path
import csv

import pandas as pd
import plotly.express as px
import streamlit as st


# ============================================================
# PANHERPES — Reads × montagem × referência
# V3 = proximidade baseada em SAM/CIGAR/NM/AS
# V4 = PCoA homóloga reads + referência + consenso
# ============================================================

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


FILES = {
    "v3_summary": DATA / "integration_summary_v3.tsv",
    "v3_reads": DATA / "read_reference_alignment_v3.tsv.gz",
    "v3_pair": DATA / "reference_pairwise_reads_v3.tsv.gz",

    "v4_pcoa": DATA / "integrated_homology_pcoa_v4.tsv.gz",
    "v4_summary": DATA / "integrated_homology_summary_v4.tsv",
    "v4_cons_ref": DATA / "integrated_consensus_reference_v4.tsv",
    "v4_errors": DATA / "integrated_homology_errors_v4.tsv",

    # opcional: permite selecionar também as demais execuções
    "diversity_R": DATA / "diversity_summary_R.tsv",
}


@st.cache_data(show_spinner=False)
def read_tsv(path):
    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        quoting=csv.QUOTE_NONE,
        engine="python",
    )


def to_num(df, cols):
    df = df.copy()

    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(
                df[col],
                errors="coerce",
            )

    return df


def fmt(x, decimals=2, suffix=""):
    if pd.isna(x):
        return "—"

    return f"{float(x):.{decimals}f}{suffix}"


def contagem(series):
    if series.empty:
        return "—"

    vc = series.value_counts()

    return " · ".join(
        f"{nome}: {n}"
        for nome, n in vc.items()
    )


# ============================================================
# VALIDAR ARQUIVOS
# ============================================================

required = [
    "v3_summary",
    "v3_reads",
    "v3_pair",
    "v4_pcoa",
    "v4_summary",
    "v4_cons_ref",
    "v4_errors",
]

missing = [
    FILES[x].name
    for x in required
    if not FILES[x].exists()
]

if missing:
    st.error(
        "Arquivos necessários ausentes em data/: "
        + ", ".join(missing)
    )
    st.stop()


# ============================================================
# CARREGAR
# ============================================================

v3_summary = read_tsv(FILES["v3_summary"])

v3_summary = to_num(
    v3_summary,
    [
        "sampled_reads_R",
        "sampled_reads_mapped",
        "sampled_reads_canonical_confirmed",
        "median_divergence_aligned_pct",
        "q25_divergence_aligned_pct",
        "q75_divergence_aligned_pct",
        "median_query_coverage_pct",
        "reads_coverage_ge_70",
        "reads_coverage_ge_80",
        "reads_coverage_ge_90",
    ],
)


v3_reads = read_tsv(
    FILES["v3_reads"]
)

v3_reads = to_num(
    v3_reads,
    [
        "MAPQ",
        "NM",
        "AS",
        "read_length",
        "query_aligned_bases",
        "reference_aligned_bases",
        "alignment_columns",
        "soft_clipped_bases",
        "inserted_bases",
        "deleted_bases",
        "divergence_aligned_pct",
        "identity_aligned_est_pct",
        "query_coverage_pct",
    ],
)


v3_pair = read_tsv(
    FILES["v3_pair"]
)

v3_pair = to_num(
    v3_pair,
    [
        "divergence_A_pct",
        "divergence_B_pct",
        "delta_B_minus_A_pct",
        "query_coverage_A_pct",
        "query_coverage_B_pct",
        "AS_A",
        "AS_B",
        "MAPQ_A",
        "MAPQ_B",
        "min_query_coverage_pct",
    ],
)


v4_pcoa = read_tsv(
    FILES["v4_pcoa"]
)

v4_pcoa = to_num(
    v4_pcoa,
    [
        "original_length",
        "PCoA1",
        "PCoA2",
        "PCoA1_var_pct",
        "PCoA2_var_pct",
    ],
)


v4_summary = read_tsv(
    FILES["v4_summary"]
)

v4_cons_ref = read_tsv(
    FILES["v4_cons_ref"]
)

v4_cons_ref = to_num(
    v4_cons_ref,
    [
        "comparable_columns",
        "differences",
        "divergence_consensus_reference_pct",
        "identity_consensus_reference_pct",
        "internal_gap_base_columns",
    ],
)


v4_errors = read_tsv(
    FILES["v4_errors"]
)


# ============================================================
# CONJUNTOS DE EXECUÇÕES
# ============================================================

exec_v3 = set(
    v3_summary["execution"]
)

exec_v4 = set(
    v4_pcoa["execution"]
)

exec_error = set(
    v4_errors["execution"]
)

exec_canonical = set(
    v3_reads["execution"]
)

exec_sem_canonical = (
    exec_v3 -
    exec_canonical
)


# incluir 221 execuções, se o resumo R estiver disponível
all_exec = set(exec_v3)

if FILES["diversity_R"].exists():

    try:

        dr = read_tsv(
            FILES["diversity_R"]
        )

        if "execution" in dr.columns:
            all_exec.update(
                dr["execution"]
            )

    except Exception:
        pass


all_exec = sorted(
    all_exec
)


# ============================================================
# CABEÇALHO
# ============================================================

st.title(
    "🧬 Reads × montagem × referência"
)

st.caption(
    "Integra a proximidade read→referência baseada no "
    "alinhamento real (V3: SAM/CIGAR/NM/AS) com a PCoA "
    "homóloga de reads + referência(s) + consenso(s) (V4). "
    "A PCoA original das reads permanece uma análise separada."
)


# ============================================================
# VISÃO GLOBAL
# ============================================================

c1, c2, c3, c4, c5 = st.columns(5)

c1.metric(
    "Execuções V3",
    len(exec_v3),
)

c2.metric(
    "PCoA V4",
    len(exec_v4),
)

c3.metric(
    "V4: <3 entidades",
    len(exec_error),
)

c4.metric(
    "Sem read V3 canônica",
    len(exec_sem_canonical),
)

c5.metric(
    "Alinhamentos V3",
    f"{len(v3_reads):,}".replace(",", "."),
)


with st.expander(
    "Cobertura e auditoria da integração"
):

    a, b, c, d = st.columns(4)

    a.metric(
        "Pontos V4",
        f"{len(v4_pcoa):,}".replace(",", "."),
    )

    b.metric(
        "Reads V4",
        int(
            (
                v4_pcoa["entity_type"]
                == "Read"
            ).sum()
        ),
    )

    c.metric(
        "Referências V4",
        int(
            (
                v4_pcoa["entity_type"]
                == "Reference"
            ).sum()
        ),
    )

    d.metric(
        "Consensos V4",
        int(
            (
                v4_pcoa["entity_type"]
                == "Consensus"
            ).sum()
        ),
    )

    st.markdown(
        "**Execuções com menos de 3 entidades para PCoA:**"
    )

    st.dataframe(
        v4_errors,
        hide_index=True,
        use_container_width=True,
    )

    st.markdown(
        "**Execuções V3 sem read canônica da amostra R para V4:**"
    )

    st.write(
        ", ".join(
            sorted(
                exec_sem_canonical
            )
        )
    )


st.divider()


# ============================================================
# EXECUÇÃO
# ============================================================

execution = st.selectbox(
    "Execução",
    all_exec,
)


# status da integração

if execution in exec_v4:

    st.success(
        "PCoA homóloga V4 disponível para esta execução."
    )

elif execution in exec_error:

    erro = v4_errors.loc[
        v4_errors["execution"] == execution,
        "error",
    ].iloc[0]

    st.warning(
        f"V3 disponível, mas a PCoA V4 não foi construída: {erro}"
    )

elif execution in exec_v3:

    st.warning(
        "Esta execução possui evidência V3, mas não possui read "
        "V3 canônica entre as reads amostradas pelo método R. "
        "A PCoA homóloga não foi forçada."
    )

else:

    st.info(
        "Esta execução não entrou na integração V3/V4. "
        "Consulte as páginas da análise original das reads."
    )

    st.stop()


summary_exec = v3_summary[
    v3_summary["execution"]
    == execution
].copy()


if summary_exec.empty:
    st.stop()


labels = sorted(
    summary_exec[
        "target_label"
    ].unique()
)


x1, x2, x3 = st.columns(3)

x1.metric(
    "Alvos V3",
    len(labels),
)

x2.metric(
    "Reads R amostradas",
    int(
        summary_exec[
            "sampled_reads_R"
        ].max()
    ),
)

x3.metric(
    "Alvos",
    " | ".join(labels),
)


# ============================================================
# PCoA V4
# ============================================================

st.subheader(
    "PCoA homóloga V4"
)

p = v4_pcoa[
    v4_pcoa["execution"]
    == execution
].copy()


if p.empty:

    st.info(
        "PCoA homóloga não disponível para esta execução. "
        "Isso não equivale a ausência de sinal viral."
    )

else:

    p["_size"] = p[
        "entity_type"
    ].map(
        {
            "Read": 7,
            "Consensus": 14,
            "Reference": 17,
        }
    )

    var1 = p[
        "PCoA1_var_pct"
    ].dropna()

    var2 = p[
        "PCoA2_var_pct"
    ].dropna()

    xlab = (
        f"PCoA1 ({var1.iloc[0]:.1f}%)"
        if len(var1)
        else "PCoA1"
    )

    ylab = (
        f"PCoA2 ({var2.iloc[0]:.1f}%)"
        if len(var2)
        else "PCoA2"
    )

    fig = px.scatter(
        p,
        x="PCoA1",
        y="PCoA2",
        color="entity_type",
        symbol="entity_type",
        symbol_map={
            "Read": "circle",
            "Consensus": "diamond",
            "Reference": "star",
        },
        size="_size",
        size_max=18,
        hover_data={
            "entity_id": True,
            "target_label": True,
            "mapped_targets": True,
            "original_length": True,
            "source": True,
            "_size": False,
        },
        labels={
            "PCoA1": xlab,
            "PCoA2": ylab,
            "entity_type": "Entidade",
        },
    )

    st.plotly_chart(
        fig,
        use_container_width=True,
    )

    st.caption(
        "● Read · ◆ Consenso · ★ Referência. "
        "A PCoA resume distâncias na região homóloga; "
        "não é, isoladamente, uma classificação taxonômica."
    )


# ============================================================
# V3 — ALVO INDIVIDUAL
# ============================================================

st.subheader(
    "Proximidade read → referência (V3)"
)

target = st.selectbox(
    "Alvo",
    labels,
)


sr = summary_exec[
    summary_exec[
        "target_label"
    ]
    == target
].iloc[0]


r = v3_reads[
    (v3_reads["execution"] == execution)
    &
    (v3_reads["target_label"] == target)
].copy()


m1, m2, m3, m4 = st.columns(4)

m1.metric(
    "Reads mapeadas",
    int(
        sr[
            "sampled_reads_mapped"
        ]
    ),
)

m2.metric(
    "Divergência mediana",
    fmt(
        sr[
            "median_divergence_aligned_pct"
        ],
        2,
        "%",
    ),
)

m3.metric(
    "Cobertura mediana",
    fmt(
        sr[
            "median_query_coverage_pct"
        ],
        2,
        "%",
    ),
)

m4.metric(
    "Cobertura ≥80%",
    int(
        sr[
            "reads_coverage_ge_80"
        ]
    ),
)


if len(r):

    g1, g2 = st.columns(2)

    with g1:

        fig = px.histogram(
            r,
            x="divergence_aligned_pct",
            nbins=30,
            labels={
                "divergence_aligned_pct":
                    "Divergência alinhada (%)"
            },
        )

        st.plotly_chart(
            fig,
            use_container_width=True,
        )

    with g2:

        fig = px.scatter(
            r,
            x="query_coverage_pct",
            y="divergence_aligned_pct",
            hover_data=[
                "read",
                "qname",
                "AS",
                "MAPQ",
                "CIGAR",
                "NM",
            ],
            labels={
                "query_coverage_pct":
                    "Cobertura da read (%)",

                "divergence_aligned_pct":
                    "Divergência alinhada (%)",
            },
        )

        st.plotly_chart(
            fig,
            use_container_width=True,
        )


    with st.expander(
        "Tabela read-level V3"
    ):

        cols = [
            "read",
            "qname",
            "rname",
            "MAPQ",
            "CIGAR",
            "NM",
            "AS",
            "read_length",
            "query_aligned_bases",
            "alignment_columns",
            "divergence_aligned_pct",
            "query_coverage_pct",
        ]

        st.dataframe(
            r[
                [
                    c
                    for c in cols
                    if c in r.columns
                ]
            ],
            hide_index=True,
            use_container_width=True,
        )


# ============================================================
# CONSENSO × REFERÊNCIA
# ============================================================

cr = v4_cons_ref[
    (v4_cons_ref["execution"] == execution)
    &
    (v4_cons_ref["target_label"] == target)
]


if len(cr):

    cr = cr.iloc[0]

    st.markdown(
        "#### Consenso × referência — V4"
    )

    z1, z2, z3, z4 = st.columns(4)

    z1.metric(
        "Colunas comparáveis",
        int(
            cr[
                "comparable_columns"
            ]
        ),
    )

    z2.metric(
        "Divergência",
        fmt(
            cr[
                "divergence_consensus_reference_pct"
            ],
            2,
            "%",
        ),
    )

    z3.metric(
        "Identidade alinhada",
        fmt(
            cr[
                "identity_consensus_reference_pct"
            ],
            2,
            "%",
        ),
    )

    z4.metric(
        "Gap↔base interno",
        int(
            cr[
                "internal_gap_base_columns"
            ]
        ),
    )


# ============================================================
# READS QUE ALINHAM A DUAS REFERÊNCIAS
# ============================================================

st.subheader(
    "Comparação direta entre referências"
)

pair = v3_pair[
    v3_pair["execution"]
    == execution
].copy()


if pair.empty:

    st.info(
        "Nenhuma read da amostra R possui alinhamento "
        "V3 válido em duas referências nesta execução."
    )

else:

    pair["pair"] = (
        pair["target_A"]
        + " × "
        + pair["target_B"]
    )

    pair_name = st.selectbox(
        "Par de referências",
        sorted(
            pair["pair"].unique()
        ),
    )

    pp = pair[
        pair["pair"]
        == pair_name
    ].copy()

    A = pp[
        "target_A"
    ].iloc[0]

    B = pp[
        "target_B"
    ].iloc[0]


    threshold = st.slider(
        "Cobertura mínima da read em ambas as referências (%)",
        0,
        100,
        80,
        5,
    )


    pp = pp[
        pp[
            "min_query_coverage_pct"
        ]
        >= threshold
    ].copy()


    q1, q2, q3 = st.columns(3)

    q1.metric(
        "Reads compartilhadas",
        len(pp),
    )

    q2.metric(
        "Menor divergência",
        contagem(
            pp[
                "lower_divergence_reference"
            ]
        ),
    )

    q3.metric(
        "Maior AS",
        contagem(
            pp[
                "higher_AS_reference"
            ]
        ),
    )


    informative = pp[
        pp[
            "lower_divergence_reference"
        ].isin([A, B])
        &
        pp[
            "higher_AS_reference"
        ].isin([A, B])
    ]


    discord = (
        informative[
            "lower_divergence_reference"
        ]
        != informative[
            "higher_AS_reference"
        ]
    ).sum()


    st.metric(
        "Divergência × AS discordantes",
        (
            f"{int(discord)}/{len(informative)}"
            if len(informative)
            else "—"
        ),
    )


    if len(pp):

        med = pd.DataFrame(
            {
                "Referência": [
                    A,
                    B,
                ],

                "Divergência mediana (%)": [
                    pp[
                        "divergence_A_pct"
                    ].median(),

                    pp[
                        "divergence_B_pct"
                    ].median(),
                ],

                "Cobertura mediana (%)": [
                    pp[
                        "query_coverage_A_pct"
                    ].median(),

                    pp[
                        "query_coverage_B_pct"
                    ].median(),
                ],

                "AS mediano": [
                    pp["AS_A"].median(),
                    pp["AS_B"].median(),
                ],
            }
        )


        st.dataframe(
            med,
            hide_index=True,
            use_container_width=True,
        )


        fig = px.scatter(
            pp,
            x="divergence_A_pct",
            y="divergence_B_pct",
            hover_data=[
                "read",
                "qname",
                "query_coverage_A_pct",
                "query_coverage_B_pct",
                "AS_A",
                "AS_B",
                "CIGAR_A",
                "CIGAR_B",
                "lower_divergence_reference",
                "higher_AS_reference",
            ],
            labels={
                "divergence_A_pct":
                    f"Divergência em {A} (%)",

                "divergence_B_pct":
                    f"Divergência em {B} (%)",
            },
        )


        values = pd.concat(
            [
                pp["divergence_A_pct"],
                pp["divergence_B_pct"],
            ]
        ).dropna()


        if len(values):

            lo = float(
                values.min()
            )

            hi = float(
                values.max()
            )

            fig.add_shape(
                type="line",
                x0=lo,
                y0=lo,
                x1=hi,
                y1=hi,
                line=dict(
                    dash="dash"
                ),
            )


        st.plotly_chart(
            fig,
            use_container_width=True,
        )


        lower = pp[
            "lower_divergence_reference"
        ].value_counts()

        higher = pp[
            "higher_AS_reference"
        ].value_counts()


        st.info(
            f"Entre as reads compartilhadas acima do limiar selecionado, "
            f"a menor divergência ocorre em {A}: {int(lower.get(A, 0))}, "
            f"{B}: {int(lower.get(B, 0))}; "
            f"o maior AS ocorre em {A}: {int(higher.get(A, 0))}, "
            f"{B}: {int(higher.get(B, 0))}. "
            f"Divergência e AS são métricas complementares e podem discordar. "
            f"Essas diferenças não constituem isoladamente identificação "
            f"taxonômica nem demonstração de duas populações virais."
        )


        with st.expander(
            "Tabela das reads compartilhadas"
        ):

            st.dataframe(
                pp.drop(
                    columns=[
                        "pair"
                    ],
                    errors="ignore",
                ),
                hide_index=True,
                use_container_width=True,
            )


# ============================================================
# AUDITORIA V4
# ============================================================

with st.expander(
    "Resumo técnico V4"
):

    ss = v4_summary[
        v4_summary[
            "execution"
        ]
        == execution
    ]

    if len(ss):

        st.dataframe(
            ss,
            hide_index=True,
            use_container_width=True,
        )

    elif execution in exec_error:

        st.dataframe(
            v4_errors[
                v4_errors[
                    "execution"
                ]
                == execution
            ],
            hide_index=True,
            use_container_width=True,
        )


st.divider()

st.caption(
    "Mapeamento, divergência, AS, consenso, PCoA e classificação "
    "taxonômica são camadas distintas. A mesma read alinhando a "
    "referências relacionadas não demonstra, por si só, coinfecção "
    "ou dois haplótipos independentes."
)
