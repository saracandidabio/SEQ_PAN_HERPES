#!/usr/bin/env python3
from pathlib import Path
import argparse
import shutil
import sys

FILES = {
    "dashboard_final/execution_evidence_final_v2.tsv": "execution_evidence.tsv",
    "dashboard_final/target_evidence_final.tsv": "target_evidence.tsv",
    "dashboard_final/amac_apal_final.tsv": "amac_apal.tsv",
    "dashboard_final/dashboard_qc_summary.tsv": "qc_summary.tsv",
    "dashboard_final/taxonomy_source_audit.tsv": "taxonomy_source_audit.tsv",
    "competitive_mapping/ambiguity_summary/execution_ambiguity_summary.tsv": "execution_ambiguity.tsv",
    "competitive_mapping/ambiguity_summary/pair_ambiguity_global.tsv": "pair_ambiguity_global.tsv",
    "dataset_valido/competitive_inputs/dedup_execution_summary.tsv": "dedup_summary.tsv",
    "competitive_mapping/combined_7_references_metadata.tsv": "references.tsv",
    "competitive_mapping/competitive_execution_summary.tsv": "competitive_execution_summary.tsv",
    "competitive_mapping/competitive_target_summary.tsv": "competitive_target_summary.tsv",
    "competitive_mapping/competitive_pairwise_overlap.tsv": "competitive_pairwise_overlap.tsv",
}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--dados",
        default="/Users/navio01/Downloads/PanHerpes_dashboard_data",
    )
    args = ap.parse_args()

    src = Path(args.dados).expanduser().resolve()
    dst = Path(__file__).resolve().parent / "data"
    dst.mkdir(exist_ok=True)

    missing = []

    for rel, name in FILES.items():
        s = src / rel
        d = dst / name

        if not s.exists():
            missing.append(rel)
            continue

        shutil.copy2(
            s,
            d,
        )

    print("PANHERPES - DADOS AGREGADOS")
    print("=" * 70)
    print(f"Copiados: {len(FILES)-len(missing)}/{len(FILES)}")
    print(f"Ausentes: {len(missing)}")

    for x in missing:
        print("MISSING:", x)

    if missing:
        sys.exit(1)

    print("\nOK.")
    print("Agora rode:")
    print(
        'Rscript preparar_reads_metodo_R.R '
        '"/Users/navio01/Downloads/PanHerpes_dashboard_data" "./data"'
    )

if __name__ == "__main__":
    main()
