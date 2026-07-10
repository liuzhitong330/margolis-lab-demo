"""
fetch_data.py — Margolis Lab Demo (Elyssa Margolis, UCSF)
Pulls ISH expression density for opioid receptor genes (MOR/KOR/DOR) across
addiction-relevant brain regions from the Allen Mouse Brain Atlas (ABA) API.

Key structures: VTA, lateral habenula, nucleus accumbens, SNc, PFC, amygdala.
The Margolis lab works on how opioid receptor heterogeneity across these circuits
shapes motivated behavior and addiction pharmacology.
"""

import urllib.request, json, csv, os

OUT = os.path.dirname(os.path.abspath(__file__))

# Brain regions relevant to the Margolis lab's work
STRUCTS = {
    "VTA":   749,   # Ventral tegmental area — lab's primary target
    "LHb":   186,   # Lateral habenula — opioid withdrawal, aversion
    "MHb":   483,   # Medial habenula — high MOR expression
    "NAc":   56,    # Nucleus accumbens — reward output
    "SNc":   374,   # Substantia nigra, compact — DA neurons
    "PFC":   972,   # Prelimbic cortex — top-down control
    "BLA":   295,   # Basolateral amygdala — aversive learning
    "vHipp": 407,   # Hippocampus CA1 — context
    "Hyp":   1097,  # Lateral hypothalamus — stress/homeostasis
    "PAG":   795,   # Periaqueductal gray — pain/stress
    "IPN":   100,   # Interpeduncular nucleus — habenula output
}

# Allen Brain Atlas ISH experiments (adult mouse, coronal)
EXPERIMENTS = {
    "Oprm1 (MOR)": 69860840,   # mu opioid receptor
    "Oprk1 (KOR)": 80517217,   # kappa opioid receptor
    "Oprd1 (DOR)": 100144385,  # delta opioid receptor
}

print("Querying Allen Brain Atlas ISH expression data …")
results = {name: {} for name in STRUCTS}
for gene, exp_id in EXPERIMENTS.items():
    print(f"  {gene} (experiment {exp_id}) …")
    for name, struct_id in STRUCTS.items():
        url = (f"https://api.brain-map.org/api/v2/data/query.json?"
               f"criteria=model::StructureUnionize,"
               f"rma::criteria,section_data_set[id$eq{exp_id}],"
               f"structure[id$eq{struct_id}]&num_rows=1")
        try:
            with urllib.request.urlopen(url, timeout=20) as r:
                d = json.loads(r.read())
            rows = d.get("msg", [])
            val = float(rows[0].get("expression_density", 0) or 0) if rows else 0.0
            results[name][gene] = round(val, 6)
        except Exception as e:
            print(f"    WARNING: {name} failed — {e}")
            results[name][gene] = 0.0

with open(os.path.join(OUT, "opioid_expression.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["region"] + list(EXPERIMENTS.keys()))
    for name in STRUCTS:
        w.writerow([name] + [results[name][g] for g in EXPERIMENTS])

print(f"\nWrote opioid_expression.tsv")
print(f"\nKey observations:")
genes = list(EXPERIMENTS.keys())
for region in ["VTA", "LHb", "MHb", "NAc"]:
    vals = [results[region][g] for g in genes]
    print(f"  {region}: MOR={vals[0]:.5f}  KOR={vals[1]:.5f}  DOR={vals[2]:.5f}")
