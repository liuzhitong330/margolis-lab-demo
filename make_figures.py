"""
make_figures.py — Margolis Lab Demo
Figure 1: MOR / KOR / DOR expression density across addiction-relevant brain regions
           (Allen Mouse Brain Atlas ISH data) — showing the receptor heterogeneity
           that defines circuit-specific opioid pharmacology
Figure 2: Schematic of the 2025 KOR signaling reversal finding —
           stress switches KOR from inhibitory to excitatory in cortically-projecting VTA DA neurons
"""

import csv, os

OUT = os.path.dirname(os.path.abspath(__file__))

rows = []
with open(os.path.join(OUT, "opioid_expression.tsv")) as f:
    reader = csv.DictReader(f, delimiter="\t")
    genes = [c for c in reader.fieldnames if c != "region"]
    for row in reader:
        rows.append(row)

regions = [r["region"] for r in rows]

# ── Figure 1: Grouped bar chart — opioid receptor expression by region ────────
FW, FH   = 720, 460
PAD_L    = 80
PAD_R    = 160
PAD_T    = 90
PAD_B    = 75
AW = FW - PAD_L - PAD_R
AH = FH - PAD_T - PAD_B

GENE_COLORS = {
    "Oprm1 (MOR)": "#c0392b",  # red
    "Oprk1 (KOR)": "#1a5c8a",  # blue
    "Oprd1 (DOR)": "#27ae60",  # green
}

# max expression for scale
all_vals = [float(r[g]) for r in rows for g in genes]
y_max = max(all_vals) * 1.15

n_reg = len(regions)
n_gene = len(genes)
group_w = AW / n_reg
bar_w = group_w * 0.22
gap = group_w * 0.06

def px(i_reg, i_gene):
    return PAD_L + i_reg * group_w + (group_w - n_gene*bar_w - (n_gene-1)*gap) / 2 + i_gene * (bar_w + gap)

def py(v):
    return PAD_T + AH - v / y_max * AH

# Axis
ax1  = (f'<line x1="{PAD_L}" y1="{PAD_T}" x2="{PAD_L}" y2="{PAD_T+AH}" stroke="#ccc" stroke-width="1.2"/>'
        f'<line x1="{PAD_L}" y1="{PAD_T+AH}" x2="{PAD_L+AW}" y2="{PAD_T+AH}" stroke="#ccc" stroke-width="1.2"/>')
yticks = ""
for v in [0.000, 0.003, 0.006, 0.009, 0.012]:
    if v > y_max: break
    ty = py(v)
    yticks += (f'<line x1="{PAD_L-4}" y1="{ty:.1f}" x2="{PAD_L}" y2="{ty:.1f}" stroke="#aaa" stroke-width="1"/>'
               f'<text x="{PAD_L-8}" y="{ty+4:.1f}" text-anchor="end" font-size="9" fill="#888">{v:.3f}</text>')
yticks += (f'<text transform="rotate(-90,16,{PAD_T+AH/2:.0f})" x="16" y="{PAD_T+AH/2:.0f}" '
           f'text-anchor="middle" font-size="9.5" fill="#555">ISH expression density</text>')

bars1 = ""
HIGHLIGHT_REGIONS = {"VTA", "LHb", "MHb"}
for i_reg, row in enumerate(rows):
    reg = row["region"]
    cx = PAD_L + (i_reg + 0.5) * group_w
    is_hi = reg in HIGHLIGHT_REGIONS
    # region label
    bars1 += (f'<text x="{cx:.1f}" y="{PAD_T+AH+18}" text-anchor="middle" font-size="10" '
              f'fill="{"#222" if is_hi else "#666"}" font-weight="{"700" if is_hi else "400"}">{reg}</text>')
    if is_hi:
        bars1 += (f'<rect x="{PAD_L + i_reg*group_w + 2}" y="{PAD_T}" '
                  f'width="{group_w-4:.1f}" height="{AH}" fill="#f5f5f0" rx="2"/>')
    for i_gene, gene in enumerate(genes):
        val = float(row[gene])
        bx = px(i_reg, i_gene)
        bh = val / y_max * AH
        by = PAD_T + AH - bh
        col = GENE_COLORS[gene]
        bars1 += (f'<rect x="{bx:.1f}" y="{by:.1f}" width="{bar_w:.1f}" height="{bh:.1f}" '
                  f'fill="{col}" opacity="{"0.95" if is_hi else "0.65"}" rx="1"/>')

# Region label subtitles for highlighted ones
REGION_LABELS = {
    "VTA": "Margolis lab\nprimary target",
    "LHb": "Margolis lab\nfocus: withdrawal",
    "MHb": "highest MOR",
}
for i_reg, row in enumerate(rows):
    reg = row["region"]
    if reg in REGION_LABELS:
        cx = PAD_L + (i_reg + 0.5) * group_w
        for j, line in enumerate(REGION_LABELS[reg].split("\n")):
            bars1 += (f'<text x="{cx:.1f}" y="{PAD_T+AH+33+j*11}" text-anchor="middle" '
                      f'font-size="7.5" fill="#888">{line}</text>')

# Legend
leg1 = f'<rect x="{PAD_L+AW+12}" y="{PAD_T+10}" width="140" height="72" rx="4" fill="white" stroke="#ddd" stroke-width="1"/>'
for i, gene in enumerate(genes):
    col = GENE_COLORS[gene]
    lx = PAD_L + AW + 24
    ly = PAD_T + 28 + i * 22
    short = gene.split(" ")[1].strip("()")
    leg1 += (f'<rect x="{lx}" y="{ly-8}" width="12" height="12" fill="{col}" opacity="0.9" rx="1"/>'
             f'<text x="{lx+17}" y="{ly+2}" font-size="10" fill="{col}" font-weight="700">{gene}</text>')

svg1 = f"""<svg viewBox="0 0 {FW} {FH}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{FW//2}" y="22" text-anchor="middle" font-size="13" font-weight="700" fill="#222">
    Opioid Receptor Expression Across Addiction-Relevant Circuits
  </text>
  <text x="{FW//2}" y="40" text-anchor="middle" font-size="10" fill="#666">
    Allen Mouse Brain Atlas ISH · Oprm1 (MOR), Oprk1 (KOR), Oprd1 (DOR) · expression density
  </text>
  <text x="{FW//2}" y="57" text-anchor="middle" font-size="10" fill="#444">
    VTA is KOR-dominant with almost no MOR — habenula (LHb, MHb) carries the highest MOR load
  </text>
  <text x="{FW//2}" y="72" text-anchor="middle" font-size="9.5" fill="#888">
    Highlighted regions are primary targets of the Margolis lab
  </text>
  {ax1}{yticks}{bars1}{leg1}
</svg>"""

with open(os.path.join(OUT, "opioid_expression.svg"), "w") as f:
    f.write(svg1)
print("Wrote opioid_expression.svg")


# ── Figure 2: KOR signaling reversal schematic (2025 biorXiv finding) ──────────
FW2, FH2 = 700, 400
MID = FW2 // 2

# Layout: two panels side by side — Naive (left) vs Stress (right)
# Each panel: VTA DA neuron → KOR → result

PANEL_W = 290
PANEL_H = 260
PY_TOP = 110

LEFT_X  = 30
RIGHT_X = FW2 - LEFT_X - PANEL_W

def panel(px_start, label, stressed):
    col_head = "#c0392b" if stressed else "#1a5c8a"
    col_kor  = "#8e44ad"
    arrow_col = "#e67e22" if stressed else "#27ae60"
    effect = "EXCITATORY ↑" if stressed else "INHIBITORY ↓"
    eff_col  = "#c0392b" if stressed else "#27ae60"
    eff_bg   = "#fff0f0" if stressed else "#f0fff4"
    eff_border = "#c0392b" if stressed else "#27ae60"

    cx = px_start + PANEL_W // 2
    # Panel box
    s  = (f'<rect x="{px_start}" y="{PY_TOP-18}" width="{PANEL_W}" height="{PANEL_H}" '
          f'rx="8" fill="{"#fff8f0" if stressed else "#f0f8ff"}" stroke="{col_head}" stroke-width="1.5"/>')
    # Condition label
    s += (f'<text x="{cx}" y="{PY_TOP+5}" text-anchor="middle" font-size="12" '
          f'font-weight="700" fill="{col_head}">{label}</text>')

    # Neuron body
    ny = PY_TOP + 70
    s += (f'<circle cx="{cx}" cy="{ny}" r="32" fill="{col_head}" opacity="0.12"/>'
          f'<circle cx="{cx}" cy="{ny}" r="32" fill="none" stroke="{col_head}" stroke-width="1.5"/>'
          f'<text x="{cx}" y="{ny-8}" text-anchor="middle" font-size="9" fill="{col_head}" font-weight="700">VTA DA</text>'
          f'<text x="{cx}" y="{ny+6}" text-anchor="middle" font-size="8.5" fill="{col_head}">neuron</text>'
          f'<text x="{cx}" y="{ny+19}" text-anchor="middle" font-size="7.5" fill="#888">(cortically-projecting)</text>')

    # KOR receptor on neuron surface
    kor_x = cx + 32; kor_y = ny - 10
    s += (f'<rect x="{kor_x-12}" y="{kor_y-10}" width="34" height="20" rx="4" '
          f'fill="{col_kor}" opacity="0.15"/>'
          f'<rect x="{kor_x-12}" y="{kor_y-10}" width="34" height="20" rx="4" '
          f'fill="none" stroke="{col_kor}" stroke-width="1"/>'
          f'<text x="{kor_x+5}" y="{kor_y+4}" text-anchor="middle" font-size="8.5" '
          f'fill="{col_kor}" font-weight="700">KOR</text>')

    # Arrow down from neuron to effect box
    arr_y1 = ny + 33; arr_y2 = ny + 70
    s += (f'<line x1="{cx}" y1="{arr_y1}" x2="{cx}" y2="{arr_y2}" '
          f'stroke="{arrow_col}" stroke-width="2" marker-end="url(#tri)"/>'
          f'<text x="{cx+6}" y="{arr_y1+25}" font-size="8.5" fill="{arrow_col}">'
          f'{"reversal" if stressed else "canonical"}</text>')

    # Effect box
    eff_y = ny + 72
    s += (f'<rect x="{px_start+30}" y="{eff_y}" width="{PANEL_W-60}" height="34" '
          f'rx="5" fill="{eff_bg}" stroke="{eff_border}" stroke-width="1.5"/>'
          f'<text x="{cx}" y="{eff_y+21}" text-anchor="middle" font-size="11" '
          f'fill="{eff_col}" font-weight="700">{effect}</text>')
    return s

svg2_body = panel(LEFT_X, "Naive state", False) + panel(RIGHT_X, "After stress", True)

# Dynorphin/KOR agonist label between panels
mid_y = PY_TOP + 70
svg2_body += (f'<text x="{MID}" y="{mid_y-15}" text-anchor="middle" font-size="10" fill="#8e44ad" font-weight="700">Dynorphin (KOR agonist)</text>'
              f'<line x1="{LEFT_X+PANEL_W+8}" y1="{mid_y}" x2="{RIGHT_X-8}" y2="{mid_y}" '
              f'stroke="#8e44ad" stroke-width="1.5" stroke-dasharray="5,3"/>'
              f'<text x="{MID}" y="{mid_y+15}" text-anchor="middle" font-size="9" fill="#888">same receptor, opposite output</text>')

svg2 = f"""<svg viewBox="0 0 {FW2} {FH2}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{FW2//2}" y="22" text-anchor="middle" font-size="13" font-weight="700" fill="#222">
    Stress Flips KOR Signaling in Cortically-Projecting VTA Dopamine Neurons
  </text>
  <text x="{FW2//2}" y="40" text-anchor="middle" font-size="10" fill="#666">
    Schematic of Margolis EB (2025) biorXiv — same KOR, opposite effect after stress
  </text>
  <text x="{FW2//2}" y="57" text-anchor="middle" font-size="10" fill="#444">
    Kappa opioid receptors normally inhibit DA neuron firing; stress reverses this to excitation
  </text>
  <text x="{FW2//2}" y="72" text-anchor="middle" font-size="9" fill="#888">
    Only in a subset of DA neurons projecting to prefrontal cortex — not in all VTA DA neurons
  </text>
  {svg2_body}
  <text x="{FW2//2}" y="{FH2-10}" text-anchor="middle" font-size="8.5" fill="#aaa">
    Margolis EB (2025) biorXiv · "A switch in kappa opioid receptor signaling from inhibitory to excitatory induced by stress"
  </text>
</svg>"""

with open(os.path.join(OUT, "kor_signaling.svg"), "w") as f:
    f.write(svg2)
print("Wrote kor_signaling.svg")
