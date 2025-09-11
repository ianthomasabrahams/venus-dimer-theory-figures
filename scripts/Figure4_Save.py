#!/usr/bin/env python3
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("figures", exist_ok=True)

# Execute the original, unmodified script in this process
g = {}
with open("scripts/Fig4_cryo_coherence_1x2_v2.0.0-arxiv.py", "r", encoding="utf-8") as fh:
    code = fh.read()
exec(compile(code, "scripts/Fig4_cryo_coherence_1x2_v2.0.0-arxiv.py", "exec"), g, g)

# Save the current active figure to disk
plt.gcf().savefig("figures/fig4.png", dpi=300, bbox_inches="tight")
print("Wrote figures/fig4.png")