# Venus Dimer Theory Figures (Exact Scripts + Non-interactive Saving)
*Excitonic Coupling and Photon Antibunching in Venus Yellow Fluorescent Protein Dimers: A Lindblad Master Equation Approach*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16887703.svg)](https://doi.org/10.5281/zenodo.16887703)

This repository contains the **exact figure scripts** (contents unchanged, only filenames cleaned) and thin, non-interactive **save-wrappers** so figures are written to disk deterministically (no GUI).

> Tested with: Python 3.11, Conda 24.x, PyMOL 2.5+

---

## Figures & scripts
- **Figure 1:**  
  `scripts/Figure1_VenusHomodimer_Dephasing_PhaseHalf_fs.py`  
   → wrapper: `scripts/Figure1_Save.py` → saves **`figures/fig1.png`**

- **Figure 2:**  
  `scripts/Figure2_VenusHeterodimer_EnergyMixedState_StokesShift_ps.py`  
   → wrapper: `scripts/Figure2_Save.py` → saves **`figures/fig2.png`**

- **Figure 3:**  
  `scripts/Figure3_PyMOL_Align1MYW_1GFL.py` (PyMOL)  
   → renderer: `scripts/Figure3_Render.pml` → saves **`figures/fig3.png`**  
   (Uses **1MYW** (Venus) and **1GFL** (avGFP) from RCSB PDB.)

- **Figure 4:**  
  `scripts/Figure4_Cryo_ThermalVsDephasing_2x2_ps.py`  
   → wrapper: `scripts/Figure4_Save.py` → saves **`figures/fig4.png`**

---

## Quick start
```bash
conda env create -f environment.yml
conda activate venus-dimer
make all
```
---

## Reproducibility
- Exact scripts are preserved verbatim.  
- Saving is handled by tiny wrapper scripts using a non-interactive Matplotlib backend (`Agg`).  
- Git commit for this archived version (v1.2.0-arxiv):  
  - Short hash: `45ada66`  
  - Full hash: `45ada6694c9288de6117479710c712b48000d939` (this 40-character hash matches the `commit:` field in `CITATION.cff`).  
- For complete reproducibility, export a lock file from your machine to overwrite `environment-explicit.txt`:  

```bash
conda list --explicit > environment-explicit.txt
```

---

## How to cite

If you use this code, please cite the archived release:

Abrahams, I. T. (2025). venus-dimer-theory-figures (v1.2.0-arxiv): Exact figure scripts.
Zenodo. https://doi.org/10.5281/zenodo.16892071

BibTeX

@software{abrahams_venus_dimer_2025,
  author    = {Abrahams, Ian T.},
  title     = {venus-dimer-theory-figures (v1.2.0-arxiv): Exact figure scripts},
  year      = {2025},
  publisher = {Zenodo},
  version   = {v1.2.0-arxiv},
  doi       = {10.5281/zenodo.16892071},
  url       = {https://doi.org/10.5281/zenodo.16892071}
}

---

## License

Code: MIT (see LICENSE)

Figures/data: CC BY 4.0 (see LICENSE-FIGURES)

