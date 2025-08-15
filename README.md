# Venus Dimer Theory Figures (Exact Scripts + Non-interactive Saving)
*Excitonic Coupling and Photon Antibunching in Venus Yellow Fluorescent Protein Dimers: A Lindblad Master Equation Approach*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16882983.svg)](https://doi.org/10.5281/zenodo.16882983)

This repository contains the **exact figure scripts** (contents unchanged, only filenames cleaned) and thin, non-interactive **save-wrappers** so figures are written to disk deterministically (no GUI).

> Tested with: Python 3.11, Conda 24.x, PyMOL 2.5+

## Figures & scripts
- **Figure 1:** `scripts/Figure1_VenusHomodimer_Dephasing_PhaseHalf_fs.py` → wrapper `scripts/Figure1_Save.py` saves **`figures/fig1.png`**
- **Figure 2:** `scripts/Figure2_VenusHeterodimer_EnergyMixedState_StokesShift_ps.py` → wrapper `scripts/Figure2_Save.py` saves **`figures/fig2.png`**
- **Figure 3:** `scripts/Figure3_PyMOL_Align1MYW_1GFL.py` (PyMOL) → renderer `scripts/Figure3_Render.pml` saves **`figures/fig3.png`**  
  Uses **1MYW** (Venus) and **1GFL** (avGFP) from RCSB PDB.
- **Figure 4:** `scripts/Figure4_Cryo_ThermalVsDephasing_2x2_ps.py` → wrapper `scripts/Figure4_Save.py` saves **`figures/fig4.png`**

## Quick start
```bash
conda env create -f environment.yml
conda activate venus-dimer
make all
```

## Reproducibility
- Exact scripts are preserved verbatim.
- Saving is handled by tiny wrapper scripts using a non-interactive Matplotlib backend (`Agg`).
- Git commit hash for this archived version (v1.0.1-arxiv): `13df46c2bf6acc57017e7925999b6ba0d8abe194` (`13df46c` short hash).
- Consider exporting a lock file from your machine and overwrite `environment-explicit.txt`:
  ```bash
  conda list --explicit > environment-explicit.txt
  ```

## How to cite
If you use this code, please cite the archived release:

Abrahams, I. T. (2025). *venus-dimer-theory-figures (v1.0.1-arxiv): Exact figure scripts*.  
Zenodo. https://doi.org/10.5281/zenodo.16882983

**BibTeX**
```bibtex
@software{abrahams_venus_dimer_2025,
  author  = {Abrahams, Ian T.},
  title   = {venus-dimer-theory-figures (v1.0.1-arxiv): Exact figure scripts},
  year    = {2025},
  publisher = {Zenodo},
  version = {v1.0.1-arxiv},
  doi     = {10.5281/zenodo.16882983},
  url     = {https://doi.org/10.5281/zenodo.16882983}
}
```

## License
- Code: MIT (see `LICENSE`)
- Figures/data: CC BY 4.0 (see `LICENSE-FIGURES`)
