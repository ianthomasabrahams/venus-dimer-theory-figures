# Venus Dimer Theory Figures (Exact Scripts + Non‑interactive Saving)
*Excitonic Coupling and Photon Antibunching in Venus Yellow Fluorescent Protein Dimers: A Lindblad Master Equation Approach*

This repository contains the **exact figure scripts** (contents unchanged, only filenames cleaned) and thin, non‑interactive **save wrappers** so figures are written to disk deterministically.

## Figures & scripts
- **Figure 1:** `scripts/Figure1_VenusHomodimer_Dephasing_PhaseHalf_fs.py` → wrapper `scripts/Figure1_Save.py` saves **`figures/fig1.png`**
- **Figure 2:** `scripts/Figure2_VenusHeterodimer_EnergyMixedState_StokesShift_ps.py` → wrapper `scripts/Figure2_Save.py` saves **`figures/fig2.png`**
- **Figure 3:** `scripts/Figure3_PyMOL_Align1MYW_1GFL.py` (PyMOL) → renderer `scripts/Figure3_Render.pml` saves **`figures/fig3.png`**
- **Figure 4:** `scripts/Figure4_Cryo_ThermalVsDephasing_2x2_ps.py` → wrapper `scripts/Figure4_Save.py` saves **`figures/fig4.png`**

## Quick start
```bash
conda env create -f environment.yml
conda activate venus-dimer
make all
```
This will create `figures/fig1.png` … `fig4.png`.

## PyMOL (Figure 3)
The alignment script fetches **1MYW** and **1GFL** from RCSB (internet required). The Makefile runs:
```bash
pymol -cq scripts/Figure3_Render.pml
```
If you prefer direct CLI:
```bash
pymol -cq -d "run scripts/Figure3_PyMOL_Align1MYW_1GFL.py; align_1myw_1gfl(); png figures/fig3.png, 300; quit"
```
For fully offline reproduction, place the CIF files locally and adjust the PyMOL script to `load` instead of `fetch`.

## Data
If any script references inputs, put them in `data/raw/` or set the DOI in a fetch helper and run `make data`.

## Reproducibility
- Exact scripts are preserved verbatim.
- Saving is handled by tiny wrapper scripts using a non‑interactive Matplotlib backend.
- Consider exporting a lock file from your machine and overwrite `environment-explicit.txt`:
  ```bash
  conda list --explicit > environment-explicit.txt
  ```

## License
- Code: MIT (see `LICENSE`)
- Figures/data: CC BY 4.0 (see `LICENSE-FIGURES`)