.PHONY: all figs clean

all: figs

figs: figures/fig1.png figures/fig2.png figures/fig3.png figures/fig4.png

figures/fig1.png:
	python scripts/Figure1_Save.py

figures/fig2.png:
	python scripts/Figure2_Save.py

figures/fig3.png:
	pymol -cq scripts/Figure3_Render.pml

figures/fig4.png:
	python scripts/Figure4_Save.py

clean:
	rm -rf figures/* data/processed/* results/*