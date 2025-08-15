# align_1myw_1gfl.py
# Usage (headless):
#   pymol -cq align_1myw_1gfl.py
#
# Or inside an interactive PyMOL session:
#   run align_1myw_1gfl.py
#   align_1myw_1gfl()
##*or run fetch 1myw, type=cif, async=0; fetch 1gfl, type=cif, async=0; remove 1myw and resn HOH; remove 1gfl and resn HOH; hide everything, all; show cartoon, 1myw; show cartoon, 1gfl; color olive, 1myw; color marine, 1gfl; orient 1gfl; align 1myw and polymer.protein, 1gfl and polymer.protein; zoom 1myw or 1gfl; save aligned_1myw_1gfl.pse


from pymol import cmd

def align_1myw_1gfl():
    # Clean start
    cmd.reinitialize()

    # Fetch structures (mmCIF format; wait for completion)
    cmd.fetch("1myw", type="cif", async_=0)
    cmd.fetch("1gfl", type="cif", async_=0)

    # Remove waters
    cmd.remove("1myw and resn HOH")
    cmd.remove("1gfl and resn HOH")

    # Set cartoon representation and colors before alignment
    cmd.hide("everything", "all")
    cmd.show("cartoon", "1myw")
    cmd.show("cartoon", "1gfl")
    cmd.color("olive", "1myw")
    cmd.color("marine", "1gfl")

    # Orient to target
    cmd.orient("1gfl")

    # Align 1MYW onto 1GFL using protein atoms only
    mob = "1myw and polymer.protein"
    tgt = "1gfl and polymer.protein"
    rmsd = cmd.align(mob, tgt)[0]  # RMSD value
    print(f"Alignment complete. RMSD (Å): {rmsd:.3f}")

    # Zoom to both structures
    cmd.zoom("1myw or 1gfl")

    # Save the aligned session
    cmd.save("aligned_1myw_1gfl.pse")
    print("Saved aligned session as aligned_1myw_1gfl.pse")

if __name__ == "__main__":
    align_1myw_1gfl()
