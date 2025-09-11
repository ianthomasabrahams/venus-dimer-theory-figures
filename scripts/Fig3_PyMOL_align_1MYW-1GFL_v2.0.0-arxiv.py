# align_1myw_1gfl_v2.0.0-arxiv.py
# Headless usage (recommended):
#   pymol -cq fig3_pymol_align_1myw-1gfl_v2.0.0-arxiv.py
#
# This script:
#  - Fetches 1MYW and 1GFL (mmCIF)
#  - Aligns 1MYW onto 1GFL using protein atoms
#  - Sets a clean cartoon style with colors (Venus: olive; GFP: marine)
#  - Prints the git commit hash (if available)
#  - Saves: aligned_1myw_1gfl_v2.0.0-arxiv.png and aligned_1myw_1gfl_v2.0.0-arxiv.pse

from pymol import cmd
import subprocess

SUFFIX = "v2.0.0-arxiv"

def _get_git_commit_short():
    try:
        h = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.STDOUT)
        return h.decode().strip()
    except Exception:
        return "unknown"

def align_1myw_1gfl():
    commit = _get_git_commit_short()
    print("Git commit:", commit)

    # Clean start
    cmd.reinitialize()

    # Appearance
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("ray_trace_mode", 1)
    cmd.set("antialias", 2)
    cmd.set("cartoon_transparency", 0.0)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("depth_cue", 0)
    cmd.set("orthoscopic", 1)
    cmd.set("ray_shadows", 0)

    # Fetch structures (mmCIF format; wait for completion)
    cmd.fetch("1myw", type="cif", async_=0)
    cmd.fetch("1gfl", type="cif", async_=0)

    # Remove waters
    cmd.remove("1myw and resn HOH")
    cmd.remove("1gfl and resn HOH")

    # Representations & colors
    cmd.hide("everything", "all")
    cmd.show("cartoon", "1myw")
    cmd.show("cartoon", "1gfl")
    cmd.color("olive", "1myw")
    cmd.color("marine", "1gfl")

    # Align 1MYW onto 1GFL using protein atoms
    mob = "1myw and polymer.protein"
    tgt = "1gfl and polymer.protein"
    rmsd = cmd.align(mob, tgt)[0]  # RMSD value
    print(f"Alignment complete. RMSD (Å): {rmsd:.3f}")

    # View
    cmd.orient("1gfl")
    cmd.zoom("1myw or 1gfl")

    # Optional: create a small dummy label with commit hash near the corner (non-intrusive)
    try:
        cmd.pseudoatom("commit_label", pos=[0,0,0])
        cmd.hide("spheres", "commit_label")
        cmd.label("commit_label", f"'commit: {commit}'")
        cmd.set("label_position", [0, -1, 0], "commit_label")
        cmd.set("label_size", 12, "commit_label")
        cmd.set("label_color", "grey50", "commit_label")
    except Exception:
        pass

    # Save image and session
    png_name = f"aligned_1myw_1gfl_{SUFFIX}.png"
    pse_name = f"aligned_1myw_1gfl_{SUFFIX}.pse"
    cmd.png(png_name, width=2200, height=1600, dpi=300, ray=1)
    cmd.save(pse_name)

    print("Saved:", png_name)
    print("Saved:", pse_name)

if __name__ == "__main__":
    align_1myw_1gfl()