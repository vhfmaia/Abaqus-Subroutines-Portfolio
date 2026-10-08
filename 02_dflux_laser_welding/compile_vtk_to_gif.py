"""
compile_vtk_to_gif.py
=====================
Multi-Viewport 3DEXPERIENCE / SIMULIA Animated Visual Report Compiler
for Circumferential Laser Welding Simulation.

Generates synchronized multi-viewport animated GIFs from exported VTK frames
(or from laser_therm_vtk.zip / laser_mech_vtk.zip).

Layout:
  1. Viewport 1 (Left)   : 3D Isometric Global View (European 1st Angle Projection)
  2. Viewport 2 (Center) : Top View (X-Y Plane along Z axis) - Laser Path & HAZ
  3. Viewport 3 (Right)  : 3D Cross-Section (R-Z Cutaway) - Penetration Depth & Melt Pool
"""

import os
import sys
import glob
import time
import shutil
import zipfile
import tempfile
import argparse
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
import matplotlib.colors as mcolors

# -----------------------------------------------------------------------------
# 1. 3DEXPERIENCE Standard Colormap Palette & Color Bands
# -----------------------------------------------------------------------------
COLORMAP_3DX = [
    "#0000ff",  # Blue (Lowest)
    "#0055ff",
    "#00aaff",
    "#00ffff",  # Cyan
    "#00ffaa",
    "#00ff55",
    "#00ff00",  # Green
    "#aaff00",
    "#ffff00",  # Yellow
    "#ffaa00",
    "#ff5500",
    "#ff0000"   # Red (Highest)
]

def parse_args():
    parser = argparse.ArgumentParser(description="Compile VTK frames into multi-viewport GIF")
    parser.add_argument("input", nargs="?", default=None, help="Input VTK file, directory, or ZIP archive")
    parser.add_argument("--output", "-o", default=None, help="Output GIF file path")
    parser.add_argument("--field", "-f", default=None, choices=["auto", "Temperature", "Von_Mises", "PEEQ"],
                        help="Scalar field to visualize (default: auto-detect)")
    parser.add_argument("--fps", type=int, default=15, help="Animation playback FPS")
    parser.add_argument("--max-val", type=float, default=None, help="Max scalar value for color normalization")
    parser.add_argument("--min-val", type=float, default=None, help="Min scalar value for color normalization")
    return parser.parse_args()

def locate_and_extract_vtk(input_path):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = []
    if input_path:
        candidates.append(input_path)
    candidates.extend([
        os.path.join(script_dir, "laser_therm_vtk.zip"),
        os.path.join(script_dir, "laser_mech_vtk.zip"),
        os.path.join(script_dir, "gear_helical_vtk.zip")
    ])

    zip_file = None
    for cand in candidates:
        if cand and os.path.exists(cand) and cand.endswith(".zip"):
            zip_file = cand
            break

    temp_dir = None
    if zip_file:
        print("[INFO] Found VTK archive: %s" % zip_file)
        temp_dir = tempfile.mkdtemp(prefix="vtk_weld_")
        with zipfile.ZipFile(zip_file, 'r') as zf:
            zf.extractall(temp_dir)
        vtk_files = sorted(glob.glob(os.path.join(temp_dir, "*.vtk")))
    else:
        search_dir = input_path if (input_path and os.path.isdir(input_path)) else script_dir
        vtk_files = sorted(glob.glob(os.path.join(search_dir, "*.vtk")))

    return vtk_files, temp_dir, zip_file

def main():
    args = parse_args()
    vtk_files, temp_dir, zip_source = locate_and_extract_vtk(args.input)

    if not vtk_files:
        print("[ERROR] No VTK files found to process. Please provide an ODB or VTK zip.")
        sys.exit(1)

    print("[INFO] Found %d VTK frames." % len(vtk_files))

    # Inspect first frame to detect active scalar field
    sample_mesh = pv.read(vtk_files[0])
    point_arrays = list(sample_mesh.point_data.keys())
    cell_arrays = list(sample_mesh.cell_data.keys())
    print("[INFO] Point fields: %s" % point_arrays)
    print("[INFO] Cell fields: %s" % cell_arrays)

    active_field = args.field
    if active_field is None or active_field == "auto":
        if "Temperature" in point_arrays:
            active_field = "Temperature"
        elif "Von_Mises" in cell_arrays or "Von_Mises" in point_arrays:
            active_field = "Von_Mises"
        elif "PEEQ" in cell_arrays or "PEEQ" in point_arrays:
            active_field = "PEEQ"
        else:
            all_fields = point_arrays + cell_arrays
            active_field = all_fields[0] if all_fields else None

    if not active_field:
        print("[ERROR] No visualizable scalar field detected in VTK.")
        sys.exit(1)

    print("[INFO] Visualizing field: '%s'" % active_field)

    # Determine scalar bounds
    min_val = args.min_val if args.min_val is not None else 25.0
    if args.max_val is not None:
        max_val = args.max_val
    elif active_field == "Temperature":
        max_val = 1600.0  # Melting / peak weld temperature
        cmap = "inferno"
    elif active_field == "Von_Mises":
        min_val = 0.0 if args.min_val is None else args.min_val
        max_val = 850.0   # 16MnCr5 ultimate tensile ceiling
        cmap = COLORMAP_3DX
    elif active_field == "PEEQ":
        min_val = 0.0 if args.min_val is None else args.min_val
        max_val = 0.15
        cmap = COLORMAP_3DX
    else:
        min_val = 0.0 if args.min_val is None else args.min_val
        max_val = 1.0
        cmap = COLORMAP_3DX

    # Output filename
    if args.output:
        out_gif = args.output
    elif zip_source:
        base = os.path.splitext(os.path.basename(zip_source))[0].replace("_vtk", "")
        out_gif = f"{base}_simulation.gif"
    else:
        out_gif = "laser_welding_simulation.gif"

    # PyVista Off-screen Multi-Viewport Renderer
    pv.set_plot_theme("document")
    plotter = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(1800, 600))
    plotter.set_background("white")

    rendered_images = []
    print(f"[INFO] Rendering {len(vtk_files)} frames to '{out_gif}'...")

    for idx, vtk_p in enumerate(vtk_files):
        mesh = pv.read(vtk_p)

        # Clear previous actors
        plotter.subplot(0, 0)
        plotter.clear_actors()
        plotter.add_mesh(mesh, scalars=active_field, cmap=cmap, clim=[min_val, max_val],
                         show_edges=False, smooth_shading=True, show_scalar_bar=False)
        plotter.camera_position = [(38.0, -45.0, 35.0), (0.0, 0.0, 10.0), (0.0, 0.0, 1.0)]
        plotter.add_text("Isometric 3D View", font_size=11, position="upper_left", color="black")

        plotter.subplot(0, 1)
        plotter.clear_actors()
        plotter.add_mesh(mesh, scalars=active_field, cmap=cmap, clim=[min_val, max_val],
                         show_edges=False, smooth_shading=True, show_scalar_bar=False)
        plotter.view_xy()
        plotter.camera.zoom(1.15)
        plotter.add_text("Top View (XY Flank)", font_size=11, position="upper_left", color="black")

        plotter.subplot(0, 2)
        plotter.clear_actors()
        # Half cutaway plane (Y >= 0)
        cut_mesh = mesh.clip(normal=[0, 1, 0], origin=[0, 0, 0], invert=True)
        plotter.add_mesh(cut_mesh, scalars=active_field, cmap=cmap, clim=[min_val, max_val],
                         show_edges=False, smooth_shading=True, show_scalar_bar=True,
                         scalar_bar_args={"title": f"{active_field}", "n_labels": 5, "vertical": True})
        plotter.camera_position = [(0.0, -50.0, 10.0), (0.0, 0.0, 10.0), (0.0, 0.0, 1.0)]
        plotter.add_text("RZ Cross-Section Cut", font_size=11, position="upper_left", color="black")

        img = plotter.screenshot(return_img=True)
        rendered_images.append(Image.fromarray(img))

    plotter.close()

    # Save animated GIF using Pillow
    if rendered_images:
        print(f"[INFO] Saving animated GIF ({len(rendered_images)} frames) at {args.fps} FPS...")
        duration_ms = int(1000.0 / args.fps)
        rendered_images[0].save(
            out_gif,
            save_all=True,
            append_images=rendered_images[1:],
            duration=duration_ms,
            loop=0,
            optimize=True
        )
        gif_size_mb = os.path.getsize(out_gif) / (1024.0 * 1024.0)
        print(f"[SUCCESS] Multi-viewport GIF generated: {out_gif} ({gif_size_mb:.2f} MB)")

    # Cleanup temp directory
    if temp_dir and os.path.exists(temp_dir):
        shutil.rmtree(temp_dir, ignore_errors=True)

if __name__ == "__main__":
    main()
