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
    parser.add_argument("--stride", "-s", type=int, default=2, help="Frame step stride (default: 2)")
    parser.add_argument("--scale", type=float, default=0.75, help="Output resolution scale factor (default: 0.75)")
    parser.add_argument("--max-val", type=float, default=None, help="Max scalar value for color normalization")
    parser.add_argument("--min-val", type=float, default=None, help="Min scalar value for color normalization")
    return parser.parse_args()

def extract_frame_metadata(vtk_path):
    time_val = 0.0
    step_name = ""
    try:
        with open(vtk_path, "r") as f:
            f.readline()
            line2 = f.readline()
            parts = line2.split()
            if "Time" in parts:
                t_idx = parts.index("Time")
                time_val = float(parts[t_idx + 1])
            if "Step" in parts:
                s_idx = parts.index("Step")
                step_name = parts[s_idx + 1]
    except Exception:
        pass
    return time_val, step_name

def locate_and_extract_vtk(input_path):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = []
    if input_path:
        candidates.append(input_path)
    candidates.extend([
        os.path.join(script_dir, "disk_heatsource_th_vtk.zip"),
        os.path.join(script_dir, "laser_therm_vtk.zip"),
        os.path.join(script_dir, "disk_heatsource_me_vtk.zip"),
        os.path.join(script_dir, "laser_mech_vtk.zip"),
        os.path.join(script_dir, "gear_helical_vtk.zip")
    ])
    candidates.extend(sorted(glob.glob(os.path.join(script_dir, "*_vtk.zip"))))

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

    total_available = len(vtk_files)
    if args.stride > 1:
        vtk_files = vtk_files[::args.stride]
        print("[INFO] Applied stride %d: using %d / %d frames." % (args.stride, len(vtk_files), total_available))
    else:
        print("[INFO] Processing all %d VTK frames." % len(vtk_files))

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

    # Font setup for HUD banner
    banner_h = 42
    try:
        font = ImageFont.truetype("arial.ttf", 16)
        font_bold = ImageFont.truetype("arialbd.ttf", 16)
    except Exception:
        font = ImageFont.load_default()
        font_bold = font

    rendered_images = []
    print(f"[INFO] Rendering {len(vtk_files)} frames to '{out_gif}'...")
    t0 = time.time()

    for idx, vtk_p in enumerate(vtk_files):
        sim_time, step_name = extract_frame_metadata(vtk_p)
        mesh = pv.read(vtk_p)

        # Detect peak position for local melt pool tracking
        scalars_arr = mesh.point_data.get(active_field, None)
        if scalars_arr is not None and len(scalars_arr) > 0:
            max_idx = np.argmax(scalars_arr)
            peak_val = float(scalars_arr[max_idx])
            peak_x, peak_y, peak_z = mesh.points[max_idx]
            theta = np.arctan2(peak_y, peak_x)
        else:
            peak_val = 0.0
            peak_x, peak_y, peak_z = 0.0, 0.0, 20.0
            theta = 0.0

        # Subplot 0: 3D Isometric Global View
        plotter.subplot(0, 0)
        plotter.clear_actors()
        plotter.add_mesh(mesh, scalars=active_field, cmap=cmap, clim=[min_val, max_val],
                         show_edges=False, smooth_shading=True, show_scalar_bar=False)
        plotter.camera_position = [(80.0, -80.0, 60.0), (0.0, 0.0, 8.0), (0.0, 0.0, 1.0)]
        plotter.add_text("3D Isometric View", font_size=11, position="upper_left", color="black")

        # Subplot 1: Top View (XY)
        plotter.subplot(0, 1)
        plotter.clear_actors()
        plotter.add_mesh(mesh, scalars=active_field, cmap=cmap, clim=[min_val, max_val],
                         show_edges=False, smooth_shading=True, show_scalar_bar=False)
        plotter.view_xy()
        plotter.camera.zoom(0.95)
        plotter.add_text("Top View (X-Y Seam)", font_size=11, position="upper_left", color="black")

        # Subplot 2: Melt Pool Zoom (Tracking Beam)
        plotter.subplot(0, 2)
        plotter.clear_actors()
        plotter.add_mesh(mesh, scalars=active_field, cmap=cmap, clim=[min_val, max_val],
                         show_edges=False, smooth_shading=True, show_scalar_bar=True,
                         scalar_bar_args={"title": "Temp (°C)" if active_field == "Temperature" else active_field,
                                          "n_labels": 5, "vertical": True, "fmt": "%.0f"})
        cam_offset = np.array([np.cos(theta + 0.3) * 16.0, np.sin(theta + 0.3) * 16.0, 12.0])
        plotter.camera_position = [(peak_x + cam_offset[0], peak_y + cam_offset[1], peak_z + cam_offset[2]),
                                   (peak_x, peak_y, peak_z - 2.0),
                                   (0.0, 0.0, 1.0)]
        plotter.add_text("Melt Pool Zoom (Local HAZ)", font_size=11, position="upper_left", color="white")

        img = plotter.screenshot(return_img=True)
        frame_im = Image.fromarray(img)

        # Composite HUD Header Banner
        banner = Image.new("RGB", (frame_im.width, frame_im.height + banner_h), (245, 247, 250))
        banner.paste(frame_im, (0, banner_h))
        draw = ImageDraw.Draw(banner)
        draw.text((18, 11), "ABAQUS / DFLUX: Circumferential Laser Welding", fill=(15, 23, 42), font=font_bold)
        unit_str = "°C" if active_field == "Temperature" else ""
        hud_info = f"Time: {sim_time:.3f} s  |  Peak: {peak_val:.1f} {unit_str}  |  Frame: {idx + 1} / {len(vtk_files)}"
        draw.text((frame_im.width - 480, 11), hud_info, fill=(30, 41, 59), font=font)
        draw.line([(0, banner_h - 1), (frame_im.width, banner_h - 1)], fill=(226, 232, 240), width=1)

        # Optional downscale for web/README optimization
        if args.scale and args.scale < 1.0:
            new_w = int(banner.width * args.scale)
            new_h = int(banner.height * args.scale)
            banner = banner.resize((new_w, new_h), Image.Resampling.LANCZOS)

        rendered_images.append(banner)

        if (idx + 1) % 10 == 0 or (idx + 1) == len(vtk_files):
            elapsed = time.time() - t0
            fps_render = (idx + 1) / elapsed if elapsed > 0 else 0
            print(f"[INFO] Rendered {idx + 1}/{len(vtk_files)} frames ({elapsed:.1f} s, {fps_render:.1f} fps)...")

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
