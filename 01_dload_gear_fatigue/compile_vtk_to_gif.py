"""
compile_vtk_to_gif.py
=====================
Compilador Multi-Viewport 3DEXPERIENCE SIMULIA:
  1. Painel 1 (Esquerda) : Vista Isometrica (3D Global - Padrao Europeu 1o Diedro)
                           Face Y a Direita, Lateral a Esquerda, Topo em Cima.
                           Escala visual balanceada com a Vista Frontal.
                           Legenda compacta com wireframe no Canto Inferior Esquerdo (livre de obstrucao).
  2. Painel 2 (Centro)   : Vista Frontal (Eixo Y) - Cara da coroa no plano X-Z (Escala balanceada).
  3. Painel 3 (Direita)  : Corte RZ 3D com Profundidade (Fatia sem malha interna + Tensao Continua/Smeared).

Gera um GIF animado 100% estavel, com paleta global indexada (cores exatas com hashcode fixo do inicio ao fim).
"""

import os
import sys
import glob
import time
import shutil
import zipfile
import tempfile
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
import matplotlib.colors as mcolors

# -------------------------------------------------------------
# 0. Localizacao de Ficheiros e Auto-Extracao do ZIP
# -------------------------------------------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))

zip_candidates = [
    os.path.join(script_dir, "3DX", "gear_helical_vtk.zip"),
    os.path.join(script_dir, "gear_helical_vtk.zip"),
    os.path.join(script_dir, "01_dload_gear_fatigue", "gear_helical_vtk.zip")
]

vtk_zip_path = None
for z in zip_candidates:
    if os.path.exists(z):
        vtk_zip_path = z
        break

temp_extract_dir = None
if vtk_zip_path:
    print(f"[INFO] Arquivo ZIP encontrado: {vtk_zip_path}")
    temp_extract_dir = tempfile.mkdtemp(prefix="vtk_extract_")
    print(f"[INFO] A extrair ficheiros para directoria temporaria: {temp_extract_dir}...")
    with zipfile.ZipFile(vtk_zip_path, 'r') as zf:
        zf.extractall(temp_extract_dir)
    vtk_files = sorted(glob.glob(os.path.join(temp_extract_dir, "*.vtk")))
else:
    vtk_files = sorted(glob.glob(os.path.join(script_dir, "*.vtk")))
    if not vtk_files:
        vtk_files = sorted(glob.glob(os.path.join(script_dir, "3DX", "*.vtk")))

if not vtk_files:
    print("[ERRO] Nenhum ficheiro VTK encontrado para processar.")
    sys.exit(1)

total_frames = len(vtk_files)
print(f"[INFO] Total de frames VTK a processar: {total_frames}")

# -------------------------------------------------------------
# 1. Analise Cinematica e Extracao de Tempos
# -------------------------------------------------------------
times = []
for fpath in vtk_files:
    t_val = 0.0
    with open(fpath, "r", encoding="utf-8", errors="ignore") as f:
        for _ in range(5):
            line = f.readline()
            if "Time" in line:
                parts = line.strip().split()
                try:
                    t_val = float(parts[parts.index("Time") + 1])
                except (ValueError, IndexError):
                    pass
                break
    times.append(t_val)

raw_thetas = []
for fpath in vtk_files:
    mesh = pv.read(fpath)
    vm = mesh.cell_data.get("Von_Mises", None)
    if vm is None:
        vm = mesh.point_data.get("Von_Mises", None)
    val_max = float(np.nanmax(vm))
    if val_max > 1.0:
        c_idx = int(np.nanargmax(vm))
        cp = mesh.cell_centers().points[c_idx]
        th = np.arctan2(cp[0], cp[2])
        raw_thetas.append(th)
    else:
        raw_thetas.append(None)

first_valid = next(th for th in raw_thetas if th is not None)
clean_thetas = [th if th is not None else first_valid for th in raw_thetas]
unwrapped = np.unwrap(clean_thetas)

t_ramp = 0.20
ramp_indices = [i for i, t in enumerate(times) if t <= t_ramp + 1e-4]
ramp_end_idx = ramp_indices[-1] if ramp_indices else 0

th_start = unwrapped[ramp_end_idx]
th_end = unwrapped[-1]
delta_total = th_end - th_start
t_start = times[ramp_end_idx]
t_end = times[-1]
t_roll = max(t_end - t_start, 1e-6)

print(f"[CINEMATICA] t_ramp = {t_start:.3f} s (frame {ramp_end_idx}) -> theta = {np.rad2deg(th_start):.2f} deg")
print(f"[CINEMATICA] t_end  = {t_end:.3f} s (frame {total_frames - 1}) -> theta = {np.rad2deg(th_end):.2f} deg")
print(f"[CINEMATICA] Rotacao total: {np.rad2deg(delta_total):.2f} deg ao longo de {t_roll:.3f} s")

smooth_thetas = []
for t in times:
    if t <= t_start:
        th = th_start
    else:
        frac = (t - t_start) / t_roll
        th = th_start + frac * delta_total
    smooth_thetas.append(th)

# -------------------------------------------------------------
# 2. Paleta 3DEXPERIENCE com Hashcodes Exatos e Imutaveis
# -------------------------------------------------------------
bg_3dx = "#e0e2e7"

# 12 cores exatas das bandas (do azul 0 MPa ao vermelho 100+ MPa)
colors_3dx = [
    "#0000ff",  # 0-10 MPa (Azul profundo)
    "#005cff",  # 10-20
    "#00b9ff",  # 20-30
    "#00ffe7",  # 30-40 (Ciano)
    "#00ff8b",  # 40-50
    "#00ff2e",  # 50-60 (Verde)
    "#2eff00",  # 60-70
    "#8bff00",  # 70-80 (Lima)
    "#e7ff00",  # 80-85 (Amarelo)
    "#ffb900",  # 85-90 (Laranja claro)
    "#ff5c00",  # 90-95 (Laranja escuro)
    "#ff0000"   # 95-100+ MPa (Vermelho)
]
cmap_3dx = mcolors.LinearSegmentedColormap.from_list("3dx_mises", colors_3dx, N=12)

# Lista ordenada de cores para as amostras da legenda (Topo: Vermelho -> Fundo: Azul)
swatch_hex_colors = list(reversed(colors_3dx))

# Cores fixas reservadas na paleta global indexada
fixed_palette_hex = [
    "#ff0000", "#ff5c00", "#ffb900", "#e7ff00", "#8bff00", "#2eff00",
    "#00ff2e", "#00ff8b", "#00ffe7", "#00b9ff", "#005cff", "#0000ff",
    "#969696", "#000000", "#ffffff", "#e0e2e7", "#444444"
]

def hex_to_rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))

# Gerador da Legenda Compacta 3DEXPERIENCE (Zero Flicker + Wireframe)
def create_compact_legend_widget():
    labels = ['>100', '90', '80', '70', '60', '50', '40', '30', '20', '10', '<0']
    
    w_card, h_card = 120, 260
    card = Image.new('RGB', (w_card, h_card), (255, 255, 255))
    draw = ImageDraw.Draw(card)
    
    # Wireframe exterior do cartao
    draw.rectangle([0, 0, w_card - 1, h_card - 1], outline='black', width=1)
    
    try:
        font_title = ImageFont.truetype('arial.ttf', 10)
        font_sub = ImageFont.truetype('arial.ttf', 9)
        font_lbl = ImageFont.truetype('arial.ttf', 9)
    except:
        font_title = font_sub = font_lbl = ImageFont.load_default()
        
    draw.text((8, 6), 'Stress components', fill='black', font=font_title)
    draw.text((8, 18), 'Von Mises (MPa)', fill='black', font=font_sub)
    draw.text((8, 29), 'Centroid', fill='#444444', font=font_sub)
    
    bx, by = 10, 44
    bw, bh = 20, 13
    
    # Amostras de cor com wireframe individual por caixa (Hashcodes exatos)
    for i, col in enumerate(swatch_hex_colors):
        y0 = by + i * bh
        y1 = y0 + bh
        draw.rectangle([bx, y0, bx + bw, y1], fill=col, outline='black', width=1)
        if i < len(labels):
            draw.text((bx + bw + 5, y1 - 5), labels[i], fill='black', font=font_lbl)
            
    # Caixa "No Value"
    y_nv = by + len(swatch_hex_colors) * bh + 8
    draw.rectangle([bx, y_nv, bx + bw, y_nv + bh], fill='#969696', outline='black', width=1)
    draw.text((bx + bw + 5, y_nv + 2), 'No Value', fill='black', font=font_lbl)
    
    return card

legend_card = create_compact_legend_widget()

try:
    font_header = ImageFont.truetype('arial.ttf', 13)
except:
    font_header = ImageFont.load_default()

# -------------------------------------------------------------
# 3. Renderizacao Sequencial dos Frames (3 Vistas Sincronizadas)
# -------------------------------------------------------------
pil_frames = []
t_proc_start = time.time()

for idx, fpath in enumerate(vtk_files):
    mesh = pv.read(fpath)
    th_frame = smooth_thetas[idx]
    th_deg = np.rad2deg(th_frame)

    # Feature edges exteriores para o modelo global 3D
    edges_global = mesh.extract_feature_edges(
        boundary_edges=True,
        feature_edges=True,
        manifold_edges=False,
        feature_angle=25.0
    )

    # Versao continua suave (nodal averaged / smeared) para a seccao em corte
    mesh_smooth = mesh.cell_data_to_point_data()

    # Transformacao para o corte 3D: rotacionar por -th_deg em torno de Y
    mesh_aligned_sm = mesh_smooth.rotate_y(-th_deg, inplace=False)

    # Corpo solido 3D retido por tras do plano de corte (X >= 0 e raio Z >= 68 mm)
    half_solid = mesh_aligned_sm.clip(normal=[1, 0, 0], origin=[0, 0, 0], invert=False)
    tooth_solid = half_solid.clip(normal=[0, 0, 1], origin=[0, 0, 68.0], invert=False)
    edges_tooth = tooth_solid.extract_feature_edges(
        boundary_edges=True,
        feature_edges=True,
        manifold_edges=False,
        feature_angle=25.0
    )

    # Fatia transversal no plano de corte X = 0 com tensao continua (SEM malha interna)
    cut_slice = mesh_aligned_sm.slice(normal=[1, 0, 0], origin=[0, 0, 0])
    cut_slice_tooth = cut_slice.clip(normal=[0, 0, 1], origin=[0, 0, 68.0], invert=False)
    cut_perim = cut_slice_tooth.extract_feature_edges(
        boundary_edges=True,
        feature_edges=False,
        manifold_edges=False
    )

    plotter = pv.Plotter(shape=(1, 3), window_size=[1920, 640], off_screen=True)
    plotter.set_background(bg_3dx)

    # ---------------------------------------------------------
    # VISTA 1: ISOMETRIC VIEW (3D Global - Padrao Europeu / 1o Diedro)
    # Face Y a Direita, Lateral a Esquerda, Topo em Cima.
    # Escala visual calibrada para igualar a Vista Frontal (~515-540 px)
    # ---------------------------------------------------------
    plotter.subplot(0, 0)
    plotter.add_mesh(
        mesh,
        scalars="Von_Mises",
        cmap=cmap_3dx,
        clim=[0.0, 100.0],
        n_colors=12,
        show_edges=False,
        show_scalar_bar=False
    )
    plotter.add_mesh(edges_global, color="black", line_width=1.0)
    plotter.camera.focal_point = ( 25.0, -15.0, -5.0)
    plotter.camera.position =    (240.0, 330.0, 170.0)
    plotter.camera.up = (0.0, 0.0, 1.0)
    plotter.camera.zoom(0.96)

    # ---------------------------------------------------------
    # VISTA 2: FRONTAL VIEW (Y-Axis)
    # Escala visual calibrada para igualar a Vista 1 (~515 px)
    # ---------------------------------------------------------
    plotter.subplot(0, 1)
    plotter.add_mesh(
        mesh,
        scalars="Von_Mises",
        cmap=cmap_3dx,
        clim=[0.0, 100.0],
        n_colors=12,
        show_edges=False,
        show_scalar_bar=False
    )
    plotter.add_mesh(edges_global, color="black", line_width=1.0)
    plotter.camera.focal_point = (0.0, -15.0, 0.0)
    plotter.camera.position = (0.0, 475.0, 0.0)
    plotter.camera.up = (0.0, 0.0, 1.0)

    # ---------------------------------------------------------
    # VISTA 3: 3D RZ SECTION (Cut & Depth - Smeared Continuous Stress)
    # ---------------------------------------------------------
    plotter.subplot(0, 2)
    plotter.add_mesh(
        tooth_solid,
        scalars="Von_Mises",
        cmap=cmap_3dx,
        clim=[0.0, 100.0],
        n_colors=12,
        show_edges=False,
        show_scalar_bar=False
    )
    plotter.add_mesh(edges_tooth, color="black", line_width=1.0)

    # Fatia no plano de corte frontal (Smeared Continuous, show_edges=False)
    plotter.add_mesh(
        cut_slice_tooth,
        scalars="Von_Mises",
        cmap=cmap_3dx,
        clim=[0.0, 100.0],
        n_colors=12,
        show_edges=False,
        show_scalar_bar=False
    )
    # Delimitacao do perimetro exterior da fatia
    plotter.add_mesh(cut_perim, color="black", line_width=1.0)
    
    # Camara 3D estatica focada no corte e profundidade
    plotter.camera.focal_point = (10.0, -15.0, 84.0)
    plotter.camera.position = (-85.0, -45.0, 105.0)
    plotter.camera.up = (0.0, 0.0, 1.0)

    raw_img = plotter.screenshot(return_img=True)
    plotter.close()

    # Pos-processamento PIL: Sobreposicao dos Titulos e Legenda Estatica
    pil_img = Image.fromarray(raw_img)
    draw = ImageDraw.Draw(pil_img)

    # Titulos dos paineis em Ingles
    draw.text((20, 15), "1. Isometric View (3D Global)", fill="black", font=font_header)
    draw.text((660, 15), "2. Frontal View (Y-Axis)", fill="black", font=font_header)
    draw.text((1300, 15), "3. 3D RZ Section (Cut & Depth)", fill="black", font=font_header)

    # Legenda 3DEXPERIENCE Compacta no Canto Inferior Esquerdo (livre de obstrucao)
    pil_img.paste(legend_card, (15, 365))

    pil_frames.append(pil_img)

    del mesh, mesh_smooth, mesh_aligned_sm, half_solid, tooth_solid, edges_tooth, cut_slice, cut_slice_tooth, cut_perim

    if (idx + 1) % 10 == 0 or (idx + 1) == total_frames:
        elapsed = time.time() - t_proc_start
        fps_rate = (idx + 1) / elapsed
        print(f"  Frame {idx + 1}/{total_frames} processado ({fps_rate:.1f} fps)...")

# -------------------------------------------------------------
# 4. Construcao da Paleta Master Imutavel & Gravacao do GIF
# -------------------------------------------------------------
fps = 5 if total_frames <= 35 else 8  # Half presentation speed (~120ms per frame)
output_gif = "gear_torque_simulation.gif"
print(f"[INFO] A construir paleta master global com hashcodes exatos...")

# Extrair RGBs fixos
fixed_rgbs = [hex_to_rgb(c) for c in fixed_palette_hex]

# Construir imagem de amostragem combinada para calibrar cores da malha
sample_img = Image.new('RGB', (1920, 640), (224, 226, 231))
sample_step = max(len(pil_frames) // 4, 1)
for i, f_idx in enumerate(range(0, len(pil_frames), sample_step)):
    # Amostrar fatias de varios frames
    box = (i * 480, 0, (i + 1) * 480, 640)
    crop_f = pil_frames[f_idx].crop(box)
    sample_img.paste(crop_f, (i * 480, 0))

# Quantizar a amostra para 256 cores base
temp_p = sample_img.quantize(colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
raw_palette = list(temp_p.getpalette())

# Bloquear as primeiras posicoes da paleta com os hashcodes exatos e imutaveis
for i, rgb in enumerate(fixed_rgbs):
    raw_palette[i * 3 : i * 3 + 3] = list(rgb)

master_palette_image = Image.new('P', (1, 1))
master_palette_image.putpalette(raw_palette)

print(f"[INFO] A quantizar todos os {len(pil_frames)} frames contra a paleta master...")
quantized_frames = [
    f.quantize(palette=master_palette_image, dither=Image.Dither.NONE)
    for f in pil_frames
]

print(f"[INFO] A gravar GIF animado perfeitamente estavel em '{output_gif}' a {fps} FPS...")
quantized_frames[0].save(
    output_gif,
    save_all=True,
    append_images=quantized_frames[1:],
    duration=int(1000.0 / fps),
    loop=0,
    optimize=True
)

# Sincronizar tambem para a pasta 3DX
dest_3dx = os.path.join(script_dir, "3DX", output_gif)
if os.path.exists(os.path.join(script_dir, "3DX")):
    shutil.copyfile(output_gif, dest_3dx)
    print(f"[INFO] Copia sincronizada em '{dest_3dx}'.")

# Sincronizar para a pasta 01_dload_gear_fatigue
dest_portfolio = os.path.join(script_dir, "01_dload_gear_fatigue", output_gif)
if os.path.exists(os.path.join(script_dir, "01_dload_gear_fatigue")):
    shutil.copyfile(output_gif, dest_portfolio)
    print(f"[INFO] Copia sincronizada em '{dest_portfolio}'.")

# -------------------------------------------------------------
# 5. Limpeza de ficheiros temporarios
# -------------------------------------------------------------
if temp_extract_dir and os.path.exists(temp_extract_dir):
    shutil.rmtree(temp_extract_dir, ignore_errors=True)
    print("[INFO] Ficheiros VTK temporarios eliminados.")

print(f"\n[SUCESSO] Animacao sincronizada gerada com sucesso: {output_gif}")