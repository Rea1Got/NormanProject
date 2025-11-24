import matplotlib.pyplot as plt 
import numpy as np
import json
from scipy import stats

with open("cfg/cfg.json", "r") as f:
    cfg = json.load(f)
    coord_file = cfg["coord_file"]
    num_molecules = cfg["num_molecules"]
    density = cfg["density"]


coord_data = np.loadtxt(coord_file)
T = len(coord_data)  # количество временных точек
trajectory = coord_data.reshape(T, num_molecules, 3)
total_volume = num_molecules / density
L = total_volume ** (1/3)
L_array = [L, L, L]

def save_to_xyz(coords, output_filename, particle_type="A", box_size=None):
    """
    Сохраняет координаты в формате XYZ для Ovito
    
    Parameters:
    - coords: массив формы (T, num_molecules, 3)
    - output_filename: имя выходного файла
    - particle_type: тип частиц (по умолчанию "A")
    - box_size: размеры ячейки [Lx, Ly, Lz], если None - вычисляется из данных
    """
    T, num_molecules, _ = coords.shape
    
    # Если размер ячейки не задан, вычисляем из максимальных координат
    if box_size is None:
        box_size = np.max(coords, axis=(0, 1)) - np.min(coords, axis=(0, 1))
    
    with open(output_filename, 'w') as f:
        for t in range(T):
            # Заголовок для каждого кадра
            f.write(f"{num_molecules}\n")
            f.write(f"Lattice=\"{box_size[0]} 0.0 0.0 0.0 {box_size[1]} 0.0 0.0 0.0 {box_size[2]}\" ")
            f.write(f"Properties=species:S:1:pos:R:3 Time={t}\n")
            
            # Координаты частиц
            for i in range(num_molecules):
                x, y, z = coords[t, i]
                f.write(f"{particle_type} {x:.6f} {y:.6f} {z:.6f}\n")

save_to_xyz(trajectory, "tmp/200_300K/trajectory.xyz", box_size=L_array)

