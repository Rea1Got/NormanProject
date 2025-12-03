import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist

with open("cfg/cfg.json", "r") as f:
    cfg = json.load(f)
    coord_file = cfg["coord_file"]
    velocity_file = cfg["velocity_file"]
    num_molecules = cfg["num_molecules"]
    total_steps = cfg["total_steps"]
    dt = cfg["dt"]  
    density = cfg["density"]

velocity_data_loaded = np.loadtxt(velocity_file)
velocity_data = velocity_data_loaded.reshape(-1, num_molecules, 3)

N = num_molecules
L = (N / density)**(1/3) 

data = np.loadtxt(coord_file)  # shape: (n_frames, 3*N)
n_frames = data.shape[0]  

def compute_g_r(coords_array, N, L, r_max=None, n_bins=250, sample_every=150):
    if r_max is None:
        r_max = L / 2      

    dr = r_max / n_bins
    bins = np.linspace(0, r_max, n_bins+1)
    r_centers = (bins[1:] + bins[:-1]) / 2
    histogram = np.zeros(n_bins)
    
    n_frames_total = coords_array.shape[0]
    
    # Используем только каждый sample_every-ный кадр
    frame_indices = range(0, n_frames_total, sample_every)
    n_frames_used = len(frame_indices)
    
    print(f"Используется {n_frames_used} из {n_frames_total} кадров")
    print(f"r_max = {r_max}, dr = {dr}")
    
    # Вычисляем все попарные индексы (треугольная матрица)
    i_indices, j_indices = np.triu_indices(N, k=1)
    n_pairs = len(i_indices)
    print(f"Количество уникальных пар частиц: {n_pairs}")
    
    for counter, frame_idx in enumerate(frame_indices):
        # Получаем координаты для этого кадра
        coords = coords_array[frame_idx].reshape(N, 3)
        
        # Векторно вычисляем все расстояния сразу
        r_vectors = coords[i_indices] - coords[j_indices]
        r_vectors = r_vectors - L * np.round(r_vectors / L)
        distances = np.sqrt(np.sum(r_vectors**2, axis=1))
        
        # Фильтруем по r_max и добавляем в гистограмму
        mask = distances < r_max
        valid_distances = distances[mask]
        
        if len(valid_distances) > 0:
            bin_indices = (valid_distances / dr).astype(int)
            # Убедимся, что индексы в пределах массива
            bin_indices = np.minimum(bin_indices, n_bins - 1)
            
            # Каждая пара дает вклад 1 (не 2!)
            np.add.at(histogram, bin_indices, 1)
        
        if (counter + 1) % 10 == 0:
            print(f"Обработано кадров: {counter + 1}/{n_frames_used}")
    
    density = N / (L**3)
    print(f"Плотность системы ρ = {density}")
    
    d_volumes = 4.0 * np.pi * r_centers**2 * dr
    
    mol_in_d_vol = density * d_volumes
    
    # 4. Общее число пар, которые мы могли бы наблюдать
    #    за все кадры = n_frames_used * N * идеальное_число_в_слое
    total_pairs_possible = n_frames_used * N * mol_in_d_vol
    
    # 5. g(r) = (наблюдаемое_число_пар) / (общее_число_пар_в_идеальном_газе)
    g_r = (2.0 * histogram) / total_pairs_possible
    
    return r_centers, g_r, histogram

# Расчет g(r)
r, g_r, hist = compute_g_r(data, N, L, r_max=L/2, n_bins=400, sample_every=50)

def calculate_velocity_squared(velocity_data):
    velocity_squared = np.sum(velocity_data**2, axis=2)  # форма: (T_v, num_molecules)
    mean_squared_velocity = np.mean(velocity_squared)
    
    return mean_squared_velocity 

rms_vel = calculate_velocity_squared(velocity_data)

plt.figure(figsize=(10, 6))
plt.plot(r, g_r, 'b-', linewidth=2, label=f'g(r), pho={density}, T={rms_vel/3:.4f}')
plt.axhline(y=1, color='r', linestyle='--', alpha=0.5, label='Идеальный газ (g(r)=1)')
plt.xlabel('Расстояние r', fontsize=12)
plt.ylabel('g(r)', fontsize=12)
plt.title('Радиальная функция распределения', fontsize=14)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

# # Дополнительно: интеграл для получения координационного числа
# dr = r[1] - r[0]
# coord_number = 4 * np.pi * density * np.cumsum(r**2 * g_r) * dr
#
# plt.figure(figsize=(10, 6))
# plt.plot(r, coord_number, 'g-', linewidth=2)
# plt.xlabel('Расстояние r', fontsize=12)
# plt.ylabel('n(r) (координационное число)', fontsize=12)
# plt.title('Интегральная координационная функция', fontsize=14)
# plt.grid(True, alpha=0.3)
# plt.tight_layout()
#
plt.show()
