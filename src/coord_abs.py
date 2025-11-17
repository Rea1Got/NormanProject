import matplotlib.pyplot as plt 
import numpy as np
import json
from scipy import stats

with open("cfg/cfg.json", "r") as f:
    cfg = json.load(f)
    coord_abs_file = cfg["coord_abs_file"]
    num_molecules = cfg["num_molecules"]
    total_steps = cfg["total_steps"]
    dt = cfg["dt"]  # шаг по времени в безразмерных единицах

# Загрузка и подготовка данных
coord_data = np.loadtxt(coord_abs_file)
T = len(coord_data)  # количество временных точек
print(f"Number of time points: {T}")
print(f"Number of molecules: {num_molecules}")

# Преобразование в траектории: [время, молекула, координаты]
trajectory = coord_data.reshape(T, num_molecules, 3)

# Параметры анализа
window_size = int(input("Window size: "))  
step_size = int(input("Step size: "))
cut_index = int(input("Cut first points: "))

# Проверка корректности параметров
if window_size > T:
    raise ValueError(f"Window size ({window_size}) cannot be larger than total steps ({T})")
if cut_index >= window_size:
    raise ValueError(f"Cut index ({cut_index}) must be less than window size ({window_size})")

print(f"\nAnalysis parameters:")
print(f"Window size: {window_size} points")
print(f"Window step: {step_size} points") 
print(f"Cut first {cut_index} points")

def compute_msd_optimized(trajectory, window_size, step_size):
    """
    Оптимизированное вычисление MSD с скользящим окном
    """
    T, N, dim = trajectory.shape
    num_windows = (T - window_size) // step_size + 1
    
    msd_results = np.zeros((num_windows, window_size))
    
    for window_idx in range(num_windows):
        start = window_idx * step_size
        end = start + window_size
        window_traj = trajectory[start:end]  # [window_size, N, 3]
        
        # Вычисляем MSD для всех временных лагов в окне
        for tau in range(window_size):
            if tau == 0:
                msd_results[window_idx, tau] = 0.0
            else:
                # Смещение за время tau
                displacements = window_traj[tau] - window_traj[0]
                # MSD = средний квадрат смещения
                msd_results[window_idx, tau] = np.mean(displacements**2)
        
        if (window_idx + 1) % max(1, num_windows // 10) == 0:
            print(f"Processed {window_idx + 1}/{num_windows} windows")
    
    return msd_results

# Вычисляем MSD
print("\nComputing MSD...")
msd_results = compute_msd_optimized(trajectory, window_size, step_size)

# Усредняем по всем окнам
average_msd = np.mean(msd_results, axis=0)

# Обрезаем начальные точки и создаем массив времен
tau = np.arange(window_size) * dt  # переводим шаги в безразмерное время
average_msd_cut = average_msd[cut_index:]
tau_cut = tau[cut_index:]

print(f"\nMSD results shape: {msd_results.shape}")
print(f"Final MSD array size: {len(average_msd_cut)}")

# Линейная регрессия для определения коэффициента диффузии
slope, intercept, r_value, p_value, std_err = stats.linregress(tau_cut, average_msd_cut)
D = slope / 6.0  # для 3D: MSD = 6Dt

print("\n" + "="*50)
print("RESULTS:")
print("="*50)
print(f"Diffusion coefficient D = {D:.6e}")
print(f"Slope = {slope:.6f}")
print(f"R² = {r_value**2:.6f}")
print(f"Standard error = {std_err:.6e}")

# Визуализация
plt.figure(figsize=(10, 6))

# Экспериментальные данные
plt.plot(tau_cut, average_msd_cut, 'bo', alpha=0.6, markersize=3, label='MSD data')

# Линейная аппроксимация
fit_line = intercept + slope * tau_cut
plt.plot(tau_cut, fit_line, 'r-', linewidth=2, 
         label=f'Linear fit: MSD = {slope:.4f}t + {intercept:.4f}')

plt.xlabel('Time (reduced units)')
plt.ylabel('MSD (reduced units)')
plt.title(f'Mean Squared Displacement\nD = {D:.4e}, R² = {r_value**2:.4f}')
plt.legend()
plt.grid(True, alpha=0.3)

# Добавляем текст с результатами в углу графика
textstr = f'D = {D:.4e}\nR² = {r_value**2:.4f}\nN = {num_molecules}'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=12,
         verticalalignment='top', bbox=props)

plt.tight_layout()
plt.savefig('msd_analysis.png', dpi=300)
plt.show()

# Дополнительная информация
print("\nAdditional information:")
print(f"Time range analyzed: {tau_cut[0]:.3f} to {tau_cut[-1]:.3f} (reduced units)")
print(f"MSD range: {average_msd_cut[0]:.3f} to {average_msd_cut[-1]:.3f} (reduced units)")
