import matplotlib.pyplot as plt 
import numpy as np
import json
from scipy.optimize import curve_fit

with open("cfg/cfg.json", "r") as f:
    cfg = json.load(f)
    coord_abs_file = cfg["coord_abs_file"]
    num_molecules = cfg["num_molecules"]
    total_steps = cfg["total_steps"]
    snapshot = cfg["snapshot"]
    dt = cfg["dt"]

coord_data = np.loadtxt(coord_abs_file)
T = len(coord_data)  # num of time points
print(f"Number of points: {T}")
trajectory = coord_data.reshape(T, num_molecules, 3)

# Функция для линейной аппроксимации
def linear_func(t, a, b):
    return a * t + b

print("Enter window size and step size:")
window_size = int(input("Window size: "))  
step_size = int(input("Step size: ")) 

print(f"Analysis of trajectory with {T} time points and {num_molecules} molecules.")
print(f"Window size: {window_size} points ({window_size * 100 / T}% out of points)")
print(f"Window step size: {step_size} points ({step_size * 100 / window_size}% out of window size)")



# Анали здля каждого скользящего окна
def compute_msd_in_window(window_trajectory):
    """
    Вычисляет MSD для временного окна траектории
    window_trajectory: массив формы (window_size, num_molecules, 3)
    Возвращает: массив MSD для каждого временного лага в окне
    """
    window_length = window_trajectory.shape[0]
    msd = np.zeros(window_length)
    
    for tau in range(window_length):
        if tau == 0:
            msd[tau] = 0.0  
        else:
            displacements = window_trajectory[tau:] - window_trajectory[:-tau]
            squared_distances = np.sum(displacements**2, axis=2)
            msd[tau] = np.mean(squared_distances)
    
    return msd

msd_results = []
for start_idx in range(0, T - window_size + 1, step_size):
    end_idx = start_idx + window_size
    
    window_traj = trajectory[start_idx:end_idx]
    window_msd = compute_msd_in_window(window_traj)
    msd_results.append(window_msd)
    
    if (start_idx // step_size) % 10 == 0:
        print(f"Processed window starting at {start_idx}")

msd_results = np.array(msd_results)
print(f"MSD results shape: {msd_results.shape}")  # (number_of_window, window_size)

average_msd_per_tau = np.mean(msd_results, axis=0)
print(f"Average MSD per tau shape: {average_msd_per_tau.shape}")

print("\nMSD for first 10 time lags (averaged over all windows):")
for tau in range(min(10, window_size)):
    print(f"tau = {tau}: MSD = {average_msd_per_tau[tau]:.6f}")

tau = [i for i in range(window_size)]
avg_D = [average_msd_per_tau[i]/(6*i) for i in range(1, window_size)]
avg_D = np.array(avg_D)
std_D = avg_D.std()
calculated_D = avg_D.mean()
print(f"Средний коэффициент диффузии: {calculated_D:.4e} ± {std_D:.4e}")
print(f"Относительная погрешность: {std_D/calculated_D*100:.5}%") 


# Визуализация результатов

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

ax1.hist(avg_D, bins=20, alpha=0.7, edgecolor='black')
ax1.axvline(calculated_D, color='r', linestyle='--', label=f'Среднее: {calculated_D:.4e}')
ax1.set_xlabel('Коэффициент диффузии')
ax1.set_ylabel('Частота')
ax1.set_title('Распределение коэффициентов диффузии по окнам')
ax1.legend()
ax1.grid(True, alpha=0.3)
############
theoretical_msd = [6 * calculated_D * i for i in range(window_size)] 

ax2.plot(tau, average_msd_per_tau, 'o', alpha=0.7)
ax2.plot(tau, theoretical_msd, 'r--', linewidth=2, 
            label=f'Рассчитанная (D={calculated_D:.2e})')
ax2.set_xlabel('Количество шагов')
ax2.set_ylabel('MSD')
ax2.set_title('Зависимость MSD от времени в пределах окна')
ax2.grid(True, alpha=0.3)
ax2.legend()
plt.tight_layout()
plt.show()
