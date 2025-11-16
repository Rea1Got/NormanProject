import matplotlib.pyplot as plt 
import numpy as np
import json
from scipy.optimize import curve_fit
from scipy import stats

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

print("Enter window size and step size:")
window_size = int(input("Window size: "))  
step_size = int(input("Step size: "))
cut_index = int(input("Cut first points: "))

print(f"Analysis of trajectory with {T} time points and {num_molecules} molecules.")
print(f"Window size: {window_size} points ({window_size * 100 / T}% out of points)")
print(f"Window step size: {step_size} points ({step_size * 100 / window_size}% out of window size)")

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
            displacements = window_trajectory[tau, :, :] - window_trajectory[0, :, :]
            msd[tau] = np.mean(displacements**2)

    return msd

msd_results = []
for start_idx in range(0, T - window_size + 1, step_size):
    end_idx = start_idx + window_size
    
    window_traj = trajectory[start_idx:end_idx]
    window_msd = compute_msd_in_window(window_traj)
    msd_results.append(window_msd)
    
    if (start_idx // (T - window_size + 1)) % 100 == 0:
        print(f"Completed: {start_idx/(T - window_size + 1) * 100:.1f}%")

msd_results = np.array(msd_results)
print(f"MSD results shape: {msd_results.shape}")

average_msd_per_tau = np.mean(msd_results, axis=0)

# Обрезаем начальные точки
average_msd_per_tau = average_msd_per_tau[cut_index:]
tau = np.arange(cut_index, window_size)

print(f"Average MSD per tau shape: {average_msd_per_tau.shape}")

print("\nMSD for first 10 time lags after cut (averaged over all windows):")
for i, tau_val in enumerate(tau[:10]):
    print(f"tau = {tau_val}: MSD = {average_msd_per_tau[i]:.6f}")

# КОРРЕКТНЫЙ РАСЧЕТ КОЭФФИЦИЕНТА ДИФФУЗИИ
# Линейная регрессия для MSD(tau)
slope, intercept, r_value, p_value, std_err = stats.linregress(tau, average_msd_per_tau)
D_from_slope = slope / 6  # Коэффициент диффузии из наклона (для 3D: MSD = 6Dt)

# Теоретическая MSD на основе линейной регрессии
theoretical_msd = intercept + slope * np.array(tau)

# Старый метод для сравнения (менее точный)
avg_D_old = average_msd_per_tau / (6 * tau)
std_D_old = avg_D_old.std()
calculated_D_old = avg_D_old.mean()

print("\n=== РЕЗУЛЬТАТЫ ===")
print(f"Коэффициент диффузии из усреднения (старый метод): {calculated_D_old:.4e} ± {std_D_old:.4e}")
print(f"Коэффициент диффузии из наклона MSD (корректный): {D_from_slope:.4e}")
print(f"Коэффициент детерминации R²: {r_value**2:.4f}")
print(f"Относительная погрешность (старый метод): {std_D_old/calculated_D_old*100:.5f}%")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# График 1: Распределение коэффициентов диффузии
weights = [x/sum(avg_D_old) for x in avg_D_old]
ax1.hist(avg_D_old, bins=20, alpha=0.7, edgecolor='black', weights=weights)
ax1.axvline(calculated_D_old, color='r', linestyle='--', label=f'Среднее: {calculated_D_old:.4e}')
ax1.axvline(D_from_slope, color='g', linestyle='-', label=f'Из наклона: {D_from_slope:.4e}')
ax1.set_xlabel('Коэффициент диффузии')
ax1.set_ylabel('Частота')
ax1.set_title('Распределение коэффициентов диффузии (старый метод)')
ax1.legend()
ax1.grid(True, alpha=0.3)

# График 2: MSD vs tau с правильным расчетом наклона
ax2.plot(tau, average_msd_per_tau, 'o', alpha=0.7, label='Экспериментальные данные')
ax2.plot(tau, theoretical_msd, 'r-', linewidth=2, 
         label=f'Линейная аппроксимация (D={D_from_slope:.2e})')
ax2.set_xlabel('Количество шагов (τ)')
ax2.set_ylabel('MSD')
ax2.set_title('Зависимость MSD от времени')
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.show()

print(f"\n=== ДЕТАЛИ АППРОКСИМАЦИИ ===")
print(f"Наклон (slope): {slope:.6f}")
print(f"Пересечение (intercept): {intercept:.6f}")
print(f"Стандартная ошибка наклона: {std_err:.6f}")
print(f"p-value: {p_value:.6f}")
