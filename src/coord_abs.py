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

print("Enter window size and step size:")
window_size = int(input("Window size: "))  
step_size = int(input("Step size: "))
cut_index = int(input("Cut first points: "))

print(f"Analysis of trajectory with {T} time points and {num_molecules} molecules.")
print(f"Window size: {window_size} points ({window_size * 100 / T}% out of points)")
print(f"Window step size: {step_size} points ({step_size * 100 / window_size}% out of window size)")

############
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
            # displacements = window_trajectory[tau:] - window_trajectory[:-tau]
            # squared_distances = np.sum(displacements**2, axis=2)
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

############
msd_results = np.array(msd_results)
print(f"MSD results shape: {msd_results.shape}")  # (number_of_window, window_size)

average_msd_per_tau = np.mean(msd_results, axis=0)

# cut_index = 200 
average_msd_per_tau = average_msd_per_tau[cut_index:]
tau = np.arange(cut_index, window_size)

print(f"Average MSD per tau shape: {average_msd_per_tau.shape}")

print("\nMSD for first 10 time lags after cut (averaged over all windows):")
for i, tau_val in enumerate(tau[:10]):
    print(f"tau = {tau_val}: MSD = {average_msd_per_tau[i]:.6f}")

avg_D = average_msd_per_tau / (6 * tau)
std_D = avg_D.std()
calculated_D = avg_D.mean()

print(f"Средний коэффициент диффузии: {calculated_D:.4e} ± {std_D:.4e}")
print(f"Относительная погрешность: {std_D/calculated_D*100:.5f}%")


############ 
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
weights = [x/sum(avg_D) for x in avg_D]
ax1.hist(avg_D, bins=20, alpha=0.7, edgecolor='black', weights = weights)
ax1.axvline(calculated_D, color='r', linestyle='--', label=f'Среднее: {calculated_D:.4e}')
ax1.set_xlabel('Коэффициент диффузии')
ax1.set_ylabel('Частота')
ax1.set_title('Распределение коэффициентов диффузии по окнам')
ax1.legend()
ax1.grid(True, alpha=0.3)
########
theoretical_msd = [6 * calculated_D * t for t in tau]

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
