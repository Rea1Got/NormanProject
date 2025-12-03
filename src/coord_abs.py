import matplotlib.pyplot as plt 
import numpy as np
import json
from scipy import stats

with open("cfg/cfg.json", "r") as f:
    cfg = json.load(f)
    coord_abs_file = cfg["coord_abs_file"]
    velocity_data = cfg["velocity_file"]
    num_molecules = cfg["num_molecules"]
    total_steps = cfg["total_steps"]
    dt = cfg["dt"]  # шаг по времени в безразмерных единицах
    density = cfg["density"]

coord_data = np.loadtxt(coord_abs_file)
velocity_data = np.loadtxt(velocity_data)
T_v = len(velocity_data)
T = len(coord_data)  # количество временных точек
print(f"Number of time points: {T}")
print(f"Number of molecules: {num_molecules}")

trajectory = coord_data.reshape(T, num_molecules, 3)

velocity_coord = velocity_data.reshape(T_v, num_molecules, 3)
speed = np.linalg.norm(velocity_coord, axis=2)  # (T_v, num_molecules)
velocity_mean = np.mean(speed)

window_size = int(input("Window size: "))  
step_size = int(input("Step size: "))
cut_index = int(input("Cut first points: "))

print(f"\nAnalysis parameters:")
print(f"Window size: {window_size} points")
print(f"Window step: {step_size} points") 
print(f"Cut first {cut_index} points")

def compute_msd(trajectory, window_size, step_size):
    """
    Вычисление MSD с скользящим окном
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
                displacements = window_traj[tau] - window_traj[0]
                msd_results[window_idx, tau] = np.mean(displacements**2)
        
        if (window_idx + 1) % max(1, num_windows // 10) == 0:
            print(f"Processed {window_idx + 1}/{num_windows} windows")
    
    return msd_results

# Вычисляем MSD
print("\nComputing MSD...")
msd_results = compute_msd(trajectory, window_size, step_size)

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
mean_free_path = 3 * D / velocity_mean

print("\n" + "="*50)
print("RESULTS:")
print("="*50)
print(f"Diffusion coefficient D = {D:.6e}")
print(f"Slope = {slope:.6f}")
print(f"Standard error = {std_err:.6e}")
print(f"Mean free path = {mean_free_path:.6f}")
print(f"Cross section = {1/(mean_free_path * density * np.sqrt(2)):.4f}")
# print(f"R² = {r_value**2:.6f}")

plt.figure(figsize=(10, 6))
plt.plot(tau_cut, average_msd_cut, 'bo', alpha=0.6, markersize=3, label='MSD data')

fit_line = intercept + slope * tau_cut
plt.plot(tau_cut, fit_line, 'r-', linewidth=2, 
         label=f'Linear fit: MSD = {slope:.4f}t + {intercept:.4f}')

plt.xlabel('Time (reduced units)')
plt.ylabel('MSD (reduced units)')
# plt.title(f'Mean Squared Displacement\nD = {D:.4e}, R² = {r_value**2:.4f}')
plt.legend()
plt.grid(True, alpha=0.3)

textstr = f'D = {D:.4e}\nR² = {r_value**2:.4f}\nN = {num_molecules}'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=12,
         verticalalignment='top', bbox=props)

plt.tight_layout()
plt.savefig('msd_analysis.png', dpi=300)
plt.show()
#############################################################################

plt.figure(figsize=(10, 6))
plt.loglog(tau[1:], average_msd[1:], 'bo', alpha=0.6, markersize=3, label='MSD data')

x_fit = np.logspace(np.log10(tau[1:].min()), np.log10(tau[1:].max()), 100)
y_fit_diffusive = 6 * D * x_fit  # MSD = 6Dt
plt.plot(x_fit, y_fit_diffusive, 'r--', linewidth=2, 
         label=f'Diffusive slope (MSD ~ t)')

first_point_scale = average_msd[1] / (tau[1]**2)
y_fit_ballistic = first_point_scale * x_fit**2
plt.plot(x_fit, y_fit_ballistic, 'g--', linewidth=2, 
         label='Ballistic slope (MSD ~ t²)')

plt.xlabel('Time (reduced units)')
plt.ylabel('MSD (reduced units)')
plt.title('Mean Squared Displacement (Log-Log Scale)')
plt.legend()
plt.grid(True, alpha=0.3)

textstr = f'D = {D:.4e}'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=12,
         verticalalignment='top', bbox=props)

t_intersect_analytic = 6 * D / first_point_scale
if t_intersect_analytic >= x_fit.min() and t_intersect_analytic <= x_fit.max():
    t_intersect = t_intersect_analytic
    msd_intersect = 6 * D * t_intersect  # или a_ballistic * t_intersect**2
    
    print(f"\n=== Point of cross ===")
    print(f"Time of intersection: {t_intersect:.6f} ")
    print(f"MSD in point of intersection: {msd_intersect:.6f} ")
    print(f"Mean free path (from intersection) = {velocity_mean * t_intersect:.6f}")
    
    plt.plot(t_intersect, msd_intersect, 'ro', markersize=8, label='Intersection point')

plt.tight_layout()
plt.savefig('msd_analysis_log_log.png', dpi=300)
plt.show()

#############################################################################

print("\nAdditional information:")
print(f"Time range analyzed: {tau_cut[0]:.3f} to {tau_cut[-1]:.3f} (reduced units)")
print(f"MSD range: {average_msd_cut[0]:.3f} to {average_msd_cut[-1]:.3f} (reduced units)")
