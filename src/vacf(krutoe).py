import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.integrate import cumtrapz
import json


with open("cfg/cfg.json", "r") as f:
    cfg = json.load(f)
    coord_abs_file = cfg["coord_abs_file"]
    velocity_file = cfg["velocity_file"]
    num_molecules = cfg["num_molecules"]
    total_steps = cfg["total_steps"]
    dt = cfg["dt"]  
    snapshot = cfg["snapshot"]
velocity_data_loaded = np.loadtxt(velocity_file)
velocities = velocity_data_loaded.reshape(-1, num_molecules, 3)

print(f"Original data shape: {velocities.shape}")

squared_speeds = np.sum(velocities**2, axis=2) 
mean_squared_speed = np.mean(squared_speeds)
print(f"Средний квадрат скорости (вручную): {mean_squared_speed:.6f}")

# velocities = velocities[250000:, :, :]
#
# print(f"Cut data shape: {velocities.shape}")
#
# T = velocities.shape[0]
# acf_val = []
# for i in range(num_molecules):
#     for j in range(3):
#         acf = [np.sum(velocities[:T-k] * velocities[k:]) for k in range(T)]
#         acf_val.append(acf)

velocities = velocities[150000:, :, :]  
print(f"Cut data shape: {velocities.shape}")

T = velocities.shape[0]
acf_val = []

print("Computing ACF...")
for i in range(num_molecules):
    if i % 20 == 0:  
        print(f"Processing molecule {i}/{num_molecules}")
    for j in range(3):
        v_data = velocities[:, i, j]  # Берем все точки времени
        
        acf = signal.correlate(v_data, v_data, mode='full', method='auto')
        acf_positive = acf[acf.size//2:]
        acf_norm = acf_positive / (len(v_data) - np.arange(len(acf_positive)))
        
        acf_val.append(acf_norm)


acf_val = np.array(acf_val)

acf_val = np.mean(acf_val, axis=0)  
print(acf_val[0])
acf_val_normal = acf_val / acf_val[0] 

t = np.arange(len(acf_val_normal))

plt.plot(t * dt, acf_val_normal, 'r-')
plt.xlabel('Время (у.е.)', fontsize=12)
plt.ylabel('АКФС', fontsize=12)
plt.ylim(-0.02, 1)
plt.title('Нормированная автокорреляционная функция скорости')
plt.grid(True)
plt.tight_layout()
plt.show()

t = np.arange(len(acf_val_normal)) * dt
y = acf_val

threshold = -1.00
cutoff_index = np.argmax(y < threshold)

if cutoff_index == 0 and y[0] >= threshold:
    cutoff_index = len(y)
    print("АКФ не опустилась ниже порога, интегрируем всю кривую")
else:
    print(f"Интегрируем АКФ до времени {t[cutoff_index]:.4f} (индекс {cutoff_index})")

t_cut = t[1:cutoff_index]
y_cut = y[1:cutoff_index]

cumulative_integral = np.zeros(len(y_cut))
for i in range(1, len(y_cut)):
    cumulative_integral[i] = cumulative_integral[i-1] + (y_cut[i-1] + y_cut[i]) * (t_cut[i] - t_cut[i-1]) / 2

D =  cumulative_integral[-1]
print(f"Коэффициент диффузии D = {D:.6f}")
print(f"<v^2> = {mean_squared_speed:.6f}")
print(f"Normal ACF integral = {cumulative_integral[-1]:.6f}")

# def exponential_fit(x, a, b, c):
#     """Экспоненциальный фит: D = a + b*exp(c*x)"""
#     # Используем ограничения, чтобы избежать переполнения
#     return a + b * np.exp(c * x)

plt.plot(1/t_cut, cumulative_integral, 'r-')
plt.xlabel('1/t', fontsize=12)
plt.ylabel('Значение интеграла', fontsize=12)
plt.title('Сходимость интеграла')
plt.grid(True)
plt.tight_layout()
plt.show()

plt.loglog(1/t_cut, cumulative_integral, 'r-')
plt.xlabel('1/t', fontsize=12)
plt.ylabel('Значение интеграла', fontsize=12)
plt.title('Сходимость интеграла')
plt.grid(True)
plt.tight_layout()
plt.show()


