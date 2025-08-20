import matplotlib.pyplot as plt 
import numpy as np
import json

NUM_BARS = 40  # number of points on plot
THERM_BALANCE = 150  # number of first sample after thermodynamic balance 
EXCLUDE_LAST = 15  # exclude last EXCLUDE_LAST points

with open("cfg/cfg.json", "r") as f:
    cfg = json.load(f)
    file_data = cfg["velocity_file"]
    epsilone_real = cfg["epsilon_real"]
    m_real = cfg["mass_real"]
    num_molecules = cfg["num_molecules"]
    total_steps = cfg["total_steps"]
    snapshot = cfg["snapshot"]
    f.close()

data_generated = np.loadtxt(file_data)
vel_coord_2 = [[], [], []]
velocity = []
kinetic_energy_per_step = []
for i in range(THERM_BALANCE, len(data_generated)):
    step_energy = 0
    for j in range(0, len(data_generated[i]), 3):
        velocity.append(np.sqrt(data_generated[i][j]**2 + data_generated[i][j+1]**2 + data_generated[i][j+2]**2))
        vel_coord_2[0].append(data_generated[i][j]**2)
        vel_coord_2[1].append(data_generated[i][j+1]**2)
        vel_coord_2[2].append(data_generated[i][j+2]**2)

velocity.sort()
for i in range(3):
    vel_coord_2[i].sort()

vel_max = velocity[-1]
vel_coord_max_2 = [vel_coord_2[0][-1], vel_coord_2[1][-1], vel_coord_2[2][-1]]

prob = np.zeros(NUM_BARS)

for i in range(len(velocity)):
    for j in range(NUM_BARS):
        if (j/NUM_BARS*vel_max <= velocity[i] < (j+1)/NUM_BARS*vel_max):
            prob[j] += 1/len(velocity)

all_squared = np.concatenate(vel_coord_2)
vel_np = np.array(velocity)
v2_mean = vel_np.mean()**2

total_max = max(vel_coord_max_2)
bin_edges = np.linspace(0, total_max, NUM_BARS + 1)
hist, _ = np.histogram(all_squared, bins=bin_edges)
prob_avg = hist / len(all_squared)

mid_bins = (bin_edges[:-1] + bin_edges[1:]) / 2

valid_mask = prob_avg > 0
x_data = mid_bins[valid_mask][1:-EXCLUDE_LAST]  
y_log = np.log(prob_avg[valid_mask])[1:-EXCLUDE_LAST]  

coefficients = np.polyfit(x_data, y_log, deg=1)
k = coefficients[0]
b = coefficients[1]

x_fit = np.linspace(0, total_max, 200)
y_fit = np.exp(b + k * x_fit)

y_pred = k * x_data + b
ss_res = np.sum((y_log - y_pred)**2)
ss_tot = np.sum((y_log - np.mean(y_log))**2)
r_squared = 1 - (ss_res / ss_tot)

########################################################################## 
print(f"Экспонента: y = {np.exp(b):.6f} * exp({k:.6f}x)")
print(f"R² = {r_squared:.6f}\n")

print(f"Средний квадрат скорости = {v2_mean:.6f}")
print(f"Температура через энергию = {v2_mean/3:.6f}")
print(f"Температура через аппроксимации = {-1/(2*k):.6f}")
# error_pct = abs(k_B_corrected - BOLTZMANN) / BOLTZMANN * 100
# print(f"Погрешность: {error_pct:.6f}%")
########################################################################## 

# kinetic_energy_per_step = []
# for i in range(len(data_generated)):
#     step_energy = 0.5 * np.sum(data_generated[i]**2)
#     kinetic_energy_per_step.append(step_energy)
#
# plt.figure(figsize=(12, 6))
# plt.plot(kinetic_energy_per_step)
# plt.xlabel('Номер шага')
# plt.ylabel('Кинетическая энергия системы')
# plt.title('Эволюция кинетической энергии')
# # plt.legend()
# plt.grid(True)
# plt.savefig('energy_evolution.png')
# plt.show()
#
########################################################################## 
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(8, 10))
plt.subplots_adjust(hspace=0.4)  

axs[0].plot(
    [i/NUM_BARS*vel_max for i in range(NUM_BARS)],
    prob,
    'o', 
    color='red',
    markersize=5,
    mec='black'
)
axs[0].set_xlabel('$v$, абсолютная скорость молекулы')
axs[0].set_ylabel('Вероятность')
axs[0].set_title('Распределение абсолютных скоростей', loc='left')
axs[0].grid(True)

axs[1].semilogy(
    mid_bins,
    prob_avg,
    'o',
    color='blue',
    markersize=5,
    mec='black',
    label='Данные'
)

axs[1].semilogy(
    x_fit,
    y_fit,
    'r-',
    linewidth=2,
    label=f'Аппроксимация: $y = {np.exp(b):.3f}e^{{{k:.3f}x}}$'
)

axs[1].set_xlabel('$v^2$, квадрат компоненты скорости')
axs[1].set_ylabel('Вероятность')
axs[1].set_title('Усредненное распределение квадрата компоненты скорости', loc='left')
axs[1].grid(True)

fig.suptitle('Распределения скоростей молекул', fontsize=16)
plt.savefig('velocity_plot.png')
plt.show()
########################################################################## 
print('Clear file ' + file_data + ': y/n?')
if (input().lower() == 'y'):
    with open(file_data, 'w') as file:
        file = ''
