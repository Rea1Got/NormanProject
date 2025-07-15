import matplotlib.pyplot as plt 
import numpy as np
import json

NUM_BARS = 30
THERM_BALANCE = 500  # number of first sample after thermodynamic balance 

with open("cfg/cfg.json", "r") as f:
    cfg = json.load(f)
    file_data = cfg["velocity_file"]
    num_molecules = cfg["num_molecules"]
    total_steps = cfg["total_steps"]
    snapshot = cfg["snapshot"]
    f.close()

data_generated = np.loadtxt(file_data)
vel_coord = [[], [], []]
velocity = []
for i in range(THERM_BALANCE, len(data_generated)):
    for j in range(0, len(data_generated[i]), 3):
        velocity.append(np.sqrt(data_generated[i][j]**2 + data_generated[i][j+1]**2 + data_generated[i][j+2]**2))
        vel_coord[0].append(data_generated[i][j]**2)
        vel_coord[1].append(data_generated[i][j+1]**2)
        vel_coord[2].append(data_generated[i][j+2]**2)
velocity.sort()
for i in range(3):
    vel_coord[i].sort()

vel_max = velocity[-1]
vel_coord_max = [vel_coord[0][-1], vel_coord[1][-1], vel_coord[2][-1]]

prob = np.zeros(NUM_BARS)
# prob_coord = [np.zeros(NUM_BARS) for _ in range(3)]

for i in range(len(velocity)):
    for j in range(NUM_BARS):
        if (j/NUM_BARS*vel_max <= velocity[i] < (j+1)/NUM_BARS*vel_max):
            prob[j] += 1/len(velocity)
        # for k in range(3):
        #     if (j/NUM_BARS*vel_coord_max[k] <= vel_coord[k][i] < (j+1)/NUM_BARS*vel_coord_max[k]):
        #         prob_coord[k][j] += 1/len(velocity)

all_squared = np.concatenate(vel_coord)
total_max = max(vel_coord_max)
bin_edges = np.linspace(0, total_max, NUM_BARS + 1)
hist, _ = np.histogram(all_squared, bins=bin_edges)
prob_avg = hist / len(all_squared)

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
    bin_edges[:-1],
    prob_avg,
    'o',
    color='blue',
    markersize=5,
    mec='black'
)
axs[1].set_xlabel('$v^2$, квадрат компоненты скорости')
axs[1].set_ylabel('Вероятность')
axs[1].set_title('Усредненное распределение квадрата компоненты скорости', loc='left')
axs[1].grid(True)

fig.suptitle('Распределения скоростей молекул', fontsize=16)
plt.show()
########################################################################## 
print('Clear file ' + file_data + ': y/n?')
if (input().lower() == 'y'):
    with open(file_data, 'w') as file:
        file = ''
