import matplotlib.pyplot as plt 
import numpy as np
import matplotlib as mpl
import json

NUM_BARS = 25
THERM_BALANCE = 100  # number of first sample after thermodynamic balance 

with open("cfg/cfg.json", "r") as f:
    cfg = json.load(f)
    num_molecules = cfg["num_molecules"]
    total_steps = cfg["total_steps"]
    snapshot = cfg["snapshot"]
    f.close()

data_generated = np.loadtxt("tmp/velocity.txt")
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
prob_coord = [np.zeros(NUM_BARS) for _ in range(3)]

for i in range(len(velocity)):
    for j in range(NUM_BARS):
        if (j/NUM_BARS*vel_max <= velocity[i] < (j+1)/NUM_BARS*vel_max):
            prob[j] += 1/len(velocity)
        for k in range(3):
            if (j/NUM_BARS*vel_coord_max[k] <= vel_coord[k][i] < (j+1)/NUM_BARS*vel_coord_max[k]):
                prob_coord[k][j] += 1/len(velocity)

########################################################################## 
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 10))
plt.subplots_adjust(wspace=0.3, hspace=0.4)  
#######################
axs[0, 0].plot(
    [i/NUM_BARS*vel_max for i in range(NUM_BARS)],
    prob,
    'o', 
    color='red',
    markersize=5,
    mec='black'
)
axs[0, 0].set_xlabel('$v$, абсолютная скорость молекулы')
axs[0, 0].set_ylabel('Вероятность')
axs[0, 0].set_title('Распределение абсолютных скоростей', loc='left')
axs[0, 0].grid(True)
#######################
titles = ['X-компонента скорости', 'Y-компонента скорости', 'Z-компонента скорости']
positions = [(0,1), (1,0), (1,1)]

for i in range(3):
    row, col = positions[i]
    
    ax_component = [j/NUM_BARS*vel_coord_max[i] for j in range(NUM_BARS)]
    
    axs[row, col].semilogy(
        ax_component,
        prob_coord[i],
        'o',
        color='blue' if i==0 else 'green' if i==1 else 'purple',
        markersize=5,
        mec='black'
    )
    axs[row, col].set_xlabel(f'$v_{"xyz"[i]}^2$, компонента скорости')
    axs[row, col].set_ylabel('Вероятность')
    axs[row, col].set_title(titles[i], loc='left')
    axs[row, col].grid(True)

fig.suptitle('Распределения скоростей молекул', fontsize=16)
plt.show()
########################################################################## 
print('Clear files? y/n')
if (input().lower() == 'y'):
    with open('tmp/velocity.txt', 'w') as file:
        file = ''

