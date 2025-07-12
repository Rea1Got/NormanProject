import matplotlib.pyplot as plt 
import numpy as np
import matplotlib as mpl
import json

num_bars = 20

with open("cfg/cfg.json", "r") as f:
    cfg = json.load(f)
    num_molecules = cfg["num_molecules"]
    total_steps = cfg["total_steps"]
    snapshot = cfg["snapshot"]
    f.close()

data_generated = np.loadtxt("tmp/velocity.txt")
velocity = []
for i in range(len(data_generated)):
    for j in range(0, len(data_generated[i]), 3):
        velocity.append(np.sqrt(data_generated[i][j]**2 + data_generated[i][j+1]**2 + data_generated[i][j+2]**2))
velocity.sort()

vel_max = velocity[-1]
prob = np.zeros(num_bars)

for i in range(len(velocity)):
    for j in range(num_bars):
        if (j/num_bars*vel_max <= velocity[i] < (j+1)/num_bars*vel_max):
            prob[j] +=  1/len(velocity)
     
##########################################################################
ax = []
for i in range(num_bars):
    ax.append(i/num_bars*vel_max)

fig = plt.figure(figsize=(12, 8))
########################
# plt.bar(ax, prob)
plt.plot(ax, prob, 'o-r', c='k', lw=0, alpha = 1, mec='r', mew=0, ms=5, label = 'Усредненное значение времени')
plt.xlabel('N, количество переменных')
plt.ylabel(r'$t, 10^{-3} с$')
plt.title('Рис. 1(a) Линейный поиск (случайное событие)', loc='left')
plt.grid(True)
plt.legend()
plt.show()
########################################################################## 


print('Clear files? y/n')
if (input().lower() == 'y'):
    with open('tmp/velocity.txt', 'w') as file:
        file = ''

