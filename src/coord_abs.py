import matplotlib.pyplot as plt 
import numpy as np
import json

with open("cfg/cfg.json", "r") as f:
    cfg = json.load(f)
    coord_abs_file = cfg["coord_abs_file"]
    num_molecules = cfg["num_molecules"]
    total_steps = cfg["total_steps"]
    snapshot = cfg["snapshot"]
    dt = cfg["dt"]
    f.close()

coord_data = np.loadtxt(coord_abs_file)
rms_values = []
for step in coord_data:
    molecules = step.reshape(num_molecules, 3)
    sq_distances = np.sum(molecules**2, axis=1)
    rms = np.sqrt(np.mean(sq_distances))
    rms_values.append(rms)

# Построение графика
steps = np.arange(len(rms_values))
steps = [i*snapshot for i in steps]
plt.figure(figsize=(10, 6))
plt.plot(steps, rms_values, 'r-', linewidth=2)
plt.xlabel('Номер шага')
plt.ylabel('Среднеквадратичное абсолютное положение')
plt.title('Зависимость RMS положения от номера шага')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

