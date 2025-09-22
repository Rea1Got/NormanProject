import matplotlib.pyplot as plt 
import numpy as np
import json
from scipy.optimize import curve_fit

WINDOW_SIZE_PERCENT = 1  # Размер окна в процентах от общего времени (рекомендуется 5-20%)
STEP_SIZE_PERCENT = 5   # Шаг сдвига окна в процентах от размера окна (рекомендуется 5-20%)
MAX_LAG_PERCENT = 50      # Максимальный лаг в процентах от размера окна (рекомендуется 10-30%)
MIN_LAG_POINTS = 5        # Минимальное количество точек для аппроксимации

# Загрузка конфигурации
with open("cfg/cfg.json", "r") as f:
    cfg = json.load(f)
    coord_abs_file = cfg["coord_abs_file"]
    num_molecules = cfg["num_molecules"]
    total_steps = cfg["total_steps"]
    snapshot = cfg["snapshot"]
    dt = cfg["dt"]

# Чтение данных координат
coord_data = np.loadtxt(coord_abs_file)
T = len(coord_data)  # Общее количество временных точек
print(f"Number of points: {T}")
# Преобразование данных в массив с размерностью (время, молекулы, координаты)
trajectory = coord_data.reshape(T, num_molecules, 3)

# Функция для линейной аппроксимации
def linear_func(t, a, b):
    return a * t + b

# Вычисление параметров анализа на основе процентных значений
print("Enter window size and step size:")
window_size = int(input())  # Размер окна
step_size = int(input())  # Шаг сдвига окна

print(f"Analysis of trajectory with {T} time points and {num_molecules} molecules.")
print(f"Window size: {window_size} points ({window_size * 100 / T}% out of points)")
print(f"Window step size: {step_size} points ({step_size * 100 / window_size}% out of window size)")

# Массивы для хранения результатов
diffusion_coeffs = []
time_values = []
msd_curves = []

# Анализ для каждого скользящего окна


# Усреднение коэффициентов диффузии
if diffusion_coeffs:
    avg_D = np.mean(diffusion_coeffs)
    std_D = np.std(diffusion_coeffs)
    print(f"Средний коэффициент диффузии: {avg_D:.4e} ± {std_D:.4e}")
    print(f"Относительная погрешность: {std_D/avg_D*100:.5}%") 
    # Визуализация результатов
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # График 1: Распределение коэффициентов диффузии
    ax1.hist(diffusion_coeffs, bins=20, alpha=0.7, edgecolor='black')
    ax1.axvline(avg_D, color='r', linestyle='--', label=f'Среднее: {avg_D:.4e}')
    ax1.set_xlabel('Коэффициент диффузии')
    ax1.set_ylabel('Частота')
    ax1.set_title('Распределение коэффициентов диффузии по скользящим окнам')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # График 2: MSD кривые для нескольких окон
    colors = plt.cm.viridis(np.linspace(0, 1, min(5, len(msd_curves))))
    for i, msd in enumerate(msd_curves[-6:-1]):
        time_lags = np.arange(len(msd)) * dt * snapshot
        ax2.plot(time_lags[1:], msd[1:], color=colors[i], alpha=0.7, 
                label=f'Окно {i+1}, D={diffusion_coeffs[i]:.2e}')
    
    # Добавляем теоретическую линию для среднего коэффициента диффузии
    theoretical_msd = 6 * avg_D * time_lags
    ax2.plot(time_lags, theoretical_msd, 'r--', linewidth=2, 
            label=f'Теоретическая (D={avg_D:.2e})')
    
    ax2.set_xlabel('Время')
    ax2.set_ylabel('MSD')
    ax2.set_title('Кривые MSD для различных окон')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    # Дополнительный анализ: зависимость коэффициента диффузии от времени
    plt.figure(figsize=(10, 6))
    plt.plot(time_values, diffusion_coeffs, 'o-', alpha=0.7)
    plt.axhline(avg_D, color='r', linestyle='--', label=f'Среднее: {avg_D:.4e}')
    plt.xlabel('Время начала окна')
    plt.ylabel('Коэффициент диффузии')
    plt.title('Зависимость коэффициента диффузии от времени начала окна')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    
else:
    print("Не удалось вычислить коэффициент диффузии")
