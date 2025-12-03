import matplotlib.pyplot as plt 
import numpy as np
import json
from scipy import stats

def calculate_rms_velocity(velocity_data):
    velocity_squared = np.sum(velocity_data**2, axis=2)  # форма: (T_v, num_molecules)
    
    mean_squared_velocity = np.mean(velocity_squared)
    
    rms_velocity = np.sqrt(mean_squared_velocity)
    
    return rms_velocity


def precompute_normalized_velocities(velocity_per_molecule):
    """
    Предварительно вычисляет нормированные векторы скоростей для всего массива
    
    Parameters:
    velocity_per_molecule: array of shape (num_molecules, T_v, 3)
    
    Returns:
    normalized_velocities: array of shape (num_molecules, T_v, 3) с нормированными векторами
    valid_mask: boolean mask indicating which vectors have non-zero norm
    """
    # Вычисляем нормы всех векторов
    norms = np.linalg.norm(velocity_per_molecule, axis=2)  # форма: (num_molecules, T_v)
    
    # Создаем маску для ненулевых норм
    valid_mask = norms > 1e-10
    
    # Создаем массив для нормированных векторов
    normalized_velocities = np.zeros_like(velocity_per_molecule)
    
    # Нормируем только валидные векторы
    normalized_velocities[valid_mask] = velocity_per_molecule[valid_mask] / norms[valid_mask, np.newaxis]
    
    return normalized_velocities, valid_mask


def calculate_rotation_time_optimized_v2(normalized_velocities, valid_mask, angle_threshold_rad, 
                                        max_time_window=1000, skip_steps_after_found=1,
                                        step_start_points=1):
    num_molecules, T_v, _ = normalized_velocities.shape
    cos_threshold = np.cos(angle_threshold_rad)
    times = []
    
    max_time_window = min(max_time_window, T_v - 1)
    
    for i in range(num_molecules):
        t_start = 0
        while t_start < T_v - 1:
            if not valid_mask[i, t_start]:
                t_start += step_start_points
                continue
                
            v0_norm = normalized_velocities[i, t_start]
            found_time = False
            
            max_dt = min(max_time_window, T_v - t_start - 1)
            
            for dt in range(1, max_dt + 1):
                t_current = t_start + dt
                
                if not valid_mask[i, t_current]:
                    continue
                    
                v1_norm = normalized_velocities[i, t_current]
                cos_angle = np.dot(v0_norm, v1_norm)
                
                if cos_angle < cos_threshold:
                    times.append(dt)
                    found_time = True
                    t_start = min(t_start + dt + skip_steps_after_found, T_v - 1)
                    break
            
            if not found_time:
                t_start += step_start_points
    
    return np.mean(times) if times else 0.0, times

def analyze_velocity_data_optimized(velocity_data, angle_threshold_degrees, 
                                   time_start=270000, time_step=5, mol_step=2,
                                   max_time_window=500, skip_steps_after_found=10,
                                   step_start_points=1):
    print(f"Original data shape: {velocity_data.shape}")
    velocity_per_molecule = velocity_data.transpose(1, 0, 2)
    velocity_optimized = velocity_per_molecule[::mol_step, time_start::time_step, :]
    print(f"Optimized shape: {velocity_optimized.shape}")
    
    print("Precomputing normalized velocities...")
    normalized_velocities, valid_mask = precompute_normalized_velocities(velocity_optimized)
    
    rms_velocity = calculate_rms_velocity(velocity_data[time_start::time_step, ::mol_step, :])
    
    angle_threshold_rad = np.radians(angle_threshold_degrees)
    
    print("Calculating rotation time...")
    average_time_dt, all_times = calculate_rotation_time_optimized_v2(
        normalized_velocities, 
        valid_mask, 
        angle_threshold_rad, 
        max_time_window,
        skip_steps_after_found,
        step_start_points
    )
    
    # Дополнительная статистика
    if all_times:
        time_std = np.std(all_times)
        time_median = np.median(all_times)
    else:
        time_std = 0
        time_median = 0
    
    results = {
        'rms_velocity': rms_velocity,
        'rotation_time_dt': average_time_dt,
        'mean_free_path': average_time_dt * rms_velocity,
        'time_statistics': {
            'std': time_std,
            'median': time_median,
            'min': min(all_times) if all_times else 0,
            'max': max(all_times) if all_times else 0,
            'count': len(all_times)
        },
        'optimization_info': {
            'time_points_used': velocity_optimized.shape[1],
            'molecules_used': velocity_optimized.shape[0],
            'reduction_factor': velocity_per_molecule.size / velocity_optimized.size,
            'skip_steps_after_found': skip_steps_after_found,
            'step_start_points': step_start_points
        }
    }
    
    return results

def main_optimized(velocity_data, angle_threshold_degrees):
    """
    Оптимизированная основная функция
    """
    print(f"Number of time points: {velocity_data.shape[0]}")
    print(f"Number of molecules: {velocity_data.shape[1]}")
    
    # Параметры оптимизации (можно настроить)
    time_start = 0  # начинаем с n-й точки
    time_step = 1        # берем каждую n-ю точку по времени
    mol_step = 1         # берем каждую n-ю молекулу
    max_time_window = 10**6 # ограничиваем окно поиска n шагами
    skip_steps_after_found = 1
    step_start_points=5

    # Анализ с оптимизацией
    results = analyze_velocity_data_optimized(
        velocity_data, 
        angle_threshold_degrees,
        time_start=time_start,
        time_step=time_step,
        mol_step=mol_step,
        max_time_window=max_time_window,
        skip_steps_after_found=skip_steps_after_found,
        step_start_points=step_start_points
    )
    
    return results

if __name__ == "__main__":
    # Параметры анализа
    angle_threshold_degrees = 90  # порог в градусах
    
    # Загрузка конфигурации
    with open("cfg/cfg.json", "r") as f:
        cfg = json.load(f)
        coord_abs_file = cfg["coord_abs_file"]
        velocity_file = cfg["velocity_file"]
        num_molecules = cfg["num_molecules"]
        total_steps = cfg["total_steps"]
        dt = cfg["dt"]  
        snapshot = cfg["snapshot"]

    # Загрузка и подготовка данных о скоростях
    velocity_data_loaded = np.loadtxt(velocity_file)
    velocity_data = velocity_data_loaded.reshape(-1, num_molecules, 3)
    
    # Анализ данных
    results = main_optimized(velocity_data, angle_threshold_degrees)
    
    # Вывод результатов
    print("\n=== РЕЗУЛЬТАТЫ АНАЛИЗА ===")
    print(f"Корень средней квадратичной скорости: {results['rms_velocity']:.4f}")
    print(f"Среднее время изменения угла на {angle_threshold_degrees}°: {results['rotation_time_dt']:.4f} dt")
    print(f"Длина свободного пробега: {results['mean_free_path']*dt:.6f}")
    print(f"Анализ. {results['time_statistics']}")
