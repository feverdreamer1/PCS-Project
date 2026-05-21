# Modelo de Ising 2D: Transiciones de Fase y Universalidad

Este repositorio contiene el código fuente, los datos y los scripts de análisis para el estudio computacional del modelo de Ising. Se analizan tres topologías diferentes mediante simulaciones Monte Carlo (Algoritmo de Metrópolis) para estudiar sus observables termodinámicos y demostrar la hipótesis de Universalidad.

## 📂 Estructura del Repositorio

### Archivos Principales
* **`plot1d.py`**: Script didáctico que compara el modelo en 1D vs 2D. Demuestra empíricamente que en 1D el ruido térmico destruye el orden a cualquier temperatura $T>0$, mientras que en 2D sí existe una fase ferromagnética estable.
* **`plot_ising.py`**: Generador de animaciones. Produce el archivo `ising_termometro.gif`, mostrando la evolución visual de los dominios magnéticos de la red a medida que varía la temperatura.
* **`comparison_1D_vs_2D.png`**: Gráfica de salida del script comparativo.
* **`LICENSE`**: Licencia pública (GNU GPLv3) del código.

### Directorio `NNN/` (Red Cuadrada: Primeros y Segundos Vecinos)
Contiene la simulación para la red cuadrada clásica (NN) y la extendida con interacciones diagonales (NNN).
* **`nnn.f90`** y **`randomnumber.f`**: Código fuente en Fortran 90 que ejecuta la simulación Monte Carlo.
* **`termodinamica.dat` / `termodinamica-next-next.dat`**: Datos termodinámicos crudos ($m$, $e$, $c_V$).
* **`correlacion.dat` / `correlacion-next-next.dat`**: Datos de la función de correlación espacial $f(r)$.
* **`plot.py`**: Script en Python que lee los datos y genera las gráficas de termodinámica, el escalamiento de tamaño finito para $T_c$, y extrae los exponentes críticos $\beta$ y $\nu$.
* **`fig_*.png`**: Gráficas de resultados ya procesadas.

### Directorio `Triangle/` (Red Triangular)
Contiene la adaptación del modelo para una topología donde cada espín interactúa con 6 vecinos.
* **`Triangulo.f90`**: Código fuente en Fortran 90 adaptado para la red triangular.
* **`termodinamica_triangular.dat` / `correlacion_triangular.dat`**: Datos crudos de la simulación.
* **`plot.py`**: Script en Python específico para procesar y graficar los resultados del modelo triangular, contrastándolos con la solución exacta de Wannier.
* **`fig_*.png`**: Gráficas de resultados ya procesadas.

## 🚀 Uso rápido

1.  **Simulaciones (Fortran):** Compila los archivos fuente usando un compilador estándar como `gfortran`.
    ```bash
    # Ejemplo para el modelo NNN
    gfortran nnn.f90 randomnumber.f -o ising_sim
    ./ising_sim
    ```
2.  **Análisis y Gráficas (Python):** Asegúrate de tener instalados `numpy`, `matplotlib` y `scipy`. Ejecuta los scripts dentro de sus respectivas carpetas para regenerar los resultados.
    ```bash
    python plot.py
    ```
