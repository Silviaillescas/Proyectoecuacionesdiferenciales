# Métodos Numéricos para Ecuaciones Diferenciales Ordinarias (EDOs)

Este repositorio contiene dos implementaciones principales para resolver ecuaciones diferenciales ordinarias (EDOs): el **método de Euler** y el **método de Runge-Kutta de cuarto orden (RK4)**. Ambos scripts están diseñados para resolver ecuaciones de primer y segundo orden, así como sistemas de ecuaciones diferenciales.

---

### 1. `Euler.py`
Este script implementa el método de **Euler** para resolver ecuaciones diferenciales de forma numérica. También incluye ejemplos y graficación de las soluciones.

#### Funcionalidades:
- Resolución de ecuaciones diferenciales de primer orden:
  - `dy/dx = -2xy`
  - `dy/dx = y^2 + 1`
  - `dy/dx = x^2 - y`
- Resolución de ecuaciones de segundo orden reescritas como sistemas:
  - `d²y/dx² + 2dy/dx + 5y = 0`
- Resolución de sistemas de ecuaciones diferenciales:
  - `dy1/dx = -3y1 + 2`
  - `dy2/dx = 4y1 + y2`

#### Uso:
- Se resuelven las ecuaciones tanto numéricamente (método de Euler) como analíticamente.
- Genera gráficos comparativos entre las soluciones analíticas y numéricas, además de representar el error absoluto.

#### Dependencias:
- `numpy`
- `matplotlib`
- `sympy`


### 2. `runge_kutta.py`
Este script implementa el método de **Runge-Kutta de cuarto orden (RK4)** para resolver ecuaciones diferenciales de forma más precisa que Euler. También incluye ejemplos y graficación de las soluciones.

#### Funcionalidades:
- Resolución de ecuaciones diferenciales de primer orden:
  - `dy/dx = y - x`
- Resolución de ecuaciones de segundo orden reescritas como sistemas:
  - `d²y/dx² + 2dy/dx + y = 0`
- Resolución de sistemas de ecuaciones diferenciales:
  - `dx/dt = x + y`
  - `dy/dt = x - y`

#### Uso:
- Resuelve las ecuaciones numéricamente utilizando el método RK4.
- Genera gráficos comparativos entre las soluciones analíticas y numéricas cuando están disponibles.
- Representa el error absoluto entre las soluciones numéricas y analíticas.

#### Dependencias:
- `numpy`
- `matplotlib`
- `math`

#### Ejecución de ambos programas:
```bash
python runge_kutta.py
python Euler.py

