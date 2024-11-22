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

#### Ejecución:
```bash
python Euler.py

---

### runge_kutta.py

```markdown
### 2. `runge_kutta.py`
Este script implementa el **método de Runge-Kutta de cuarto orden (RK4)** para resolver ecuaciones diferenciales de forma más precisa que Euler.

#### Funcionalidades:
- Resolución de ecuaciones diferenciales de primer orden:
  - `dy/dx = y - x`
- Resolución de ecuaciones de segundo orden reescritas como sistemas:
  - `y'' + 2y' + y = 0`
- Resolución de sistemas de ecuaciones diferenciales:
  - `dx/dt = x + y`
  - `dy/dt = x - y`

#### Uso:
- Resuelve las ecuaciones numéricamente utilizando RK4 y, en algunos casos, compara los resultados con las soluciones analíticas.
- Genera gráficos que representan las soluciones numéricas y, cuando es posible, la solución analítica.

#### Dependencias:
- `numpy`
- `matplotlib`
- `math`

#### Ejecución:
```bash
python runge_kutta.py
