import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, Function, Eq, dsolve, lambdify

# Definición de variables para SymPy
x = symbols('x')
y = Function('y')(x)
y1 = Function('y1')(x)
y2 = Function('y2')(x)

# Definición de las ecuaciones diferenciales
eq1 = Eq(y.diff(x), -2 * x * y)  # Primera ecuación de primer orden
eq2 = Eq(y.diff(x), y**2 + 1)  # Segunda ecuación de primer orden
eq3 = Eq(y.diff(x), x**2 - y)  # Tercera ecuación de primer orden
eq4 = Eq(y.diff(x, x) + 2 * y.diff(x) + 5 * y, 0)  # Ecuación de segundo orden

# Definición del nuevo sistema de ecuaciones diferenciales
eq_y1 = Eq(y1.diff(x), -3 * y1 + 2)  # dy1/dx = -3y1 + 2
eq_y2 = Eq(y2.diff(x), 4 * y1 + y2)  # dy2/dx = 4y1 + y2

# Resolver analíticamente las ecuaciones 1, 2, 3, 4
sol_eq1 = dsolve(eq1, y)
sol_eq2 = dsolve(eq2, y)
sol_eq3 = dsolve(eq3, y)
sol_eq4 = dsolve(eq4, y)

# Resolver analíticamente el nuevo sistema de ecuaciones diferenciales
sol_system = dsolve([eq_y1, eq_y2])

# Crear funciones numéricas para las soluciones analíticas
y1_analytic_system = lambdify(x, sol_system[0].rhs.subs('C1', 1).subs('C2', 0))
y2_analytic_system = lambdify(x, sol_system[1].rhs.subs('C1', 1).subs('C2', 0))

# Método de Euler para sistemas de ecuaciones
def euler_method_system(f, x0, Y0, h, n_steps):
    x_values = [x0]
    y1_values = [Y0[0]]
    y2_values = [Y0[1]]
    Y = np.array(Y0)

    for i in range(n_steps):
        Y_next = Y + h * f(x_values[-1], Y)
        x_next = x_values[-1] + h
        x_values.append(x_next)
        y1_values.append(Y_next[0])
        y2_values.append(Y_next[1])
        Y = Y_next

    return np.array(x_values), np.array(y1_values), np.array(y2_values)

# Método de Euler para una sola ecuación
def euler_method(f, x0, y0, h, n_steps):
    x_values = [x0]
    y_values = [y0]

    for i in range(n_steps):
        y_next = y_values[-1] + h * f(x_values[-1], y_values[-1])
        x_next = x_values[-1] + h
        x_values.append(x_next)
        y_values.append(y_next)

    return np.array(x_values), np.array(y_values)

# Definir las ecuaciones para Euler
def equation_1(x, y):
    return -2 * x * y

def equation_2(x, y):
    return y**2 + 1

def equation_3(x, y):
    return x**2 - y

def system_equation_4(x, Y):
    y1, y2 = Y  # Aquí y1 es y, y2 es y'
    dy1_dx = y2
    dy2_dx = -2 * y2 - 5 * y1  # De la ecuación d²y/dx² + 2dy/dx + 5y = 0
    return np.array([dy1_dx, dy2_dx])

# Parámetros ajustados para Euler
h_eq2 = 0.001  # Tamaño del paso reducido para la ecuación 2
h_eq3 = 0.001  # Tamaño del paso reducido para la ecuación 3
h_eq4 = 0.001  # Tamaño del paso reducido para la ecuación 4
h_system = 0.001  # Tamaño del paso reducido para el sistema de ecuaciones

# Número de pasos ajustados
n_steps_eq2 = int(1.4 / h_eq2)  # Ajustado al rango de la ecuación 2
n_steps_eq3 = 5000  # Manteniendo un rango grande con el nuevo h
n_steps_eq4 = 5000  # Ajustado para la ecuación 4
n_steps_system = 5000  # Ajustado para el sistema de ecuaciones

# Resolver numéricamente las ecuaciones individuales con los nuevos pasos
x_vals_euler_1, y_vals_euler_1 = euler_method(equation_1, 0, 1, h_eq3, n_steps_eq3)
x_vals_euler_2, y_vals_euler_2 = euler_method(equation_2, 0, 0, h_eq2, n_steps_eq2)
x_vals_euler_3, y_vals_euler_3 = euler_method(equation_3, 0, 1, h_eq3, n_steps_eq3)

# Resolver el sistema de ecuaciones
x_vals_system, y1_vals_system, y2_vals_system = euler_method_system(system_equation_4, 0, [1, 0], h_system, n_steps_system)

# Crear valores para soluciones analíticas
x_vals = np.linspace(0, 5, 100)
y_vals_analytic_1 = lambdify(x, sol_eq1.rhs.subs('C1', 1))(x_vals)
y_vals_analytic_2 = lambdify(x, sol_eq2.rhs.subs('C1', 0))(x_vals[:n_steps_eq2])
y_vals_analytic_3 = lambdify(x, sol_eq3.rhs.subs('C1', 1))(x_vals)

# Imprimir las soluciones analíticas con descripción de cada ecuación
print("===== Soluciones Analíticas =====\n")

print("Solución analítica de la ecuación 1 (dy/dx = -2xy):")
print(f"y(x) = {sol_eq1.rhs}\n")

print("Solución analítica de la ecuación 2 (dy/dx = y^2 + 1):")
print(f"y(x) = {sol_eq2.rhs}\n")

print("Solución analítica de la ecuación 3 (dy/dx = x^2 - y):")
print(f"y(x) = {sol_eq3.rhs}\n")

print("Solución analítica de la ecuación 4 (d²y/dx² + 2dy/dx + 5y = 0):")
print(f"y(x) = {sol_eq4.rhs}\n")

print("Solución analítica del sistema de ecuaciones (dy1/dx = -3y1 + 2, dy2/dx = 4y1 + y2):")
print(f"y1(x) = {sol_system[0].rhs}")
print(f"y2(x) = {sol_system[1].rhs}\n")

print("=================================")


# Graficar todas las soluciones en una figura con subgráficas
plt.figure(figsize=(16, 14))

# Ecuación 1
plt.subplot(3, 2, 1)
plt.plot(x_vals, y_vals_analytic_1, label='Solución Analítica', linestyle='dashed', color='red')
plt.plot(x_vals_euler_1, y_vals_euler_1, label='Método de Euler', color='blue')
plt.title("Ecuación 1: dy/dx = -2xy")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid()

# Ecuación 2
plt.subplot(3, 2, 2)
plt.plot(x_vals[:n_steps_eq2], y_vals_analytic_2, label='Solución Analítica', linestyle='dashed', color='red')
plt.plot(x_vals_euler_2, y_vals_euler_2, label='Método de Euler', color='blue')
plt.title("Ecuación 2: dy/dx = y^2 + 1 (Rango Ajustado)")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid()

# Ecuación 3
plt.subplot(3, 2, 3)
plt.plot(x_vals, y_vals_analytic_3, label='Solución Analítica', linestyle='dashed', color='red')
plt.plot(x_vals_euler_3, y_vals_euler_3, label='Método de Euler', color='blue')
plt.title("Ecuación 3: dy/dx = x^2 - y")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid()

# Sistema de ecuaciones diferenciales
plt.subplot(3, 2, 4)
plt.plot(x_vals_system, y1_vals_system, label='y1 (Euler)', color='blue')
plt.plot(x_vals_system, y2_vals_system, label='y2 (Euler)', color='orange')
plt.plot(x_vals, y1_analytic_system(x_vals), label='y1 (Analítica)', linestyle='dashed', color='blue')
plt.plot(x_vals, y2_analytic_system(x_vals), label='y2 (Analítica)', linestyle='dashed', color='orange')
plt.title("Sistema de Ecuaciones: dy1/dx = -3y1 + 2, dy2/dx = 4y1 + y2")
plt.xlabel("x")
plt.ylabel("y1, y2")
plt.legend()
plt.grid()

# Ecuación de segundo orden (convertida en sistema)
x_vals_euler_4, y1_vals_euler_4, y2_vals_euler_4 = euler_method_system(system_equation_4, 0, [1, 0], h_eq4, n_steps_eq4)
y_vals_analytic_4 = lambdify(x, sol_eq4.rhs.subs('C1', 1).subs('C2', 0))(x_vals)

plt.subplot(3, 2, 5)
plt.plot(x_vals, y_vals_analytic_4, label='Solución Analítica', linestyle='dashed', color='red')
plt.plot(x_vals_euler_4, y1_vals_euler_4, label='Método de Euler', color='blue')
plt.title("Ecuación 4: d²y/dx² + 2dy/dx + 5y = 0")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid()

# Ajustar el diseño de las subgráficas
plt.tight_layout()
plt.subplots_adjust(hspace=0.5, wspace=0.4)  # Más espacio entre las gráficas
plt.show()
