import numpy as np
import matplotlib.pyplot as plt
import math

class EDOSolver:
    def __init__(self):
        pass  # Constructor vacío, no es necesario inicializar variables.

    # Método para realizar un paso del método de Runge-Kutta de cuarto orden (RK4).
    def rk4_paso(self, f, t, y, h):
        # Verifica si `y` es una lista, tupla o arreglo de NumPy (caso de sistemas de ecuaciones diferenciales).
        if isinstance(y, (list, tuple, np.ndarray)):
            # Calcula k1 para cada ecuación en el sistema. 
            # `fi` es una función que evalúa cada ecuación del sistema.
            k1 = [fi(*([t] + list(y))) for fi in f]
            # Calcula k2, evaluando las funciones en el punto intermedio (t + h/2) con y ajustado por k1.
            k2 = [fi(t + h/2, *(y + np.array(k1)*h/2)) for fi in f]
            # Calcula k3, también evaluado en el punto intermedio (t + h/2) con y ajustado por k2.
            k3 = [fi(t + h/2, *(y + np.array(k2)*h/2)) for fi in f]
            # Calcula k4, evaluando las funciones en el siguiente paso (t + h) con y ajustado por k3.
            k4 = [fi(t + h, *(y + np.array(k3)*h)) for fi in f]
            
            # Retorna el valor actualizado de `y` utilizando la fórmula de RK4:
            # y_n+1 = y_n + h/6 * (k1 + 2*k2 + 2*k3 + k4).
            return y + h/6 * (np.array(k1) + 2*np.array(k2) + 2*np.array(k3) + np.array(k4))
        else:
            # Caso para una ecuación diferencial simple (no es un sistema).
            # Calcula k1, evaluando la función en el punto actual (t, y).
            k1 = f(t, y)
            # Calcula k2 en el punto intermedio (t + h/2) ajustando y con k1.
            k2 = f(t + h/2, y + h*k1/2)
            # Calcula k3 también en el punto intermedio con y ajustado por k2.
            k3 = f(t + h/2, y + h*k2/2)
            # Calcula k4 en el siguiente paso (t + h) con y ajustado por k3.
            k4 = f(t + h, y + h*k3)
            
            # Retorna el valor actualizado de `y` usando la fórmula de RK4.
            return y + h/6 * (k1 + 2*k2 + 2*k3 + k4)
    
    def resolver(self, edo_tipo, f, condiciones_iniciales, t0, tf, h, 
                solucion_analitica=None):
        # Calcula el número de pasos dividiendo el intervalo de tiempo [t0, tf] entre el tamaño del paso h.
        n = int((tf - t0) / h)
        # Genera una lista de valores de tiempo uniformemente espaciados entre t0 y tf.
        t = np.linspace(t0, tf, n+1)
        
        # Caso 1: Resolver una EDO de primer orden.
        if edo_tipo == 'primer_orden':
            # Inicializa el arreglo para almacenar las soluciones de la variable dependiente.
            y = np.zeros(n+1)
            # Asigna la condición inicial al primer elemento de la solución.
            y[0] = condiciones_iniciales[0]
            
            # Itera a través de los pasos para calcular la solución numérica usando el método RK4.
            for i in range(n):
                y[i+1] = self.rk4_paso(f, t[i], y[i], h)
            
            # Si se proporciona una solución analítica, la evalúa en los puntos de tiempo.
            if solucion_analitica:
                y_analitica = np.array([solucion_analitica(ti) for ti in t])
                # Devuelve el tiempo, la solución numérica y la solución analítica.
                return t, y, y_analitica
            
            # Si no hay solución analítica, solo devuelve tiempo y solución numérica.
            return t, y
            
        # Caso 2: Resolver una EDO de segundo orden.
        elif edo_tipo == 'segundo_orden':
            # Inicializa una matriz para almacenar la solución (y y su derivada).
            y = np.zeros((2, n+1))
            # Asigna las condiciones iniciales a la primera columna (t = t0).
            y[:, 0] = condiciones_iniciales
            
            # Itera a través de los pasos para calcular las soluciones numéricas usando RK4.
            for i in range(n):
                y[:, i+1] = self.rk4_paso(f, t[i], y[:, i], h)
            
            # Si se proporcionan soluciones analíticas (y y su derivada), las evalúa en los puntos de tiempo.
            if isinstance(solucion_analitica, (list, tuple)) and len(solucion_analitica) == 2:
                y_analitica = [
                    np.array([solucion_analitica[0](ti) for ti in t]),  # Solución para y.
                    np.array([solucion_analitica[1](ti) for ti in t]) if solucion_analitica[1] else None  # Solución para y'.
                ]
                # Devuelve tiempo, solución numérica (y y su derivada) y soluciones analíticas.
                return t, y[0, :], y[1, :], y_analitica
            
            # Si no hay solución analítica, devuelve tiempo y las soluciones numéricas (y y su derivada).
            return t, y[0, :], y[1, :]
            
        # Caso 3: Resolver un sistema de ecuaciones diferenciales.
        elif edo_tipo == 'sistema':
            # Inicializa una matriz para almacenar las soluciones de las dos variables del sistema.
            y = np.zeros((2, n+1))
            # Asigna las condiciones iniciales al primer punto (t = t0).
            y[:, 0] = condiciones_iniciales
            
            # Itera a través de los pasos para calcular las soluciones numéricas usando RK4.
            for i in range(n):
                y[:, i+1] = self.rk4_paso(f, t[i], y[:, i], h)
            
            # Devuelve el tiempo y las soluciones del sistema (para las dos variables).
            return t, y[0, :], y[1, :]


    def graficar_solucion(self, t, *soluciones, titulo, labels=None, 
                         soluciones_analiticas=None):
        # Crea una nueva figura para los gráficos con un tamaño de 12x6 pulgadas.
        plt.figure(figsize=(12, 6))
        
        # Si no se proporcionan etiquetas para las soluciones, se generan etiquetas predeterminadas.
        if labels is None:
            labels = [f'y{i+1}' for i in range(len(soluciones))]
        
        # Si no se proporcionan soluciones analíticas, inicializa con valores None.
        if soluciones_analiticas is None:
            soluciones_analiticas = [None] * len(soluciones)
            
        # Itera sobre las soluciones y sus respectivas etiquetas.
        for i, (sol, label) in enumerate(zip(soluciones, labels)):
            # Grafica la solución numérica con una línea sólida.
            plt.plot(t, sol, label=f'Numérica {label}', linestyle='-')
            
            # Si hay una solución analítica correspondiente, se grafica con una línea discontinua.
            if soluciones_analiticas[i] is not None:
                plt.plot(t, soluciones_analiticas[i], 
                         label=f'Analítica {label}', 
                         linestyle='--')
        
        # Configura la cuadrícula, la leyenda, las etiquetas de los ejes y el título del gráfico.
        plt.grid(True)
        plt.legend()
        plt.xlabel('t')  # Etiqueta del eje x (tiempo).
        plt.ylabel('y')  # Etiqueta del eje y (solución).
        plt.title(titulo)  # Título del gráfico.
        plt.show()  # Muestra el gráfico.
        
        # Crea otra figura para graficar el error entre las soluciones numéricas y analíticas.
        plt.figure(figsize=(12, 6))
        error_calculado = False  # Bandera para verificar si se calcularon errores.
        
        # Itera sobre las soluciones numéricas y analíticas.
        for i, (sol, label) in enumerate(zip(soluciones, labels)):
            if soluciones_analiticas[i] is not None:
                # Calcula el error absoluto entre la solución numérica y la analítica.
                error = np.abs(sol - soluciones_analiticas[i])
                # Grafica el error para esta solución.
                plt.plot(t, error, label=f'Error {label}')
                error_calculado = True
        
        # Si se calcularon errores, configura el gráfico de error.
        if error_calculado:
            plt.title('Error entre solución numérica y analítica')  # Título del gráfico de error.
            plt.xlabel('t')  # Etiqueta del eje x.
            plt.ylabel('Error absoluto')  # Etiqueta del eje y.
            plt.legend()  # Muestra la leyenda para los errores.
            plt.grid(True)  # Activa la cuadrícula.
            plt.show()  # Muestra el gráfico de error.

def main():
    """Función principal que resuelve tres ecuaciones específicas: 
    una de primer orden, una de segundo orden y un sistema de ecuaciones."""
    
    # Crea una instancia de la clase EDOSolver, que resuelve ecuaciones diferenciales.
    solver = EDOSolver()
    
    # Define el intervalo de tiempo [t0, tf] y el tamaño del paso h.
    t0, tf = 0.0, 10.0
    h = 0.01  # Paso de integración.

    # Define la EDO de primer orden: dy/dx = y - x.
    def edo1(x, y):
        return y - x  # Representa la derivada dy/dx.

    # Define la solución analítica de la EDO: y(t) = t + 1 + e^t.
    def edo1_analitica(x):
        return x + 1 + math.exp(x)  # Solución exacta de la ecuación.

    # Mensajes informativos sobre la ecuación a resolver.
    print("Resolviendo dy/dx = y - x...")
    print("\nEcuación de la solución analítica:")
    print("y(t) = t + 1 + e^t")
    print("Derivada: dy/dt = 1 + e^t")
    print("Condición inicial: y(0) = 1")
    
    # Resuelve la EDO numéricamente utilizando el método RK4.
    t, y, y_analitica = solver.resolver('primer_orden', edo1, [1.0], t0, tf, h, edo1_analitica)
    
    # Grafica la solución numérica y analítica.
    solver.graficar_solucion(t, y, titulo='Solución de dy/dx = y - x', 
                            labels=['y'], 
                            soluciones_analiticas=[y_analitica])

    # Define las ecuaciones para una EDO de segundo orden reescrita como un sistema:
    # y'' + 2y' + y = 0.
    def edo2_sistema(t, y, z):
        return z  # z = dy/dt.

    def edo2_sistema_z(t, y, z):
        return -2*z - y  # Representa y'' en términos de y y z.

    # Solución analítica para la EDO: y(t) = e^(-t).
    def edo2_analitica(t):
        return (1 + 0*t) * math.exp(-t)  # Solución analítica.

    # Mensajes informativos sobre la ecuación a resolver.
    print("\nResolviendo y'' + 2y' + y = 0...")
    print("\nEcuación de la solución analítica:")
    print("y(t) = e^(-t)")
    print("Primera derivada: y'(t) = -e^(-t)")
    print("Segunda derivada: y''(t) = e^(-t)")
    print("Condiciones iniciales: y(0) = 1, y'(0) = 0")
    
    # Resuelve la EDO de segundo orden utilizando el método RK4.
    t, y, dy, y_analitica = solver.resolver('segundo_orden', 
                                           (edo2_sistema, edo2_sistema_z), 
                                           [1.0, 0.0], t0, tf, h, 
                                           [edo2_analitica, None])
    
    # Grafica la solución numérica (y, y') y la solución analítica de y.
    solver.graficar_solucion(t, y, dy, 
                             titulo="Solución de y'' + 2y' + y = 0", 
                             labels=['y', "y'"],
                             soluciones_analiticas=[y_analitica[0], None])

    # Define un sistema de dos ecuaciones diferenciales acopladas:
    # dx/dt = x + y
    def sistema_dx(t, x, y):
        return x + y  # Derivada de x respecto al tiempo.

    # dy/dt = x - y
    def sistema_dy(t, x, y):
        return x - y  # Derivada de y respecto al tiempo.

    # Mensajes informativos sobre el sistema de ecuaciones.
    print("\nResolviendo el sistema dx/dt = x + y, dy/dt = x - y...")
    print("\nNo se tiene solución analítica para este sistema de ecuaciones.")
    
    # Resuelve el sistema de ecuaciones numéricamente utilizando el método RK4.
    t, x, y = solver.resolver('sistema', (sistema_dx, sistema_dy), 
                              [1.0, 1.0], t0, tf, h)
    
    # Grafica las soluciones numéricas para x y y.
    solver.graficar_solucion(t, x, y, 
                             titulo='Solución del sistema de ecuaciones', 
                             labels=['x', 'y'])

# Llama a la función principal para ejecutar las simulaciones.
main()
