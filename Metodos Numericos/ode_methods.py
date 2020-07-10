""" 
    Funciones para cada metodo de resolucion de EDOs 
    Ejemplos de como usarlas en 'Guia 6 - EDO.ipynb'
"""


import numpy as np

""" 
    Euler
    @param to: tiempo inicial
    @param h: paso
    @param n: numero de iteraciones a realizar (implicitamente esto define el tf)
    @param yo: condiciones iniciales, su tipo y tamaño en caso de ser un arreglo definen 
                si se trata de una sola ecuacion o un sistema
    @param f: callback a la funcion de carga, que debe tener el siguiente prototipo:
                f(tn, yn) -> devuelve un elemento del mismo tipo y tamaño que yn (escalar o arreglo en caso de un sistema)
"""
def euler(to, h, n, f, yo):
    if (type(yo) is int) or (type(yo) is float):
        yo = [yo]
    y = np.zeros((n+1, len(yo)))
    y[0] = yo
    t = to
    for k in range(1, n+1):
        y[k] = y[k-1] + h*f(t,y[k-1])
        t = t + h
    return y

""" 
    Heun
    @param to: tiempo inicial
    @param h: paso
    @param n: numero de iteraciones a realizar (implicitamente esto define el tf)
    @param yo: condiciones iniciales, su tipo y tamaño en caso de ser un arreglo definen 
                si se trata de una sola ecuacion o un sistema
    @param f: callback a la funcion de carga, que debe tener el siguiente prototipo:
                f(tn, yn) -> devuelve un elemento del mismo tipo y tamaño que yn (escalar o arreglo en caso de un sistema)
"""
def heun(to, h, n, f, yo):
    if (type(yo) is int) or (type(yo) is float):
        yo = [yo]
    y = np.zeros((n+1, len(yo)))
    y[0] = yo
    K1 = 0
    K2 = 0
    t = to
    for k in range(1, n+1):
        K1 = f(t, y[k-1])
        t = t + h
        K2 = f(t, y[k-1] + h*K1)
        y[k] = y[k-1] + h/2*(K1 + K2)
    return y

""" 
    Taylor
    @param order: orden del metodo requerido
    @param to: tiempo inicial
    @param h: paso
    @param n: numero de iteraciones a realizar (implicitamente esto define el tf)
    @param yo: condiciones iniciales, su tipo y tamaño en caso de ser un arreglo definen 
                si se trata de una sola ecuacion o un sistema
    @param yderiv: callback a las derivadas de y. Debe haber 'order' callbacks, que deben tener el siguiente prototipo:
                yderiv[i](tn, yn) -> devuelve un elemento del mismo tipo y tamaño que yn (escalar o arreglo en caso de un sistema)
"""
def taylor(order, to, h, n, yderiv, yo):
    if (type(yo) is int) or (type(yo) is float):
        yo = [yo]
    y = np.zeros((n+1, len(yo)))
    y[0] = yo
    t = to
    for k in range(1, n+1):
        y[k] = y[k-1]
        for i in range(order):
            y[k] += h**(i+1)*(yderiv[i])(t, y[k-1])/np.math.factorial(i+1)
        t = t + h
    return y
   
"""
    Richardson: calcula el valor de la constante de proporcionalidad C del error tal que err_global <= C.h^order
    @param y1: valor de y(tn) utilizando un paso h1
    @param h1: paso utilizado para calcular y1
    @param y2: valor de y(tn) utilizando un paso h2
    @param h2: paso utilizado para calcular y2
    @param order: orden del metodo
"""
def richardson(y1, h1, y2, h2, order):
    return (y1-y2)/(h1**order - h2**order)

"""
    wanted_step: calcula el paso requerido para un metodo de orden 'order' con constante 'c', un determinado 'error'.
"""
def wanted_step(c, error, order):
    return np.power(error/np.abs(c),1.0/order)

    
""" Funcion simple para mostrar el resultado y opcionalmente el error si el parametros 'real' no es None """
def print_result(method, n, res, real):
    print(f'{method} con {n} pasos')
    print(res)
    if real is not None:
        print(f'Error: {np.abs(res[-1][-1] - real):.5f}')
    print(' ')
    
