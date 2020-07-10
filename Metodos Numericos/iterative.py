""" Implementacion de metodos iterativos para resolucion de ecuaciones no lineales """

import numpy as np

""" 
    bisec: Metodo de la biseccion en [a, b] con tolerancia 'tol', recibe callback a la funcion para evaluarla
"""
def bisec(a, b, tol, f):
    a_new = a
    b_new = b
    n = int(np.ceil(np.log2((b_new - a_new)/tol)) - 1)
    print(f'Se precisan {n} iteraciones')
    for i in range(n):
        c = a_new + (b_new - a_new)/2
        if (f(c)*f(a_new)) >= 0:
            a_new = c
        else:
            b_new = c
        print(f'Raiz en iteracion {i + 1}: {(a_new + b_new)/2}')
    return (a_new + b_new) / 2

""" 
    fixed_point: Metodo de punto fijo arranca en xo con tolerancia 'tol', recibe callback a la funcion para evaluarla
"""
def fixed_point(xo, tol, f):
    x0 = 0
    x1 = xo
    i = 0
    while np.abs(x1 - x0) > tol or not i:
        x0 = x1
        x1 = f(x0)
        i += 1
        print(f'Raiz en iteracion {i}: x{i} = {x1}')
    return x1

""" 
    newton_raphson: Metodo de Newton en [a, b] con tolerancia 'tol', recibe callback a la funcion y su derivada para evaluarlas
"""
def newton_raphson(a, b, tol, f, fprime):
    x0 = (a + b) / 2 + 1000
    x1 = (a + b) / 2
    i = 0
    while np.abs(x1 - x0) > tol:
        x0 = x1
        x1 = x0 - f(x0)/fprime(x0)
        i += 1
        print(f'Raiz en iteracion {i}: x{i} = {x1}')
    return x1

