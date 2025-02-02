import math
import numpy as np

def line_eq(x1,y1,x0,y0,x,y):
    #Equacao da reta
    L = (y1-y0) * x - (x1 - x0) * y + (x1 - x0) * y0 - (y1 - y0) * x0
    return L >= 0
    

def inside(x_v, y_v, x, y):
    #Funcao que checa se um ponto esta dentro. 
    # Se o ponto checado em l1, l2 e l3 por positivo, ele esta dentro do triangulo e deve ser pintado

    for i in range(0,len(x_v),3):
        l1 = line_eq(x_v[i+1],y_v[i+1],x_v[i],y_v[i], x,y)
        l2 = line_eq(x_v[i+2],y_v[i+2],x_v[i+1],y_v[i+1], x,y)
        l3 = line_eq(x_v[i],y_v[i],x_v[i+2],y_v[i+2], x,y)                                     
        if (l1 and l2 and l3):
            return True
    return False
    

def baricentro(x, y, xA, xB, xC, yA, yB, yC):
    alfa = (- (x - xB)*(yC - yB) + (y - yB)*(xC - xB)) / (- (xA - xB)*(yC - yB) + (yA - yB)*(xC - xB))    
    beta = (- (x - xC)*(yA - yC) + (y - yC)*(xA - xC)) / (- (xB - xC)*(yA - yC) + (yB - yC)*(xA - xC))
    gama = 1 - alfa - beta

    return alfa, beta, gama


def getRotationMatrix(rotation):

    #Matrix de Rotações
    if rotation[0] != 0:
        rotation_matrix = np.matrix([
            [1, 0, 0, 0],
            [0, math.cos(rotation[3]), -math.sin(rotation[3]), 0],
            [0, math.sin(rotation[3]), math.cos(rotation[3]), 0],
            [0, 0, 0, 1]
        ]) 

    elif rotation[1] != 0:
        rotation_matrix = np.matrix([
            [math.cos(rotation[3]), 0, math.sin(rotation[3]), 0],
            [0, 1, 0, 0],
            [-math.sin(rotation[3]), 0, math.cos(rotation[3]), 0],
            [0, 0, 0, 1]
        ])

    else:
        rotation_matrix = np.matrix([
            [math.cos(rotation[3]), -math.sin(rotation[3]), 0, 0],
            [math.sin(rotation[3]), math.cos(rotation[3]), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ])
    return rotation_matrix

def calc_T(a, b):
    T_minus_1_sub = np.subtract(a,b)
    T_minus_1 = [x / 2 for x in T_minus_1_sub]
    return T_minus_1