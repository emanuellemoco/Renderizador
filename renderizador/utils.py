import math
import numpy as np


# #Funcao sign ideia extraida: https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
# def sign(p1x, p1y, p2x, p2y, p3x, p3y):
#     #Equacao da reta
#     return (p1x - p3x) * (p2y - p3y) - (p2x - p3x) * (p1y - p3y);
    

# def inside(x_v, y_v, x, y):
#     #Funcao que checa se um ponto esta dentro. 
#     # Se o ponto checado em l1, l2 e l3 for positivo, ele esta dentro do triangulo e deve ser pintado
#     l1 = sign(x, y, x_v[0], y_v[0],  x_v[1], y_v[1] )
#     l2 = sign(x, y, x_v[1], y_v[1],  x_v[2], y_v[2] )
#     l3 = sign(x, y, x_v[2], y_v[2],  x_v[0], y_v[0] )

#     has_neg = (l1 < 0) or (l2 < 0) or (l3 < 0);
#     has_pos = (l1 > 0) or (l2 > 0) or (l3 > 0);

#     return not (has_neg and has_pos)

def line_eq(x1,y1,x0,y0,x,y):
    #Equacao da reta
    L = (y1-y0) * x - (x1 - x0) * y + (x1 - x0) * y0 - (y1 - y0) * x0
    return L >= 0
    

def inside(x_v, y_v, x, y, par=True):
    #Funcao que checa se um ponto esta dentro. 
    # Se o ponto checado em l1, l2 e l3 por positivo, ele esta dentro do triangulo e deve ser pintado
    for i in range(0,len(x_v),3):
        if par:
            l1 = line_eq(x_v[i+1],y_v[i+1],x_v[i],y_v[i], x,y)
            l2 = line_eq(x_v[i+2],y_v[i+2],x_v[i+1],y_v[i+1], x,y)
            l3 = line_eq(x_v[i],y_v[i],x_v[i+2],y_v[i+2], x,y)                                     
        if (l1 and l2 and l3):
            return True
    return False
    




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

