#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
Biblioteca Gráfica / Graphics Library.

Desenvolvido por: Emanuelle e Leonardo 
Disciplina: Computação Gráfica
Data:
"""
import numpy as np
import gpu          # Simula os recursos de uma GPU
import math
from utils import ( 
    inside,
    getRotationMatrix,
    baricentro
)
class GL:
    """Classe que representa a biblioteca gráfica (Graphics Library)."""

    width = 800   # largura da tela
    height = 600  # altura da tela
    near = 0.01   # plano de corte próximo
    far = 1000    # plano de corte distante

    @staticmethod
    def setup(width, height, near=0.01, far=1000):
        """Definr parametros para câmera de razão de aspecto, plano próximo e distante."""
        GL.width = width
        GL.height = height
        GL.near = near
        GL.far = far
        GL.screen_matrix = np.matrix([
            [GL.width/2, 0, 0, GL.width/2],
            [0, -GL.height/2, 0, GL.height/2],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ]) 

    @staticmethod
    def get2DCoord(point):
        #Para cada triangulo
        points_2d = []
        for i in range(0, len(point)-2, +3):
 
            point3d = np.matrix([[point[i]], [point[i+1]], [point[i+2]], [1]])
            
            world_point = GL.world.dot(point3d)
            #Identificando as coordenadas ortonormais (matrix Look-At)
            # print("world_point: ",world_point)

            #matriz camera translação
            camera_translation = np.matrix([
                [1, 0, 0, GL.position[0,0]],
                [0, 1, 0, GL.position[0,1]],
                [0, 0, 1, GL.position[0,2]],
                [0, 0, 0, 1]
            ])

            camera_rotation = getRotationMatrix(GL.orientation)

            #Transformação de Look-at           (T.R)⁻1 = R⁻1 . T⁼-1
            x_matrix = np.linalg.inv(camera_rotation).dot(np.linalg.inv(camera_translation))
            #Coordenadas do Mundo -> Coordenadas Camera
            camera_point = x_matrix.dot(world_point)

            #Transformações Perspectivas
            perspective = GL.p_matrix.dot(camera_point)

            #Divisão homogenea para fazer a normalização:
            perspective_normalized = perspective/perspective[3][0]
            # perspective_normalized = np.divide(perspective, perspective[3][0])


            #Tensformação para coordenadas da tela
            screen_point = GL.screen_matrix.dot(perspective_normalized)
            # screen_point = screen_point/2 #somente para teste
            #gpu.GPU.draw_pixels([int(screen_point[0,0]), int(screen_point[0,1])], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel
            points_2d.append(screen_point[0,0])
            points_2d.append(screen_point[1,0])
        return points_2d

    @staticmethod
    def triangleSet(point, colors):
        """Função usada para renderizar TriangleSet."""
        # Nessa função você receberá pontos no parâmetro point, esses pontos são uma lista
        # de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x do
        # primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
        # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e
        # assim por diante.
        # No TriangleSet os triângulos são informados individualmente, assim os três
        # primeiros pontos definem um triângulo, os três próximos pontos definem um novo
        # triângulo, e assim por diante.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o TriangleSet
        # você pode assumir o desenho das linhas com a cor emissiva (emissiveColor).

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        # print("TriangleSet : pontos = {0}".format(point)) # imprime no terminal pontos
        # print("TriangleSet : colors = {0}".format(colors)) # imprime no terminal as cores

        # Exemplo de desenho de um pixel branco na coordenada 10, 10
        # gpu.GPU.draw_pixels([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel
        triangle_2d_coord = GL.get2DCoord(point)

        GL.fill_triangle(triangle_2d_coord, colors)
    
        
    @staticmethod
    def fill_triangle(vertices, colors, reverse=False, vertexColor=False):
        """Função para rasterizar os pixels dentro do triângulo """
        # print("TriangleSet2D : vertices = {0}".format(vertices)) # imprime no terminal
        # print("TriangleSet2D : colors = {0}".format(colors)) # imprime no terminal as cores
        # # Exemplo:
        # gpu.GPU.set_pixel(24, 8, 255, 255, 0) # altera um pixel da imagem (u, v, r, g, b)

        x_list = vertices[::2]
        y_list = vertices[1::2]
        if reverse:
            x_list = x_list[::-1]
            y_list = y_list[::-1]

        print("colors : ", colors)
        # print("VERTICES fill_triangle - " ,vertices)
        for j in range(int(min(x_list)), int(max(x_list))):
            for i in range(int(min(y_list)), int(max(y_list))): 
                if inside(x_list, y_list, j,i):
                    #caso tenha vertexColor, tem que calcular o baricentro e usar os valores de alfa, beta e gama para multiplicar pela cor
                    if vertexColor:
                        alfa, beta, gama = baricentro(j, i, x_list[0], x_list[1], x_list[2], y_list[0], y_list[1], y_list[2] )
                        gpu.GPU.draw_pixels([j,i], gpu.GPU.RGB8,  [colors[0][0]*alfa*255+colors[1][0]*beta*255+colors[2][0]*gama*255, colors[0][1]*alfa*255+colors[1][1]*beta*255+colors[2][1]*gama*255, colors[0][2]*alfa*255+colors[1][2]*beta*255+colors[2][2]*gama*255])

                    else:
                    #print("i, j: ", i, j)
                      gpu.GPU.draw_pixels([j,i], gpu.GPU.RGB8,  [colors["diffuseColor"][0]*255, colors["diffuseColor"][1]*255, colors["diffuseColor"][2]*255])




    @staticmethod
    def viewpoint(position, orientation, fieldOfView):
        """Função usada para renderizar (na verdade coletar os dados) de Viewpoint."""
        # Na função de viewpoint você receberá a posição, e campo de visão da
        # câmera virtual. Use esses dados para poder calcular e criar a matriz de projeção
        # perspectiva para poder aplicar nos pontos dos objetos geométricos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        # print("Viewpoint : ", end='')
        # print("position = {0} ".format(position), end='')
        # print("orientation = {0} ".format(orientation), end='')
        # print("fieldOfView = {0} ".format(fieldOfView))
        # ####### remover 4 linhas a cima 

        fovy = 2 * math.atan(math.tan(fieldOfView / 2) * GL.height / math.sqrt(GL.height ** 2 + GL.width ** 2))

        GL.position = np.matrix([
            [position[0], position[1], position[2], 1]
        ]) 
        GL.orientation = orientation
        
        top = GL.near * math.tan(fovy) 
        aspect = GL.width / GL.height
        right  = top * aspect
        #matrix de perspectiva
        GL.p_matrix = np.matrix([
            [GL.near/right, 0, 0, 0],
            [0, GL.near/top, 0,0], 
            [0,0, -(GL.far+GL.near)/(GL.far-GL.near), (-2 * GL.far * GL.near)/(GL.far-GL.near)], 
            [0, 0, -1, 0]
            ])


    @staticmethod
    def transform_in(translation, scale, rotation):
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_in será chamada quando se entrar em um nó X3D do tipo Transform
        # do grafo de cena. Os valores passados são a escala em um vetor [x, y, z]
        # indicando a escala em cada direção, a translação [x, y, z] nas respectivas
        # coordenadas e finalmente a rotação por [x, y, z, t] sendo definida pela rotação
        # do objeto ao redor do eixo x, y, z por t radianos, seguindo a regra da mão direita.
        # Quando se entrar em um nó transform se deverá salvar a matriz de transformação dos
        # modelos do mundo em alguma estrutura de pilha.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        # print("Transform : ", end='')
        # if translation:
        #     print("translation = {0} ".format(translation), end='') # imprime no terminal
        # if scale:
        #     print("scale = {0} ".format(scale), end='') # imprime no terminal
        # if rotation:
        #     print("rotation = {0} ".format(rotation), end='') # imprime no terminal
        # print("")
    
        # O Braga nos explicou como fazer essa parte
        scale_matrix = np.matrix([
            [scale[0], 0, 0 ,0],
            [0, scale[1], 0, 0],
            [0, 0, scale[2], 0],
            [0, 0, 0, 1]
        ])
        
        translation_matrix = np.matrix([
            [1, 0, 0, translation[0]],
            [0, 1, 0, translation[1]],
            [0, 0, 1, translation[2]],
            [0, 0, 0, 1]
        ])
        
        rotation_matrix = getRotationMatrix(rotation)
        GL.world = translation_matrix.dot(rotation_matrix).dot(scale_matrix)

        print("ENTREI AQUI")

        

        # GL.world = scale.dot(translation)  # PRECISA DISSO? ??? SE NÃO TIVER ROTAÇÃO E MULTIPLICAR POR ELA MESMA DÁ ERRO? BRAGA HElp


    @staticmethod
    def transform_out():
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_out será chamada quando se sair em um nó X3D do tipo Transform do
        # grafo de cena. Não são passados valores, porém quando se sai de um nó transform se
        # deverá recuperar a matriz de transformação dos modelos do mundo da estrutura de
        # pilha implementada.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        # print("Saindo de Transform")

    @staticmethod
    def triangleStripSet(point, stripCount, colors):
        """Função usada para renderizar TriangleStripSet."""
        # A função triangleStripSet é usada para desenhar tiras de triângulos interconectados,
        # você receberá as coordenadas dos pontos stripCount parâmetro point, esses pontos são uma
        # lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x
        # do primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
        # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e assim
        # por diante. No TriangleStripSet a quantidade de vértices a serem usados é informado
        # em uma lista chamada stripCount (perceba que é uma lista). Ligue os vértices na ordem,
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante. Cuidado com a orientação dos vértices, ou seja,
        # todos no sentido horário ou todos no sentido anti-horário, conforme especificado.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        # print("TriangleStripSet : pontos = {0} ".format(point), end='')
        # for i, strip in enumerate(stripCount):
        #     print("TriangleStripSet : colors = {0}".format(colors)) # imprime no terminal as cores
        #     print("strip[{0}] = {1} ".format(i, strip), end='')
        # print("")
        
        points_2d = GL.get2DCoord(point)
        orientation = 1
        for i in range(0,stripCount[0]-1):
            lista_vertice = []
            for j in range(i, i+3):
                lista_vertice.append(points_2d[2*(j-1)])
                lista_vertice.append(points_2d[2*(j-1)+1])
            if orientation % 2 ==0:
                GL.fill_triangle(lista_vertice, colors)
            else:
                GL.fill_triangle(lista_vertice, colors, reverse=True)
            orientation +=1
    @staticmethod
    def indexedTriangleStripSet(point, index, colors):
        """Função usada para renderizar IndexedTriangleStripSet."""
        # A função indexedTriangleStripSet é usada para desenhar tiras de triângulos
        # interconectados, você receberá as coordenadas dos pontos no parâmetro point, esses
        # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor
        # da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto, point[2]
        # o valor z da coordenada z do primeiro ponto. Já point[3] é a coordenada x do
        # segundo ponto e assim por diante. No IndexedTriangleStripSet uma lista informando
        # como conectar os vértices é informada em index, o valor -1 indica que a lista
        # acabou. A ordem de conexão será de 3 em 3 pulando um índice. Por exemplo: o
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante. Cuidado com a orientação dos vértices, ou seja,
        # todos no sentido horário ou todos no sentido anti-horário, conforme especificado.

        # # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        # print("IndexedTriangleStripSet : pontos = {0}, index = {1}".format(point, index))
        # print("IndexedTriangleStripSet : colors = {0}".format(colors)) # imprime as cores

        
        points_2d = GL.get2DCoord(point)
        orientation = 2
        for i in range(len(index)-3):
            lista_vertice = []
            for j in range(index[i], index[i]+3):
                lista_vertice.append(points_2d[j*2])
                lista_vertice.append(points_2d[j*2+1])
            if orientation % 2 ==0:
                GL.fill_triangle(lista_vertice, colors)
            else:
                GL.fill_triangle(lista_vertice, colors, reverse=True)
            orientation +=1 



    @staticmethod
    def box(size, colors):
        """Função usada para renderizar Boxes."""
        # A função box é usada para desenhar paralelepípedos na cena. O Box é centrada no
        # (0, 0, 0) no sistema de coordenadas local e alinhado com os eixos de coordenadas
        # locais. O argumento size especifica as extensões da caixa ao longo dos eixos X, Y
        # e Z, respectivamente, e cada valor do tamanho deve ser maior que zero. Para desenha
        # essa caixa você vai provavelmente querer tesselar ela em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        # print("Box : size = {0}".format(size)) # imprime no terminal pontos
        # print("Box : colors = {0}".format(colors)) # imprime no terminal as cores

        points_3d = []

        p1 = [-size[0], -size[1],  size[2]] # P1 = (-1, -1, 1)
        p2 = [-size[0],  size[1],  size[2]] # P2 = (-1, 1, 1)
        p3 = [size[0], -size[1],  size[2]]  # P3 = (1, -1, 1)
        p4 = [-size[0],  size[1], -size[2]] # P4 = (-1, 1, -1)
        p5 = [size[0],  size[1],  size[2]]  # P5 = (1, 1, 1)
        p6 = [size[0], -size[1], -size[2]]  # P6 = (1, -1, -1)
        p7 = [-size[0], -size[1], -size[2]] # P7 = (-1, -1, -1)
        p8 = [size[0],  size[1], -size[2]]  # P8 = (1, 1, -1)
        points_3d = p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8

        points_2d = GL.get2DCoord(points_3d)
    
    # quadrado 1 7 4 2 
        # 1 7 4
        vertice = points_2d[0:2] + points_2d[12:14] + points_2d[6:8]
        GL.fill_triangle(vertice, colors)

        # 4 2 1 
        vertice = points_2d[6:8] + points_2d[2:4] + points_2d[0:2]
        GL.fill_triangle(vertice, colors)

    # quadrado 3 6 8 5 
        # 3 6 8
        vertice = points_2d[4:6] + points_2d[10:12] + points_2d[14:16]
        GL.fill_triangle(vertice, colors)

        # 8 5 3 
        vertice = points_2d[14:16] + points_2d[8:10] + points_2d[4:6]
        GL.fill_triangle(vertice, colors)

    # ligacoes
        # 1 6 7 
        vertice = points_2d[0:2] + points_2d[10:12] + points_2d[12:14]
        GL.fill_triangle(vertice, colors)

        # 2 4 8 
        vertice = points_2d[2:4] + points_2d[6:8] + points_2d[14:16]
        GL.fill_triangle(vertice, colors)



    @staticmethod
    def indexedFaceSet(coord, coordIndex, colorPerVertex, color, colorIndex,
                       texCoord, texCoordIndex, colors, current_texture):
        """Função usada para renderizar IndexedFaceSet."""
        # A função indexedFaceSet é usada para desenhar malhas de triângulos. Ela funciona de
        # forma muito simular a IndexedTriangleStripSet porém com mais recursos.
        # Você receberá as coordenadas dos pontos no parâmetro cord, esses
        # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim coord[0] é o valor
        # da coordenada x do primeiro ponto, coord[1] o valor y do primeiro ponto, coord[2]
        # o valor z da coordenada z do primeiro ponto. Já coord[3] é a coordenada x do
        # segundo ponto e assim por diante. No IndexedFaceSet uma lista de vértices é informada
        # em coordIndex, o valor -1 indica que a lista acabou.
        # A ordem de conexão será de 3 em 3 pulando um índice. Por exemplo: o
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante.
        # Adicionalmente essa implementação do IndexedFace aceita cores por vértices, assim
        # se a flag colorPerVertex estiver habilitada, os vértices também possuirão cores
        # que servem para definir a cor interna dos poligonos, para isso faça um cálculo
        # baricêntrico de que cor deverá ter aquela posição. Da mesma forma se pode definir uma
        # textura para o poligono, para isso, use as coordenadas de textura e depois aplique a
        # cor da textura conforme a posição do mapeamento. Dentro da classe GPU já está
        # implementadado um método para a leitura de imagens.

        # Os prints abaixo são só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("IndexedFaceSet : ")
        if coord:
            print("\tpontos(x, y, z) = {0}, coordIndex = {1}".format(coord, coordIndex))
        print("colorPerVertex = {0}".format(colorPerVertex))
        if colorPerVertex and color and colorIndex:
            print("\tcores(r, g, b) = {0}, colorIndex = {1}".format(color, colorIndex))
        if texCoord and texCoordIndex:
            print("\tpontos(u, v) = {0}, texCoordIndex = {1}".format(texCoord, texCoordIndex))
        if current_texture:
            image = gpu.GPU.load_texture(current_texture[0])
            print("\t Matriz com image = {0}".format(image))
            print("\t Dimensões da image = {0}".format(image.shape))
        print("IndexedFaceSet : colors = {0}".format(colors))  # imprime no terminal as cores
        
        
        points_2d = GL.get2DCoord(coord)
        
        #criando a lista de vertice (x, y) com sua respectiva cor para passar para funcao fill_triangle
        i = 0
        while i < (len(coordIndex)):
            lista_vertice = []
            colors_list = []
            while coordIndex[i] != -1:
                lista_vertice.append(points_2d[2 * (coordIndex[i])])
                lista_vertice.append(points_2d[2 * (coordIndex[i]) + 1])
                colors_list.append(color[3*colorIndex[i]:3*colorIndex[i]+3])
                i+=1
            # print("lista_vertice: ", lista_vertice)
            # print("colors_list: ", colors_list)
            GL.fill_triangle(lista_vertice, colors_list, vertexColor=True)
            i+=1


# [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0], 
#  -------------   -----------     -----------   ------------
# colorIndex = [0, 1, 2, -1, 2, 3, 0, -1]        
        
        # Exemplo de desenho de um pixel branco na coordenada 10, 10
        gpu.GPU.draw_pixels([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel
        

        


    # Para o futuro (Não para versão atual do projeto.)
    def vertex_shader(self, shader):
        """Para no futuro implementar um vertex shader."""

    def fragment_shader(self, shader):
        """Para no futuro implementar um fragment shader."""
