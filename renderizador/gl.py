#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# pylint: disable=invalid-name

"""
Biblioteca Gráfica / Graphics Library.

Desenvolvido por: Emanuelle e Leonardo 
Disciplina: Computação Gráfica
Data:
"""
import numpy as np
import gpu  # Simula os recursos de uma GPU
import math
from utils import inside, getRotationMatrix, baricentro

import time         # Para operações com tempo

import gpu          # Simula os recursos de uma GPU

class GL:
    """Classe que representa a biblioteca gráfica (Graphics Library)."""

    width = 800  # largura da tela
    height = 600  # altura da tela
    near = 0.01  # plano de corte próximo
    far = 1000  # plano de corte distante

    @staticmethod
    def setup(width, height, near=0.01, far=1000):
        """Definr parametros para câmera de razão de aspecto, plano próximo e distante."""
        GL.width = width
        GL.height = height
        GL.near = near
        GL.far = far
        GL.screen_matrix = np.matrix(
            [
                [GL.width / 2, 0, 0, GL.width / 2],
                [0, -GL.height / 2, 0, GL.height / 2],
                [0, 0, 1, 0],
                [0, 0, 0, 1],
            ]
        )
        GL.matrix_list = []
        GL.zbuffer = np.matrix(np.ones((GL.height, GL.width)) * np.inf)

        GL.light = False

    @staticmethod
    def get2DCoord(point):
        
        # Para cada triangulo
        points_2d = []
        z_list = []
        for i in range(0, len(point) - 2, +3):

            point3d = np.matrix([[point[i]], [point[i + 1]], [point[i + 2]], [1]])

            world_point = GL.world.dot(point3d)
            # Identificando as coordenadas ortonormais (matrix Look-At)
            # print("world_point: ",world_point)

            # matriz camera translação
            camera_translation = np.matrix(
                [
                    [1, 0, 0, GL.position[0, 0]],
                    [0, 1, 0, GL.position[0, 1]],
                    [0, 0, 1, GL.position[0, 2]],
                    [0, 0, 0, 1],
                ]
            )

            camera_rotation = getRotationMatrix(GL.orientation)

            # Transformação de Look-at           (T.R)⁻1 = R⁻1 . T⁼-1
            x_matrix = np.linalg.inv(camera_rotation).dot(
                np.linalg.inv(camera_translation)
            )
            # Coordenadas do Mundo -> Coordenadas Camera
            camera_point = x_matrix.dot(world_point)

            # Transformações Perspectivas
            perspective = GL.p_matrix.dot(camera_point)

            # Divisão homogenea para fazer a normalização:
            perspective_normalized = perspective / perspective[3][0]
            # perspective_normalized = np.divide(perspective, perspective[3][0])

            z_list.append(perspective[2][0])

            # Tensformação para coordenadas da tela
            screen_point = GL.screen_matrix.dot(perspective_normalized)
            # screen_point = screen_point/2 #somente para teste
            # gpu.GPU.draw_pixels([int(screen_point[0,0]), int(screen_point[0,1])], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel
            points_2d.append(screen_point[0, 0])
            points_2d.append(screen_point[1, 0])

        return points_2d, z_list

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
        triangle_2d_coord, z_list = GL.get2DCoord(point)

        GL.fill_triangle(triangle_2d_coord, colors, z_list)

    @staticmethod
    def fill_triangle(vertices, colors, z, texture=None, uv=None, reverse=False, vertexColor=False, hasTexture=False):
        """Função para rasterizar os pixels dentro do triângulo"""
        # print("TriangleSet2D : vertices = {0}".format(vertices)) # imprime no terminal
        # print("TriangleSet2D : colors = {0}".format(colors)) # imprime no terminal as cores
        # # Exemplo:
        # gpu.GPU.set_pixel(24, 8, 255, 255, 0) # altera um pixel da imagem (u, v, r, g, b)
        x_list = vertices[::2]
        y_list = vertices[1::2]
        if reverse:
            x_list = x_list[::-1]
            y_list = y_list[::-1]
            z = z[::-1]

        if hasTexture:
            u_list = uv[::2]
            v_list = uv[1::2] 

        for j in range(int(np.floor(min(x_list))), int(np.ceil(max(x_list)))):
            for i in range(int(np.floor(min(y_list))), int(np.ceil(max(y_list)))): ### KD O BRAGA KD O BRAGA KD O BRAGA KD O BRAGA KD O BRAGA KD O BRAGA KD O BRAGA 
                if inside(x_list, y_list, j, i):
                    # caso tenha vertexColor, tem que calcular o baricentro e usar os valores de alfa, beta e gama para multiplicar pela cor
                    alfa, beta, gama = baricentro(j, i, x_list[0], x_list[1], x_list[2], y_list[0], y_list[1], y_list[2],)
                    
                    current_z = 1 / (alfa * (1/z[0]) + beta * (1/z[1]) + gama * (1/z[2]))
                    if current_z < GL.zbuffer[i, j]:
                        GL.zbuffer[i, j] = current_z
                        if vertexColor:
                            gpu.GPU.draw_pixels(
                                [j, i],
                                gpu.GPU.RGB8,
                                [
                                    # colors[0][0] * alfa * 255 + colors[1][0] * beta * 255 + colors[2][0] * gama * 255,
                                    # colors[0][1] * alfa * 255 + colors[1][1] * beta * 255 + colors[2][1] * gama * 255,
                                    # colors[0][2] * alfa * 255 + colors[1][2] * beta * 255 + colors[2][2] * gama * 255,
                                    current_z * (colors[0][0] * alfa * 255 / z[0] + colors[1][0] * beta * 255 / z[1] + colors[2][0] * gama * 255/ z[2]),
                                    current_z * (colors[0][1] * alfa * 255 / z[0] + colors[1][1] * beta * 255 / z[1] + colors[2][1] * gama * 255/ z[2]),
                                    current_z * (colors[0][2] * alfa * 255 / z[0] + colors[1][2] * beta * 255 / z[1] + colors[2][2] * gama * 255/ z[2]),
                                ],
                            )
                            # print("R: {0} -> {1} ".format(colors[0][0] * alfa * 255 + colors[1][0] * beta * 255 + colors[2][0] * gama * 255, current_z * (colors[0][0] * alfa * 255 / z[0] + colors[1][0] * beta * 255 / z[1] + colors[2][0] * gama * 255/ z[2]) ) )
                            # print("G: {0} -> {1} ".format(colors[0][1] * alfa * 255 + colors[1][1] * beta * 255 + colors[2][1] * gama * 255, current_z * (colors[0][1] * alfa * 255 / z[0] + colors[1][1] * beta * 255 / z[1] + colors[2][1] * gama * 255/ z[2]) ))
                            # print("B: {0} -> {1} ".format(colors[0][2] * alfa * 255 + colors[1][2] * beta * 255 + colors[2][2] * gama * 255, current_z * (colors[0][2] * alfa * 255 / z[0] + colors[1][2] * beta * 255 / z[1] + colors[2][2] * gama * 255/ z[2])))
                            # print("-" *10)

                        elif hasTexture:
                            # print(texture.shape)
                            # print("->>> ",texture)
                            x_tex = int((alfa * u_list[0] + beta * u_list[1] + gama * u_list[2]) * texture.shape[1])
                            y_tex = int((alfa * v_list[0] + beta * v_list[1] + gama * v_list[2]) * (-texture.shape[0]))

                            gpu.GPU.draw_pixels(
                                [j, i],
                                gpu.GPU.RGB8,
                                [
                                    texture[y_tex][x_tex][0],
                                    texture[y_tex][x_tex][1],
                                    texture[y_tex][x_tex][2],
                                ],
                            )
                        else:
                            # print("i, j: ", i, j)
                            gpu.GPU.draw_pixels(
                                [j, i],
                                gpu.GPU.RGB8,
                                [
                                    colors["diffuseColor"][0] * 255,
                                    colors["diffuseColor"][1] * 255,
                                    colors["diffuseColor"][2] * 255,
                                ],
                            )

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

        fovy = 2 * math.atan(
            math.tan(fieldOfView / 2)
            * GL.height
            / math.sqrt(GL.height ** 2 + GL.width ** 2)
        )

        GL.position = np.matrix([[position[0], position[1], position[2], 1]])
        GL.orientation = orientation

        top = GL.near * math.tan(fovy)
        aspect = GL.width / GL.height
        right = top * aspect
        # matrix de perspectiva
        GL.p_matrix = np.matrix(
            [
                [GL.near / right, 0, 0, 0],
                [0, GL.near / top, 0, 0],
                [
                    0,
                    0,
                    -(GL.far + GL.near) / (GL.far - GL.near),
                    (-2 * GL.far * GL.near) / (GL.far - GL.near),
                ],
                [0, 0, -1, 0],
            ]
        )

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
        scale_matrix = np.matrix(
            [
                [scale[0], 0, 0, 0],
                [0, scale[1], 0, 0],
                [0, 0, scale[2], 0],
                [0, 0, 0, 1],
            ]
        )

        translation_matrix = np.matrix(
            [
                [1, 0, 0, translation[0]],
                [0, 1, 0, translation[1]],
                [0, 0, 1, translation[2]],
                [0, 0, 0, 1],
            ]
        )

        rotation_matrix = getRotationMatrix(rotation)
        if len(GL.matrix_list) > 0:
            GL.world = (GL.matrix_list[-1]).dot(
                translation_matrix.dot(rotation_matrix).dot(scale_matrix)
            )
        else:
            GL.world = translation_matrix.dot(rotation_matrix).dot(scale_matrix)
        GL.matrix_list.append(GL.world)

        # GL.world = scale.dot(translation)  # PRECISA DISSO? ??? SE NÃO TIVER ROTAÇÃO E MULTIPLICAR POR ELA MESMA DÁ ERRO? BRAGA HElp

    @staticmethod
    def transform_out():
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_out será chamada quando se sair em um nó X3D do tipo Transform do
        # grafo de cena. Não são passados valores, porém quando se sai de um nó transform se
        # deverá recuperar a matriz de transformação dos modelos do mundo da estrutura de
        # pilha implementada.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        GL.matrix_list.pop()
        if len(GL.matrix_list) > 0:
            GL.world = GL.matrix_list[-1]
        else:
            GL.world = None

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

        points_2d, z_list = GL.get2DCoord(point)
        orientation = 1
        for i in range(0, stripCount[0] - 1):
            lista_vertice = []
            for j in range(i, i + 3):
                lista_vertice.append(points_2d[2 * (j - 1)])
                lista_vertice.append(points_2d[2 * (j - 1) + 1])
            if orientation % 2 == 0:
                GL.fill_triangle(lista_vertice, colors, z_list)
            else:
                GL.fill_triangle(lista_vertice, colors, z_list, reverse=True)
            orientation += 1

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

        points_2d, z_list = GL.get2DCoord(point)
        orientation = 2
        for i in range(len(index) - 3):
            lista_vertice = []
            for j in range(index[i], index[i] + 3):
                lista_vertice.append(points_2d[j * 2])
                lista_vertice.append(points_2d[j * 2 + 1])
            if orientation % 2 == 0:
                GL.fill_triangle(lista_vertice, colors, z_list)
            else:
                GL.fill_triangle(lista_vertice, colors, z_list, reverse=True)
            orientation += 1

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

        p1 = [-size[0], -size[1], size[2]]  # P1 = (-1, -1, 1)
        p2 = [-size[0], size[1], size[2]]  # P2 = (-1, 1, 1)
        p3 = [size[0], -size[1], size[2]]  # P3 = (1, -1, 1)
        p4 = [-size[0], size[1], -size[2]]  # P4 = (-1, 1, -1)
        p5 = [size[0], size[1], size[2]]  # P5 = (1, 1, 1)
        p6 = [size[0], -size[1], -size[2]]  # P6 = (1, -1, -1)
        p7 = [-size[0], -size[1], -size[2]]  # P7 = (-1, -1, -1)
        p8 = [size[0], size[1], -size[2]]  # P8 = (1, 1, -1)
        points_3d = p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8

        points_2d, z_list = GL.get2DCoord(points_3d)

        # quadrado 1 7 4 2
        # 1 7 4
        vertice = points_2d[0:2] + points_2d[12:14] + points_2d[6:8]
        GL.fill_triangle(vertice, colors, z_list)

        # 4 2 1
        vertice = points_2d[6:8] + points_2d[2:4] + points_2d[0:2]
        GL.fill_triangle(vertice, colors, z_list)

        # quadrado 3 6 8 5
        # 3 6 8
        vertice = points_2d[4:6] + points_2d[10:12] + points_2d[14:16]
        GL.fill_triangle(vertice, colors, z_list)

        # 8 5 3
        vertice = points_2d[14:16] + points_2d[8:10] + points_2d[4:6]
        GL.fill_triangle(vertice, colors, z_list)

        # ligacoes
        # 1 6 7
        vertice = points_2d[0:2] + points_2d[10:12] + points_2d[12:14]
        GL.fill_triangle(vertice, colors, z_list)

        # 2 4 8
        vertice = points_2d[2:4] + points_2d[6:8] + points_2d[14:16]
        GL.fill_triangle(vertice, colors, z_list)

    @staticmethod
    def indexedFaceSet(
        coord,
        coordIndex,
        colorPerVertex,
        color,
        colorIndex,
        texCoord,
        texCoordIndex,
        colors,
        current_texture,
    ):
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

        # # Os prints abaixo são só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        # print("IndexedFaceSet : ")
        # if coord:
        #     print("\tpontos(x, y, z) = {0}, coordIndex = {1}".format(coord, coordIndex))
        # print("colorPerVertex = {0}".format(colorPerVertex))
        # if colorPerVertex and color and colorIndex:
        #     print("\tcores(r, g, b) = {0}, colorIndex = {1}".format(color, colorIndex))
        # if texCoord and texCoordIndex:
        #     print("\tpontos(u, v) = {0}, texCoordIndex = {1}".format(texCoord, texCoordIndex))
        # if current_texture:
        #     image = gpu.GPU.load_texture(current_texture[0])
        #     print("\t Matriz com image = {0}".format(image))
        #     print("\t Dimensões da image = {0}".format(image.shape))
        # print("IndexedFaceSet : colors = {0}".format(colors))  # imprime no terminal as cores

        points_2d, z_list = GL.get2DCoord(coord)

        # criando a lista de vertice (x, y) com sua respectiva cor para passar para funcao fill_triangle
        vertex_color = True
        texture = False
        colors_list = []
        image_texture = None
        uv = []

        if color is None:
            vertex_color = False
            colors_list = colors
        if texCoord is not None:
            vertex_color = False
            texture = True
            image_texture = gpu.GPU.load_texture(current_texture[0])

        i = 0
        while i < (len(coordIndex)):
            lista_vertice = []
            if texture:
                uv  = []
            if vertex_color:
                colors_list = []
            while coordIndex[i] != -1:
                lista_vertice.append(points_2d[2 * (coordIndex[i])])
                lista_vertice.append(points_2d[2 * (coordIndex[i]) + 1])
                if vertex_color:
                    colors_list.append(color[3 * colorIndex[i] : 3 * colorIndex[i] + 3])
                if texture:
                    uv.append(texCoord[2*texCoordIndex[i]])
                    uv.append(texCoord[2*texCoordIndex[i]+1])

                i += 1

            GL.fill_triangle(lista_vertice, colors_list, z_list,texture=image_texture, uv=uv, vertexColor=vertex_color, hasTexture=texture)
            i += 1

    # [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0],
    #  -------------   -----------     -----------   ------------
    # colorIndex = [0, 1, 2, -1, 2, 3, 0, -1]

    @staticmethod
    def sphere(radius, colors):
        """Função usada para renderizar Esferas."""
        # A função sphere é usada para desenhar esferas na cena. O esfera é centrada no
        # (0, 0, 0) no sistema de coordenadas local. O argumento radius especifica o
        # raio da esfera que está sendo criada. Para desenha essa esfera você vai
        # precisar tesselar ela em triângulos, para isso encontre os vértices e defina
        # os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        # print("Sphere : radius = {0}".format(radius)) # imprime no terminal o raio da esfera
        print("Sphere : colors = {0}".format(colors)) # imprime no terminal as cores

        #dividir em 5 
        #desenhar por tiras 
        #no topo define o ponto central z
        #phi - muda de curva
        #theta - gira em volta de uma curva
        # de 0 a 2pi 

        if GL.light:
            normal = [0, 0, 1]
            NL = normal * GL.L  
            ambient_i  = GL.Iia * colors["diffuseColor"] * colors["shiness"]
            diffuse_i  = GL.Ii * colors["diffuseColor"] * NL
            specular_i = GL.Ii * colors["specularColor"] * (N * ((GL.L + GL.v)/np.abs(GL.L + GL.v)))**(colors["shiness"]*128)
            
        
        pontos = 12 # número de pontos 
        aneis = 12 #número de aneis
        coord_3d = []

        angulo_ponto = 2*math.pi/pontos
        angulo_anel = 2*math.pi/aneis

        for i in range(aneis):
            phi = math.pi / 2 - i * angulo_anel
            for j in range(pontos):#np.arange(0,360,(360/pontos)): 
                # print("THETA: {0} PHI: {0} ".format(theta,phi))
                theta = j * angulo_ponto
                x = radius * math.sin(phi) * math.cos(theta)
                y = radius * math.sin(phi) * math.sin(theta)
                z = radius * math.cos(phi)
                coord_3d.append(x)
                coord_3d.append(y)
                coord_3d.append(z)
                # print(f"aaaaaaaaa: ({x}, {y}, {z})")

        

        coord_2d, z_list = GL.get2DCoord(coord_3d)

        # pontos de cima e de baixo
        # lista_pontos = [0, -radius, 0, 0, radius, 0] 
        lista_pontos = [0, 0, -radius, 0, 0, radius] 

        lista_pontos2d, z = GL.get2DCoord(lista_pontos)


        x_2d = coord_2d[::2]
        y_2d = coord_2d[1::2]
        #print("x -> : ", x_2d)
        #print("y -> : ", y_2d)
        #print("x2d_len")
        for i in range(len(x_2d)-pontos):
            vertices = [] 
            z_list_atual = []
            #primeiro trangulo do retangulo
            if i % pontos == (pontos -1): #5 0 6
                #11 % 6 = 5
                vertices.append(x_2d[i])         # i
                vertices.append(y_2d[i])         # i
                z_list_atual.append(z_list[i])

                vertices.append(x_2d[i+1-pontos])       # 0 (primeiro linha)
                vertices.append(y_2d[i+1-pontos])       # 0  primeiro
                z_list_atual.append(z_list[i+1-pontos])

                vertices.append(x_2d[i+1])  # PROX
                vertices.append(y_2d[i+1])  # PROX
                z_list_atual.append(z_list[i+1])

                GL.fill_triangle(vertices, colors, z_list_atual)

                vertices = []
                z_list_atual = []
                #5 11 6
                vertices.append(x_2d[i])          # 0
                vertices.append(y_2d[i])          # 0
                z_list_atual.append(z_list[i])
                vertices.append(x_2d[i+pontos])   # linha de cima
                vertices.append(y_2d[i+pontos])   # linha de cima
                z_list_atual.append(z_list[i+pontos])
                vertices.append(x_2d[i+1])   # prox
                vertices.append(y_2d[i+1])   # prox
                z_list_atual.append(z_list[i+1])
                GL.fill_triangle(vertices, colors, z_list_atual, reverse=True)




            else:
                z_list_atual = []
                vertices.append(x_2d[i])         # 0
                vertices.append(y_2d[i])         # 0
                z_list_atual.append(z_list[i])
                vertices.append(x_2d[i+1])       # 1
                vertices.append(y_2d[i+1])       # 1
                z_list_atual.append(z_list[i+1])
                vertices.append(x_2d[i+pontos])  # 7
                vertices.append(y_2d[i+pontos])  # 7
                z_list_atual.append(z_list[i+pontos])
                GL.fill_triangle(vertices, colors, z_list_atual)

                vertices = []
                z_list_atual = []
                vertices.append(x_2d[i])          # 0
                vertices.append(y_2d[i])          # 0
                z_list_atual.append(z_list[i])
                vertices.append(x_2d[i+pontos-1]) # 6
                vertices.append(y_2d[i+pontos-1]) # 6
                z_list_atual.append(z_list[i+pontos-1])
                vertices.append(x_2d[i+pontos])   # 7
                vertices.append(y_2d[i+pontos])   # 7
                z_list_atual.append(z_list[i+pontos])
                GL.fill_triangle(vertices, colors, z_list_atual, reverse=True)
                
        
        # PARA FAZER OS TOPOS DA ESFERA
        #print("coord_2d pontos", coord_2d)
        #print("lista pontos", lista_pontos2d)
        i = len(coord_2d) - 1 - pontos
        #print("len(coord_2d) : ", len(coord_2d) )

        while i < len(x_2d) :
           
           vertices = []
           if i % pontos != (pontos -1): 
                vertices.append(x_2d[i])
                vertices.append(y_2d[i])
                vertices.append(x_2d[i+1])
                vertices.append(y_2d[i+1])
                vertices.append(lista_pontos2d[2])
                vertices.append(lista_pontos2d[3])
                GL.fill_triangle(vertices, colors, z_list)
           else:
                vertices.append(x_2d[i])
                vertices.append(y_2d[i])
                vertices.append(x_2d[i-pontos])
                vertices.append(y_2d[i-pontos])
                vertices.append(lista_pontos2d[2])
                vertices.append(lista_pontos2d[3])
                GL.fill_triangle(vertices, colors, z_list)
           i+=1

        i = 0
        while i <pontos :
           
           vertices = []
           #print("i", i)
           if i % pontos != (pontos -1): 
                vertices.append(x_2d[i])
                vertices.append(y_2d[i])
                vertices.append(x_2d[i+1])
                vertices.append(y_2d[i+1])
                vertices.append(lista_pontos2d[0])
                vertices.append(lista_pontos2d[1])
                GL.fill_triangle(vertices, colors, z_list)

           else:
                vertices.append(x_2d[i])
                vertices.append(y_2d[i])
                vertices.append(x_2d[i-pontos])
                vertices.append(y_2d[i-pontos])
                vertices.append(lista_pontos2d[0])
                vertices.append(lista_pontos2d[1])
                GL.fill_triangle(vertices, colors, z_list)

        
           i+=1
    
            # 12    13    18           
            # 13    14    18
            # 14    15    18    
            # 15    16    18
            # 16    17    18

            #12  17    18
            
        
        
        

        



    #i    i+1     i+7          #i      i+7     i+6  
    #0      1       7           0       6       7   
    #1      2       8           1       7       8
    #2      3       9           2       8       9
    #3      4       10          3       9       10     
    #4      5       11          4       10       11
    
    #5      6       12          5       11       12                       
    #5      0       6           5        1  


    @staticmethod
    def navigationInfo(headlight):
        """Características físicas do avatar do visualizador e do modelo de visualização."""
        # O campo do headlight especifica se um navegador deve acender um luz direcional que
        # sempre aponta na direção que o usuário está olhando. Definir este campo como TRUE
        # faz com que o visualizador forneça sempre uma luz do ponto de vista do usuário.
        # A luz headlight deve ser direcional, ter intensidade = 1, cor = (1 1 1),
        # ambientIntensity = 0,0 e direção = (0 0 −1).

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("NavigationInfo : headlight = {0}".format(headlight)) # imprime no terminal

    @staticmethod
    def directionalLight(ambientIntensity, color, intensity, direction):
        """Luz direcional ou paralela."""
        # Define uma fonte de luz direcional que ilumina ao longo de raios paralelos
        # em um determinado vetor tridimensional. Possui os campos básicos ambientIntensity,
        # cor, intensidade. O campo de direção especifica o vetor de direção da iluminação
        # que emana da fonte de luz no sistema de coordenadas local. A luz é emitida ao
        # longo de raios paralelos de uma distância infinita.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("DirectionalLight : ambientIntensity = {0}".format(ambientIntensity))
        print("DirectionalLight : color = {0}".format(color)) # imprime no terminal
        print("DirectionalLight : intensity = {0}".format(intensity)) # imprime no terminal
        print("DirectionalLight : direction = {0}".format(direction)) # imprime no terminal

        GL.Ilrgb = color
        GL.Ii = intensity
        GL.Iia = ambientIntensity
        GL.L = direction*(-1)
        GL.light = True
        GL.v = np.matrix([0,0,1])

        


    @staticmethod
    def pointLight(ambientIntensity, color, intensity, location):
        """Luz pontual."""  
        # Fonte de luz pontual em um local 3D no sistema de coordenadas local. Uma fonte
        # de luz pontual emite luz igualmente em todas as direções; ou seja, é omnidirecional.
        # Possui os campos básicos ambientIntensity, cor, intensidade. Um nó PointLight ilumina
        # a geometria em um raio de sua localização. O campo do raio deve ser maior ou igual a
        # zero. A iluminação do nó PointLight diminui com a distância especificada.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("PointLight : ambientIntensity = {0}".format(ambientIntensity))
        print("PointLight : color = {0}".format(color)) # imprime no terminal
        print("PointLight : intensity = {0}".format(intensity)) # imprime no terminal
        print("PointLight : location = {0}".format(location)) # imprime no terminal

    @staticmethod
    def fog(visibilityRange, color):
        """Névoa."""
        # O nó Fog fornece uma maneira de simular efeitos atmosféricos combinando objetos
        # com a cor especificada pelo campo de cores com base nas distâncias dos
        # vários objetos ao visualizador. A visibilidadeRange especifica a distância no
        # sistema de coordenadas local na qual os objetos são totalmente obscurecidos
        # pela névoa. Os objetos localizados fora de visibilityRange do visualizador são
        # desenhados com uma cor de cor constante. Objetos muito próximos do visualizador
        # são muito pouco misturados com a cor do nevoeiro.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Fog : color = {0}".format(color)) # imprime no terminal
        print("Fog : visibilityRange = {0}".format(visibilityRange))

    @staticmethod
    def timeSensor(cycleInterval, loop):
        """Gera eventos conforme o tempo passa."""
        # Os nós TimeSensor podem ser usados para muitas finalidades, incluindo:
        # Condução de simulações e animações contínuas; Controlar atividades periódicas;
        # iniciar eventos de ocorrência única, como um despertador;
        # Se, no final de um ciclo, o valor do loop for FALSE, a execução é encerrada.
        # Por outro lado, se o loop for TRUE no final de um ciclo, um nó dependente do
        # tempo continua a execução no próximo ciclo. O ciclo de um nódirection TimeSensor dura
        # cycleInterval segundos. O valor de cycleInterval deve ser maior que zero.

        # Deve retornar a fração de tempo passada em fraction_changed

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("TimeSensor : cycleInterval = {0}".format(cycleInterval)) # imprime no terminal
        print("TimeSensor : loop = {0}".format(loop))

        # Esse método já está implementado para os alunos como exemplo
        epoch = time.time()  # time in seconds since the epoch as a floating point number.
        fraction_changed = (epoch % cycleInterval) / cycleInterval

        return fraction_changed

    @staticmethod
    def splinePositionInterpolator(set_fraction, key, keyValue, closed):
        """Interpola não linearmente entre uma lista de vetores 3D."""
        # Interpola não linearmente entre uma lista de vetores 3D. O campo keyValue possui
        # uma lista com os valores a serem interpolados, key possui uma lista respectiva de chaves
        # dos valores em keyValue, a fração a ser interpolada vem de set_fraction que varia de
        # zeroa a um. O campo keyValue deve conter exatamente tantos vetores 3D quanto os
        # quadros-chave no key. O campo closed especifica se o interpolador deve tratar a malha
        # como fechada, com uma transições da última chave para a primeira chave. Se os keyValues
        # na primeira e na última chave não forem idênticos, o campo closed será ignorado.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("SplinePositionInterpolator : set_fraction = {0}".format(set_fraction))
        print("SplinePositionInterpolator : key = {0}".format(key)) # imprime no terminal
        print("SplinePositionInterpolator : keyValue = {0}".format(keyValue))
        print("SplinePositionInterpolator : closed = {0}".format(closed))

        # Abaixo está só um exemplo de como os dados podem ser calculados e transferidos
        value_changed = [0.0, 0.0, 0.0]
        
        return value_changed

    @staticmethod
    def orientationInterpolator(set_fraction, key, keyValue):
        """Interpola entre uma lista de valores de rotação especificos."""
        # Interpola rotações são absolutas no espaço do objeto e, portanto, não são cumulativas.
        # Uma orientação representa a posição final de um objeto após a aplicação de uma rotação.
        # Um OrientationInterpolator interpola entre duas orientações calculando o caminho mais
        # curto na esfera unitária entre as duas orientações. A interpolação é linear em
        # comprimento de arco ao longo deste caminho. Os resultados são indefinidos se as duas
        # orientações forem diagonalmente opostas. O campo keyValue possui uma lista com os
        # valores a serem interpolados, key possui uma lista respectiva de chaves
        # dos valores em keyValue, a fração a ser interpolada vem de set_fraction que varia de
        # zeroa a um. O campo keyValue deve conter exatamente tantas rotações 3D quanto os
        # quadros-chave no key.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("OrientationInterpolator : set_fraction = {0}".format(set_fraction))
        print("OrientationInterpolator : key = {0}".format(key)) # imprime no terminal
        print("OrientationInterpolator : keyValue = {0}".format(keyValue))

        # Abaixo está só um exemplo de como os dados podem ser calculados e transferidos
        value_changed = [0, 0, 1, 0]

        return value_changed

    # Para o futuro (Não para versão atual do projeto.)
    def vertex_shader(self, shader):
        """Para no futuro implementar um vertex shader."""
        print("AAAAAAAAAAAAAAA")

    def fragment_shader(self, shader):
        """Para no futuro implementar um fragment shader."""
        print("OIIIIIIIIIIIIIIIIII")
