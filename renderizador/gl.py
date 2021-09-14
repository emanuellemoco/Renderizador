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
        # p_matrix 4,4

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
        print("TriangleSet : pontos = {0}".format(point)) # imprime no terminal pontos
        print("TriangleSet : colors = {0}".format(colors)) # imprime no terminal as cores

        # Exemplo de desenho de um pixel branco na coordenada 10, 10
        gpu.GPU.draw_pixels([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel

        # Para cada triangulo
        for i in range(0,len(point),9):
            if i+8>len(point):
                break
            for j in range(i, i+3):
                point3d = np.matrix([
                    [point[j],point[j+1],point[j+2],1]
                ])  ##==> GL.world * isso
                world_point = point3d.dot(GL.world)
                #Identificando as coordenadas ortonormais (matrix Look-At)

                

                # ignorando o quarto valor
                eye = np.matrix([
                    [GL.position[0,0], GL.position[0,1], GL.position[0,2]]
                ])
                world_point3d = np.matrix([
                    [world_point[0,0], world_point[0,1], world_point[0,2]]
                ])       
                # print("world_point3d: ", world_point3d)
                # print("position: ", position)
                orientation = np.matrix([
                [GL.orientation[0], GL.orientation[1], GL.orientation[2]]
                ]) 
                # orientation = np.matrix([[0, 1, 0]]) 
                # print("orientation ", orientation)

                w = np.divide(np.subtract(world_point3d, eye) , np.abs(np.subtract(world_point3d, eye))) #
            
                u = np.divide(np.cross(w, orientation)  , np.abs(np.cross(w, orientation))) 
                print("u ", u)
                v = np.divide(np.cross(u, w)  , np.abs(np.cross(u, w)) ) 
                # print("v ", v)

                #ARRUMAR AQUI
                #Transformação de Look-At
                lookAt = np.matrix([
                    [u[0,0], v[0,0], - w[0,0], eye[0,0]],
                    [u[0,1], v[0,1], -w[0,1],  eye[0,1]],
                    [u[0,2], v[0,2], -w[0,2],  eye[0,2]],
                    [0, 0, 0, 1],
                ])
                #A transformação de Look-At é o inverso da matriz Look-At
                lookAt_tranform = np.linalg.inv(lookAt)
                # print("lookAt_tranform: ", lookAt_tranform)

                #Coordenadas do Mundo -> Coordenadas Camera
                camera_point = world_point.dot(lookAt_tranform)

                #Transformações Perspectivas
                perspective = camera_point.dot(GL.p_matrix )

                #Divisão homogenea para fazer a normalização:
                perspective_normalized = perspective/perspective[0,3]

                # print("perspective_normalized: ", perspective_normalized)
            
            

        
            #Translacao 

            #World (multiplicar a matriz )
            

            #gpu.GPU.draw_pixels([int(point_1[0]), int(point_1[1])], gpu.GPU.RGB8,  [colors["diffuseColor"][0]*255, colors["diffuseColor"][1]*255, colors["diffuseColor"][2]*255])
            #gpu.GPU.draw_pixels([int(point_2[0]), int(point_2[1])], gpu.GPU.RGB8,  [colors["diffuseColor"][0]*255, colors["diffuseColor"][1]*255, colors["diffuseColor"][2]*255])
            #gpu.GPU.draw_pixels([int(point_3[0]), int(point_3[1])], gpu.GPU.RGB8,  [colors["diffuseColor"][0]*255, colors["diffuseColor"][1]*255, colors["diffuseColor"][2]*255])


    


    @staticmethod
    def viewpoint(position, orientation, fieldOfView):
        """Função usada para renderizar (na verdade coletar os dados) de Viewpoint."""
        # Na função de viewpoint você receberá a posição, orientação e campo de visão da
        # câmera virtual. Use esses dados para poder calcular e criar a matriz de projeção
        # perspectiva para poder aplicar nos pontos dos objetos geométricos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Viewpoint : ", end='')
        print("position = {0} ".format(position), end='')
        print("orientation = {0} ".format(orientation), end='')
        print("fieldOfView = {0} ".format(fieldOfView))
        ####### remover 4 linhas a cima 

        fovy = 2 * math.atan(math.tan(fieldOfView / 2) * GL.height / math.sqrt(GL.height ** 2 + GL.width ** 2))

        GL.position = np.matrix([
            [position[0], position[1], position[2], 1]
        ]) 
        GL.orientation = orientation

        # GL.orientation = np.matrix([
        #     [orientation[0], orientation[1], orientation[2], orientation[3]]
        # ]) 

    
        top = GL.near * fovy 
        GL.p_matrix = np.matrix([[GL.near/GL.width, 0, 0, 0],[0, GL.near/top, 0,0], [0,0, -(GL.far+GL.near)/(GL.far-GL.near), (-2 * GL.far * GL.near)/(GL.far-GL.near)], [0, 0, -1, 0]])


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
        print("Transform : ", end='')
        if translation:
            print("translation = {0} ".format(translation), end='') # imprime no terminal
        if scale:
            print("scale = {0} ".format(scale), end='') # imprime no terminal
        if rotation:
            print("rotation = {0} ".format(rotation), end='') # imprime no terminal
        print("")
    
        # O Braga nos explicou como fazer essa parte

        scale = np.matrix([
            [scale[0], 0, 0 ,0],
            [0, scale[1], 0, 0],
            [0, 0, scale[2], 0],
            [0, 0, 0, 1]
        ])
        


        translation = np.matrix([
            [1, 0, 0, translation[0]],
            [0, 1, 0, translation[1]],
            [1, 0, 1, translation[2]],
            [0, 0, 0, 1]
        ])

        #Matrix de Rotações: https://en.wikipedia.org/wiki/Rotation_matrix
        rotationX = np.matrix([
            [1, 0, 0, 0],
            [0, math.cos(rotation[3]), -math.sin(rotation[3]), 0],
            [0, math.sin(rotation[3]), -math.cos(rotation[3]), 0],
            [0, 0, 0, 1]
        ]) 

        rotationY = np.matrix([
            [math.cos(rotation[3]), 0, math.sin(rotation[3]), 0],
            [0, 1, 0, 0],
            [-math.sin(rotation[3]), 0, math.sin(rotation[3]), 0],
            [0, 0, 0, 1]
        ])
        
        rotationZ = np.matrix([
            [math.cos(rotation[3]), -math.sin(rotation[3]), 0, 0],
            [math.sin(rotation[3]), math.cos(rotation[3]), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ])

        if rotation[0] != 0:
            GL.world = scale.dot(rotationX).dot(translation)
        elif rotation[1] != 0:
            GL.world = scale.dot(rotationY).dot(translation)
        elif rotation[2] != 0:
            GL.world = scale.dot(rotationZ).dot(translation)
        else:
            GL.world = scale.dot(translation)
        


    @staticmethod
    def transform_out():
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_out será chamada quando se sair em um nó X3D do tipo Transform do
        # grafo de cena. Não são passados valores, porém quando se sai de um nó transform se
        # deverá recuperar a matriz de transformação dos modelos do mundo da estrutura de
        # pilha implementada.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Saindo de Transform")

    @staticmethod
    def triangleStripSet(point, stripCount, colors):
        """Função usada para renderizar TriangleStripSet."""
        # A função triangleStripSet é usada para desenhar tiras de triângulos interconectados,
        # você receberá as coordenadas dos pontos no parâmetro point, esses pontos são uma
        # lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x
        # do primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
        # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e assim
        # por diante. No TriangleStripSet a quantidade de vértices a serem usados é informado
        # em uma lista chamada stripCount (perceba que é uma lista). Ligue os vértices na ordem,
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante. Cuidado com a orientação dos vértices, ou seja,
        # todos no sentido horário ou todos no sentido anti-horário, conforme especificado.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("TriangleStripSet : pontos = {0} ".format(point), end='')
        for i, strip in enumerate(stripCount):
            print("strip[{0}] = {1} ".format(i, strip), end='')
        print("")
        print("TriangleStripSet : colors = {0}".format(colors)) # imprime no terminal as cores

        # Exemplo de desenho de um pixel branco na coordenada 10, 10
        gpu.GPU.draw_pixels([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel

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

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("IndexedTriangleStripSet : pontos = {0}, index = {1}".format(point, index))
        print("IndexedTriangleStripSet : colors = {0}".format(colors)) # imprime as cores

        # Exemplo de desenho de um pixel branco na coordenada 10, 10
        gpu.GPU.draw_pixels([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel

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
        print("Box : size = {0}".format(size)) # imprime no terminal pontos
        print("Box : colors = {0}".format(colors)) # imprime no terminal as cores

        # Exemplo de desenho de um pixel branco na coordenada 10, 10
        gpu.GPU.draw_pixels([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel

    @staticmethod
    def indexedFaceSet(coord, coordIndex, colorPerVertex, color, colorIndex,
                       texCoord, texCoordIndex, colors, current_texture):
        """Função usada para renderizar IndexedFaceSet."""
        # A função indexedFaceSet é usada para desenhar malhas de triângulos. Ela funciona de
        # forma muito simular a IndexedTriangleStripSet porém com mais recursos.
        # Você receberá as coordenadas dos pontos no parâmetro cord, esses
        # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor
        # da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto, point[2]
        # o valor z da coordenada z do primeiro ponto. Já point[3] é a coordenada x do
        # segundo ponto e assim por diante. No IndexedFaceSet uma lista informando
        # como conectar os vértices é informada em coordIndex, o valor -1 indica que a lista
        # acabou. A ordem de conexão será de 3 em 3 pulando um índice. Por exemplo: o
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante.
        # Adicionalmente essa implementação do IndexedFace suport cores por vértices, assim
        # a se a flag colorPerVertex estiver habilidades, os vértices também possuirão cores
        # que servem para definir a cor interna dos poligonos, para isso faça um cálculo
        # baricêntrico de que cor deverá ter aquela posição. Da mesma forma se pode definir uma
        # textura para o poligono, para isso, use as coordenadas de textura e depois aplique a
        # cor da textura conforme a posição do mapeamento. Dentro da classe GPU já está
        # implementadado um método para a leitura de imagens.

        # Os prints abaixo são só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("IndexedFaceSet : ")
        if coord:
            print("\tpontos(x, y, z) = {0}, coordIndex = {1}".format(coord, coordIndex))
        if colorPerVertex:
            print("\tcores(r, g, b) = {0}, colorIndex = {1}".format(color, colorIndex))
        if texCoord:
            print("\tpontos(u, v) = {0}, texCoordIndex = {1}".format(texCoord, texCoordIndex))
        if current_texture:
            image = gpu.GPU.load_texture(current_texture[0])
            print("\t Matriz com image = {0}".format(image))
        print("IndexedFaceSet : colors = {0}".format(colors))  # imprime no terminal as cores

        # Exemplo de desenho de um pixel branco na coordenada 10, 10
        gpu.GPU.draw_pixels([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel
