# -*- coding: UTF-8 -*-

"""
Rotinas de operação de nós X3D.

Desenvolvido por:
Disciplina: Computação Gráfica
Data:
"""

import gpu          # Simula os recursos de uma GPU

#################################################################################
# NÃO USAR MAIS ESSE ARQUIVO. AS ROTINAS DEVEM SER IMPLEMENTADAS AGORA NO gl.GL #
#################################################################################

# web3d.org/documents/specifications/19775-1/V3.0/Part01/components/geometry2D.html#Polypoint2D
def polypoint2D(point, colors):
    # Nessa função você receberá pontos no parâmetro point, esses pontos são uma lista
    # de pontos x, y sempre na ordem. Assim point[0] é o valor da coordenada x do
    # primeiro ponto, point[1] o valor y do primeiro ponto. Já point[2] é a
    # coordenada x do segundo ponto e assim por diante. Assuma a quantidade de pontos
    # pelo tamanho da lista e assuma que sempre vira uma quantidade par de valores.
    # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Polypoint2D
    # você pode assumir o desenho dos pontos com a cor emissiva (emissiveColor).

    # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
    print("Polypoint2D : pontos = {0}".format(point)) # imprime no terminal pontos
    print("Polypoint2D : colors = {0}".format(colors)) # imprime no terminal as cores

    # Exemplo:
    # gpu.GPU.set_pixel(3, 1, 0, 255, 255) # altera um pixel da imagem (u, v, r, g, b)
    # cuidado com as cores, o X3D especifica de (0,1) e o Framebuffer de (0,255) 
    
    # Para cada ponto, pinta o pixel de acordo (convertendo as coordenadas para int) 
    for i in range(0,len(point),2):
        gpu.GPU.set_pixel(int(point[i]), int(point[i+1]), colors["emissiveColor"][0]*255, colors["emissiveColor"][1]*255, colors["emissiveColor"][2]*255)   
    

def drawLine(lineSegments, colors):
    # Ideia do Código  de: http://rosettacode.org/wiki/Bitmap/Bresenham%27s_line_algorithm#Python

    #cria uma lista com as coordenadas X e outra para as coordenadas Y
    u_list = (lineSegments[::2])
    v_list = (lineSegments[1::2])

    u_list = list(map(int, u_list))
    v_list = list(map(int, v_list))

    sx = -1 if u_list[0] > u_list[1] else 1       # if x0 > x1 
    sy = -1 if v_list[0] > v_list[1] else 1       # if y0 > y1 
    
    u, v = u_list[0], v_list[0]

    # if u_list[0] > u_list[1]:    
    #     v = v_list[1]
    # else:
    #     v = v_list[0]
        
    # if v_list[0] > v_list[1]:    
    #     u = u_list[1]
    # else:
    #     u = u_list[0]
    #calcula o coeficiente 
    dx = abs(u_list[1] - u_list[0])
    dy = abs(v_list[1] - v_list[0])
    if dx > dy: 
        err = dx / 2.0
        while u != u_list[1]:
            gpu.GPU.set_pixel(u, v, colors["emissiveColor"][0]*255, colors["emissiveColor"][1]*255, colors["emissiveColor"][2]*255)   
            err -= dy
            if err < 0:
                v += sy
                err += dx
            u += sx

    else:
        err = dy / 2.0
        
        while v != v_list[1]:
            gpu.GPU.set_pixel(u, v, colors["emissiveColor"][0]*255, colors["emissiveColor"][1]*255, colors["emissiveColor"][2]*255)   
            err -= dx
            if err < 0:
                u += sx
                err += dy
            v += sy

    # if dx > dy: 
    #     for i in range(int(min(u_list)),int(max(u_list))):
    #         #Arredonda V
    #         gpu.GPU.set_pixel(i, v, colors["emissiveColor"][0]*255, colors["emissiveColor"][1]*255, colors["emissiveColor"][2]*255)   
        
    #         #Adiciona o coef angular em V para pintar o prox ponto
    #         v+= m
    # else:
    #     for i in range(int(min(v_list)),int(max(v_list))):
    #         print(u)
    #         #Arredonda U
    #         gpu.GPU.set_pixel(int(u), i, colors["emissiveColor"][0]*255, colors["emissiveColor"][1]*255, colors["emissiveColor"][2]*255)   
        
    #         #Adiciona o coef angular em V para pintar o prox ponto
    #         u+= m

# web3d.org/documents/specifications/19775-1/V3.0/Part01/components/geometry2D.html#Polyline2D
def polyline2D(lineSegments, colors):
    """Função usada para renderizar Polyline2D."""
    # Nessa função você receberá os pontos de uma linha no parâmetro lineSegments, esses
    # pontos são uma lista de pontos x, y sempre na ordem. Assim point[0] é o valor da
    # coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto. Já point[2] é
    # a coordenada x do segundo ponto e assim por diante. Assuma a quantidade de pontos
    # pelo tamanho da lista. A quantidade mínima de pontos são 2 (4 valores), porém a
    # função pode receber mais pontos para desenhar vários segmentos. Assuma que sempre
    # vira uma quantidade par de valores.
    # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Polyline2D
    # você pode assumir o desenho das linhas com a cor emissiva (emissiveColor).
    
    # Faz a função de desenhar a linha para cada linha no lineSegments    
    for i in range(0,int(len(lineSegments)-2),2):
        drawLine(lineSegments[i:i+4], colors)
    
    print("Polyline2D : lineSegments = {0}".format(lineSegments)) # imprime no terminal
    print("Polyline2D : colors = {0}".format(colors)) # imprime no terminal as cores
    # Exemplo:
    pos_x = gpu.GPU.width//2
    pos_y = gpu.GPU.height//2
    # gpu.GPU.set_pixel(pos_x, pos_y, 255, 0, 0) # altera um pixel da imagem (u, v, r, g, b)


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
        if l1 and l2 and l3:
            return True
    return False
    

# web3d.org/documents/specifications/19775-1/V3.0/Part01/components/geometry2D.html#TriangleSet2D
def triangleSet2D(vertices, colors):
    """Função usada para renderizar TriangleSet2D."""
    # Nessa função você receberá os vertices de um triângulo no parâmetro vertices,
    # esses pontos são uma lista de pontos x, y sempre na ordem. Assim point[0] é o
    # valor da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto.
    # Já point[2] é a coordenada x do segundo ponto e assim por diante. Assuma que a
    # quantidade de pontos é sempre multiplo de 3, ou seja, 6 valores ou 12 valores, etc.
    # O parâmetro colors é um dicionário com os tipos cores possíveis, para o TriangleSet2D
    # você pode assumir o desenho das linhas com a cor emissiva (emissiveColor).
    print("TriangleSet2D : vertices = {0}".format(vertices)) # imprime no terminal
    print("TriangleSet2D : colors = {0}".format(colors)) # imprime no terminal as cores
    # # Exemplo:
    # gpu.GPU.set_pixel(24, 8, 255, 255, 0) # altera um pixel da imagem (u, v, r, g, b)

    x_list = vertices[::2]
    y_list = vertices[1::2]

    for i in range(int(min(x_list)), int(max(x_list))):
        for j in range(int(min(y_list)), int(max(y_list))):
            if inside(x_list, y_list, i+0.5,j+0.5):
                gpu.GPU.set_pixel(i, j, colors["emissiveColor"][0]*255, colors["emissiveColor"][1]*255, colors["emissiveColor"][2]*255)



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

def transform_out():
    """Função usada para renderizar (na verdade coletar os dados) de Transform."""
    # A função transform_out será chamada quando se sair em um nó X3D do tipo Transform do
    # grafo de cena. Não são passados valores, porém quando se sai de um nó transform se
    # deverá recuperar a matriz de transformação dos modelos do mundo da estrutura de
    # pilha implementada.

    # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
    print("Saindo de Transform")

def triangleStripSet(point, stripCount, colors):
    """Função usada para renderizar TriangleStripSet."""
    # A função triangleStripSet é usada para desenhar tiras de triângulos interconectados,
    # você receberá as coordenadas dos pontos no parâmetro point, esses pontos são uma
    # lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x
    # do primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
    # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e assim
    # por diante. No TriangleStripSet a quantidade de vértices a serem usados é informado
    # em uma lista chamada stripCount (perceba que é uma lista).

    # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
    print("TriangleStripSet : pontos = {0} ".format(point), end='')
    for i, strip in enumerate(stripCount):
        print("strip[{0}] = {1} ".format(i, strip), end='')
    print("")
    print("TriangleStripSet : colors = {0}".format(colors)) # imprime no terminal as cores

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
    # depois 2, 3 e 4, e assim por diante.

    # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
    print("IndexedTriangleStripSet : pontos = {0}, index = {1}".format(point, index))
    print("IndexedTriangleStripSet : colors = {0}".format(colors)) # imprime no terminal as cores

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
    print("IndexedFaceSet : colors = {0}".format(colors)) # imprime no terminal as cores
