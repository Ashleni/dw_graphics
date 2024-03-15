from display import *
from matrix import *
from gmath import *

def add_mesh(polygons, file):
    f = open(file + '.obj', "r").read().split("\n")
    vertices = []
    
    for x in f:
        line = x.split()
        if len(line) != 0:
            if line[0] == "v":
                vertice = []
                for coord in line[1:4]:
                    vertice.append(float(coord))
                vertices.append(vertice)
            if line[0] == "f":
                vertice = []
                for points in line[1:]:
                    point = points.split("/")
                    vertice.append(int(point[0]))
                    #vertice.append(int(point))
                p0 = vertices[vertice[0]-1]
                p1 = vertices[vertice[1]-1]
                p2 = vertices[vertice[2]-1]
                if len(vertice) == 3 or len(vertice) == 4:
                    add_polygon(polygons, p0[0], p0[1], p0[2], 
                                          p1[0], p1[1], p1[2], 
                                          p2[0], p2[1], p2[2])
                    if len(vertice) == 4:
                        p3 = vertices[vertice[3]-1]
                        add_polygon(polygons, p0[0], p0[1], p0[2], 
                                              p2[0], p2[1], p2[2], 
                                              p3[0], p3[1], p3[2])

def draw_scanline(x0, z0, x1, z1, y, screen, zbuffer, color):
    if x0 > x1:
        tx = x0
        tz = z0
        x0 = x1
        z0 = z1
        x1 = tx
        z1 = tz

    x = x0
    z = z0
    delta_z = (z1 - z0) / (x1 - x0 + 1) if (x1 - x0 + 1) != 0 else 0

    while x <= x1:
        plot(screen, zbuffer, color, x, y, z)
        x+= 1
        z+= delta_z

def scanline_convert(polygons, i, screen, zbuffer, color):
    flip = False
    BOT = 0
    TOP = 2
    MID = 1

    points = [ (polygons[i][0], polygons[i][1], polygons[i][2]),
               (polygons[i+1][0], polygons[i+1][1], polygons[i+1][2]),
               (polygons[i+2][0], polygons[i+2][1], polygons[i+2][2]) ]

    # alas random color, we hardly knew ye
    #color = [0,0,0]
    #color[RED] = (23*(i/3)) %256
    #color[GREEN] = (109*(i/3)) %256
    #color[BLUE] = (227*(i/3)) %256

    points.sort(key = lambda x: x[1])
    if (points[BOT][0] > points[MID][0] and points[BOT][1] == points[MID][1]):
        BOT = 1
        MID = 0
    x0 = points[BOT][0]
    z0 = points[BOT][2]
    x1 = points[BOT][0]
    z1 = points[BOT][2]
    y = int(points[BOT][1])
    

    distance0 = int(points[TOP][1]) - y * 1.0 + 1
    distance1 = int(points[MID][1]) - y * 1.0 + 1
    distance2 = int(points[TOP][1]) - int(points[MID][1]) * 1.0 + 1

    dx0 = (points[TOP][0] - points[BOT][0]) / distance0 if distance0 != 0 else 0
    dz0 = (points[TOP][2] - points[BOT][2]) / distance0 if distance0 != 0 else 0
    dx1 = (points[MID][0] - points[BOT][0]) / distance1 if distance1 != 0 else 0
    dz1 = (points[MID][2] - points[BOT][2]) / distance1 if distance1 != 0 else 0

    while y <= int(points[TOP][1]):
        if ( not flip and y >= int(points[MID][1])):
            flip = True

            dx1 = (points[TOP][0] - points[MID][0]) / distance2 if distance2 != 0 else 0
            dz1 = (points[TOP][2] - points[MID][2]) / distance2 if distance2 != 0 else 0
            x1 = points[MID][0]
            z1 = points[MID][2]

        #draw_line(int(x0), y, z0, int(x1), y, z1, screen, zbuffer, color)
        draw_scanline(int(x0), z0, int(x1), z1, y, screen, zbuffer, color)
        x0+= dx0
        z0+= dz0
        x1+= dx1
        z1+= dz1
        y+= 1
        


def add_polygon( polygons, x0, y0, z0, x1, y1, z1, x2, y2, z2 ):
    add_point(polygons, x0, y0, z0)
    add_point(polygons, x1, y1, z1)
    add_point(polygons, x2, y2, z2)

def draw_polygons_normal( polygons, screen, zbuffer, view, ambient, light, symbols, reflect):
    ### GOURAUD BATTLE PLAN:
    # 1. Loop through polygons and calculate the surface normal, then add the three vertecies to a dict (vertex : normal)
    #       - if the vertex is already in the dict, then we update the normal (old_normal + new_normal), this inevitably causes each vertex to be a vector of all connected polygon normals.
    #       - idk I think you loop through the dict of vertices and normalize everything again
    # 2. Calculate the color at each vertex normal (I = IA + ID + IS)
    # 3. Interpolate the color as we scanline (so we need delta_y of color and delta_x of color (I think for each R,G,B value))
    if len(polygons) < 2:
        print('Need at least 3 points to draw')
        return

    point = 0

    vertecii = {}

    # Add vertex(s) and their normals to a dictionary 
    while point < len(polygons) - 2:
        normal = calculate_normal(polygons, point)[:]
        if normal[2] > 0:
            for i in range(3): #loop through each vertex POINT OF ERROR 1
                # print(f"...attempting to create vertex {i}")
                tuplex = (polygons[point + i][0], polygons[point + i][1], polygons[point + i][2])
                tuplex_str = ", ".join(map(str, tuplex))
                # print(f"completed vertex {i}: {tuplex_str}!")
                if tuplex in vertecii:
                    # print(f"...attempting to update vertex {i}")
                    vertecii[tuplex][0] += normal[0]
                    vertecii[tuplex][1] += normal[1]
                    vertecii[tuplex][2] += normal[2]
                    # print("completed update!")
                else:
                    # print("...attempting to add vertex to dictionary")
                    vertecii[tuplex] = normal
                    # print("completed vertex adding!")
        point+= 3
    print("completed dictionary!")

    # Attempt to loop through the dictionary and normalize everything and then calculate the color 
    print("...attempting to normalize vectors in dictionary and convert to colors")
    for key in vertecii:
        normalize(vertecii[key])
        vertecii[key] = get_lighting(vertecii[key], view, ambient, light, symbols, reflect)
    print("completed normalization and obtained colors!")

    point = 0

    while point < len(polygons) - 2:

        normal = calculate_normal(polygons, point)[:]

        #print normal
        if normal[2] > 0:

            color = get_lighting(normal, view, ambient, light, symbols, reflect )
            scanline_convert(polygons, point, screen, zbuffer, color)

            # draw_line( int(polygons[point][0]),
            #            int(polygons[point][1]),
            #            polygons[point][2],
            #            int(polygons[point+1][0]),
            #            int(polygons[point+1][1]),
            #            polygons[point+1][2],
            #            screen, zbuffer, color)
            # draw_line( int(polygons[point+2][0]),
            #            int(polygons[point+2][1]),
            #            polygons[point+2][2],
            #            int(polygons[point+1][0]),
            #            int(polygons[point+1][1]),
            #            polygons[point+1][2],
            #            screen, zbuffer, color)
            # draw_line( int(polygons[point][0]),
            #            int(polygons[point][1]),
            #            polygons[point][2],
            #            int(polygons[point+2][0]),
            #            int(polygons[point+2][1]),
            #            polygons[point+2][2],
            #            screen, zbuffer, color)
        point+= 3

def draw_scanline_gouraud(x0, z0, x1, z1, y, screen, zbuffer, x0_color, x1_color):
    if x0 > x1:
        tx = x0
        tz = z0
        x0 = x1
        z0 = z1
        x1 = tx
        z1 = tz
        tx0_color = x0_color
        x0_color = x1_color
        x1_color = tx0_color
    
    x = x0
    z = z0
    delta_z = (z1 - z0) / (x1 - x0 + 1) if (x1 - x0 + 1) != 0 else 0

    while x <= x1:
        # print(x0_color)
        if ((x1 - x0) == 0):
            color = [int(x1_color[0]),int(x1_color[1]),int(x1_color[2])]
            # color =[0,0,0]
        # elif x0 - x == 0:
        #     color = [int(x0_color[0]),int(x0_color[1]),int(x0_color[2])]
        #     # color =[0,0,0]
        else:
            a = ((x1 - x)/(x1 - x0))
            b = ((x - x0)/(x1 - x0))
            var_x0_color = [x0_color[0] * a, x0_color[1] * a,x0_color[2] * a]
            var_x1_color = [x1_color[0] * b, x1_color[1] * b,x1_color[2] * b]
            color = [int(a + b) for a, b in zip(var_x0_color, var_x1_color)]
            # if color[1] >220 and color[2] >220:
            #     print(f"{x0_color} {x1_color} {color} {var_x0_color} {var_x1_color} {((x1 - x)/(x1 - x0))} {((x - x0)/(x1 - x0))}")
            #     color = [0,0,0]
        # if x1_color == x0_color:
        #     color = [0,0,0]
        plot(screen, zbuffer, color, x, y, z)
        x+= 1
        z+= delta_z

def scanline_convert_gouraud(polygons, i, screen, zbuffer, vertexes):
    flip = False
    BOT = 0
    TOP = 2
    MID = 1

    points = [ (polygons[i][0], polygons[i][1], polygons[i][2]),
               (polygons[i+1][0], polygons[i+1][1], polygons[i+1][2]),
               (polygons[i+2][0], polygons[i+2][1], polygons[i+2][2]) ]

    # alas random color, we hardly knew ye
    #color = [0,0,0]
    #color[RED] = (23*(i/3)) %256
    #color[GREEN] = (109*(i/3)) %256
    #color[BLUE] = (227*(i/3)) %256

    points.sort(key = lambda x: x[1])

    x0 = points[BOT][0]
    z0 = points[BOT][2]
    x1 = points[BOT][0] 
    z1 = points[BOT][2]
    y = int(points[BOT][1])
    y_double = points[BOT][1]
    dy = (points[TOP][1] - points[BOT][1])/(int(points[TOP][1]) - y) if (int(points[TOP][1]) - y) != 0 else 0

    distance0 = int(points[TOP][1]) - y * 1.0 + 1
    distance1 = int(points[MID][1]) - y * 1.0 + 1
    distance2 = int(points[TOP][1]) - int(points[MID][1]) * 1.0 + 1

    # delta y of color change 

    dx0 = (points[TOP][0] - points[BOT][0]) / distance0 if distance0 != 0 else 0
    dz0 = (points[TOP][2] - points[BOT][2]) / distance0 if distance0 != 0 else 0
    dx1 = (points[MID][0] - points[BOT][0]) / distance1 if distance1 != 0 else 0
    dz1 = (points[MID][2] - points[BOT][2]) / distance1 if distance1 != 0 else 0

    

    # change in color from bottom to top 
    while y <= int(points[TOP][1]):
        if ( not flip and y >= int(points[MID][1])):
            flip = True

            dx1 = (points[TOP][0] - points[MID][0]) / distance2 if distance2 != 0 else 0
            dz1 = (points[TOP][2] - points[MID][2]) / distance2 if distance2 != 0 else 0
            x1 = points[MID][0]
            z1 = points[MID][2]

        #draw_line(int(x0), y, z0, int(x1), y, z1, screen, zbuffer, color)
        x0_color = get_intensity(x0, y_double, points, 0.2, vertexes)
        x1_color = get_intensity(x1, y_double, points, 0.2, vertexes)
        draw_scanline_gouraud(int(x0), z0, int(x1), z1, y, screen, zbuffer, x0_color, x1_color)
        x0+= dx0
        z0+= dz0
        x1+= dx1
        z1+= dz1
        y+= 1
        y_double+=dy

def get_intensity(x, y, vertecies, threshold, vertex_color):
    vertecies = get_vertecies(x, y, vertecies, threshold)
    if (not vertecies[0]):
        print ("Error! Error! Error!")
    else:
        BOT = 0
        TOP = 1
        vertecies.pop(0)
        if vertecies[0][1] > vertecies[1][1]:
            BOT = 1
            TOP = 0
        intensity_TOP = vertex_color[vertecies[TOP]]
        intensity_BOT = vertex_color[vertecies[BOT]]
        yB = vertecies[BOT][1]
        yT = vertecies[TOP][1]
        if intensity_BOT == intensity_TOP:
            return intensity_BOT
        if yB == yT:
            if abs(x - vertecies[TOP][0]) < 0.1 and abs(y - vertecies[TOP][1]) < 0.1:
                # return [0,0,0]
                return intensity_TOP
            elif abs(x - vertecies[BOT][0]) < 0.1 and abs(y - vertecies[BOT][1]) < 0.1:
                return intensity_BOT
            else:
                # print(f"{x} {y} {intensity_TOP} {intensity_BOT} {vertecies[BOT][1]} {vertecies[TOP][1]} {vertecies[BOT]} {vertecies[TOP]} {vertecies[2]}")
                # print("PPPPPPP")
                return intensity_TOP
        if yT - yB < 1:
            if TOP == 0:
                # return[0,0,0]
                return vertex_color[vertecies[1]]
            else:
                return vertex_color[vertecies[1]]
        slope = (yT - yB)/(vertecies[TOP][0]/vertecies[BOT][0])
        r = 1 * ((slope * (x - vertecies[TOP][0])) + vertecies[TOP][1])
        a = ((y - yB)/(yT - yB))
        b = ((yT - y)/(yT - yB))
        # normalizer = 1/(a + b)
        # a *= normalizer
        # b *= normalizer
        intensity_TOP = [intensity_TOP[0] * a,intensity_TOP[1] * a,intensity_TOP[2] * a]
        intensity_BOT = [intensity_BOT[0] * b,intensity_BOT[1] * b,intensity_BOT[2] * b]
        intensity = [intensity_TOP[0] + intensity_BOT[0],intensity_TOP[1] + intensity_BOT[1],intensity_TOP[2] + intensity_BOT[2]]
        if(intensity[0] > 255 or intensity[1] > 255 or intensity[2] > 255):

            # print(vertecies)
            # print(f' WAAAA {x}, {y}, {r}, YYYY: {yB} {yT} {vertecies[BOT]}, {vertecies[TOP]}, {intensity_TOP}, {intensity_BOT}, {intensity}, {a}, {b}')
            intensity_TOP = vertex_color[vertecies[TOP]]
            intensity_BOT = vertex_color[vertecies[BOT]]
            r = 1 * ((slope * (x - vertecies[TOP][0])) + vertecies[TOP][1])
            a = ((r - yB)/(yT - yB))
            b = ((yT - r)/(yT - yB))
            intensity_TOP = [intensity_TOP[0] * a,intensity_TOP[1] * a,intensity_TOP[2] * a]
            intensity_BOT = [intensity_BOT[0] * b,intensity_BOT[1] * b,intensity_BOT[2] * b]
            intensity = [intensity_TOP[0] + intensity_BOT[0],intensity_TOP[1] + intensity_BOT[1],intensity_TOP[2] + intensity_BOT[2]]
            # if (y-yB) > (yT - y):
            #     return vertex_color[vertecies[TOP]]
            # else:
            #     return vertex_color[vertecies[BOT]]
        return intensity 

def get_vertecies(x, y, vertecies, threshold):
    for vertex_0 in vertecies:
        for vertex_1 in vertecies:
            if abs(math.dist([x, y], [vertex_0[0], vertex_0[1]]) + math.dist([x, y], [vertex_1[0], vertex_1[1]]) - math.dist([vertex_0[0], vertex_0[1]], [vertex_1[0], vertex_1[1]])) < threshold:
                return [True, vertex_0, vertex_1, vertecies]
    # print(str(x) + ", " + str(y))
    # print(vertecies)
    for vertex_0 in vertecies:
        if abs(y - vertex_0[1]) < 2 and abs(x - vertex_0[0]) < 2:
            # print(str(x) + ", " + str(y))
            # print(vertecies)
            return [True, vertex_0, vertex_0, vertecies]
    return [True, vertecies[0], vertecies[1], vertecies]

def draw_polygons_gouraud( polygons, screen, zbuffer, view, ambient, light, symbols, reflect):
    ### GOURAUD BATTLE PLAN:
    # 1. Loop through polygons and calculate the surface normal, then add the three vertecies to a dict (vertex : normal)
    #       - if the vertex is already in the dict, then we update the normal (old_normal + new_normal), this inevitably causes each vertex to be a vector of all connected polygon normals.
    #       - idk I think you loop through the dict of vertices and normalize everything again
    # 2. Calculate the color at each vertex normal (I = IA + ID + IS)
    # 3. Interpolate the color as we scanline (so we need delta_y of color and delta_x of color (I think for each R,G,B value))
    if len(polygons) < 2:
        print('Need at least 3 points to draw')
        return

    point = 0

    vertecii = {}

    # Add vertex(s) and their normals to a dictionary 
    while point < len(polygons) - 2:
        normal = calculate_normal(polygons, point)[:]
        if normal[2] > 0:
            for i in range(3): #loop through each vertex POINT OF ERROR 1
                # print(f"...attempting to create vertex {i}")
                tuplex = (polygons[point + i][0], polygons[point + i][1], polygons[point + i][2])
                tuplex_str = ", ".join(map(str, tuplex))
                # print(f"completed vertex {i}: {tuplex_str}!")
                if tuplex in vertecii:
                    # print(f"...attempting to update vertex {i}")
                    vertecii[tuplex][0] += normal[0]
                    vertecii[tuplex][1] += normal[1]
                    vertecii[tuplex][2] += normal[2]
                    # print("completed update!")
                else:
                    # print("...attempting to add vertex to dictionary")
                    vertecii[tuplex] = normal
                    # print("completed vertex adding!")
        point+= 3
    print("completed dictionary!")

    # Attempt to loop through the dictionary and normalize everything and then calculate the color 
    print("...attempting to normalize vectors in dictionary and convert to colors")
    for key in vertecii:
        #normalize(vertecii[key])
        vertecii[key] = get_lighting(vertecii[key], view, ambient, light, symbols, reflect)
    print("completed normalization and obtained colors!")

    point = 0

    while point < len(polygons) - 2:

        normal = calculate_normal(polygons, point)[:]

        #print normal
        if normal[2] > 0:

            scanline_convert_gouraud(polygons, point, screen, zbuffer, vertecii)

            # draw_line( int(polygons[point][0]),
            #            int(polygons[point][1]),
            #            polygons[point][2],
            #            int(polygons[point+1][0]),
            #            int(polygons[point+1][1]),
            #            polygons[point+1][2],
            #            screen, zbuffer, color)
            # draw_line( int(polygons[point+2][0]),
            #            int(polygons[point+2][1]),
            #            polygons[point+2][2],
            #            int(polygons[point+1][0]),
            #            int(polygons[point+1][1]),
            #            polygons[point+1][2],
            #            screen, zbuffer, color)
            # draw_line( int(polygons[point][0]),
            #            int(polygons[point][1]),
            #            polygons[point][2],
            #            int(polygons[point+2][0]),
            #            int(polygons[point+2][1]),
            #            polygons[point+2][2],
            #            screen, zbuffer, color)
        point+= 3

def draw_polygons( polygons, screen, zbuffer, view, ambient, light, symbols, reflect, shade_type):
    if shade_type == "gouraud":
        draw_polygons_gouraud( polygons, screen, zbuffer, view, ambient, light, symbols, reflect)
    elif shade_type == "flat":
        draw_polygons_normal( polygons, screen, zbuffer, view, ambient, light, symbols, reflect)
    else:   
        draw_polygons_gouraud( polygons, screen, zbuffer, view, ambient, light, symbols, reflect)

def add_box( polygons, x, y, z, width, height, depth ):
    x1 = x + width
    y1 = y - height
    z1 = z - depth

    #front
    add_polygon(polygons, x, y, z, x1, y1, z, x1, y, z)
    add_polygon(polygons, x, y, z, x, y1, z, x1, y1, z)

    #back
    add_polygon(polygons, x1, y, z1, x, y1, z1, x, y, z1)
    add_polygon(polygons, x1, y, z1, x1, y1, z1, x, y1, z1)

    #right side
    add_polygon(polygons, x1, y, z, x1, y1, z1, x1, y, z1)
    add_polygon(polygons, x1, y, z, x1, y1, z, x1, y1, z1)
    #left side
    add_polygon(polygons, x, y, z1, x, y1, z, x, y, z)
    add_polygon(polygons, x, y, z1, x, y1, z1, x, y1, z)

    #top
    add_polygon(polygons, x, y, z1, x1, y, z, x1, y, z1)
    add_polygon(polygons, x, y, z1, x, y, z, x1, y, z)
    #bottom
    add_polygon(polygons, x, y1, z, x1, y1, z1, x1, y1, z)
    add_polygon(polygons, x, y1, z, x, y1, z1, x1, y1, z1)

def add_sphere(polygons, cx, cy, cz, r, step ):
    points = generate_sphere(cx, cy, cz, r, step)

    lat_start = 0
    lat_stop = step
    longt_start = 0
    longt_stop = step

    step+= 1
    for lat in range(lat_start, lat_stop):
        for longt in range(longt_start, longt_stop):

            p0 = lat * step + longt
            p1 = p0+1
            p2 = (p1+step) % (step * (step-1))
            p3 = (p0+step) % (step * (step-1))

            if longt != step - 2:
                add_polygon( polygons, points[p0][0],
                             points[p0][1],
                             points[p0][2],
                             points[p1][0],
                             points[p1][1],
                             points[p1][2],
                             points[p2][0],
                             points[p2][1],
                             points[p2][2])
            if longt != 0:
                add_polygon( polygons, points[p0][0],
                             points[p0][1],
                             points[p0][2],
                             points[p2][0],
                             points[p2][1],
                             points[p2][2],
                             points[p3][0],
                             points[p3][1],
                             points[p3][2])


def generate_sphere( cx, cy, cz, r, step ):
    points = []

    rot_start = 0
    rot_stop = step
    circ_start = 0
    circ_stop = step

    for rotation in range(rot_start, rot_stop):
        rot = rotation/float(step)
        for circle in range(circ_start, circ_stop+1):
            circ = circle/float(step)

            x = r * math.cos(math.pi * circ) + cx
            y = r * math.sin(math.pi * circ) * math.cos(2*math.pi * rot) + cy
            z = r * math.sin(math.pi * circ) * math.sin(2*math.pi * rot) + cz

            points.append([x, y, z])
            #print 'rotation: %d\tcircle%d'%(rotation, circle)
    return points

def add_torus(polygons, cx, cy, cz, r0, r1, step ):
    points = generate_torus(cx, cy, cz, r0, r1, step)

    lat_start = 0
    lat_stop = step
    longt_start = 0
    longt_stop = step

    for lat in range(lat_start, lat_stop):
        for longt in range(longt_start, longt_stop):

            p0 = lat * step + longt;
            if (longt == (step - 1)):
                p1 = p0 - longt;
            else:
                p1 = p0 + 1;
            p2 = (p1 + step) % (step * step);
            p3 = (p0 + step) % (step * step);

            add_polygon(polygons,
                        points[p0][0],
                        points[p0][1],
                        points[p0][2],
                        points[p3][0],
                        points[p3][1],
                        points[p3][2],
                        points[p2][0],
                        points[p2][1],
                        points[p2][2] )
            add_polygon(polygons,
                        points[p0][0],
                        points[p0][1],
                        points[p0][2],
                        points[p2][0],
                        points[p2][1],
                        points[p2][2],
                        points[p1][0],
                        points[p1][1],
                        points[p1][2] )


def generate_torus( cx, cy, cz, r0, r1, step ):
    points = []
    rot_start = 0
    rot_stop = step
    circ_start = 0
    circ_stop = step

    for rotation in range(rot_start, rot_stop):
        rot = rotation/float(step)
        for circle in range(circ_start, circ_stop):
            circ = circle/float(step)

            x = math.cos(2*math.pi * rot) * (r0 * math.cos(2*math.pi * circ) + r1) + cx;
            y = r0 * math.sin(2*math.pi * circ) + cy;
            z = -1*math.sin(2*math.pi * rot) * (r0 * math.cos(2*math.pi * circ) + r1) + cz;

            points.append([x, y, z])
    return points


def add_circle( points, cx, cy, cz, r, step ):
    x0 = r + cx
    y0 = cy
    i = 1

    while i <= step:
        t = float(i)/step
        x1 = r * math.cos(2*math.pi * t) + cx;
        y1 = r * math.sin(2*math.pi * t) + cy;

        add_edge(points, x0, y0, cz, x1, y1, cz)
        x0 = x1
        y0 = y1
        i+= 1

def add_curve( points, x0, y0, x1, y1, x2, y2, x3, y3, step, curve_type ):

    xcoefs = generate_curve_coefs(x0, x1, x2, x3, curve_type)[0]
    ycoefs = generate_curve_coefs(y0, y1, y2, y3, curve_type)[0]

    i = 1
    while i <= step:
        t = float(i)/step
        x = t * (t * (xcoefs[0] * t + xcoefs[1]) + xcoefs[2]) + xcoefs[3]
        y = t * (t * (ycoefs[0] * t + ycoefs[1]) + ycoefs[2]) + ycoefs[3]
        #x = xcoefs[0] * t*t*t + xcoefs[1] * t*t + xcoefs[2] * t + xcoefs[3]
        #y = ycoefs[0] * t*t*t + ycoefs[1] * t*t + ycoefs[2] * t + ycoefs[3]

        add_edge(points, x0, y0, 0, x, y, 0)
        x0 = x
        y0 = y
        i+= 1


def draw_lines( matrix, screen, zbuffer, color ):
    if len(matrix) < 2:
        print('Need at least 2 points to draw')
        return

    point = 0
    while point < len(matrix) - 1:
        draw_line( int(matrix[point][0]),
                   int(matrix[point][1]),
                   matrix[point][2],
                   int(matrix[point+1][0]),
                   int(matrix[point+1][1]),
                   matrix[point+1][2],
                   screen, zbuffer, color)
        point+= 2

def add_edge( matrix, x0, y0, z0, x1, y1, z1 ):
    add_point(matrix, x0, y0, z0)
    add_point(matrix, x1, y1, z1)

def add_point( matrix, x, y, z=0 ):
    matrix.append( [x, y, z, 1] )



def draw_line( x0, y0, z0, x1, y1, z1, screen, zbuffer, color ):

    #swap points if going right -> left
    if x0 > x1:
        xt = x0
        yt = y0
        zt = z0
        x0 = x1
        y0 = y1
        z0 = z1
        x1 = xt
        y1 = yt
        z1 = zt

    x = x0
    y = y0
    z = z0
    A = 2 * (y1 - y0)
    B = -2 * (x1 - x0)
    wide = False
    tall = False

    if ( abs(x1-x0) >= abs(y1 - y0) ): #octants 1/8
        wide = True
        loop_start = x
        loop_end = x1
        dx_east = dx_northeast = 1
        dy_east = 0
        d_east = A
        distance = x1 - x + 1
        if ( A > 0 ): #octant 1
            d = A + B/2
            dy_northeast = 1
            d_northeast = A + B
        else: #octant 8
            d = A - B/2
            dy_northeast = -1
            d_northeast = A - B

    else: #octants 2/7
        tall = True
        dx_east = 0
        dx_northeast = 1
        distance = abs(y1 - y) + 1
        if ( A > 0 ): #octant 2
            d = A/2 + B
            dy_east = dy_northeast = 1
            d_northeast = A + B
            d_east = B
            loop_start = y
            loop_end = y1
        else: #octant 7
            d = A/2 - B
            dy_east = dy_northeast = -1
            d_northeast = A - B
            d_east = -1 * B
            loop_start = y1
            loop_end = y

    dz = (z1 - z0) / distance if distance != 0 else 0

    while ( loop_start < loop_end ):
        plot( screen, zbuffer, color, x, y, z )
        if ( (wide and ((A > 0 and d > 0) or (A < 0 and d < 0))) or
             (tall and ((A > 0 and d < 0) or (A < 0 and d > 0 )))):

            x+= dx_northeast
            y+= dy_northeast
            d+= d_northeast
        else:
            x+= dx_east
            y+= dy_east
            d+= d_east
        z+= dz
        loop_start+= 1
    plot( screen, zbuffer, color, x, y, z )
