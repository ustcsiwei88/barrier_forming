import os
import sys
import random
import numpy as np
from shapely.geometry import Polygon, Point

def run(set_num, each):
    of = open('input.txt', 'w')
    of.write('0.0 0.0 1.0 1.0\n')

    of.write(str(set_num) + '\n')
    ob = []
    sets = [[] for i in range(set_num+1)]

    cnt = {3:0, 4:0, 5:0, 6:0}
    polygons = []

    for i in range(80):
        j = random.randint(3, 6)
        f = open('data/poly_' + str(j) + '_' + str(cnt[j]) + '.txt')
        cnt[j] = cnt[j] + 1
        f.readline()
        f.readline()
        f.readline()
        p = []
        dx = random.uniform(0, 0.6)
        dy = random.uniform(0, 0.6)
        for _ in range(j):
            pt = list(map(float, f.readline().split(' ')))
            
            p.append((pt[0] * .4 + dx, pt[1] * .4 + dy))
        p = np.array(p)

        poly = Polygon(p)
        if np.sum([p[val][0] * p[(val + 1)%len(p)][1] - p[(val + 1)%len(p)][0] * p[val][1] for val in range(len(p))]) < 0:
            p = np.flip(p, axis = 1)
            poly = Polygon(p)
        inter = False
        for po in polygons:
            if po.intersects(poly):
                inter = True
        if not inter:
            polygons.append(poly)
        if len(polygons) >= each:
            break


    for i in range(1):
        # st = each * i
        of.write(str(each) + '\n')
        for j in range(each):
            poly = polygons[j]
            coords = poly.exterior.coords
            of.write(str(len(coords) - 1) + '\n')
            for ind, it in enumerate(coords):
                if ind == len(coords) - 1:
                    break
                of.write("%.4f %.4f "%(it[0], it[1]))
            of.write('\n')

    for i in range(set_num):
        of.write(str(each) + '\n')
        for j in range(each):
            found = False
            while not found:
                it = np.random.rand(2) * 0.9 + 0.05
                p = Polygon([it-(0.01,0.01), it+(0.01,-0.01), it+(0,0.01)])
                found = True
                for poly in polygons:
                    if p.intersects(poly):
                        found = False
                        break
                if found:
                    polygons.append(p)
            of.write('1\n%.4f %.4f\n'%(it[0], it[1]))

run(3, 3)
os.system("./min_seg_multi_sets")
