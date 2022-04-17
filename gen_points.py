import sys
import random
import numpy as np
import os
from shapely.geometry import Polygon


def run(set_num, each):
    of = open('input.txt', 'w')
    of.write('0.0 0.0 1.0 1.0\n')


    of.write(str(set_num) + '\n')
    ob = []
    sets = [[] for i in range(set_num+1)]
    polys = []
    of.write('0\n')

    for i in range(set_num):
        of.write(str(each) + '\n')
        for j in range(each):
            inter = True
            while inter:
                it = np.random.rand(2) * 0.9 + 0.05
                inter = False
                poly = Polygon([it-(0.01,0.01), it+(0.01,-0.01), it+(0,0.01)])
                for p in polys:
                    if p.intersects(poly):
                        inter = True
                        break
            polys.append(poly)
            of.write("1\n%.4f %.4f\n"%(it[0], it[1]))

run(3, 3)
os.system("./min_seg_multi_sets")