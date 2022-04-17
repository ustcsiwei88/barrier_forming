import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
import sys
# from shapely.geometry.polygon import Polygon

color = ['darkorange', 'lime', 'darkmagenta', 'steelblue']

if __name__ == '__main__':
    f = open('input.txt')
    bounds = [float(v) for v in f.readline().rstrip(' \n').split(' ')]
    # print(bounds)
    fig, ax = plt.subplots()

    num_sets = int(f.readline())
    
    # obstacle
    num_obj = int(f.readline())
    for i in range(num_obj):
        f.readline()
        # l = f.readline()
        arr = [float(v) for v in f.readline().rstrip(' \n').split(' ')]
        print(arr)
        p = []
        for j in range(len(arr) // 2):
            p.append([arr[j*2], arr[j*2+1]])
        p = Polygon(p, color = 'grey')
        ax.add_patch(p)

    for set_id in range(num_sets):
        num_obj = int(f.readline())
        for i in range(num_obj):
            poly_pts = int(f.readline())
            
            arr = [float(v) for v in f.readline().rstrip(' \n').split(' ')]
            if poly_pts == 1:
                ax.plot(arr[0], arr[1], '.', color = color[set_id], lw = 2, ms = 10)
                continue
            # print(arr)
            p = []
            for j in range(len(arr) // 2):
                p.append([arr[j*2], arr[j*2+1]])
            p = Polygon(p, color = color[set_id])
            ax.add_patch(p)


    if 0:
        f_seg = open('segments.txt')
        for l in f_seg.readlines():
            seg = list(map(float, l.split(' ')))
            # print(seg)
            ax.plot([seg[0], seg[2]], [seg[1], seg[3]], 'k.--', lw = 1)
    if 1:
        # f_cut = open('cuts.txt')
        f_cut = open('min_segment_cuts.txt')
        for l in f_cut.readlines():
            seg = list(map(float, l.split(' ')))
            ax.plot([seg[0], seg[2]], [seg[1], seg[3]], 'r.-', lw = 2)
        
    # fig.add_subplot(111)
    # ax = plt.axis()
    ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], 'k-')
    ax.axis('off')
    ax.set_xlim(bounds[0], bounds[2])
    ax.set_ylim(bounds[1], bounds[3])
    ax.set_aspect(1)
    plt.savefig('instance.png', pad_inches = 0.01)
    plt.show()