import os 
import sys
import random

def run(set_num, each):
    f = open('input.txt', 'w')
    f.write('0.0 0.0 1.0 1.0\n')

    # set_num = 3

    f.write(str(set_num) + '\n')
    ob = []
    sets = [[] for i in range(set_num+1)]

    roll = []
    for i in range(7):
        for j in range(7):
            roll.append((i,j))

    random.shuffle(roll)
    random.shuffle(roll)

    for i in range(set_num + 1):
        for j in range(i*each, i*each + each):
            sets[i].append(roll[j])

                

    for s in sets:
        f.write(str(len(s)) + '\n')
        for it in s:
            f.write("4\n")
            for dx, dy in [(-1, -1), (1, -1), (1, 1), (-1, 1)]:
                f.write("%.4f %.4f"%(it[0] * 1./7 + 1./14 + dx * 1./28,  it[1] * 1./7 + 1./14 + dy * 1./28) + ' ')
            f.write('\n')


run(3, 3)
os.system("./min_seg_multi_sets")