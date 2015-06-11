#!/usr/bin/env python3

import sys


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage:\n ./cleanData.py fileName')
        exit(1)
    key = 0
    str2num = {}
    l_max = 0
    r_max = 0
    with open(sys.argv[1], "r") as fin:
        for line in fin:
            ls = line.split()
            # print str2num[ls[0]], str2num[ls[1]]
            for x in ls:
                if x not in str2num:
                    str2num[x] = key
                    key += 1
            l_max = max(l_max, str2num[ls[0]])
            r_max = max(r_max, str2num[ls[1]])

    print(l_max, r_max, file=sys.stderr)        
    
            

    with open(sys.argv[1], "r") as fin:
        for line in fin:
            ls = line.split()
            print(str2num[ls[0]], str2num[ls[1]])

     # N = 10
     # for i in range(N):
     #     for j in range(r_max):
     #         print(i, j)
     #         print(j, i)
            
            

    
