#!/usr/bin/env python

import sys
from sets import Set


def q_gram(s, q=2):
    ls = [s[i:-q + i + 1] for i in xrange(q-1)]
    ls.append(s[q-1:])

    return Set(["".join(x) for x in apply(zip, ls)])


if __name__ == "__main__":
    q = 2
    if len(sys.argv) > 2:
        q = int(sys.argv[2])

    str2num = {}

    id = 0
    with open(sys.argv[1]) as fin:
        for line in fin:
            ls = line.split(',')
            ls = q_gram(ls[1], q)
            for x in ls:
                if x not in str2num:
                    str2num[x] = id
                    id += 1


    with open(sys.argv[1]) as fin:
        for line in fin:
            ls = line.split(',')
            node = int(ls[0]) - 1
            ls = q_gram(ls[1], q)
            for x in ls:
                print node, str2num[x]






                    



