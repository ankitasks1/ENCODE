# Read in known complexes (filtered for the dataset)

# Go through each predicted complex and find complexes with shared members
# How accurate is the complex in terms of sensitivity and specificity

import sys

hits = sys.argv[1]
comps = sys.argv[2]

complexes = dict()

# Read in complexes as sets

with open(comps) as c:
    for x in c.readlines():
        x = x.rstrip()
        v = x.split('\t')
        prots = v[1].split(',')
        # print(v[0],set(prots))
        complexes[v[0]] = set(prots)

# print(complexes)
# print(len(complexes))

with open(hits) as h:
    for x in h.readlines():
        x = x.rstrip()
        prots = x.split('\t')
        protscap = prots[:2]
        # print(protscap)
        names = list()
        for i in protscap:
            # print(i)
            names.append(i)
        # print(names)
        #convert names(list) to set so that later on 'intersection' function can be applied which works with set

        prot_set = set(names)
        # print(prot_set)
        for y in complexes:
            inter = complexes[y].intersection(prot_set)
            # print(complexes[y], prot_set)
            # print(inter)
            if len(inter) > 1:
                print ("{}\t{}\t{}\t{}".format(complexes[y], prot_set, inter, 'known'))

