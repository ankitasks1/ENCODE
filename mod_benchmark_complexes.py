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


with open(hits) as h:
    for x in h.readlines():
        x = x.rstrip()
        prots = x.split('\t')
        # print(prots)
        names = list()
        for i in prots:
            # print(i)
            names.append(i)

        prot_set = set(names)
        for y in complexes:
            inter = complexes[y].intersection(prot_set)
            # print(complexes[y], prot_set)
            if len(inter) > 1:
                #print ("{}\t{}".format(y, inter))
                # sens = number of known complex members found / total complex members 
                sens = len(inter) / len(complexes[y])
                # spec = number of correct complex members / size of predicted complex
                spec = len(inter) / len(prot_set)
                print("{}\t{}\t{}\t{}\t{}\t{}".format(y, complexes[y], prot_set, inter, sens, spec))