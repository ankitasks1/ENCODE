# Determine the number of linkses from the database which have at least two proteins from the ENCODE data
# Print out only the members of each links which are found in the ENCODE data

import os,sys


info_file = sys.argv[1]
links_file = sys.argv[2]

info_id = dict()

with open(info_file) as i:
    for m in i.readlines():
        m = m.rstrip()
        n = m.split('\t')
        # print(n)
        n_id = n[0].replace('9606.', '')
        info_id[n_id] = n[1]

# print(len(info_id))

if os.path.exists('protein_id_links_id_assigned.txt'):
    os.remove('protein_id_links_id_assigned.txt')

with open('protein_id_links_id_assigned.txt', 'a') as pl:
    with open(links_file) as c:
        for x in c.readlines()[1:]:
            x = x.rstrip()
            v = x.split(' ')
            # print(v)
            sideA = v[0].replace('9606.', '')
            sideB = v[1].replace('9606.', '')
            score = v[2]
            # print(info_id[sideA], info_id[sideB], score)
            pl.write(''.join(info_id[sideA] + '\t' + info_id[sideB] + '\t' +  score + '\n'))
    pl.close()


