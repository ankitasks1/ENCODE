# Determine the number of linkses from the database which have at least two proteins from the ENCODE data
# Print out only the members of each links which are found in the ENCODE data

import os,sys


info_file = sys.argv[1]
links_file = sys.argv[2]
encode_file = sys.argv[3]

encode = set()

with open(encode_file) as e:
    for x in e.readlines():
        x = x.rstrip()
        v = x.split('\t')
        z = v[0].split('_')
        # print(z[0])
        encode.add(z[0])
# print(encode)

info_id = dict()

with open(info_file) as i:
    for m in i.readlines():
        m = m.rstrip()
        n = m.split('\t')
        # print(n)
        n_id = n[0].replace('9606.', '')
        info_id[n_id] = n[1]

# print(len(info_id))

if os.path.exists('protein_symbol_id_links.txt'):
    os.remove('protein_symbol_id_links.txt')

with open('protein_symbol_id_links.txt', 'a') as pl:
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
            prots = info_id[sideA], info_id[sideB]
            # print(prots)
            found = list()
            for p in prots:
                # print(p)
                if p in encode:
                    # print('yes',p)
                    # found is a list which will store elements of prots each time a loop runs, 
                    # so sometimes both elements will be present in encode data and sometimes only one. 
                    # Next we take only those elements where both sideA and sideB is present in encode data.
                    found.append(p)
            # print(found)
            if len(found) >= 2:
                l = "\t".join(found)
                # Pair of CAPs which are present in encode data but may not be predicted to interact by mcl
                print(l + '\t' + score)
    pl.close()



