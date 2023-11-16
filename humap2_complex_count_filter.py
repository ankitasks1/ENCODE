# Determine the number of complexes from the database which have at least two proteins from the ENCODE data
# Print out only the members of each complex which are found in the ENCODE data

import sys

complex_file = sys.argv[1]
encode_file = sys.argv[2]

encode = set()

with open(encode_file) as e:
    for x in e.readlines():
        x = x.rstrip()
        v = x.split('\t')
        z = v[0].split('_')
        # print(z[0])
        encode.add(z[0])

# print(encode)
with open(complex_file) as c:
    for x in c.readlines():
        x = x.rstrip()
        v = x.split(',')
        # print(''.join(v[1] + '\t' +  v[18]))
        prots = v[3].split(' ')
        # print((prots))
        found = list()
        for p in prots:
            # print(p)
            if p in encode:
                # print('yes',p)
                found.append(p)
        # print(found)
        if len(found) >= 2:
            l = ",".join(found)
            print(v[0]+'\t'+l)
            #print(v[1], found, v[18]) 
