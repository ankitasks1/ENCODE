import os,sys

# prepare databases
print('corum')
os.system('python3 ./scripts/corum_complex_count_filter.py humanComplexes.txt ' + sys.argv[1] + ' > humanComplexes_filtered.txt')

print('humap2')
os.system('python3 ./scripts/humap2_complex_count_filter.py humap2_complexes_20200809.txt ' + sys.argv[1] + '  > humap2_complexes_20200809_filtered.txt')


# Run mcxload

os.system('mcxload -abc ' + sys.argv[1] + ' --stream-mirror -write-tab data.tab -o data.mci')

# Run mcl
i_val = [1.4, 2, 4, 8, 16, 20, 22]

for val in i_val:
    # print(val)
    if isinstance(val, float):
        print(str(int(val * 10)))
        os.system('mcl data.mci -I ' + str(val))
        os.system('mcxdump -icl out.data.mci.I' + str(int(val * 10)) + ' -tabr data.tab -o dump.data.mci.I' +  str(int(val * 10)))
        os.system('python3 ./scripts/mod_benchmark_complexes.py dump.data.mci.I' + str(int(val * 10)) + ' humanComplexes_filtered.txt > interactions.mci.CAPs_K562_corum_' + str(int(val * 10)) + '.txt')
        os.system('python3 ./scripts/mod_benchmark_complexes.py dump.data.mci.I' + str(int(val * 10)) + ' humap2_complexes_20200809_filtered.txt > interactions.mci.CAPs_K562_humap2_' + str(int(val * 10)) + '.txt')
    else:
        print(str(val * 10))
        os.system('mcl data.mci -I ' + str(val))
        os.system('mcxdump -icl out.data.mci.I' +  str(val * 10) + ' -tabr data.tab -o dump.data.mci.I' +  str(val * 10))
        os.system('python3 ./scripts/mod_benchmark_complexes.py dump.data.mci.I' + str(val * 10) + ' humanComplexes_filtered.txt > interactions.mci.CAPs_K562_corum_' + str(val * 10) + '.txt')
        os.system('python3 ./scripts/mod_benchmark_complexes.py dump.data.mci.I' + str(val * 10) + ' humap2_complexes_20200809_filtered.txt > interactions.mci.CAPs_K562_humap2_' + str(val * 10) + '.txt')
