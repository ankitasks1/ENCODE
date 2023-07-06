import os
import sys
import argparse
import glob

parser = argparse.ArgumentParser(prog='cominecounter', description='''
usage example: python combinecounts.py  --currentpath /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/pairwise/ --processed_file_1 /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/pairwise/allouts_1/ --processed_file_2 /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/pairwise/allouts_2/ --processed_merged_2 /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/pairwise/merged_out/ --output testfinalout.txt
''')
parser.add_argument('--currentpath', help='current path of analysis',  required=True)
parser.add_argument('--processed_file_1', help='path of temporary files generated from samplesheet 1',   required=True)
parser.add_argument('--processed_file_2', help='path of temporary files generated from samplesheet 2',  required=True)
parser.add_argument('--processed_merged_2', help='path of temporary files generated from merging sample of samplesheet 2',  required=True)
parser.add_argument('--output', help='path of output files',  required=True)


# Print help message if no arguments
if len(sys.argv)==1:
   parser.print_help(sys.stderr)
   sys.exit(1)

args = parser.parse_args()

Intersect = {}
Original = {}
Cleaned = {}

for files in glob.glob(''.join(args.processed_merged_2 + '*.count.txt')):
    with open(files) as myfile:
        # myfile_id = myfile.name.replace(args.processed_merged_2, "").replace('.count.txt', "")
        # myfile_id = myfile_id.split('_')
        myfile = myfile.readlines()
        pairs = myfile[0].strip().split('\t')
        secondpair = {}
        secondpair[pairs[1]] = pairs[2]
        Intersect[pairs[0]] = secondpair

print(Intersect)
with open(''.join(args.currentpath + 'original_1.count1.txt')) as original:
    original = original.readlines()
    original = set(original)
    for elements in original:
        elements = elements.strip().split('\t')
        Original[elements[0].replace('_peaks_id.bed', '')] = elements[1]

with open(''.join(args.currentpath + 'original_2.count2.txt')) as original:
    original = original.readlines()
    original = set(original)
    for elements in original:
        elements = elements.strip().split('\t')
        Original[elements[0].replace('_peaks_id.bed', '')] = elements[1]

with open(''.join(args.processed_file_1 + 'cleaned_1.count.txt')) as cleaned:
    cleaned = cleaned.readlines()
    cleaned = set(cleaned)
    for elements in cleaned:
        elements = elements.strip().split('\t')
        Cleaned[elements[0].replace('_peaks_id.bed_filt_freeblack_1.bed', '')] = elements[1]

with open(''.join(args.processed_file_2 + 'cleaned_2.count.txt')) as cleaned:
    cleaned = cleaned.readlines()
    cleaned = set(cleaned)
    for elements in cleaned:
        elements = elements.strip().split('\t')
        Cleaned[elements[0].replace('_peaks_id.bed_filt_freeblack_2.bed', '')] = elements[1]

for values in Intersect:
    for values2 in Intersect[values]:
        # print(values, values2, Original[values], Cleaned[values], Original[values2], Cleaned[values2], Intersect[values][values2])
        with open(''.join(str(args.processed_merged_2) + args.output), 'a') as myoutfile:
            myoutfile.write(''.join(values + '\t' + values2 + '\t' + Original[values] + '\t' + Original[values2] + '\t' + Cleaned[values] + '\t' + Cleaned[values2] + '\t' + Intersect[values][values2] + '\n'))
            myoutfile.close()