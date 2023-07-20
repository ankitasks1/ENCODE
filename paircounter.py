import os
import sys
import argparse

parser = argparse.ArgumentParser(prog='Paircounter', description='''
usage example: python3 paircounter.py --inputfile ADNP_ENCFF739AJO_peaks_id.bed_filt_freeblack_1.bed_int_merged_bedfiles_2.bed --currentpath /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/pairwise/ --processed_file_1 /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/pairwise/allouts_1/ --processed_file_2 /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/pairwise/allouts_2/ --processed_merged_2 /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/pairwise/merged_out/ --outputfile intersected.counts.txt
''')
parser.add_argument('--inputfile', help='Name of the input file ', required=True)
parser.add_argument('--outputbase', help='Name of the output file ', required=False)
parser.add_argument('--currentpath', help='current path of analysis',  required=True)
parser.add_argument('--processed_file_1', help='path of temporary files generated from samplesheet 1',   required=True)
parser.add_argument('--processed_file_2', help='path of temporary files generated from samplesheet 2',  required=True)
parser.add_argument('--processed_merged_2', help='path of temporary files generated from merging sample of samplesheet 2',  required=True)

# Print help message if no arguments
if len(sys.argv)==1:
   parser.print_help(sys.stderr)
   sys.exit(1)

args = parser.parse_args()

with open(''.join(args.processed_merged_2 + args.inputfile)) as myfile:
    CAP_1 = {}
    CAPs_int = {}
    for lines in myfile.readlines():
        lines = lines.strip().split('\t')
        CAP_1[lines[11]] = {}
        if lines[23] not in CAPs_int:
            CAPs_int[lines[23]] = set() #create a set to store multiple peaks
        CAPs_int[lines[23]].add(lines[3])
        
        CAP_1[lines[11]] = CAPs_int

# if os.path.exists(''.join(str(args.processed_merged_2) + 'intersected.counts.txt')):
#     os.remove(''.join(str(args.processed_merged_2)  + 'intersected.counts.txt'))

for mainCAP in CAP_1:
    print(mainCAP)
    with open(''.join(str(args.processed_merged_2) + str(mainCAP) + '_' + args.outputbase), 'a') as myoutfile:
        for overlappedCAP in CAP_1[mainCAP]:
            print(mainCAP, overlappedCAP, len(CAP_1[mainCAP][overlappedCAP]))
            myoutfile.write(''.join(str(mainCAP)  + '\t' + str(overlappedCAP)  + '\t' + str(len(CAP_1[mainCAP][overlappedCAP]))+ '\n'))
        myoutfile.close()
    