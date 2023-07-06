import os
import sys
import argparse

parser = argparse.ArgumentParser(prog='Paircounter', description='''
usage example: python paircounter.py --inputfile ADNP_ENCFF739AJO_peaks_id.bed_filt_freeblack_1.bed_int_merged_bedfiles_2.bed --currentpath /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/pairwise/ --processed_file_1 /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/pairwise/allouts_1/ --processed_file_2 /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/pairwise/allouts_2/ --processed_merged_2 /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/pairwise/merged_out/
''')
parser.add_argument('--inputfile', help='Name of the input file ', required=True)
parser.add_argument('--currentpath', help='current path of analysis',  required=True)
parser.add_argument('--processed_file_1', help='path of temporary files generated from samplesheet 1',   required=True)
parser.add_argument('--processed_file_2', help='path of temporary files generated from samplesheet 2',  required=True)
parser.add_argument('--processed_merged_2', help='path of temporary files generated from merging sample of samplesheet 2',  required=True)

def count_items_with_substring(lst, substring):
    count = 0
    for item in lst:
        if substring in item:
            count += 1
    return count

# Print help message if no arguments
if len(sys.argv)==1:
   parser.print_help(sys.stderr)
   sys.exit(1)

args = parser.parse_args()

CAP_1 = {}
CAPs_int = {}
elements2 = []
with open(''.join(args.processed_merged_2 + args.inputfile)) as myfile:

    myfile = myfile.readlines()
    myfile_count = len(myfile)
    for lines in myfile:
        lines = lines.strip().split('\t')
        elements2.append(''.join(lines[3] + "_" + lines[23]))
        elements2 = list(set(elements2))
        countresult = count_items_with_substring(elements2, lines[23])
        CAPs_int[lines[23]] = countresult
        CAP_1[lines[11]] = CAPs_int


for CAPf in CAP_1:
    for intCAPf in CAP_1[CAPf]:
        print(''.join(str(CAPf) + '\t' + str(intCAPf) + '\t' + str(CAP_1[CAPf][intCAPf])))
        # remove if output file already exist
        if os.path.exists(''.join(str(args.processed_merged_2) + CAPf + '_' + intCAPf + '.count.txt')):
            os.remove(''.join(str(args.processed_merged_2) + CAPf + '_' + intCAPf + '.count.txt'))
        with open(''.join(str(args.processed_merged_2) + CAPf + '_' + intCAPf + '.count.txt'), 'a') as myoutfile:
            myoutfile.write(''.join(str(CAPf) + '\t' + str(intCAPf) + '\t' + str(CAP_1[CAPf][intCAPf])+ '\n'))
            myoutfile.close()
    
