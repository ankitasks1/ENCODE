import os
import sys
import argparse
from pybedtools import BedTool

parser = argparse.ArgumentParser(prog='Paircounter', description='''
usage example: python3 summit_distancecounter.py --inputfile_1 ZNF589_ENCFF345IHK_peaks_id.bed_summit_freeblackhot_sorted_merged_1.bed --inputfile_2 ZNF592_ENCFF847QJI_peaks_id.bed_summit_freeblackhot_sorted_merged_2.bed  --currentpath /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/summit_pairwise/ --processed_file_1 /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/summit_pairwise/allouts_1/ --processed_file_2 /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/summit_pairwise/allouts_2/ --processed_merged_2 /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/summit_pairwise/merged_out/ --outputbase intersected.counts.txt
''')

parser.add_argument('--inputfile_1', help='Name of the input file 1', required=True)
parser.add_argument('--inputfile_2', help='Name of the input file 2', required=True)
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

with open(''.join(args.inputfile_1), 'r') as mysummitfile_1:
#    print(mysummitfile_1)
   peakfile_1_name = os.path.basename(mysummitfile_1.name).split('_')[0]
#    print(peakfile_1_name)
   peaksbed_1 = BedTool(mysummitfile_1)
   with open(''.join(args.inputfile_2), 'r') as mysummitfile_2:
    # print(mysummitfile_2)
    peakfile_2_name = os.path.basename(mysummitfile_2.name).split('_')[0]
    # print(peakfile_2_name)
    peaksbed_2 = BedTool(mysummitfile_2)

    # #Distance between the peakfile_1 and peakfile_2 summits
    peaksbed_1_and_2_dis = peaksbed_1.closest(peaksbed_2, d=True)
    # print(peaksbed_1_and_2_dis)
    CAP_1 = {}
    CAPs_int = {}
    for lines in peaksbed_1_and_2_dis:
    #   print(lines[0])
      values = float(lines[8]) 
      # As recommended the CAPs will potentially bind together if they have summit within 150 bp
      if 0.0 <= values <= 150.0:
            # print(lines[0:9])
            peakID = '_'.join(lines[0:3])
            # print(peakID)
            CAP_ID = lines[3].split('_')[0]
            CAP_1[CAP_ID] = {}
            CAP_int_ID = lines[7].split('_')[0]
            if CAP_int_ID not in CAPs_int:
                CAPs_int[CAP_int_ID] = set() #create a set to store multiple peaks
            CAPs_int[CAP_int_ID].add(peakID)
            CAP_1[CAP_ID] = CAPs_int
# print(CAP_1)


# if os.path.exists(''.join(str(args.processed_merged_2) + 'intersected.counts.txt')):
#     os.remove(''.join(str(args.processed_merged_2)  + 'intersected.counts.txt'))

for mainCAP in CAP_1:
    print(mainCAP)
    with open(''.join(str(args.processed_merged_2) + str(mainCAP) + '_' + args.outputbase), 'a') as myoutfile:
        for overlappedCAP in CAP_1[mainCAP]:
            print(mainCAP, overlappedCAP, len(CAP_1[mainCAP][overlappedCAP]))
            myoutfile.write(''.join(str(mainCAP)  + '\t' + str(overlappedCAP)  + '\t' + str(len(CAP_1[mainCAP][overlappedCAP]))+ '\n'))
        myoutfile.close()
