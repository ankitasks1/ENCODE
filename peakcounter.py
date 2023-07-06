import os
import sys
import argparse

parser = argparse.ArgumentParser(prog='Peakcounter', description='''
usage example: python peakcounter.py --inputfile ADNP_ENCFF739AJO_peaks_id.bed_filt_freeblack_1.bed --currentpath /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/pairwise/ --processed_file /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/pairwise/allouts_1/ --type cleaned
''')
parser.add_argument('--inputfile', help='Name of the input file ', required=True)
parser.add_argument('--currentpath', help='current path of analysis',  required=False)
parser.add_argument('--processed_file', help='path of temporary files generated from samplesheet',   required=False)
parser.add_argument('--location', help='location file ', required=True)
parser.add_argument('--type', help='Name of the type of file', required=True)

# Print help message if no arguments
if len(sys.argv)==1:
   parser.print_help(sys.stderr)
   sys.exit(1)

args = parser.parse_args()

CAPs_int = {}
if args.location == 'current_1':
    with open(''.join(args.currentpath + args.inputfile)) as myfile:
        myfile = myfile.readlines()
        myfile_count = len(myfile)
        CAPs_int[args.inputfile] = myfile_count
    
    # print(CAPs_int)
    for lines in CAPs_int:
        print(lines, CAPs_int[lines])
        with open(''.join(str(args.currentpath) + str(args.type) + '.count1.txt'), 'a') as myoutfile:
            myoutfile.write(''.join(str(lines) + '\t' + str(CAPs_int[lines]) + '\n'))
            myoutfile.close()

elif args.location == 'current_2':
    with open(''.join(args.currentpath + args.inputfile)) as myfile:
        myfile = myfile.readlines()
        myfile_count = len(myfile)
        CAPs_int[args.inputfile] = myfile_count

    # print(CAPs_int)
    for lines in CAPs_int:
        print(lines, CAPs_int[lines])
        
        with open(''.join(str(args.currentpath) + str(args.type) + '.count2.txt'), 'a') as myoutfile:
            myoutfile.write(''.join(str(lines) + '\t' + str(CAPs_int[lines]) + '\n'))
            myoutfile.close()

elif args.location == 'other':
    with open(''.join(args.processed_file + args.inputfile)) as myfile:
        myfile = myfile.readlines()
        myfile_count = len(myfile)
        CAPs_int[args.inputfile] = myfile_count

    # print(CAPs_int)
    for lines in CAPs_int:
        print(lines, CAPs_int[lines])
        
        with open(''.join(str(args.processed_file) + str(args.type) + '.count.txt'), 'a') as myoutfile:
            myoutfile.write(''.join(str(lines) + '\t' + str(CAPs_int[lines]) + '\n'))
            myoutfile.close()
