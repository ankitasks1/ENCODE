import os
import sys
import argparse

parser = argparse.ArgumentParser(prog='Peakcounter', description='''
usage example: python3 get_peaks_summit_bed.py --inputfile ADNP_ENCFF739AJO_peaks_id.bed_filt_freeblack_1.bed --currentpath /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/summit_pairwise/ --processed_file /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/summit_pairwise/allouts_1/ --location processed_1
''')
parser.add_argument('--inputfile', help='Name of the input file ', required=True)
parser.add_argument('--currentpath', help='current path of analysis',  required=False)
parser.add_argument('--processed_file', help='path of temporary files generated from samplesheet',   required=False)
parser.add_argument('--location', help='location file ', required=True)

# Print help message if no arguments
if len(sys.argv)==1:
   parser.print_help(sys.stderr)
   sys.exit(1)

args = parser.parse_args()
if args.location == 'processed_1':
    bedfilename = args.inputfile.split('_filt_freeblack_1.bed')[0]

    if os.path.exists(''.join(bedfilename + '_summit_freeblack_1.bed')):
        os.remove(''.join(bedfilename + '_summit_freeblack_1.bed'))
    with open(''.join(args.processed_file + args.inputfile)) as myfile:
        with open(''.join(str(args.processed_file) + bedfilename + '_summit_freeblack_1.bed'), 'a') as mysummitbed:   
            for lines in myfile.readlines():
                lines = lines.strip().split('\t')
                start = str(int(float(lines[1])+((int(lines[2])-int(lines[1]))/2)))
                end = str(int((float(lines[2])-((int(lines[2])-int(lines[1]))/2)) + float(1)))
                summits = ''.join(lines[0] + '\t' + start + '\t' + end + '\t' + lines[3] + '\t' + lines[4] + '\t' + lines[5] + '\t' + lines[6] + '\t' + lines[7] + '\t' + lines[8] + '\t' + lines[9] + '\t' + lines[10] + '\t' + lines[11] + '\n')
                # print(summits)
                mysummitbed.write(summits)
            mysummitbed.close()
