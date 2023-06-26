import os
import sys
import pandas as pd
from sortedcontainers import SortedDict
from pybedtools import BedTool
import argparse
import time

print('''
----------------------------------------------------------------------------------------------------
|                                        Peak Pairwiser                                            |
|                                      Author: Ankit Verma                                         |
|                                  Contact: ankitverma9079@gmail.com                               |
| Note: This code will work only with MACS2 generated peaks which are modified  with 12 columns    |
----------------------------------------------------------------------------------------------------
''')

print('\nAnalysis started')

filePath = os.getcwd()
print('\nOutput will be written to: ' + filePath)


parser = argparse.ArgumentParser(prog='peakpairwiser', description='''
Description of usage: Overlap ChIP-Seq peaks with genomic bins. The peaks should be obtained from MACS2 peak caller.
Peak file should have chr at col1 , start at col2, end at col3, MACS score at col7. 
This format you generally obtain from MACS2 callpeak.
usage example: python makebins_from_peaks_q_pairwise.py --input1 s_samples_list.txt --input2 s_samples_list.txt --blacklist dummy_blacklist.bed --output pairiwse_intersection_out.txt
''')
parser.add_argument('--input1', help='Name of the input file 1 containing all the *peaks.bed', required=True)
parser.add_argument('--input2', help='Name of the input file 2 containing all the *peaks.bed', required=True)
parser.add_argument('--blacklist', help='Name of the blacklist regions file', required=True)
parser.add_argument('--output', help='BaseName of the output  file', default='_bin_out.txt', required=True)
parser.add_argument('--evalue', help='E-value for peak called, column 9', type=float, default=1.3, required=False)
parser.add_argument('--peakscore', help='MACS score, column 7, mandatory --scorebased along with this option', type=float, default=0, required=False)
parser.add_argument('--scorebased', help='Required to filter based on MACS score, column 7', action='store_true', required=False)
parser.add_argument('--path', help='current path of analysis',  default=filePath, required=False)


# Print help message if no arguments
if len(sys.argv)==1:
   parser.print_help(sys.stderr)
   sys.exit(1)

args = parser.parse_args()

# read peak file list
def readpeaklist_1(inputpeaklist_1):
    samplelista_1 = open(''.join(filePath + '/' + inputpeaklist_1), "r")
    return list(samplelista_1)

def readpeaklist_2(inputpeaklist_2):
    samplelista_2 = open(''.join(filePath + '/' + inputpeaklist_2), "r")
    return list(samplelista_2)


# read blacklist regions using pybedtools
def blacklistregions(blacklist):
    blacklista = open(''.join(filePath + '/' +  blacklist).strip())
    blacklista = blacklista.read().strip().split('\n')
    blacklist_name = "blacklist"
    blacklist_bed = BedTool(blacklista)
    return blacklist_name, blacklist_bed

# read bedfile inside peaklist using pybedtools
def peakspairwiser(peaklist_1, peaklist_2, blacklist_bed):
    serialnum = 0
    for list_1_peaks in peaklist_1:
        list_1_peaks = list_1_peaks.split('\n')[0]
        for list_2_peaks in peaklist_2:
            list_2_peaks = list_2_peaks.split('\n')[0]
            # print(''.join(list_1_peaks + '   v/s    ' +list_2_peaks))
            # Create dictionary to store peakfile1 and peakfile2
            peakfiledict_1 = {}
            peakfiledict_2 = {}
            overlaps = {}
            counts = {}
            """ Read the Peakfile_1"""
            peakfile_1 = open(list_1_peaks, 'r')
            peakfile_1_name = os.path.basename(peakfile_1.name).replace("_peaks_id.bed", "")
            peaksbed_1 = BedTool(list_1_peaks)
            # print(peaksbed_1)
            # filter the object by column 9 (8 index python) (q-value > 1.3 = q value 0.05
            # get the number of intervals in the filtered object
            if args.scorebased:
                peaksbed_1_filtered = peaksbed_1.filter(lambda x: float(x[7]) > args.peakscore)
            else:
                peaksbed_1_filtered = peaksbed_1.filter(lambda x: float(x[9]) > args.evalue)
            """ Intersect the peaksbed_1_filtered with the blacklisted regions and obtain non overlapping peaks """
            peaksbed_1_free_of_blacklist = peaksbed_1_filtered.intersect(blacklist_bed, v=True)
            peaksbed_1_filtered_count = peaksbed_1_free_of_blacklist.count()
            peakfiledict_1[peakfile_1_name] = peaksbed_1_free_of_blacklist, peaksbed_1_filtered_count
            # print(peakfile_1_name)
            counts[peakfile_1_name] = peaksbed_1_filtered_count
            """ Read the Peakfile_2"""
            peakfile_2 = open(list_2_peaks, 'r')
            peakfile_2_name = os.path.basename(peakfile_2.name).replace("_peaks_id.bed", "")
            peaksbed_2 = BedTool(list_2_peaks)
            # print(peakfile_2_name)
            # filter the object by column 9 (8 index python) (q-value > 1.3 = q value 0.05
            # get the number of intervals in the filtered object
            if args.scorebased:
                peaksbed_2_filtered = peaksbed_2.filter(lambda x: float(x[7]) > args.peakscore)
            else:
                peaksbed_2_filtered = peaksbed_2.filter(lambda x: float(x[9]) > args.evalue)
            """ Intersect the peaksbed_2_filtered with the blacklisted regions and obtain non overlapping peaks """
            peaksbed_2_free_of_blacklist = peaksbed_2_filtered.intersect(blacklist_bed, v=True)
            peaksbed_2_filtered_count = peaksbed_2_free_of_blacklist.count()
            peakfiledict_2[peakfile_2_name] = peaksbed_2_free_of_blacklist, peaksbed_2_filtered_count
            # print(peakfile_2_name)
            counts[peakfile_2_name] = peaksbed_2_filtered_count

            #Intersect the peakfile_1 and peakfile_2
            intersect_file = peaksbed_1_free_of_blacklist.intersect(peaksbed_2_free_of_blacklist, wa=True, wb=True)
            # print(intersect_file)
            """ Create Pairwise dictionaries """
            for intcont in intersect_file:
                # first peak file ids dictionaries
                if intcont[11] not in overlaps:
                    overlaps[intcont[11]] = {} #create an empty dictionary for CAPs from file 1
                # dictionaries of second peaks ids intersected with first peaks
                if intcont[23] not in overlaps[intcont[11]]:
                    overlaps[intcont[11]][intcont[23]] = set() #create a set to store multiple coordinates
                # add peak number of peaks from first peak file which are intersected with second peak file
                overlaps[intcont[11]][intcont[23]].add(intcont[3])

                # second peaks ids dictionaries
                if intcont[23] not in overlaps:
                    overlaps[intcont[23]] = {} #create an empty dictionary for CAPs from file 2
                # dictionaries of first peaks ids intersected with second peaks
                if intcont[11] not in overlaps[intcont[23]]:
                    overlaps[intcont[23]][intcont[11]] = set() #create a set to store multiple coordinates
                overlaps[intcont[23]][intcont[11]].add(intcont[15])
            # print(len(overlaps))
            # Note CAPs which have multiple replicate can be deduplicated later on by untagging ENCODE ID.

            if len(overlaps) == 0:
                serialnum += 1
                print(''.join(str(serialnum) + "\t" + peakfile_1_name + "\t" + peakfile_2_name + "\t" + str(len(peaksbed_1_free_of_blacklist)) + "\t" + str(len(peaksbed_2_free_of_blacklist)) + "\t" + "0" + "\t" + "0"))
                serialnum += 1
                print(''.join(str(serialnum) + "\t" + peakfile_2_name + "\t" + peakfile_1_name + "\t" + str(len(peaksbed_2_free_of_blacklist)) + "\t" + str(len(peaksbed_1_free_of_blacklist)) + "\t" + "0" + "\t" + "0"))

                with open(args.output, "a") as myfinalout:
                    myfinalout.write(''.join(peakfile_1_name + "\t" + peakfile_2_name + "\t" + str(len(peaksbed_1_free_of_blacklist)) + "\t" + str(len(peaksbed_2_free_of_blacklist)) + "\t" + "0" + "\t" + "0" + "\n" + peakfile_2_name + "\t" + peakfile_1_name + "\t" + str(len(peaksbed_2_free_of_blacklist)) + "\t" + str(len(peaksbed_1_free_of_blacklist)) + "\t" + "0" + "\t" + "0" + "\n"))
                    myfinalout.close()
            else:
                for i in overlaps:
                    for j in overlaps[i]:
                        serialnum += 1
                        print(''.join(str(serialnum) + "\t" + i + "\t" + j + "\t" + str(counts[i]) + "\t" + str(counts[j]) + "\t" + str(len(overlaps[i][j])) + "\t" + str(float(int(len(overlaps[i][j]))) / float(int(counts[i])))))
                        with open(args.output, "a") as myfinalout:
                            myfinalout.write(''.join(i + "\t" + j + "\t" + str(counts[i]) + "\t" + str(counts[j]) + "\t" + str(len(overlaps[i][j])) + "\t" + str(float(int(len(overlaps[i][j]))) / float(int(counts[i]))) + "\n"))
                            myfinalout.close()


# remove if output file already exist
if os.path.exists(args.output):
    os.remove(args.output)
print('\n----> Reading peaklists,\n Peaklists = ' + str(args.input1) + '   ' + str(args.input2) + '\n')
peaklistout_1 = readpeaklist_1(args.input1)
peaklistout_2 = readpeaklist_2(args.input2)
print(f"Input lists are: {peaklistout_1 + peaklistout_2}")

print('\n----> Reading blacklist regions,\n Blacklist = ' + args.blacklist + '\n')
outputblacklist = blacklistregions(args.blacklist)
# print(outputblacklist)

# print(peakspairwiser(peaklistout_1, peaklistout_2, outputblacklist[1]))
if args.scorebased:
    print('\n----> Cleaning and filtering peaksfile\n Score-based = ' + str(args.peakscore) + '\n')
    print('\n----> Doing pairwise intersection for peakfiles\n')
    pybedtoolsout = peakspairwiser(peaklistout_1, peaklistout_2, outputblacklist[1])
else:
    print('\n----> Cleaning and filtering peaksfile\n E-value-based = ' + str(args.evalue) + '\n')
    print('\n----> Doing pairwise intersection for peakfiles\n')
    pybedtoolsout = peakspairwiser(peaklistout_1, peaklistout_2, outputblacklist[1])

print(f"\nAnalysis finished\n")





