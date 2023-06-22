import os
import sys
import pandas as pd
from sortedcontainers import SortedDict
from pybedtools import BedTool
import argparse
import time

print('''
---------------------------------------------------------------------------------
|                           Peak and Genomic-sites marker                        |
|                                Author: Ankit Verma                             |
|                        Contact: ankitverma9079@gmail.com                       |
| Note: This code will work only with MACS2 generated peak files with 10 columns |
---------------------------------------------------------------------------------
''')

print('\nAnalysis started')

filePath = os.getcwd()

parser = argparse.ArgumentParser(prog='peakmarker', description='''
Description of usage: Overlap ChIP-Seq peaks with genomic bins. The peaks should be obtained from MACS2 peak caller.
Peak file should have chr at col1 , start at col2, end at col3, MACS score at col7. 
This format you generally obtain from MACS2 callpeak.
''')
parser.add_argument('--input', help='Name of the input file containing all the *peaks.bed')
parser.add_argument('--blacklist', help='Name of the blacklist regions file')
parser.add_argument('--output', help='BaseName of the output  file', default='_bin_out.txt')
parser.add_argument('--gbinsize', help='Genomic bin size',  default=1000)
parser.add_argument('--evalue', help='E-value for peak called, column 9', type=float, default=1.3)
parser.add_argument('--peakscore', help='MACS score, column 7, mandatory --scorebased along with this option', type=float, default=0)
parser.add_argument('--scorebased', help='Required to filter based on MACS score, column 7', action='store_true')
parser.add_argument('--filepath', help='current path of analysis',  default=filePath)

# Print help message if no arguments
if len(sys.argv)==1:
   parser.print_help(sys.stderr)
   sys.exit(1)

args = parser.parse_args()

print('\nOutput will be written to: ' + args.filepath)

# read peak file list
def readpeaklist(inputpeaklist):
    samplelista = open(''.join(args.filepath + '/' + inputpeaklist).strip())
    samplelista = samplelista.read().strip().split('\n')
    return samplelista

# read blacklist regions using pybedtools
def blacklistregions(blacklist):
    blacklista = open(''.join(args.filepath + '/' +  blacklist).strip())
    blacklista = blacklista.read().strip().split('\n')
    blacklist_name = "blacklist"
    blacklist_bed = BedTool(blacklista)
    return blacklist_name, blacklist_bed


# read bedfile inside peaklist using pybedtools
def readpeaksfile(eachpeakfiles, blacklist_bed):
    pybedtooldict = {}
    for samples in eachpeakfiles:
        print(''.join(args.filepath +  samples))
        peakfile = open(''.join(args.filepath +  samples), 'r')
        peakfile_name = os.path.basename(peakfile.name).split("_")[0]
        peaks = BedTool(''.join(args.filepath + samples))
        # print(peaks)
        # filter the object by column 9 (8 index python) (q-value > 1.3 = q value 0.05
        # get the number of intervals in the filtered object
        if args.scorebased:
            peaks_filtered = peaks.filter(lambda x: float(x[6]) > args.peakscore)
        else:
            peaks_filtered = peaks.filter(lambda x: float(x[8]) > args.evalue)
        # Intersect the filtered peaks with the blacklisted regions and obtain non overlapping peaks
        peaks_free_of_blacklist = peaks_filtered.intersect(blacklist_bed, v=True)
        peaks_filtered_count = peaks_free_of_blacklist.count()
        # print(samples, peaks_filtered_count)
        if peakfile_name not in pybedtooldict:
            pybedtooldict[peakfile_name] = peaks_free_of_blacklist, peaks_filtered_count
            # print(peakfile_name, "unique name")
            # print(pybedtooldict[peakfile_name][0])

        elif peakfile_name in pybedtooldict:
            tempdict = {}
            tempdict[peakfile_name] = peaks_free_of_blacklist
            # print(tempdict[peakfile_name])
            # print(peakfile_name, "repeated name found")
            # Concatenate peaks if the same CAps is present
            # print(pybedtooldict[peakfile_name][0].cat(tempdict[peakfile_name]))
            concatenated_peaks_free_of_blacklist = pybedtooldict[peakfile_name][0].cat(tempdict[peakfile_name])
            concatenated_peaks_filtered_count = concatenated_peaks_free_of_blacklist.count()
            pybedtooldict[peakfile_name] = concatenated_peaks_free_of_blacklist, concatenated_peaks_filtered_count

    return pybedtooldict

# create bins as per the bin size and store in dictionary
def quantifybins(pybedtooldict):
    fulldict = {}
    for peakfile_name in pybedtooldict:
        print("Performing analysis for: ", peakfile_name)
        # print(pybedtooldict[peakfile_name][0])
        bedfiles_cleaned = pybedtooldict[peakfile_name][0]
        peaks_filtered_count = pybedtooldict[peakfile_name][1]
        # print(peaks_filtered_count)
        # Create and empty dictionary to store CAPs
        bindict = {}
        if peakfile_name not in bindict:
            # Create and empty dictionary within the dictionary to store peaks for each CAPs
            bindict[peakfile_name] = {}
            # bindict[peakfile_name][peakfile_name] = 1
            # print(bindict)
            # Apply condition if the file is not empty after intersection
            if peaks_filtered_count > 0:
                # print("Filtered Bedtools object is not empty")
                # Create Bins and store in dictionary with key value 1
                for peak_lines in bedfiles_cleaned:
                    # print(peak_lines)
                    chrom = peak_lines[0]
                    # Divide by binsize
                    newstart = (str(int(peak_lines[1]) / int(args.gbinsize))).split('.')[0]
                    # print(newstart)
                    newend = (str(int(peak_lines[2]) / int(args.gbinsize))).split('.')[0]
                    # print(chrom, newstart, newend)
                    # If bin is equal, (not short:No) that is bin/1kb Start = End (both start and end fall in the same genomic bin)
                    if int(newstart) == int(newend):
                        # print(chrom,"Start:",peak_lines[1], "End:",peak_lines[2], "BinStart:",newstart, "BinEnd: ",newend, "No")
                        # Store Start
                        binkey = ''.join(chrom + '_' + newstart)
                        bindict[peakfile_name][binkey] = 1
                        # Store End
                        binkey = ''.join(chrom + '_' + newend)
                        bindict[peakfile_name][binkey] = 1
                        # print(bindict)
            #         # If bin is short (:Yes) that is bin/1kb Start < End
                    elif (int(newend) - int(newstart)) >= 1:
                        # print(peakfile_name, "Start:",lines[1], "End:",lines[2], "BinStart:",newstart, "BinEnd: ",newend, "Yes")
                        # print(int(newend) - int(newstart))
                        for hiddenbins in range(int(newstart), int(newend), 1):
                            # print(hiddenbins)
                            # Store Start
                            binkey = ''.join(chrom + '_' + newstart)
                            bindict[peakfile_name][binkey] = 1
                            # Store Bins in-between start and end
                            binkey = ''.join(chrom + '_' + str(hiddenbins))
                            bindict[peakfile_name][binkey] = 1
                            # Store End
                            binkey = ''.join(chrom + '_' + newend)
                            bindict[peakfile_name][binkey] = 1
                            # print(bindict)
            # just testing bindict[peakfile_name][peakfile_name] = 1

        # print(bindict)
        # print(bindict)
        # The below script will open args.output with every CAPs and write the output into it (preferred way)
        for i in SortedDict(bindict):
            for j in SortedDict(bindict[i]):
                # print(i, j)
                print(''.join(i + '\t' + j + '\t' + str(bindict[i][j])))
                with open(args.output, 'a') as myfile:
                    myfile.write(''.join(i + '\t' + j + '\t' + str(bindict[i][j])) + '\n')
                    myfile.close()
        fulldict.update(bindict)

    # print(fulldict)
    # The below script will open args.output only one time and write the output of fulldict into it (memory intensive way, not sure if it will work)
    # for i in SortedDict(fulldict):
    #     for j in SortedDict(fulldict[i]):
    #         # print(i, j)
    #         print(''.join(i + '\t' + j + '\t' + str(fulldict[i][j])))
    #         with open(args.output, 'a') as myfile:
    #             myfile.write(''.join(i + '\t' + j + '\t' + str(fulldict[i][j])) + '\n')
    #             myfile.close()


# remove if output file already exist
if os.path.exists(args.output):
    os.remove(args.output)
print('\n----> Reading peaklist,\n Peaklist = ' + args.input + '\n')
peaklistout = readpeaklist(args.input)

print('\n----> Reading blacklist regions,\n Blacklist = ' + args.blacklist + '\n')
outputblacklist = blacklistregions(args.blacklist)
# print(readpeaksfile(peaklistout, outputblacklist[1]))

if args.scorebased:
    print('\n----> Cleaning and filtering peaksfile\n Score-based = ' + str(args.peakscore) + '\n')
    pybedtoolsout = readpeaksfile(peaklistout, outputblacklist[1])
else:
    print('\n----> Cleaning and filtering peaksfile\n E-value-based = ' + str(args.evalue) + '\n')
    pybedtoolsout = readpeaksfile(peaklistout, outputblacklist[1])
print(pybedtoolsout)


print('\n----> Creating bins for each peak files as per the given binsize, \nBinsize = ' + args.gbinsize + '\n')
quantifybins(pybedtoolsout)
