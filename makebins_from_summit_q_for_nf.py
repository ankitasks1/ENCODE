import os
import sys
from sortedcontainers import SortedDict
from pybedtools import BedTool
import argparse

print('''
-------------------------------------------------------------------------------------
|                           Summit and Genomic-sites marker                         |
|                                Author: Ankit Verma                                |
|                        Contact: ankitverma9079@gmail.com                          |
| Note: This code will work only with any peak files with 3 columnsi.e chr-start-end|
-------------------------------------------------------------------------------------
''')

print('\nAnalysis started')


parser = argparse.ArgumentParser(prog='summitmarker', description='''
Description of usage: Make genomic bins from ChIP-Seq peaks summit. The peaks can be obtained from any peak caller.
Peak file should have chr at col1 , start at col2, end at col3
This format you generally chr   start   end
''')
                                 
parser.add_argument('--inputfile', help='Name of the input file containing all the *peaks.bed')
parser.add_argument('--blacklist', help='Name of the blacklist regions file')
parser.add_argument('--output', help='BaseName of the output  file', default='_bin_out.txt')
parser.add_argument('--gbinsize', help='Genomic bin size',  default=1000)
parser.add_argument('--path', help='path to the summit bed files location')

# Print help message if no arguments
if len(sys.argv)==1:
   parser.print_help(sys.stderr)
   sys.exit(1)

args = parser.parse_args()


print('\nOutput will be written to: ' + args.path)

# read peak file list
def readpeaklist(path, inputpeaklist):
    samplelista = open(''.join(path + inputpeaklist).strip())
    samplelista = samplelista.read().strip().split('\n')
    mybedlist = []
    for bedfiles in samplelista:
        bedfiles = bedfiles.strip().split('\t')
        mybedlist.append(bedfiles[0])
    return mybedlist

# read bedfile inside peaklist using pybedtools
def readpeaksfile(path, eachpeakfiles):
    pybedtooldict = {}
    for samples in eachpeakfiles:
        # print(samples)
        peakfile = open(''.join(path + samples), 'r')
        peakfile_name = os.path.basename(peakfile.name).split("_")[0]
        # print(peakfile_name)
        peaks = BedTool(''.join(path + samples))
        # print(peaks)
        # No filtering as summit does not have column for filtering and it is already the output of filtered data
        # get the number of intervals in the peak object
        peaks_count = peaks.count()
        # print(samples, peaks_count)
        if peakfile_name not in pybedtooldict:
            pybedtooldict[peakfile_name] = peaks, peaks_count
            # print(peakfile_name, "unique name")
            # print(pybedtooldict[peakfile_name][0])

        elif peakfile_name in pybedtooldict:
            tempdict = {}
            tempdict[peakfile_name] = peaks
            # print(tempdict[peakfile_name])
            # print(peakfile_name, "repeated name found")
            # Concatenate peaks if the same CAps is present
            # print(pybedtooldict[peakfile_name][0].cat(tempdict[peakfile_name]))
            concatenated_peaks = pybedtooldict[peakfile_name][0].cat(tempdict[peakfile_name])
            concatenated_peaks_count = concatenated_peaks.count()
            pybedtooldict[peakfile_name] = concatenated_peaks, concatenated_peaks_count

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
        with open(args.output, 'a') as myfile:
            for i in SortedDict(bindict):
                for j in SortedDict(bindict[i]):
                # print(i, j)
                    # print(''.join(i + '\t' + j + '\t' + str(bindict[i][j])))
                    myfile.write(''.join(i + '\t' + j + '\t' + str(bindict[i][j])) + '\n')
            myfile.close()
        fulldict.update(bindict)


# remove if output file already exist
# if os.path.exists(args.output):
#     os.remove(args.output)
print('\n----> Reading peaklist,\n Peaklist = ' + args.inputfile + '\n')
peaklistout = readpeaklist(args.path, args.inputfile)
# print(peaklistout)

# read peakfile   
# print(readpeaksfile(args.path, peaklistout))                                                                                         
pybedtoolsout = readpeaksfile(args.path, peaklistout)

# create bins
print('\n----> Creating bins for each peak files as per the given binsize, \nBinsize = ' + args.gbinsize + '\n')
quantifybins(pybedtoolsout)


