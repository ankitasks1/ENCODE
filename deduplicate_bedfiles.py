import os,sys
# usage: python deduplicate_bedfiles.py  bedfiles.txt samplelist_1.txt
bedfiles = {}
newlist = {}
# Make a dictionary of CAPs and respective bed files
with open(sys.argv[1]) as mybedfiles:
    for lines in mybedfiles.readlines():
        lines = lines.rstrip().split('\n')
        CAP = lines[0].rstrip().split('_')
        # print(CAP)
        if CAP[0] in bedfiles:
            bedfiles[CAP[0]].add(lines[0])
        else:
            bedfiles[CAP[0]] = {lines[0]}

# Collect CAPs with >1 bed files and concatenate
for CAPs in  bedfiles:
    if len(bedfiles[CAPs]) > 1:
        # print(CAPs, bedfiles[CAPs])
        newlist[CAPs] = ''.join(CAPs + '_CONCAT_' + 'peaks_id.bed')
        # Iterate through the files and concatenate their contents
        combined_content = ""
        for filename in bedfiles[CAPs]:
            with open(filename, 'r') as file:
                file_content = file.read()
                combined_content += file_content
        
        updatedbed = combined_content.strip().split('\n')
        if os.path.exists(''.join(CAPs + '_CONCAT_' + 'peaks_id.bed')):
            os.remove(''.join(CAPs + '_CONCAT_' + 'peaks_id.bed'))
        with open(''.join(CAPs + '_CONCAT_' + 'peaks_id.bed'), 'w') as myoutputfile:
            for rows in updatedbed:
                    rows = rows.strip().split('\t')
                    myoutputfile.write(''.join('\t'.join(rows) + '\n'))
            myoutputfile.close()
        
    else:
        newlist[CAPs] = bedfiles[CAPs].pop()

if os.path.exists(sys.argv[2]):
    os.remove(sys.argv[2])

with open(sys.argv[2], 'a') as mynewlist:
     mynewlist.write(''.join('beds' + '\n'))
     mynewlist.close()

with open(sys.argv[2], 'a') as mynewlist:
    for outbeds in newlist:
        mynewlist.write(''.join(newlist[outbeds] + '\n'))
    mynewlist.close()

