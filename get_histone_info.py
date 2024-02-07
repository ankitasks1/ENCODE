
#  how to run: python get_histone_info.py files_H1.txt experiment_report_2024_2_6_15h_22m_hESC_H1.tsv
# this will copy the encode file to new with histone and encode id
import os, re
import sys

def read_accession(ID):
    ids = ID.read().strip().split('\n')
    del ids[0]
    accession_list = []
    for links in ids:
        links = links.split('/')
        accession = links[6].split('.')
        accession_list.append(accession)
    return accession_list

def read_exp_sheet(sheet):
    experiment_sheet = sheet.read().strip().split('\n')
    del experiment_sheet[0]
    content_list = []
    for content in experiment_sheet:
        content = content.split('\t')
        newcontent = ''.join(content[6] + '\t' +content[15])
        content_list.append(newcontent)
    return content_list


def find_accession_in_sheet(acc, cont):
    accession_in_sheet_list = []
    count = 0
    for info in acc:
        protein = info[0]
        extension = info[1]
        for lines in cont:
            # print(code, lines)
            # use regular expression
            for matchpos in re.finditer(protein, lines):
                count += 1
                # print(matchpos.group() +'\t' + lines)
                accession_in_sheet_list.append(''.join(str(count) + '\t' + matchpos.group()  + '\t' + protein + '\t' + extension  + '\t' + lines))
    return accession_in_sheet_list

def rename_and_adjust(code_and_sheet):
    if os.path.exists('new_script_to_add_id.sh'):
        os.remove('new_script_to_add_id.sh')
    else:
        print('Reading the BED files')
        for i in code_and_sheet:
            i = i.split('\t')
            print(i[1] + '.'  + i[3] + '\t' + i[4] + '_' + i[1] + '.'  + i[3])
            os.system('cp ' + i[1] + '.'  + i[3] + ' ' + i[4] + '_' + i[1] + '.'  + i[3])

ID = open(sys.argv[1], 'r')
sheet = open(sys.argv[2], 'r')


myacc = read_accession(ID)
mysheet = read_exp_sheet(sheet)
match_acc_sheet = find_accession_in_sheet(myacc, mysheet)

# print(myacc)
# print(mysheet)
# print(match_acc_sheet)
rename_and_adjust(match_acc_sheet)
