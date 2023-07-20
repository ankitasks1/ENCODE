import os, re


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
        newcontent = ''.join(content[5] + '\t' +content[14] + '\t' + content[44])
        content_list.append(newcontent)
    return content_list


def find_accession_in_sheet(acc, cont):
    accession_in_sheet_list = []
    count = 0
    for code in acc:
        code = code[0]
        for lines in cont:
            #print(code, lines)

            for matchpos in re.finditer(code, lines):
                count += 1
                accession_in_sheet_list.append(''.join(str(count) + '\t' + matchpos.group()  + '\t' + code  + '\t' + lines))
    return accession_in_sheet_list


def rename_and_adjust(code_and_sheet):
    if os.path.exists('new_script_to_add_id.sh'):
        os.remove('new_script_to_add_id.sh')
        print('Reading the BED files')
        for i in code_and_sheet:
            i = i.split('\t')
            os.system('cp ' + i[1] + '.bed.gz ' + i[3] + '_' + i[1] + '_peaks.bed.gz')
            os.system('gunzip ' + i[3] + '_' + i[1] + '_peaks.bed.gz')
            script_to_run = ''.join('awk ' + """'{print $1"\\t"$2"\\t"$3"\\t""peaks_"NR"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10"\\t""" + '''""''' + i[3] + '_' + i[1] + '"' + """}' """  + i[3] + '_' + i[1] + '_peaks.bed > ' + i[3] + '_' + i[1] + '_peaks_id.bed')
            with open('new_script_to_add_id.sh', 'a') as myscript:
                myscript.write(script_to_run + '\n')
                myscript.close()
        print('Finished unzipping the data')
        print('\n.\n.\n')
        print('Adding peaks IDs to BED files')
        print('\n.\n.\n')
        print('Now run \n 1) chmod +x ./new_script_to_add_id.sh \n 2)  ./new_script_to_add_id.sh')


ID = open("files.txt", 'r')
sheet = open('experiment_report_2022_11_22_12h_26m.tsv', 'r')


myacc = read_accession(ID)
mysheet = read_exp_sheet(sheet)
match_acc_sheet = find_accession_in_sheet(myacc, mysheet)
rename_and_adjust(match_acc_sheet)

