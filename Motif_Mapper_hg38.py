import glob,datetime

def reverse(motif):
    motif=motif.replace('A',' t ')
    motif=motif.replace('T',' a ')
    motif=motif.replace('C',' g ')
    motif=motif.replace('G',' c ')
    motif=motif.upper()   
    motif=motif.split()
    motif=motif[::-1]
    motif=''.join(motif)   
    return motif

def find_motif(motif,file_cromosoma):
    file_cromosoma = open( file_cromosoma ,'r')
    file_cromosoma.close
    sequenza_cromosoma = file_cromosoma.read().split()
    nome_sequenza_cromosoma = sequenza_cromosoma[0][1:]
    del sequenza_cromosoma[0]
    print('Chromosome '+nome_sequenza_cromosoma[3:])
    sequenza_cromosoma = ''.join(sequenza_cromosoma)
    lunghezza_sequenza_cromosoma = len(sequenza_cromosoma)
    motif = motif
    reverse_motif = reverse(motif)
    lunghezza_motif = len(motif)
    numero_finestre=lunghezza_sequenza_cromosoma-(lunghezza_motif-1)




    file_bed=open('motif['+motif+'].hg38.bed','a')
    contatore = 0
    while contatore != numero_finestre:
        if  sequenza_cromosoma[contatore:(contatore+lunghezza_motif)].upper()== motif:
            file_bed.write(nome_sequenza_cromosoma+'\t'+str(contatore)+'\t'+str(contatore+(lunghezza_motif))+'\t'+motif+'\t'+str(contatore)+'\t'+'+'+'\n')
        elif sequenza_cromosoma[contatore:(contatore+lunghezza_motif)].upper()== reverse_motif:
            file_bed.write(nome_sequenza_cromosoma+'\t'+str(contatore)+'\t'+str(contatore+(lunghezza_motif))+'\t'+reverse_motif+'\t'+str(contatore)+'\t'+'-'+'\n')            
         
        contatore = contatore + 1    
    file_bed.close()
    
print("""_________________________________________________""")
print("""                                                 """)
print("""Mark Motif Mapper (Created by Cammisa Marco 2012)""")
print("""_________________________________________________""")
print("""                                                 """)
motif = raw_input("Inserisci la sequenza del motif dal 5' ---> 3' = ")

motif = motif.strip()
motif = motif.upper()
print ("Il motif cercato e = 5'- "+motif+" -3'")
print ('Inizio analisi...')

data_ed_ora = str(datetime.datetime.now())
data_ed_ora=data_ed_ora.replace(' ','_')
data_ed_ora=data_ed_ora.replace(':','-')
data_ed_ora=data_ed_ora.replace('.','-')

file_bed=open('motif['+motif+'].hg38.bed','w')
file_bed.write('')
file_bed.close()


files_genoma = glob.glob("hg38/*.fa")
for file_cromosoma in files_genoma:
    find_motif(motif,file_cromosoma)

print('Il risultato si trova nel file : '+'motif['+motif+'].hg38.bed')
fine=raw_input("""Premere invio per terminare il programma""")
    
