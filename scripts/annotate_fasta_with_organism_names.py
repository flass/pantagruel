fin = open('/home/flassall/AgroRhizo2018/02.gene_alignments/full_cdsfam_alignments_species_code/AGRHIZOC004685.codes.aln', 'r')
fout = open('/home/flassall/Dropbox/Frontier/ML-tree-Florent_437genomes_AgrobacteriumRhizobium_group/AGRHIZOC004685_telA_CDS.aln', 'w')

dcodeorga = {}
with open('/home/flassall/AgroRhizo2018/03.database/organism_codes.tab') as fcodeorga:
    for line in fcodeorga:
        lsp = line.rstrip('\n').split('\t')
        code = lsp[0]
        orga = lsp[1]
        stra = lsp[2]
        if not (stra in orga):
            orga += ' str. '+stra
        elif stra:
            if ' str. ' not in orga:
                orga = orga.replace(stra, 'str. '+stra)
        dcodeorga[code] = orga.replace('(', '').replace(')', '')

for line in fin:
    if line.startswith('>'):
        cdscode = line.strip('>\n')
        code = cdscode.split('_')[0]
        fout.write(">%s %s\n"%(cdscode, dcodeorga[code]))
    else:
        fout.write(line)

fin.close()
fout.close()
