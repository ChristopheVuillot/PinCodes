LISTSTOCLEAN = ['list_n_k_syst43_dim3.txt',
                'list_n_k_syst33_dim2.txt',
                'list_n_k_syst34_dim2.txt',
                'list_n_k_syst33.txt']

for FILE in LISTSTOCLEAN:
    with open(FILE, 'r') as filetoclean:
        lines = filetoclean.readlines()
        uniquelines = set(lines[1:])
    with open(FILE, 'w') as filecleaned:
        filecleaned.write('n k d\n')
        filecleaned.writelines(list(uniquelines))
