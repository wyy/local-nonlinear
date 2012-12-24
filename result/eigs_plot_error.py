# plot
fin = open('log-eigs.txt')
fout = open('eigs-plot.txt','w')
str_header = []
eigs_origin = []
for line in fin:
    header, eigs = line.split(':')
    str_header.append(header)
    eigs_origin.append(eigs.split())
num = len(eigs_origin[0])
eigs_out = [[row[i] for row in eigs_origin] for i in range(num)]
fout.write('|'.join(str_header) + '\n')
for row in eigs_out:
    fout.write(' '.join(row) + '\n')
fin.close()
fout.close()
# error
error_lr = []
error_irs = []
for i in range(num):
    exact = float(eigs_origin[0][i])
    lr = float(eigs_origin[1][i])
    irs = float(eigs_origin[-1][i])
    error_lr.append(abs(lr-exact)/exact)
    error_irs.append(abs(irs-exact)/exact)
fout = open('eigs-error.txt','w')
fout.write(str_header[1] + ':')
for i in range(num):
    fout.write(' {:.2%}'.format(error_lr[i]))
fout.write('\n')
fout.write(str_header[-1] + ':')
for i in range(num):
    fout.write(' {:.2%}'.format(error_irs[i]))
fout.write('\n')
fout.close()
