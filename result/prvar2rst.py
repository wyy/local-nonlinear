fr = open('PRVAR.lis')
fw = open('ansys-rst.txt','w')

i = 0
for line in fr:
    i += 1
    if i <= 6:
        continue
    elif i <= 26:
        fw.write(line)
    else:
        i = 1

fw.close()
fr.close()
