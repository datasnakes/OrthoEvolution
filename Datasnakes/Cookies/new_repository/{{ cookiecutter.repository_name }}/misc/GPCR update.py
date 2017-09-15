import os
import csv
os.chdir('C:\\Users\\rgilmore\\Desktop\\db\\Accession Numbers\\Combining Accession Files\\Combo and Additions on 10-25')
count = 0
old_genes = []
old_acc = []
old_dict = {}
new_dict = {}
# with open('Master_Accession_File_combo.csv', 'r+') as csvfile1:
#     re = csv.reader(csvfile1, delimiter =',')
#     wr = csv.writer(csvfile1)
#     with open('Gene_Update.csv', 'r') as csvfile2:
#         update = csv.reader(csvfile2, delimiter = ',')
#         for line1 in re:
#             for line2 in update:
#                 line2[count] = str(line2[count])
#                 if line1[0] in line2:
#                     line1[0] = line2[0]
#                     wr.writerow(line1)
#                 count += 1
with open('Master_RNA_Accession_File_version_split.csv', 'r') as csvfile1:
    r1 = csv.reader(csvfile1)
    lines1 = [l for l in r1]
    with open('Gene_Update.csv', 'r') as csvfile2:
        r2 = csv.reader(open('Gene_Update.csv'))
        lines2 = [li for li in r2]
        count = 0
        with open('Gene_Update_Isoforms2.csv', 'r') as csvfile6:
            reader = csv.reader(csvfile6)
            l_c = 0
            isoform_list = []
            iso_acc_list = []
            head_iso = []
            for _Line6 in reader:
                l_c += 1
                isoform = _Line6[0]
                accession = _Line6[1]

                head, sep, tail = _Line6[0].partition('.')

                isoform_list.append(isoform)
                iso_acc_list.append(accession)
                if head in new_dict.keys():
                    if type(head) is not dict is True:
                        input('Keys are crazy right now')

                if head not in new_dict.keys():
                    new_dict[head] = {}
                    head_iso.append(head)

                new_dict[head][isoform] = [accession, '', '', '', '']
            print('head_iso ',head_iso)
            print('isoform_list ', isoform_list)
            print('iso_acc_list', iso_acc_list)
            print(new_dict.keys())
            print(new_dict.values())
            input('and then?')

        if count == 0:
            GPCR_list = []
            Acc_list = []
            for line2 in lines2:
                if line2[0] == 'HGNC symbol':
                    continue
                GPCR_list.append(line2[0])
                Acc_list.append(line2[1])

            print(GPCR_list)
        for line1 in lines1:
            count += 1
            if count == 1:
                with open('New_Organism_List.csv', 'r') as csvfile4:
                    r3 = csv.reader(csvfile4)
                    for line3 in r3:
                        for cell in line3:
                            if cell not in line1:
                                line1.append(cell)

            for line2 in lines2:
                x = str(line1[2])
                y = line1[1]
                if x.upper() in line2:
                    print('Gene name before: ', line1[1])

                    line1[1] = line2[0]
                    print('Gene name after: ', line1[1])
                    print('')
                    old_acc.append(line1[2])
                    old_genes.append(line1[1])
                    old_dict[line1[2]] = line1[1]
                    if x.upper() in iso_acc_list:
                        for _isoform in isoform_list:
                            if y in _isoform:
                                if x.upper() in new_dict[y][_isoform]:
                                    print('Gene name before: ', line1[1])
                                    line1[1] = _isoform
                                    print('Gene name after: ', line1[1])
                                    print('')
                                    input('isoform ok?')

                elif line1[1] in line2:
                    print('Gene name before: ', line1[1])
                    line1[1] = line2[0]
                    print('Gene name after: ', line1[1])
                    print('')
                    old_acc.append(line1[2])
                    old_genes.append(line1[1])
                    old_dict[line1[2]] = line1[1]
with open('NewMaster_Accession_File_Template.csv', 'r') as csvfile5:
    reader = csv.reader(csvfile5)
    l_c = 0
    for _Line5 in reader:
        l_c += 1
        if l_c == 1:
            continue
        temp = _Line5[0:5]
        t_emp = temp.pop(0)
        temp.insert(1, '')
        new_dict[t_emp] = temp






with open('Master_Accession_File_update.csv', 'w', newline= '\n', encoding="utf-8") as csvfile3:
    writer = csv.writer(csvfile3)
    writer.writerows(lines1) #lines1 contains updated gene names
    for item in GPCR_list:
        item, sep, tail = item.partition('.')
        if item in old_genes:
            continue
        elif item in new_dict.keys():
            print(item, 'in new_dict.keys()')
            if isinstance(new_dict[item], dict):
                print('is a dict')
                for _x in new_dict[item].keys():
                    temp = []
                    temp2 = []
                    temp2 = ['None', _x]
                    temp2 = temp2 + new_dict[item][_x]
                    temp.append(temp2)
                    writer.writerows(temp)
            else:
                print('is not a dict')
                temp = []
                temp2 = []
                temp2 = ['None', item]
                temp2 = temp2 + new_dict[item]
                temp.append(temp2)
                writer.writerows(temp)
        else:
            print(item)
            input('what happened?')


    print(temp)
    print(new_dict.keys())
    print(new_dict.values())
