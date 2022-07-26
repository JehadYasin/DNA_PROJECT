with open("fASTA.txt", "r") as file:
    fasta = [item.strip() for item in file.readlines()]
def gc_content(seq):
    return float(round(((seq.count("G") + seq.count("C")) / len(seq)) * 100, 6))

fastadict = {}
fastalabel = ""

for line in fasta:
    if ">" in line:
        fastalabel = line
        fastadict[fastalabel] = ""
    else:
        fastadict[fastalabel] += line
vli = []
names = []
for key, item in fastadict.items():
    names.append(f'{key}\n')
    vli.append(gc_content(item))
merged_list = [(names[i], vli[i]) for i in range(0, len(vli))]
final = [''.join(map(str, x)) for x in merged_list]
finaldic = {}
for i in final:
    finaldic[i[1:15]] = i[15:]
g = max(finaldic, key=finaldic.get)
print(f'{g}{finaldic[g]}')



