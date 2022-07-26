import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import typing

style.use("fivethirtyeight")
fig = plt.figure()


import collections
from utils import colored
complement = {"A": "T", "G": "C", "T": "A", "C": "G"}
rna_complement = {"A": "U", "T": "A", "G": "C", "C": "G"}
DNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}


def reverse(seq):
    '''Providing the complementary DNA strand'''
    #return "".join([colored(complement[nuc]) for nuc in seq])
    mapping = str.maketrans("ACGT", "TGCA")
    return colored(seq.translate(mapping))

def length(seq):
    length = len(seq)
    NucLength = dict(collections.Counter(seq))
    return f'''
Length: {length} nucleotides
Nucleotide occurance:
Adenine: {NucLength["A"]}
Cytosine: {NucLength["C"]}
Guanine: {NucLength["G"]}
Thymine: {NucLength["T"]}'''

def rna(seq):
    print(" ")
    return "RNA sequence: " + "5' " + "-".join([colored(rna_complement[nuc]) for nuc in seq]) + " 3'"

def helix(seq):
    print('''
    
Double Helix
    ''')
    return f'''
5' {colored(seq)} 3'
   {len(seq) * "I"}
3' {reverse(seq)} 5'
'''
def GC_Content(seq):
    return round(((seq.count("G") + seq.count("C")) / len(seq)) * 100)
def GC_kmer(seq, k):
    res = ["GC-Content"]
    gc_li = []
    n = 0
    if len(seq) - int(len(seq)) > 0:
        for l in range(int(round(len(seq) / k)) - 1):
            gc_li.clear()
            for i in range(k):
                gc_li.append(seq[n])
                n += 1
            percent = ((gc_li.count("G") + gc_li.count("C")) / len(gc_li)) * 100
            res.append(round(percent, 1))
        return res
    else:
        for l in range(int(round(len(seq) / k))):
            gc_li.clear()
            for i in range(k):
                try:
                    gc_li.append(seq[n])
                    n += 1
                except IndexError:
                    pass

            percent = ((gc_li.count("G") + gc_li.count("C")) / len(gc_li)) * 100
            res.append(round(percent, 1))
        return res

def AT_Content(seq):
    return round(((seq.count("A") + seq.count("T")) / len(seq)) * 100)
def AT_kmer(seq, k):
    res1 = ["AT-Content"]
    gc_li = []
    n = 0
    if len(seq) - int(len(seq)) > 0:
        for l in range(int(round(len(seq) / k)) - 1):
            gc_li.clear()
            for i in range(k):
                gc_li.append(seq[n])
                n += 1
            percent = ((gc_li.count("A") + gc_li.count("T")) / len(gc_li)) * 100
            res1.append(round(percent, 1))
        return res1
    else:
        for l in range(int(round(len(seq) / k))):
            gc_li.clear()
            for i in range(k):
                try:
                    gc_li.append(seq[n])
                    n += 1
                except IndexError:
                    pass
            percent = ((gc_li.count("A") + gc_li.count("T")) / len(gc_li)) * 100
            res1.append(round(percent, 1))
        return res1

def matlab_plot(seq, k):
    y1 = GC_kmer(seq, k)
    y1.remove("GC-Content")
    y2 = AT_kmer(seq, k)
    y2.remove("AT-Content")
    x = [i for i in range(len(y1))]
    plt.plot(x, y1, label="GC content", marker = "D")
    plt.plot(x, y2, label="AT content", marker = "D", linestyle= "dashed")
    plt.title("Nucleotide Content through DNA sequence")
    plt.xlabel(F"{k}-Mer Order (Sequence Position)")
    plt.ylabel(f"Percentage of Nucleotides in {k}-mer (%)")
    plt.legend()
    plt.show()

def translation(seq, init: int):
    return [DNA_Codons[seq[pos: pos + 3]] for pos in range(init, len(seq) - 2, 3)]








