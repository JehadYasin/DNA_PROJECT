import DNAtoolkit


DNASEQ = input("Enter the nucleotide sequence of DNA template strand (5' - 3'): ").upper()

def full(seq):
    print(f'\n[1] -- {DNAtoolkit.rna(seq)}\n')
    print(f"[2] -- {DNAtoolkit.length(seq)} ")
    print(f'[3] -- {DNAtoolkit.helix(seq)}')
    print(f"[4] -- 5' {DNAtoolkit.reverse(seq)} 3' [reverse complement]")
    print(f"\n[5] -- GC-Content: {DNAtoolkit.GC_Content(seq)} %")
    print(f"\n[5] -- AT-Content: {DNAtoolkit.AT_Content(seq)} %")
    kinput = int(input("Enter a value of K for K-mer  content calculation: "))
    print(f"[6] -- {kinput}mer GC content 5'-3' (%):  {DNAtoolkit.GC_kmer(seq, kinput)}")
    print(f"[7] -- {kinput}mer AT content 5'-3' (%):  {DNAtoolkit.AT_kmer(seq, kinput)}")
    initinput = input("Enter the translation starting point: ")
    print(F' Amino Acid sequence: {DNAtoolkit.translation(seq, int(initinput))}')
    ask = input("Do you want to show a quick graphical display (Yes/No): ")
    if ask.upper() == "YES":
        DNAtoolkit.matlab_plot(seq, kinput)

    with open("graph.csv", "w") as file:
        li = DNAtoolkit.GC_kmer(seq, kinput)
        li1 = DNAtoolkit.AT_kmer(seq, kinput)
        merged_list = [(li[i], li1[i]) for i in range(0, len(li))]
        final = [', '.join(map(str, x)) for x in merged_list]

        file.writelines([F'{item}\n' for item in final])



full(DNASEQ)
#print(DNAtoolkit.transcribe("ACTGCTACGAGCTCAAAGATTTCAGATCCAGGTTTGTTGCAGTAGGACCACCTAACCTGAGCCCCTGCGCACAGAGTTTTCGAGGGTACACCGTCAACCTAAGCGTCGGTCGTTGCCCATTTAAACTCGTCGTCGTCACAGTAGAGGGTTGGTACAGTACGGCTATTTTGGCCTGGAGTGATAGCGTAGCGCCGGCTCAATCATCGCCGTCTGTGTCTTCTTAGTCCTTAGTTGATCACAAATCTCTGAACCGCTTCTCCCCCGAGTACGTGTCTGGCGGGACGCTGTGCTGTTCCACGTAGTACGTCAGCACCCTTCATCATGACGCTTTAAGGCAAGCACTTAGTCTACGTCGGGAACTGACGGCCCCTGATGAACCATCGTCTCTGGTTTGGGTTTATACACGGCCCAACTTTCACCCGCTTACCAATCATGAACTCTCCAGCCACCATTGTACGAACGCGGCAGGGAAAGGAGCAAATCGTAGAAAGCCCTTGTGTTGTAAAACGTAGCTCCATAGGGCGCATCTCATGGTAAATAAGGATCGTCGCGGCATGCTGGATTATAGATGTGTAGTGACCGTTTGCGAATTCATGCCGGGGCGGGCGCGTTTCGGATTCATCAAAACATTAGTAAAACGCATTGTCACGCCTGTACTCAGACGTCATTACCCATCTGTGTGAAACCTGCGCCCAGATCGGGCCATGAGGTTCATTGTTCAGCTAGCATCACGGGAGTGAATCGACCTGCGTGAAAGGGGCTTCAAAAACCATTCTAGTTAAGTAGGTCTAAGAGAAACTGGGGTACGCTAAACCGAAGTATTTCCGCAATCGCTGCGTTAGTTAGGATGAGGTTCGTGAATTGGGGGATGATGCCCGGCTGAATGAATTTGCTCAGTGTCATTAAATTACGCACAAGAAATCATCGTAAGACTGGCGCATGTCA"))















