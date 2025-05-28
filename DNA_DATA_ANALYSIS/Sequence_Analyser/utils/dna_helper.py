from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight

#task -1 len of dna 




def dna_len(seq):                                                 #dna len
    return len(seq)

def get_atgc(seq):                                                #Calculate ATGC content percentage
    length = len(seq)        
    A = round(seq.count("A")/length*100,2)
    T = round(seq.count("T")/length*100,2)
    G = round(seq.count("G")/length*100,2)
    C = round(seq.count("C")/length*100,2)
    return A,T,G,C

def dna_stability(seq):                                            # stability of dna       
    length = len(seq)    
    GC = round((seq.count("G") + seq.count("C"))/length*100)

    if GC < 40 :                          # chack condition 
        return "low" 
    elif GC <=60 and GC > 40 :
        return "Moderate"  
    else:
        return "high"

def melt_temp(seq):                                                  # melting temperature 
    return round(mt.Tm_NN(seq,strict=False),2)

def get_GC(seq):                                                         # function for GC %
    length = len(seq)   
    return round((seq.count("G") + seq.count("C"))/length*100,2)

#

# def molecular_weight(seq):                                           #function for molecular weight
#     mol_weight = 0
#     for i in range(len(seq)):
#         if seq[i] == "A":
#             mol_weight = mol_weight + 331.2 
#         elif seq[i] == "T":
#             mol_weight = mol_weight + 322.2
#         elif seq[i] == "G":
#             mol_weight = mol_weight + 347.2
#         elif seq[i] == "C":
#             mol_weight = mol_weight + 307.2
#     return round(mol_weight/1000,2)




def Molecular_weight(seq): 
    mol_weight = molecular_weight(seq)
    mol_weight_kda = mol_weight / 1000
    return round(mol_weight_kda, 2)
    

def complement_seq(seq):                                       # Reverse Complement of gene (dna)
    return Seq(seq).complement()

def reverse_complement_seq(seq):               #reverse complement sequence
    return Seq(seq).reverse_complement()
