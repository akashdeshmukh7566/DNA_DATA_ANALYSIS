from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import IsoelectricPoint as IP
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import io
import base64

def mrna_len(mRNA):                          #1.mRNA mRNA
    return len(mRNA)

def amino_acid_mRNA(mRNA):                               #2function for list of Amino Acid mRNAuence
    stop_codon = ["UAA", "UAG", "UGA"]
    codon = ""        
    i = 0 
    for i in range(i,len(mRNA),3):
        code = mRNA[i:i+3]
        if code in stop_codon:
            break
        else:
            codon = codon + code
    mrna_seq = Seq(codon)                             #biopython mRNA method object convert myrna string to mRNAuence 
    return  mrna_seq.translate()                                        #print amino acid mRNAuence by translate method 

def amino_acid_len(mRNA):                                #3 amino_acid_len
    return len(amino_acid_mRNA(mRNA))

def codon_lists(mRNA):                               #5.function for list of codon
    stop_codon = ["UAA", "UAG", "UGA"]
    codon = []          
    i = 0
    for i in range(i,len(mRNA),3):
        code = mRNA[i:i+3]
        if code in stop_codon:
            break
        else:
            codon.append(code)         
    return codon

def codon_len(mRNA):                            #6 codon len   
    return len(codon_lists(mRNA))

def amino_acid_rows(mRNA):
    print("amino_acid_rows_called--------mRNA------->>",mRNA)
    #1
    phe1 = codon_lists(mRNA).count("UUU")/codon_len(mRNA)*100
    phe2 = codon_lists(mRNA).count("UUC")/codon_len(mRNA)*100
    leu1 = codon_lists(mRNA).count("UUA")/codon_len(mRNA)*100
    leu2 = codon_lists(mRNA).count("UUG")/codon_len(mRNA)*100

    #2
    leu3 = codon_lists(mRNA).count("CUU")/codon_len(mRNA)*100
    leu4 = codon_lists(mRNA).count("CUC")/codon_len(mRNA)*100
    leu5 = codon_lists(mRNA).count("CUA")/codon_len(mRNA)*100
    leu6 = codon_lists(mRNA).count("CUG")/codon_len(mRNA)*100

    #3
    ile1 = codon_lists(mRNA).count("AUU")/codon_len(mRNA)*100
    ile2 = codon_lists(mRNA).count("AUC")/codon_len(mRNA)*100
    ile3 = codon_lists(mRNA).count("AUA")/codon_len(mRNA)*100

    met = codon_lists(mRNA).count("AUG")/codon_len(mRNA)*100

    #4
    val1 = codon_lists(mRNA).count("GUU")/codon_len(mRNA)*100
    val2 = codon_lists(mRNA).count("GUC")/codon_len(mRNA)*100
    val3 = codon_lists(mRNA).count("GUA")/codon_len(mRNA)*100
    val4 = codon_lists(mRNA).count("GUG")/codon_len(mRNA)*100

    #5
    ser1 = codon_lists(mRNA).count("UCU")/codon_len(mRNA)*100
    ser2 = codon_lists(mRNA).count("UCC")/codon_len(mRNA)*100
    ser3 = codon_lists(mRNA).count("UCA")/codon_len(mRNA)*100
    ser4 = codon_lists(mRNA).count("UCG")/codon_len(mRNA)*100

    #6
    pro1 = codon_lists(mRNA).count("CCU")/codon_len(mRNA)*100
    pro2 = codon_lists(mRNA).count("CCC")/codon_len(mRNA)*100
    pro3 = codon_lists(mRNA).count("CCA")/codon_len(mRNA)*100
    pro4 = codon_lists(mRNA).count("CCG")/codon_len(mRNA)*100

    #7
    thr1 = codon_lists(mRNA).count("ACU")/codon_len(mRNA)*100
    thr2 = codon_lists(mRNA).count("ACC")/codon_len(mRNA)*100
    thr3 = codon_lists(mRNA).count("ACA")/codon_len(mRNA)*100
    thr4 = codon_lists(mRNA).count("ACG")/codon_len(mRNA)*100

    #8
    ala1 = codon_lists(mRNA).count("GCU")/codon_len(mRNA)*100
    ala2 = codon_lists(mRNA).count("GCC")/codon_len(mRNA)*100
    ala3 = codon_lists(mRNA).count("GCA")/codon_len(mRNA)*100
    ala4 = codon_lists(mRNA).count("GCG")/codon_len(mRNA)*100

    #9
    tyr1 = codon_lists(mRNA).count("UAU")/codon_len(mRNA)*100
    tyr2 = codon_lists(mRNA).count("UAC")/codon_len(mRNA)*100

    #10
    his1 = codon_lists(mRNA).count("CAU")/codon_len(mRNA)*100
    his2 = codon_lists(mRNA).count("CAC")/codon_len(mRNA)*100

    gln1 = codon_lists(mRNA).count("CAA")/codon_len(mRNA)*100
    gln2 = codon_lists(mRNA).count("CAG")/codon_len(mRNA)*100

    #11
    asn1 = codon_lists(mRNA).count("AAU")/codon_len(mRNA)*100
    asn2 = codon_lists(mRNA).count("AAC")/codon_len(mRNA)*100

    lys1 = codon_lists(mRNA).count("AAA")/codon_len(mRNA)*100
    lys2 = codon_lists(mRNA).count("AAG")/codon_len(mRNA)*100

    #12
    asp1 = codon_lists(mRNA).count("GAU")/codon_len(mRNA)*100
    asp2 = codon_lists(mRNA).count("GAC")/codon_len(mRNA)*100

    glu1 = codon_lists(mRNA).count("GAA")/codon_len(mRNA)*100
    glu2 = codon_lists(mRNA).count("GAG")/codon_len(mRNA)*100

    #13
    cys1 = codon_lists(mRNA).count("UGU")/codon_len(mRNA)*100
    cys2 = codon_lists(mRNA).count("UGC")/codon_len(mRNA)*100

    trp = codon_lists(mRNA).count("UGG")/codon_len(mRNA)*100
    #14
    arg1 = codon_lists(mRNA).count("CGU")/codon_len(mRNA)*100
    arg2 = codon_lists(mRNA).count("CGC")/codon_len(mRNA)*100
    arg3 = codon_lists(mRNA).count("CGA")/codon_len(mRNA)*100
    arg4 = codon_lists(mRNA).count("CGG")/codon_len(mRNA)*100

    #15
    ser5 = codon_lists(mRNA).count("AGU")/codon_len(mRNA)*100
    ser6 = codon_lists(mRNA).count("AGC")/codon_len(mRNA)*100

    arg5 = codon_lists(mRNA).count("AGA")/codon_len(mRNA)*100
    arg6 = codon_lists(mRNA).count("AGG")/codon_len(mRNA)*100

    #16
    gly1 = codon_lists(mRNA).count("GGU")/codon_len(mRNA)*100
    gly2 = codon_lists(mRNA).count("GGC")/codon_len(mRNA)*100
    gly3 = codon_lists(mRNA).count("GGA")/codon_len(mRNA)*100
    gly4 = codon_lists(mRNA).count("GGG")/codon_len(mRNA)*100

    amino_acids =["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

    percentage = [
        round(ala1+ala2+ala3+ala4, 2),
        round(arg1+arg2+arg3+arg4+arg5+arg6, 2),
        round(asn1+asn2, 2),
        round(asp1+asp2, 2),
        round(cys1+cys2, 2),
        round(gln1+gln2, 2),
        round(glu1+glu2, 2),
        round(gly1+gly2+gly3+gly4, 2),
        round(his1+his2, 2),
        round(ile1+ile2+ile3, 2),
        round(leu1+leu2+leu3+leu4+leu5+leu6, 2),
        round(lys1+lys2, 2),
        round(met, 2),
        round(phe1+phe2, 2),
        round(pro1+pro2+pro3+pro4, 2),
        round(ser1+ser2+ser3+ser4+ser5+ser6, 2),
        round(thr1+thr2+thr3+thr4, 2),
        round(trp, 2),
        round(tyr1+tyr2, 2),
        round(val1+val2+val3+val4, 2)
    ]
    return amino_acids, percentage



def protein_bar_chart(amino_acids,percentage):
    bar_colors = "#FF6700"
    fig,ax = plt.subplots(figsize=(9.7,3.6))
    bars = ax.bar(amino_acids,percentage,color=bar_colors)

    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x()+bar.get_width()/2,
            height,
            f'{round(height,2)}%',
            ha = 'center',
            va='bottom',
            fontsize=7,
            color='black'        
    )
    
    plt.xlabel('amino acid')
    plt.ylabel('percentage(%)')     

    plt.xticks(rotation=45,ha='right')
    plt.tight_layout()
    
    img_buf = io.BytesIO()
    plt.savefig(img_buf, format='png')
    img_buf.seek(0)
    img_base64 = base64.b64encode(img_buf.getvalue()).decode('utf-8')

    return img_base64   
    
# print(f"bar chart protein- {protein_bar_chart(amino_acids,percentage)}")
# print(protein_bar_chart(amino_acids,percentage))

#molecular weight
amino_acid_weight = {
    
    "A" : 71.08,
    "G" : 57.05,
    "S" : 87.08,
    "P" : 97.12,
    "V" : 99.13,
    "T" : 101.11,
    "C" : 103.15,
    "I" : 113.16,
    "L" : 113.16,
    "N" : 114.11,
    "D" : 115.09,
    "Q" : 128.13,
    "K" : 128.18,
    "E" : 129.12,
    "M" : 131.20,
    "H" : 137.14,
    "F" : 147.18,
    "R" : 156.19,
    "Y" : 163.18,
    "W" : 186.22,
    "O" : 113.11,
    "U" : 121.09
    }



#MOlecular Weight -
# AA_mRNA = amino_acid_mRNA(mRNA)
def mol_weight(AA_mRNA):      #function for molecular weight of amino acid 
    
    total_weight = 0
    for AA in AA_mRNA:
        if AA in amino_acid_weight:
            total_weight = total_weight + amino_acid_weight[AA]           
    return round(total_weight/1000,2)

# print(f"M_W - {mol_weight(AA_mRNA)}")


#isoElectric_point
def isoelectric_point(AA_mRNA):
    ip = IP.IsoelectricPoint(AA_mRNA)
    pI_value = ip.pi()
    return round(pI_value,2)

# print(f"PI - {isoelectric_point(AA_mRNA)}")



amino_acid_composition = {
    'A': {'C': 3, 'H':5,  'N':1, 'O':1},
    'R': {'C': 6, 'H': 12,'N':4, 'O':1},
    'N': {'C': 4, 'H':6,  'N':2, 'O':2},
    'D': {'C': 4, 'H': 5, 'N': 1, 'O': 3},  
    'C': {'C': 3, 'H': 5, 'N': 1, 'O': 1, 'S': 1}, 
    'E': {'C': 5, 'H': 7, 'N': 1, 'O': 3},  
    'Q': {'C': 5, 'H': 8, 'N': 2, 'O': 2},  
    'G': {'C': 2, 'H': 3, 'N': 1, 'O': 1}, 
    'H': {'C': 6, 'H': 7, 'N': 3, 'O': 1},  
    'I': {'C': 6, 'H': 11, 'N': 1, 'O': 1},  
    'L': {'C': 6, 'H': 11, 'N': 1, 'O': 1}, 
    'K': {'C': 6, 'H': 12, 'N': 2, 'O': 1},  
    'M': {'C': 5, 'H': 9, 'N': 1, 'O': 1, 'S': 1},  
    'F': {'C': 9, 'H': 9, 'N': 1, 'O': 1}, 
    'P': {'C': 5, 'H': 7, 'N': 1, 'O': 1},  
    'U': {'C': 5, 'H': 5, 'N': 1, 'O': 2},
    'S': {'C': 3, 'H': 5, 'N': 1, 'O': 2},  
    'T': {'C': 4, 'H': 7, 'N': 1, 'O': 2},
    'W': {'C': 11, 'H': 10, 'N': 2, 'O':1},  
    'Y': {'C': 9, 'H': 9, 'N': 1, 'O': 2}, 
    'V': {'C': 5, 'H': 9, 'N': 1, 'O': 1},  
    'O': {'C': 5, 'H': 7, 'N': 1, 'O': 2}
}   


def calculate_atomic_composition(AA_mRNA):
    atom_composition = {'C':0, 'H':0, 'N':0, 'O':0, 'S':0}
    
    for aa in AA_mRNA:
        if aa in amino_acid_composition:
            for atom,count in amino_acid_composition[aa].items():
                atom_composition[atom] = atom_composition[atom] + count
    
    return atom_composition   
# print(f"atomic_composition{calculate_atomic_composition(AA_mRNA)}") 

# atomic_composition = calculate_atomic_composition(AA_mRNA)
     
     
# print("Atomic composition - ")

# C = atomic_composition['C']
# H = atomic_composition['H']
# N = atomic_composition['N']
# O = atomic_composition['O']
# S = atomic_composition['S']
 
 
# def molecular_formula(c,H,N,O,S):           
#     return c,H,N,O,S
# print(f"{molecular_formula(c,H,N,O,S)}")

# ATOM_NUMBER = C+H+N+O+S

def atomic_percent(Molecular_Formula,ATOM_NUMBER):

    Carbon = round((Molecular_Formula[0]/ATOM_NUMBER)*100,2)
    Hydrogen = round((Molecular_Formula[1]/ATOM_NUMBER)*100,2)
    Nitrogen = round((Molecular_Formula[2]/ATOM_NUMBER)*100,2)
    Oxygen = round((Molecular_Formula[3]/ATOM_NUMBER)*100,2)
    Sulfur = round((Molecular_Formula[4]/ATOM_NUMBER)*100,2)
    print("Carbon,Hydrogen,Nitrogen,Oxygen,Sulfur-------->>",Carbon,Hydrogen,Nitrogen,Oxygen,Sulfur)
    return Carbon,Hydrogen,Nitrogen,Oxygen,Sulfur

# print(f"atomic % {atomic_percent(ATOM_NUMBER)}")
# atom_percent = atomic_percent(ATOM_NUMBER)    

def atom_number(C,H,N,O,S):                 #ATOMIC NUMBER OF PROTEIN
    return C+H+N+O+S

# print(f"atomic number : {atom_number(c,H,N,O,S)}")

def calculate_charge_residues(AA_mRNA):
    negatively_charged = ['D','E']          #aspartic acid(D) ,glutamic acid(E)
    positively_charged = ['R','K','H']      #arginine(R) ,lysine(K),hitidin(H)
    
    negative_count = 0
    positive_count = 0
    for aa in AA_mRNA:
        if aa in negatively_charged:
            negative_count = negative_count + 1
            
        elif aa in positively_charged:
             positive_count = negative_count + 1   

    return negative_count,positive_count
# print(f" charge + , - : {calculate_charge_residues(AA_mRNA)}")

hydropathy_index = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8, 
    'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5, 
    'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
}


def calculate_gravy(AA_mRNA):
    total_hydropathy = 0
    for aa in AA_mRNA:
        total_hydropathy = total_hydropathy +  hydropathy_index[aa]  
    gravy = total_hydropathy / len(AA_mRNA)                          
    return round(gravy,2)
# print(f"gravy : {calculate_gravy(AA_mRNA)}")

def calculate_net_charge(AA_mRNA):                  #calculate net chareg of protein
    analysis = ProteinAnalysis(AA_mRNA)
    net_charge = analysis.charge_at_pH(7)
    return round(net_charge, 2)
# print(f"calcu.net_charge : {calculate_net_charge(AA_mRNA)}")

