# from Bio.SeqUtils.ProtParam import ProteinAnalysis
# from Bio.Seq import Seq
# from Bio.SeqUtils import IsoelectricPoint as IP 
import matplotlib.pyplot as plt
import io
import base64
from Bio.SeqUtils import molecular_weight

def mrna_len(mRNA):                          #1.mRNA len 
    return len(mRNA)

def mrna_augc(mRNA):                       #2.mRNA AUGC
    mRNA_len = len(mRNA)
    
    A1 = round(mRNA.count("A")/mRNA_len*100,2)
    U1 = round(mRNA.count("U")/mRNA_len*100,2)
    G1 = round(mRNA.count("G")/mRNA_len*100,2)
    C1 = round(mRNA.count("C")/mRNA_len*100,2)
    return A1,U1,G1,C1

def mrna_molecular_weight(mRNA):                   #3.function for molecular weight
    mol_weight = molecular_weight(mRNA,seq_type='RNA')
    mol_weight_kda = mol_weight / 1000
    return round(mol_weight_kda, 2)

def stability(mRNA):                                #4 stability
    mRNA_len = len(mRNA)
    GC = (mRNA.count("G") + mRNA.count("C"))/mRNA_len*100  
    
    if GC < 40:     
        return "low"
    elif GC <=60 and GC > 40:               
        return "Moderate"  
    else:                                    
        return "high"

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


#codon expression-

def codon_bar_chart(codon_list):    
    
    codons = [
        "UUU", "UUC", "UUA", "UUG", 
        "CUU", "CUC", "CUA", "CUG", 
        "AUU", "AUC", "AUG", "AUA", 
        "GUU", "GUC", "GUA", "GUG", 
        "UCU", "UCC", "UCA", "UCG", 
        "CCU", "CCC", "CCA", "CCG", 
        "ACU", "ACC", "ACA", "ACG", 
        "GCU", "GCC", "GCA", "GCG", 
        "UAU", "UAC", "UAA", "UAG", 
        "CAU", "CAC", "CAA", "CAG", 
        "AAU", "AAC", "AAA", "AAG", 
        "GAU", "GAC", "GAA", "GAG", 
        "UGU", "UGC", "UGA", "UGG", 
        "CGU", "CGC", "CGA", "CGG", 
        "AGU", "AGC", "AGA", "AGG", 
        "GGU", "GGC", "GGA", "GGG"
    ]

    codon_names = []
    codon_percentage = [] 

    total_codons = len(codon_list)
    for codon in codons:
        count = codon_list.count(codon)
        percentage = (count/total_codons)*100
        codon_names.append(codon)
        codon_percentage.append(round(percentage,2))

    bar_color = "#8F87F1"
    fig,ax = plt.subplots(figsize=(17,4))
    bars = ax.bar(codon_names,codon_percentage,color=bar_color)
    

    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width()/2,
            height,
            f'{round(height,1)}%',
            ha='center',
            va='bottom',
            fontsize=6,
            color='black' 
            
        )       
    plt.xlabel('Codons')
    plt.ylabel('Percentage (%)')

    plt.xticks(rotation=90)
    plt.tight_layout()
    
    img = io.BytesIO()
    plt.savefig(img, format='jpg',dpi=120)
    img.seek(0)
    img_base64 = base64.b64encode(img.getvalue()).decode('utf-8')
    return img_base64



#try
