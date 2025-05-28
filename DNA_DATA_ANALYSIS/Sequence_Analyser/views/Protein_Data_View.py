from django.template import loader
from django.http import HttpResponse
# from django.shortcuts import render
from ..utils.protein_helper import *

def protein(request):
    template = loader.get_template('protein_data.html')
    Dna_Seq = request.session.get('Dna_Seq')
    mRNA = Dna_Seq.replace("T", "U")
    
    mRNA_Len = mrna_len(mRNA)
    Amino_Acid_Len = amino_acid_len(mRNA)
    # codon_lists(mRNA)
    # codon_len(mRNA)
    
    amino_acids, percentage = amino_acid_rows(mRNA)
    print("amino_acids------in-view--called---->>",amino_acids)
    print("percentage------in-view--called---->>",percentage)
    
    graph = protein_bar_chart(amino_acids, percentage)

    AA_mRNA = amino_acid_mRNA(mRNA)
    Molecular_Weight = mol_weight(AA_mRNA)
    Isoelectric_Point = isoelectric_point(AA_mRNA)
    atomic_composition = calculate_atomic_composition(AA_mRNA)
    C = atomic_composition['C']
    H = atomic_composition['H']
    N = atomic_composition['N']
    O = atomic_composition['O']
    S = atomic_composition['S']

    Molecular_Formula = C,H,N,O,S

    total_atom = C+H+N+O+S

    Atomic_Percent = atomic_percent(Molecular_Formula, total_atom)
    
    negative_count,positive_count = calculate_charge_residues(AA_mRNA)
    Calculate_Gravy = calculate_gravy(AA_mRNA)
    Calculate_Net_Charge = calculate_net_charge(AA_mRNA)
    Isoelectric_Point = isoelectric_point(AA_mRNA)

    context={
        "mRNA_LENGTH"      : mRNA_Len,
        "PROTEIN_LENGTH"   : Amino_Acid_Len,
        "Molecular_Formula": Molecular_Formula,
        "MOLECULAR_WEIGHT" : Molecular_Weight,
        "TOTAL_ATOM"       : total_atom,
        "Atomic_Percent"   : Atomic_Percent,
        "Amino_Acid_seq"   : AA_mRNA,
        "negative_count"   : negative_count,
        "positive_count"   : positive_count,
        "gravy_score"      : Calculate_Gravy,
        "net_charge"       : Calculate_Net_Charge,
        "PI"               : Isoelectric_Point,
        "Graph"            : graph,

        "show_sidebar"     : True
         
    }
    return HttpResponse(template.render(context, request))
    # return render(request, 'protein_data.html')
