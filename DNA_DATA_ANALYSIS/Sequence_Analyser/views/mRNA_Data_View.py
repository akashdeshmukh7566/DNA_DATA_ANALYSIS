from django.template import loader
from django.http import HttpResponse
from ..utils.mRNA_helper import *
from Bio.SeqUtils import molecular_weight

def mRNA(request):
    Dna_Seq = request.session.get('Dna_Seq')
    mRNA = Dna_Seq.replace("T", "U")
    print("mRNA--------------->>",mRNA)
    template = loader.get_template('mRNA_data.html')
    mRNA_len = mrna_len(mRNA)
    mRNA_sequence = mRNA
    codon_length = codon_len(mRNA)
    A1,U1,G1,C1 = mrna_augc(mRNA)
    MOLECULAR_WEIGHT = mrna_molecular_weight(mRNA)
    mRNA_STABILITY = stability(mRNA)
    codon_list = codon_lists(mRNA)
    CODON_EXPRESSION = codon_bar_chart(codon_list)

    context = {
        "mRNA_len"         : mRNA_len,
        "mRNA_sequence"    : mRNA_sequence,
        "codon_len"        : codon_length,
        "A1"               : A1,
        "U1"               : U1,
        "G1"               : G1,
        "C1"               : C1,
        "MOLECULAR_WEIGHT" : MOLECULAR_WEIGHT,
        "mRNA_STABILITY"   : mRNA_STABILITY,
        "codon_list"       : codon_list,
        "Codon_Expression" : CODON_EXPRESSION,

        "show_sidebar"     : True
        }
    return HttpResponse(template.render(context, request))