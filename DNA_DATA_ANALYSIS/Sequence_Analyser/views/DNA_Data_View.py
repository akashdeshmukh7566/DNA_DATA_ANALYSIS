from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils import IsoelectricPoint as IP
from django.template import loader
from django.http import HttpResponse
from ..utils.dna_helper import *

def dna(request):
    # breakpoint()
    Dna_Seq = request.session.get('Dna_Seq')
    print("Dna_Seq_from_session----------->>",Dna_Seq)
    template = loader.get_template('dna_data.html')

    dna_length = dna_len(Dna_Seq)
    A,T,G,C = get_atgc(Dna_Seq)
    MOLECULAR_WEIGHt = Molecular_weight(Dna_Seq)
    MELTING_TEMPERATURE = melt_temp(Dna_Seq)
    Stability = dna_stability(Dna_Seq)
    GC_percentage = get_GC(Dna_Seq)
    COMPLEMENT_SEQUENCE = complement_seq(Dna_Seq)
    REVERSE_COMPLEMENT = reverse_complement_seq(Dna_Seq)

    context = {
        "Dna_len" : dna_length,        
        "DNA_seq" : Dna_Seq,
        "A"       : A,
        "T"       : T,
        "G"       : G,
        "C"       : C,
        "MOLECULAR_WEIGHt"    : MOLECULAR_WEIGHt,
        "MELTING_TEMPERATURE" : MELTING_TEMPERATURE,
        "Stability"           : Stability,
        "GC_percentage"       : GC_percentage,
        "COMPLEMENT_SEQUENCE" : COMPLEMENT_SEQUENCE,    
        "REVERSE_COMPLEMENT"  : REVERSE_COMPLEMENT,  
        "show_sidebar"        : True

    }      
    return HttpResponse(template.render(context,request))