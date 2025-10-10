from django.contrib import admin
from django.urls import path
from .views.DNA_Upload_View import upload_file
from .views.DNA_Data_View import dna
from .views.mRNA_Data_View import mRNA
from .views.Protein_Data_View import protein
from django.conf import settings
# from django.conf.urls import static

urlpatterns = [
    path('admin/', admin.site.urls),
    path("", upload_file, name="upload_file"),
    path("dna_data/", dna, name="dna_data"),
    path("mRNA_data/", mRNA, name="mRNA_data"),
    path("protein_data/", protein, name="protein_data"),
]

# + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

#hello