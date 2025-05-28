from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render
from ..forms.upload_file_form import UploadFileForm
import pyautogui

# Get screen resolution
screen_width, screen_height = pyautogui.size()
print(f"Screen resolution:---------->>> {screen_width} x {screen_height}")

def read_txt(DNA_file):
    dna_seq = DNA_file.read().decode('utf-8').upper()
    return dna_seq

def upload_file(request):
    if request.method == "POST":
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            uploaded_file = request.FILES.get('file_field')
            if uploaded_file:
                Dna_Seq = read_txt(uploaded_file)
                print("DNA_Seq---------->>", Dna_Seq)
            request.session["Dna_Seq"] = Dna_Seq
            return HttpResponseRedirect("/dna_data/")
    else:
        form = UploadFileForm()
    return render(request, "index.html", {"form": form, "show_sidebar": False, "Screen_Width": screen_width}) 
