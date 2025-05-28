from django import forms 

class UploadFileForm(forms.Form): 
	file_field = forms.FileField() 


# class DNAForm(forms.Form):
#     file_field = forms.FileField()

#     def clean_file_field(self):
#         file = self.cleaned_data.get('file_field')

#         if not file:
#             raise forms.ValidationError("No file uploaded!")

#         return file
