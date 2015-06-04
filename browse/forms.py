from django import forms
from django.forms import ModelForm
from browse.models import Sequence, Features

class SearchForm(ModelForm):
    class Meta:
        model = Sequence
        exclude = ('reviewed',)

class FilterForm(ModelForm):
    class Meta:
        model = Sequence
        exclude = ('reviewed', 'sequence')

class FeatureForm(ModelForm):
    class Meta:
        model = Features
        exclude = ('reviewed',)

SEARCH_CHOICES = (("blastp", "Search HistoneDB (blastp)"), ("hmmer", "Classify your sequence(s) (hmmer)"))
class UploadFileForm(forms.Form):
    title = forms.CharField(max_length=50)
    type = forms.ChoiceField(widget=forms.RadioSelect, choices=SEARCH_CHOICES)
    sequences = forms.CharField(widget=forms.Textarea)
    file  = forms.FileField()

