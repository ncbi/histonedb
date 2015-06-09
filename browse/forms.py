from django import forms
from django.forms import ModelForm
from browse.models import Histone, Variant, Sequence, Features
from browse.search import search_types

class AdvancedFilterForm(ModelForm):
    core_histone = forms.ModelChoiceField(Histone.objects.all())
    taxonomy     = forms.CharField(max_length=50)
    sequence     = forms.CharField(max_length=50)
    score        = forms.FloatField()
    evalue       = forms.FloatField()

    class Meta:
        model = Sequence
        fields = ["id", "core_histone", "variant", "gene", "splice", "taxonomy", "header", "sequence"]

    def __init__(self, *args, **kwargs):
        super(AdvancedFilterForm, self).__init__(*args, **kwargs)
        self.fields['id'].label = "GI"
        self.fields['core_histone'].label = "Histone"
        self.fields['evalue'].label = "E-value"
        self.fields['variant'].help_text = "Structurally distinct monophyletic clade of a histone family. Enter new or old variant names. Supports new dot (.) syntax."
        self.fields['taxonomy'].help_text = "Supports all NCBI Taxonomy names and IDs"
        self.fields['gene'].help_text = "Gene number, or more specially phylogenetic branch point order"
        self.fields['splice'].help_text = "Splice isoform index, or paralogous sequence number"

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

