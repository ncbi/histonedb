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
    unique       = forms.BooleanField()

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
        self.fields['unique'].help_text = "Only show unique sequences where no organism has multiple sequences that are identical"

class FeatureForm(ModelForm):
    class Meta:
        model = Features
        exclude = ('reviewed',)

SEARCH_CHOICES = (("blastp", "Search HistoneDB (blastp)"), ("hmmer", "Classify your sequence(s) (hmmer)"))
class AnalyzeFileForm(forms.Form):
    type = forms.ChoiceField(widget=forms.RadioSelect, choices=SEARCH_CHOICES, initial=SEARCH_CHOICES[0][0])
    sequences = forms.CharField(widget=forms.Textarea)
    file  = forms.FileField()

    def __init__(self, *args, **kwargs):
        super(AnalyzeFileForm, self).__init__(*args, **kwargs)
        self.fields['sequences'].help_text = "Max 50 Sequences. However, running > 1 sequence in blastp will not be meaningful."
