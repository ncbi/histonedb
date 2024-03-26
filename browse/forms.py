from django import forms
from django.forms import ModelForm
from browse.models import Histone, Sequence


# from autocomplete_light import ModelForm as AutoModelForm

class AdvancedFilterForm(ModelForm):
    hist_type = forms.ModelChoiceField(Histone.objects.all())
    taxonomy = forms.CharField(max_length=50)
    sequence = forms.CharField(max_length=50)
    score = forms.FloatField()
    evalue = forms.FloatField()
    refseq = forms.BooleanField()
    unique = forms.BooleanField()

    tabs = {"Basic": ["hist_type", "variant", "taxonomy"]}
    tabs_order = ["Basic", "Advanced"]

    class Meta:
        model = Sequence
        fields = ["id", "hist_type", "variant", "taxonomy", "header", "sequence"]  # "gene", "splice"],

    def __init__(self, *args, **kwargs):
        super(AdvancedFilterForm, self).__init__(*args, **kwargs)
        self.fields['id'].label = "Accession"
        self.fields['id'].tab = "Advanced"
        self.fields['hist_type'].label = "Histone"
        self.fields['hist_type'].tab = "Basic"
        self.fields[
            'variant'].help_text = "Structurally distinct monophyletic clade of a histone family. Enter new or old variant names. Supports new dot (.) syntax."
        self.fields['variant'].tab = "Basic"
        self.fields['taxonomy'].help_text = "Supports all NCBI Taxonomy names and IDs"
        self.fields['taxonomy'].tab = "Basic"
        # self.fields['gene'].help_text = "Gene number, or more specially phylogenetic branch point order"
        # self.fields['gene'].tab = "Advanced"
        # self.fields['splice'].help_text = "Splice isoform index, or paralogous sequence number"
        # self.fields['splice'].tab = "Advanced"
        self.fields['header'].tab = "Advanced"
        self.fields['header'].help_text = "Current Genbank Description"
        self.fields['sequence'].help_text = "Search for sequence motifs"
        self.fields['sequence'].tab = "Advanced"
        self.fields['score'].tab = "Advanced"
        self.fields['score'].help_text = "Bitscore from HMM used for classification"
        self.fields['evalue'].label = "E-value"
        self.fields['evalue'].tab = "Advanced"
        self.fields['refseq'].help_text = "Only show sequences from RefSeq"
        self.fields['refseq'].tab = "Advanced"
        self.fields[
            'unique'].help_text = "Only show unique sequences where no organism has multiple sequences that are identical"
        self.fields['unique'].tab = "Advanced"

        self.tabs["Advanced"] = [field for field in self.fields if field not in self.tabs["Basic"]]


class AnalyzeFileForm(forms.Form):
    sequence = forms.CharField(widget=forms.Textarea, required=False)
    file = forms.FileField(required=False)

    def __init__(self, *args, **kwargs):
        super(AnalyzeFileForm, self).__init__(*args, **kwargs)
        self.fields['sequence'].help_text = "Input one protein sequence in FASTA format."
        self.fields['file'].label = "Or, upload file"
