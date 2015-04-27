from django import forms
from django.forms import ModelForm
from server.models import Sequence, Features

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