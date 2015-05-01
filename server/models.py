from django.db import models
from server.phylocore_models import * #includes Taxonomy

class Histone(models.Model):
	id             = models.CharField(max_length=25, primary_key=True)
	taxonomic_span = models.CharField(max_length=25)
	description    = models.CharField(max_length=255)

class Variant(models.Model):
	"""Most variants map to 
	H2A.X -> multiple species, same varaint
	H2A.10 -> one species, different varaint that are species speficific
	"""
	id            = models.CharField(max_length=25, primary_key=True)
	core_type     = models.ForeignKey(Histone, related_name="variants")
	taxonmic_span = models.CharField(max_length=25)
	description   = models.CharField(max_length=255)

class Sequence(models.Model):
	id       = models.CharField(max_length=25, primary_key=True) #GI
	variant  = models.ForeignKey(Variant, related_name="sequences")
	gene     = models.IntegerField(null=True, validators=[MaxValueValidator(15),MinValueValidator(1)])
	splice   = models.IntegerField(null=True, validators=[MaxValueValidator(15),MinValueValidator(1)]) 
	taxonomy = models.ForeignKey(Taxonomy)
	header   = models.CharField(max_length=255)
	sequence = models.TextField()
	reviewed = models.BooleanField()

class Score(models.Model):
	sequence    = models.ForeignKey(Sequence, related_name="scores")
	variant     = models.ForeignKey(Variant, related_name="+")
	score       = models.FloatField() #Multiple scores for different models? Have the variant defined by the top scoring?
	evalue      = models.FloatField()
	program     = models.CharField(max_length=25)
	specificity = models.IntegerField(null=True, validators=[MaxValueValidator(100),MinValueValidator(1)])
	hmmStart    = models.IntegerField()
	hmmEnd      = models.IntegerField()
	seqStart    = models.IntegerField()
	seqEnd      = models.IntegerField()

class Features(models.Model):
	sequence             = models.OneToOneField(Sequence, primary_key=True, related_name="features") 
	alphaN_start         = models.IntegerField()
	alphaN_end           = models.IntegerField()
	alpha1_start         = models.IntegerField()
	alpha1_end           = models.IntegerField()
	alpha1ext_start      = models.IntegerField()
	alpha1ext_end        = models.IntegerField()
	alpha2_start         = models.IntegerField()
	alpha2_end           = models.IntegerField()
	alpha3_start         = models.IntegerField()
	alpha3_end           = models.IntegerField()
	alpha3ext_start      = models.IntegerField()
	alpha3ext_end        = models.IntegerField()
	alphaC_start         = models.IntegerField()
	alphaC_end           = models.IntegerField()
	beta1_start          = models.IntegerField()
	beta1_end            = models.IntegerField()
	beta2_start          = models.IntegerField()
	beta2_end            = models.IntegerField()
	loopL1_start         = models.IntegerField()
	loopL1_end           = models.IntegerField()
	loopL2_start         = models.IntegerField()
	loopL2_end           = models.IntegerField()
	mgarg1_start         = models.IntegerField()
	mgarg1_end           = models.IntegerField()
	mgarg2_start         = models.IntegerField()
	mgarg2_end           = models.IntegerField()
	mgarg3_start         = models.IntegerField()
	mgarg3_end           = models.IntegerField()
	docking_domain_start = models.IntegerField()
	docking_domain_end   = models.IntegerField()
	core                 = models.FloatField()
