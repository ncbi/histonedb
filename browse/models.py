from django.core.validators import MaxValueValidator, MinValueValidator

from django.db import models
from djangophylocore.models import Taxonomy

class Histone(models.Model):
	id             = models.CharField(max_length=25, primary_key=True)
	taxonomic_span = models.CharField(max_length=25)
	description    = models.CharField(max_length=255)

	def __unicode__(self):
		return self.id

class Variant(models.Model):
	"""Most variants map to 
	H2A.X -> multiple species, same varaint
	H2A.10 -> one species, different varaint that are species speficific
	"""
	id            = models.CharField(max_length=25, primary_key=True)
	core_type     = models.ForeignKey(Histone, related_name="variants")
	taxonmic_span = models.CharField(max_length=25) #models.ForeignKey(Taxonomy)?
	description   = models.CharField(max_length=255)
	hmmthreshold  = models.FloatField(null=True)
	aucroc        = models.IntegerField(null=True)

	def __unicode__(self):
		return self.id

class OldStyleVariant(models.Model):
	updated_variant = models.ForeignKey(Variant, related_name="old_names")
	name            = models.CharField(max_length=255, primary_key=True)
	gene            = models.IntegerField(null=True, validators=[MaxValueValidator(15),MinValueValidator(1)])
	splice          = models.IntegerField(null=True, validators=[MaxValueValidator(15),MinValueValidator(1)])
	taxonomy        = models.ForeignKey(Taxonomy, related_name="+")

	def __unicode__(self):
		return "{} (now called {})".format(self.name, self.updated_variant.id)

class Sequence(models.Model):
	id       = models.CharField(max_length=25, primary_key=True) #GI
	variant  = models.ForeignKey(Variant, related_name="sequences")
	gene     = models.IntegerField(null=True, validators=[MaxValueValidator(15),MinValueValidator(1)])
	splice   = models.IntegerField(null=True, validators=[MaxValueValidator(15),MinValueValidator(1)]) 
	taxonomy = models.ForeignKey(Taxonomy)
	header   = models.CharField(max_length=255)
	sequence = models.TextField()
	reviewed = models.BooleanField()

	def __unicode__(self):
		return "{} [Varaint={}; Organism={}]".format(self.id, self.full_variant_name, self.taxonomy.name)

	@property
	def full_variant_name(self):
		name = self.variant.id
		if self.gene:
			name += ".{}".format(self.gene)
		if self.splice:
			name += ".s{}".format(self.splice)
		return name

class Score(models.Model):
	id                     = models.IntegerField(primary_key=True)
	sequence               = models.ForeignKey(Sequence, related_name="all_model_scores")
	variant                = models.ForeignKey(Variant, related_name="+")
	above_threshold        = models.BooleanField()
	score                  = models.FloatField()
	evalue                 = models.FloatField()
	hmmStart               = models.IntegerField()
	hmmEnd                 = models.IntegerField()
	seqStart               = models.IntegerField()
	seqEnd                 = models.IntegerField()
	used_for_classifiation = models.BooleanField()

class Features(models.Model):
	sequence             = models.OneToOneField(Sequence, primary_key=True, related_name="features") 
	alphaN_start         = models.IntegerField(null=True)
	alphaN_end           = models.IntegerField(null=True)
	alpha1_start         = models.IntegerField(null=True)
	alpha1_end           = models.IntegerField(null=True)
	alpha1ext_start      = models.IntegerField(null=True)
	alpha1ext_end        = models.IntegerField(null=True)
	alpha2_start         = models.IntegerField(null=True)
	alpha2_end           = models.IntegerField(null=True)
	alpha3_start         = models.IntegerField(null=True)
	alpha3_end           = models.IntegerField(null=True)
	alpha3ext_start      = models.IntegerField(null=True)
	alpha3ext_end        = models.IntegerField(null=True)
	alphaC_start         = models.IntegerField(null=True)
	alphaC_end           = models.IntegerField(null=True)
	beta1_start          = models.IntegerField(null=True)
	beta1_end            = models.IntegerField(null=True)
	beta2_start          = models.IntegerField(null=True)
	beta2_end            = models.IntegerField(null=True)
	loopL1_start         = models.IntegerField(null=True)
	loopL1_end           = models.IntegerField(null=True)
	loopL2_start         = models.IntegerField(null=True)
	loopL2_end           = models.IntegerField(null=True)
	mgarg1_start         = models.IntegerField(null=True)
	mgarg1_end           = models.IntegerField(null=True)
	mgarg2_start         = models.IntegerField(null=True)
	mgarg2_end           = models.IntegerField(null=True)
	mgarg3_start         = models.IntegerField(null=True)
	mgarg3_end           = models.IntegerField(null=True)
	docking_domain_start = models.IntegerField(null=True)
	docking_domain_end   = models.IntegerField(null=True)
	core                 = models.FloatField()

class Structure(models.Model):
	sequence = models.OneToOneField(Sequence, primary_key=True, related_name="structures")
	pdb      = models.CharField(max_length=25)
	mmdb     = models.CharField(max_length=25)
	chain    = models.CharField(max_length=25)

class Publication(models.Model):
	id       = models.IntegerField(primary_key=True) #PubmedID
	variants = models.ManyToManyField(Variant)
	cited    = models.BooleanField() 
	
