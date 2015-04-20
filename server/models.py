from django.db import models

class Sequence(models.Model):
	variant              = models.ForeignKey("Variant")
	gene                 = models.IntegerField()
	splice               = models.IntegerField() 
	GI                   = models.CharField(max_length=25)
	species              = models.ForeignKey("Taxon")
	score                = models.FloatField() 
	evalue               = models.FloatField()
	header               = models.CharField(max_length=255)
	program              = models.CharField(max_length=25)
	sequence             = models.TextField() 
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
	reviewed             = models.BooleanField()
 
class Taxon(models.Model):
	id       = models.IntegerField(primary_key=True)
	species  = models.CharField(max_length=255)
	genus    = models.CharField(max_length=255)
	family   = models.CharField(max_length=255)
	phylum   = models.CharField(max_length=255)
	domain   = models.CharField(max_length=255)
	distance = models.FloatField()

class Histone(models.Model):
	core_type     = models.CharField(max_length=25)
	taxonomic_span = models.CharField(max_length=25)
	description   = models.CharField(max_length=255)

class Variant(models.Model):
	"""Most variants map to 
	H2A.X -> multiple species, same varaint
	H2A.10 -> one species, different varaint that are species speficific
	"""
	core_type     = models.ForeignKey("Histone")
	taxonmic_span = models.CharField(max_length=25)
	description   = models.CharField(max_length=255)  