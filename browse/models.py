from django.core.validators import MaxValueValidator, MinValueValidator

from django.db import models
from djangophylocore.models import Taxonomy

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class Histone(models.Model):
    id             = models.CharField(max_length=25, primary_key=True)
    taxonomic_span = models.CharField(max_length=25)
    description    = models.CharField(max_length=255)

    def __unicode__(self):
        return self.id

    def get_absolute_url(self):
        from django.core.urlresolvers import reverse
        return reverse('browse.views.browse_variants', args=[str(self.id)])

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

    def get_absolute_url(self):
        from django.core.urlresolvers import reverse
        return reverse('browse.views.browse_variant', args=[str(self.core_type.id), str(self.id)])

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

    @property
    def description(self):
        return "{}|{}|{}".format(self.id, self.taxonomy.name.split(" ")[0], self.full_variant_name)

    def to_dict(self):
        return {"id":self.description, "seq":self.sequence}

    def to_biopython(self, ungap=False):
        seq = Seq(self.sequence)
        try:
            score_desc = self.all_model_scores.fiter(used_for_classifiation=True).first().description()
        except:
            score_desc = ""
        if ungap:
            seq = seq.ungap("-")
        return SeqRecord(
            seq, 
            id=self.description,
            description=score_desc,
            )

    def format(self, format="fasta", ungap=False):
        return self.to_biopython(ungap=ungap).format(format)
    

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

    def __unicode__(self):
        return "<{} variant={}; score={}; above_threshold={}; used_for_classifiation={} >".format(self.sequence.id, self.variant.id, self.score, self.above_threshold, self.used_for_classifiation)

    def description(self):
        return "[Score: {}; Evalue:{}]"

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

    def __unicode__(self):
        """Returns Jalview GFF format"""
        
        features = [
            ("alphaN", "", "helix"),
            ("alpha1", "", "helix"),
            ("alpha1ext", "", "helix"),
            ("alpha2", "", "helix"),
            ("alpha3", "", "helix"),
            ("alpha3ext", "", "helix"),
            ("alphaC", "", "helix"),
            ("beta1", "", "stand"),
            ("beta2", "", "stand"),
            ("beta3", "", "stand"),
            ("loopL1", "", "loop"),
            ("loopL2", "", "loop"),
            ("mgarg1", "Minor Groov Arg 1", "residue"),
            ("mgarg2", "Minor Groov Arg 2", "residue"),
            ("mgarg3", "Minor Groov Arg 3", "residue"),
            ("docking_domain", "Docking Domain", "domain")
        ]
        outl = ""
        for feature, description, type in features:
            start = str(getattr(self, "{}_start".format(feature), -1))
            end = str(getattr(self, "{}_end".format(feature), -1))

            if "-1" in (start, end):
                continue

            description = description or feature

            id = self.sequence.description
            outl += "\t".join((description, id, "-1", start, end, type))
            outl += "\n"
        return outl

    @classmethod
    def gff_colors(cls):
        return """domain\tred
chain\t225,105,0
residue\t105,225,35
helix\tff0000
strand\t00ff00
loop\tcccccc
"""

class Structure(models.Model):
    sequence = models.OneToOneField(Sequence, primary_key=True, related_name="structures")
    pdb      = models.CharField(max_length=25)
    mmdb     = models.CharField(max_length=25)
    chain    = models.CharField(max_length=25)

class Publication(models.Model):
    id       = models.IntegerField(primary_key=True) #PubmedID
    variants = models.ManyToManyField(Variant)
    cited    = models.BooleanField() 
    
