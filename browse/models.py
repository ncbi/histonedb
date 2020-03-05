from django.core.validators import MaxValueValidator, MinValueValidator
from django.conf import settings
from django.db import models
from djangophylocore.models import Taxonomy

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from collections import defaultdict
import os

class Histone(models.Model):
    id             = models.CharField(max_length=25, primary_key=True)
    taxonomic_span = models.CharField(max_length=100)
    description    = models.CharField(max_length=1000)

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
    hist_type     = models.ForeignKey(Histone, related_name="variants")
    taxonomic_span = models.CharField(max_length=100) #models.ForeignKey(Taxonomy)?
    description   = models.CharField(max_length=1000)
    hmmthreshold  = models.FloatField(null=True) # parameter used in hmmersearch during sequence annotation
    aucroc        = models.IntegerField(null=True) # another parameter - these paramters are calculated during testing phase of manage.py buildvariants

    def __unicode__(self):
        return self.id

    def get_absolute_url(self):
        from django.core.urlresolvers import reverse
        return reverse('browse.views.browse_variant', args=[str(self.hist_type.id), str(self.id)])

#This is to handle other names for the same variants.like cenH3, CENPA, etc.
class OldStyleVariant(models.Model):
    updated_variant = models.ForeignKey(Variant, related_name="old_names")
    name            = models.CharField(max_length=255, primary_key=True)
    gene            = models.IntegerField(null=True, validators=[MaxValueValidator(15),MinValueValidator(1)])
    splice          = models.IntegerField(null=True, validators=[MaxValueValidator(15),MinValueValidator(1)])
    taxonomy        = models.ForeignKey(Taxonomy, related_name="+")

    def __unicode__(self):
        return "{} (now called {})".format(self.name, self.updated_variant.id)

class TemplateSequence(models.Model):
    variant  = models.CharField(max_length=255) #Not a foreign key; Maybe it is "General". It is just used to specify path
    taxonomy = models.ForeignKey(Taxonomy)

    def __unicode__(self):
        return "{}_{}".format(self.variant, self.taxonomy.name)

    def path(self):
        return os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast", "{}.fasta".format(str(self)))

    def get_sequence(self):
        return SeqIO.parse(self.path(), "fasta").next()


class Sequence(models.Model):
    id       = models.CharField(max_length=255, primary_key=True) #GI, superseeded by ACCESSION
    # id       = models.CharField(max_length=255, primary_key=True,db_index=True) #GI, superseeded by ACCESSION
    variant  = models.ForeignKey(Variant, related_name="sequences")
    gene     = models.IntegerField(null=True, validators=[MaxValueValidator(15),MinValueValidator(1)])
    splice   = models.IntegerField(null=True, validators=[MaxValueValidator(15),MinValueValidator(1)]) 
    taxonomy = models.ForeignKey(Taxonomy)
    header   = models.CharField(max_length=255)
    sequence = models.TextField()
    reviewed = models.BooleanField() 

    # class Meta:
    #     ordering=['taxonomy__name']

    def __eq__(self, other):
        if isinstance(other, Sequence):
            return (self.id == other.id)
        else:
            return False

    def __unicode__(self):
        return self.format() #"{} [Varaint={}; Organism={}]".format(self.id, self.full_variant_name, self.taxonomy.name)

    @property
    def gi(self):
        return self.id

    @property
    def full_variant_name(self):
        try:
            name = self.variant.id
        except:
            name = ""
        if self.gene:
            name += ".{}".format(self.gene)
        if self.splice:
            name += ".s{}".format(self.splice)
        return name

    @property
    def description(self):
        desc = self.id
        try:
            desc += "|{}".format(self.taxonomy.name.split(" ")[0])
        except:
            pass

        if self.full_variant_name:
            desc += "|{}".format(self.full_variant_name)

        return desc

    @property
    def short_description(self):
        return self.long_to_short_description(self.description)

    @staticmethod
    def long_to_short_description(desc):
        try:
            gi,tax,var=desc.replace("canonical","ca").split('|')
            return "{0}..{1}|{2:<.10}..|{3}".format(gi[0:2],gi[-2:],tax,var)
        except:
            return desc

    def to_dict(self, id=False, ref=False):
        return {"name":self.description if not id else self.id, "seq":self.sequence, "ref":ref}

    def to_biopython(self, ungap=False):
        seq = Seq(self.sequence)
        try:
            score_desc = self.all_model_scores.fiter(used_for_classification=True).first().description()
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
    """
    The score class, assigns a bunch of score entries to the sequence. For each variant a score.
    """
    # id                      = models.IntegerField(primary_key=True)
    # id = models.AutoField(primary_key=True)
    sequence                = models.ForeignKey(Sequence, related_name="all_model_scores")
    variant                 = models.ForeignKey(Variant, related_name="+")
    above_threshold         = models.BooleanField()
    score                   = models.FloatField()
    evalue                  = models.FloatField()
    hmmStart                = models.IntegerField()
    hmmEnd                  = models.IntegerField()
    seqStart                = models.IntegerField()
    seqEnd                  = models.IntegerField()
    used_for_classification = models.BooleanField()
    regex                   = models.BooleanField()

    def __unicode__(self):
        return "<{} variant={}; score={}; above_threshold={}; used_for_classification={} >".format(self.sequence.id, self.variant.id, self.score, self.above_threshold, self.used_for_classification)

    def description(self):
        return "[Score: {}; Evalue:{}]"

class FeatureManager(models.Manager):
    def from_dict(self, template, features, save=False):
        """Create model from secondary structure dictionary

        Parameters:
        -----------
        sequence : Sequence
        ss_dict : dict
            Created from tools.hist_ss
        """
        objs = [Feature(
            template=template,
            start=start,
            end=end,
            name=name,
            description="",
            color="") for name, (start, end) in features]
        
        if save:
            [o.save() for o in objs]

        return objs

    def gff(self, sequence_label="Consensus", features=None):
        """#Old colors:
domain\t990099
chain\t225,105,0
residue\t105,225,35
helix\tff0000
strand\t00ff00
loop\tcccccc
extension\tffff66
"""
        #assert isinstance(sequence_label, str) or , "Sequence label must be a string, not {}".format(str(type(sequence_label)))
        colors = {}
        gff_features = ""
        if features is None:
            features =  self.all()
        for feature in features:
            try:
                colorName = colors[feature.color[1:]]
            except KeyError:
                colorName = "color{}".format(len(colors))
                colors[feature.color[1:]] = colorName
            gff_features += feature.gff(str(sequence_label), colorName)
        gff_colors = "\n".join(["{}\t{}".format(name, color) for color, name in colors.iteritems()])

        return "{}\n{}".format(gff_colors, gff_features)

    def to_dict(self, features=None):
        if features is None:
            features =  self.all()
        
        return {feature.name:(feature.start, feature.end) for feature in features}

class Feature(models.Model):
    id          = models.CharField(max_length=255, primary_key=True)
    template    = models.ForeignKey(TemplateSequence, null=True)
    start       = models.IntegerField()
    end         = models.IntegerField()
    name        = models.CharField(max_length=600)
    description = models.CharField(max_length=600)
    color       = models.CharField(max_length=25)
    objects     = FeatureManager()

    class Meta:
        ordering = ["start"]

    def __unicode__(self):
        """Returns Jalview GFF format"""
        return self.gff(str(self.template))

    def gff(self, sequence_label=None, featureType="{}"):
        tmp = ""
        if sequence_label is None:
            tmp += "color1\t{}\n".format(self.color)
            sequence_label = str(self.template)
            featureType = "color1"

        tmp += "\t".join((self.name, sequence_label, "-1", str(self.start), str(self.end), featureType))
        tmp += "\n"
        return tmp

class Publication(models.Model):
     id       = models.IntegerField(primary_key=True) #PubmedID
     variants = models.ManyToManyField(Variant)
     cited    = models.BooleanField()

