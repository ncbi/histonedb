#-*- coding: utf-8 -*-

#Originally from Taxomanie's DjangoPhyloCore (https://code.google.com/p/taxomanie/)
#License: GNU GPL v3

from django.db import models
from django.conf import settings

import os

##################################################
#             TaxonomyReference                  #
##################################################

class TaxonomyReference( object ):
    """
    This object provides some useful methods in order to deal with taxonomy
    """

    #
    # Name related
    #

    def is_valid_name( self, name ):
        """
        return True if name is in the taxonomy
        return False otherwise
        """
        #if Taxonomy.objects.filter( name = name ):
        if name in TAXONOMY_TOC:
            return True
        return False

    def is_scientific_name( self, name ):
        """
        return True if name is a scientific name
        return False otherwise
        """
        if Taxonomy.objects.filter( name = name,
          type_name = 'scientific name' ).count():
            return True
        return False

    def is_homonym( self, name ):
        """
        return True if name is an homonym
        return False otherwise
        """
        if Taxonomy.objects.filter( 
          name = name, type_name = 'homonym' ).count():
            return True
        return False

    def is_synonym( self, name ):
        """
        return True if name is a synonym
        return False otherwise
        """
        if Taxonomy.objects.filter(
          name = name, type_name = 'synonym' ).count():
            return True
        return False

    def is_common( self, name ):
        """
        return True if name is a common name
        return False otherwise
        """
        if Taxonomy.objects.filter( 
          name = name, type_name = 'common').count():
            return True
        return False

    def is_bad_taxon( self, name ):
        """
        return True if name is a bad taxon name
        return False otherwise
        """
        if BadTaxa.objects.filter(name = name).count():
            return True
        elif not self.is_valid_name( name = name ):
            return True
        else:
            return False

    def get_name_from_common( self, common_name ):
        """
        return all taxa which have common_name
        """
        return Taxonomy.objects.filter( taxa_from_common__common__name = common_name )

    def get_name_from_synonym( self, synonym ):
        """
        return all taxa which have synonym
        """
        return Taxonomy.objects.filter( taxa_from_synonym__synonym__name = synonym )

    def get_name_from_homonym( self, homonym ):
        """
        return all taxa which have homonym
        """
        return Taxonomy.objects.filter( taxa_from_homonym__homonym__name = homonym )

    def get_object_from_name( self, taxon_name ):
        """
        return the django object from taxon name.

        Name is checked in this order :
          * scientific name
          * synonym
          * homonym
          * common
          * bad taxa
        if the name is not found, a RuntimeError is raised as it must never be
        append
        """
        if not self.is_valid_name( taxon_name ):
            if BadTaxa.objects.filter( name = taxon_name ).count():
                return BadTaxa.objects.get( name = taxon_name )
            raise ValueError, '%s not found in the database' % taxon_name
        else:
            return Taxonomy.objects.filter( name = taxon_name )[0]
        # We never should be here
        raise RuntimeError, 'Something very wrong appened'

    def strip_taxon_name( self, taxon_name, delimiter='_' ):
        """
        Strip Taxon name in order to keep the root name and to remove
        all user staff.
        """
        for delim in delimiter:
            taxon_name = taxon_name.replace( delim, ' ' )
        while not self.is_valid_name( taxon_name ):
            new_taxon_name = ' '.join( taxon_name.split()[:-1] )
            if new_taxon_name:
                taxon_name = new_taxon_name
            else:
                break
        return taxon_name

    def correct( self, name, guess = False ):
        # TODO TEST
        """
        return all scientific names associated to name as a query_set

        The name will be checked in the following order:
        - if name is a scientific name then return None
        - check if name is a homonym
        - check if name is a synonym
        - check if name is a common name
        - check if name is a misspell name
        if no result is found, an empty list is returned

        if guess is True, the system will try to correct the name 
        """
        if self.is_scientific_name( name ):
            return None
        else:
            homonyms_list = self.get_name_from_homonym( name )
            if homonyms_list:
                return homonyms_list
            else:
                synonym_list = self.get_name_from_synonym( name )
                if synonym_list:
                    return synonym_list
                else:
                    common_list = self.get_name_from_common( name )
                    if common_list:
                        return common_list
                    elif guess:
                        #from phylocore.spellcheck import SpellCheck
                        #splchk = SpellCheck( self.TAXONOMY.iterkeys() )
                        #return splchk.correct( name )
                        #from lib.spell import Suggestion
                        #sugg = Suggestion( settings.TAXONOMY_ENGINE )
                        #return list(sugg.correct( name ))
                        import difflib
                        possibilities = list( TAXONOMY_TOC )
                        return difflib.get_close_matches( name, possibilities, n=10, cutoff=0.7 )
                    else:
                        return [0]

    #
    # Taxa related
    #

    def __sort( self, x, y ):
        return int(y.parents.count() - x.parents.count())

    def get_common_parents( self, taxa_list ):
        # XXX a refactoriser
        """
        select * from djangophylocore_parentsrelation where parent_id IN (select parent_id from djangophylocore_parentsrelation where taxa_id = 10114 except select parent_id from djangophylocore_parentsrelation where taxa_id = 9989 ) and taxa_id = 10114;
        """
        """
        return a list of all common parents beetween 2 taxa.
        Note that taxa1 and taxa2 must be Taxa objects
        """
        if taxa_list:
            if len( taxa_list ) == 1:
                return taxa_list[0].parents
            first_taxon, second_taxon = taxa_list[:2]
            intersection = set( first_taxon.parents ).intersection( set( second_taxon.parents) ) 
            for taxon in taxa_list[2:]:
                intersection.intersection_update( set( taxon.parents) ) 
            return  sorted( list( intersection ), self.__sort  )
        return []

    def get_first_common_parent( self, taxa_list ):
        """
        return the first common parent of taxa_list
        if root is passed to the list, None is returned as root has no parents
        """
        common_parents = self.get_common_parents( taxa_list )
        if common_parents:
            return self.get_common_parents( taxa_list )[0]
        else:
            return None

    def get_interval_parents( self, taxon1, taxon2 ):
        """
        return parents list beetween 2 taxa.
        Note that taxon1 must be a parent of taxon2
        """
        assert taxon1 in taxon2.parents, "%s is not a parent of %s" % (taxon1, taxon2) 
        difference = set( taxon2.parents ).difference( set( taxon1.parents ) )
        return sorted( list( difference ), self.__sort )[:-1]

    def get_reference_graph( self, taxa_list ):
        """
        Take a taxa list, search in reference all parents names and
        return a networkx.DiGraph tree.
        """
        import networkx as NX
        tree = NX.DiGraph() 
        already_done = set([])
        for taxon in taxa_list:
            while taxon.name != 'root' and taxon not in already_done:
                tree.add_edge( taxon.parent, taxon )
                already_done.add( taxon )
                taxon = taxon.parent
        return tree

##################################################
#                 Rank                           #
##################################################

class Rank( models.Model ):
    """
    Species, genus, kingdom...
    """
    name = models.CharField( max_length = 80 )
    class Meta:
        ordering = ['name']

    def __unicode__( self ):
        return "%s" % self.name


##################################################
#                 Taxonomy                       #
##################################################

class Taxonomy( models.Model ):
    """
    A Taxonomy object represent an element in the taxonomy. It can be a
    scientific name or a synonym, an homonym or a common name. All taxa wich
    are present of the taxonomy reference (itis, ncbi...) is a Taxonomy
    object.

    This object herites from TaxonomyReference

    This object is a Django model. See the django documentation for more
    details.

    >>> mus = Taxonomy.objects.get( name = "mus" )

    # Get all scientific names
    >>> scientific_names = Taxonomy.objects.filter( type_name = "scientific name" )
    """
    name = models.CharField( max_length = 200 )#, unique = True )
    type_name = models.CharField( max_length = 50 )
    rank = models.ForeignKey( Rank, related_name = 'taxa', null = True )
    parent = models.ForeignKey( 'self', related_name = 'direct_children', null = True )
    _parents = models.ManyToManyField( 'self', through = 'ParentsRelation', related_name = 'children', symmetrical=False )
    class Meta:
        ordering = ['name']
        unique_together = ['name', 'type_name']

    def __unicode__( self ):
        return "%s (%s)" % ( self.name, self.type_name )

    def get_homonyms( self ):
        return Taxonomy.objects.filter( homonym_from_taxa__taxon = self )
    homonyms = property( get_homonyms, None, None,
        """
        return a queryset of all homonyms related to the Taxonomy object

        >>> taxon = Taxonomy.objects.get( name = 'echinops <plant>' )
        >>> taxon.homonyms.all()
        [<Taxonomy: echinops (homonym)>]
        """)

    def get_synonyms( self ):
        return Taxonomy.objects.filter( synonym_from_taxa__taxon = self )
    synonyms = property( get_synonyms, None, None,
        """
        return a queryset of all synonyms related to the Taxonomy object

        >>> taxon = Taxonomy.objects.get( name = 'mus' )
        >>> taxon.synonyms.all()
        [<Taxonomy: nannomys (synonym)>]
        """)

    def get_commons( self ):
        return Taxonomy.objects.filter( common_from_taxa__taxon = self )
    commons = property( get_commons, None, None,
        """
        return a queryset of all commons related to the Taxonomy object

        >>> taxon = Taxonomy.objects.get( name = 'mus musculus' )
        >>> taxon.commons.all()
        [<Taxonomy: house mouse (common)>, <Taxonomy: mouse (common)>]
        """)

    def get_scientific_names( self ):
        if self.type_name == 'homonym':
            return Taxonomy.objects.filter( taxa_from_homonym__homonym = self )
        elif self.type_name == 'synonym':
            return Taxonomy.objects.filter( taxa_from_synonym__synonym = self )
        elif self.type_name == 'common':
            return Taxonomy.objects.filter( taxa_from_common__common = self )
        else:
            return Taxonomy.objects.none()
    scientifics = property( get_scientific_names, None, None,
        """
        return a queryset of all scientific names related to the Taxonomy object

        >>> taxon = Taxonomy.objects.get( name = 'echinops' )
        >>> taxon.scientifics.all()
        [<Taxonomy: echinops <mammal> (scientific name)>, <Taxonomy: echinops <plant> (scientific name)>]
        """)

    def get_parents( self, regenerate = False ):
        if regenerate:
            self.regenerate_parents()    
        if settings.DATABASE_ENGINE == "sqlite3":
            parents_list =  Taxonomy.objects.extra( 
              tables = ['server_parentsrelation'],
              where = ["server_taxonomy.id = server_parentsrelation.parent_id and server_parentsrelation.taxon_id = %s"],
              params = [self.id],
              order_by = ['server_parentsrelation."index"'] )
        else:
            parents_list =  Taxonomy.objects.extra( 
              tables = ['server_parentsrelation'],
              where = ["server_taxonomy.id = server_parentsrelation.parent_id and server_parentsrelation.taxon_id = %s"],
              params = [self.id],
              order_by = ['server_parentsrelation.index'] )
        return parents_list
        #return [i.parent for i in self.parents_relation_taxas.all()]
    parents = property( get_parents, None, None,
        """
        return a queryset of all parents related to the Taxonomy object

        >>> taxon = Taxonomy.objects.get( name = 'mus musculus' )
        >>> taxon.parent
        <Taxonomy: mus (scientific name)>
        """)

    def get_children( self, regenerate = False ):
        if regenerate:
            return [i.taxon for i in self.parents_relation_parents.all()]
    children = property( get_children, None, None,
        """
        get all children related list to the Taxonomy object

        >>> taxon = Taxonomy.objects.get( name = 'mus' )
        >>> taxon.children.all()
        [<Taxonomy: mus musculus (scientific name)>]
        """)

#    def get_id_in_source( self, source_name ):
#        return self.fromsource_set.get( source__name = source_name ).taxon_id_in_source
     
    def regenerate_parents( self ):
        """
        Regenerate parents list of the taxa. This method is useful if we add
        taxa by hand to the taxonomy.
        """
        ParentsRelation.objects.filter( taxon = self ).delete() 
        if self.name != 'root':
            parent = self.parent
            index = 0
            while parent.name != 'root':
                ParentsRelation.objects.create( 
                    taxon = self,
                    parent = Taxa.objects.get( name = parent.name ),
                    index = index )
                parent = parent.parent
                index += 1
            ParentsRelation.objects.create( 
                taxon = self,
                parent = Taxa.objects.get( name = 'root'),
                index = index )


##################################################
#           ParentsRelation                      #
##################################################

class ParentsRelation( models.Model ):
    taxon = models.ForeignKey( Taxonomy, related_name='parents_relation_taxa' )
    parent = models.ForeignKey( Taxonomy, related_name = 'parents_relation_parents' )
    index = models.IntegerField()
    class Meta:
        unique_together = ('taxon', 'parent' )
        ordering = ['index']

    def __unicode__( self ):
        return "%s > %s (%s)" % (self.parent, self.taxon, self.index )

##################################################
#               Common staffs                    #
##################################################

class RelCommonTaxa( models.Model ):
    common = models.ForeignKey( Taxonomy, related_name='common_from_taxa' )
    taxon = models.ForeignKey( Taxonomy, related_name = 'taxa_from_common' )
    language = models.CharField( max_length = 80, null = True )
    class Meta:
        ordering = ['taxon','common']
        unique_together = ('common', 'taxon')

    def __unicode__( self ):
        return "%s -> (%s)" % ( self.common, self.taxon )

##################################################
#               Synonym staffs                   #
##################################################

class RelSynonymTaxa( models.Model ):
    synonym = models.ForeignKey( Taxonomy, related_name='synonym_from_taxa' )
    taxon = models.ForeignKey( Taxonomy, related_name = 'taxa_from_synonym' )
    class Meta:
        ordering = ['taxon','synonym']
        unique_together = ('synonym', 'taxon')

    def __unicode__( self ):
        return "%s -> (%s)" % ( self.synonym, self.taxon )

##################################################
#               Homonym staffs                   #
##################################################

class RelHomonymTaxa( models.Model ):
    homonym = models.ForeignKey( Taxonomy, related_name='homonym_from_taxa' )
    taxon = models.ForeignKey( Taxonomy, related_name = 'taxa_from_homonym' )
    class Meta:
        ordering = ['taxon','homonym']
        unique_together = ('homonym', 'taxon')

    def __unicode__( self ):
        return "%s -> (%s)" % ( self.homonym, self.taxon )

##################################################
#                   BadTaxa                      #
##################################################

class BadTaxa( models.Model ):
    name = models.CharField( max_length = 200, unique = True )
    nb_occurence = models.IntegerField( default = 0 )
    class Meta:
        ordering = ['name']

    def __unicode__( self ):
        return "%s (%s)" % ( self.name, self.nb_occurence )


#############################################
#                Signals                    #
#############################################

TAXONOMY_TOC = None
BADTAXA_TOC = None

def get_taxonomy_toc( test = False ):
    try:
        import cPickle as pickle
    except:
        import pickle
    global TAXONOMY_TOC
    localDir = os.path.dirname(__file__)
    absDir = os.path.join(os.getcwd(), localDir)
    if not TAXONOMY_TOC:
        if test:
            TAXONOMY_TOC = pickle.load( open( os.path.join( absDir, 'taxonomy_toc_test') ) )
        else:
            TAXONOMY_TOC = pickle.load( open( os.path.join( absDir, 'taxonomy_toc_%s' % settings.TAXONOMY_ENGINE ) ) )
    else:
        TAXONOMY_TOC = globals()['TAXONOMY_TOC']
    return TAXONOMY_TOC

