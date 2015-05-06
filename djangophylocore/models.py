#-*- coding: utf-8 -*-

from django.db import models
from django.db.models import signals, Q
from django.conf import settings
from django.db import transaction
from django.db import connection

from lib.phylogelib import NewickParser
from lib.nexus import Nexus, remove_nexus_comments
from pyparsing import ParseException
import datetime, re, sys, os, codecs

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
              tables = ['djangophylocore_parentsrelation'],
              where = ["djangophylocore_taxonomy.id = djangophylocore_parentsrelation.parent_id and djangophylocore_parentsrelation.taxon_id = %s"],
              params = [self.id],
              order_by = ['djangophylocore_parentsrelation."index"'] )
        else:
            parents_list =  Taxonomy.objects.extra( 
              tables = ['djangophylocore_parentsrelation'],
              where = ["djangophylocore_taxonomy.id = djangophylocore_parentsrelation.parent_id and djangophylocore_parentsrelation.taxon_id = %s"],
              params = [self.id],
              order_by = ['djangophylocore_parentsrelation.index'] )
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

##################################################
#           Phylogenetic Tree                    #
##################################################

class TaxonomyTreeOccurence( models.Model ):
    taxon = models.ForeignKey( Taxonomy, related_name = 'taxonomy_occurences' )
    tree = models.ForeignKey( 'Tree', related_name = 'taxonomy_occurences' )
    user_taxon_name = models.CharField( max_length = 200, null = True )
    nb_occurence = models.IntegerField( default = 0 )
    class Meta:
        unique_together = ( 'taxon', 'tree', 'user_taxon_name' )

    def __unicode__( self ):
        return u'%s (%s) %s' % ( self.taxon, self.nb_occurence, self.tree )

class Tree( models.Model, TaxonomyReference ):
    """
    This object represents a phylogenetic tree. It can be valid (well formed)
    or not. If not, the column_error attribute will indicate the index of the
    failure in the source.

    the _from_collection attribute indicates if the Tree was build from a
    TreeCollection or not.

    This object herites from TaxonomyReference

    This object is a Django model. See the django documentation for more
    details.

    >>> nwk_tree = "( echinops, (rattus, ( mus, azerty, black rat ), nannomys ))"
    >>> tree = Tree.objects.create( source = nwk_tree, name = "2")
    >>> tree.is_valid
    True
    >>> tree.scientifics.all()
    [<Taxonomy: mus (scientific name)>, <Taxonomy: rattus (scientific name)>]
    >>> tree.bad_taxa.all()
    [<BadTaxa: azerty (0)>]
    >>> tree.synonyms.all()
    [<Taxonomy: nannomys (synonym)>]
    >>> tree.homonyms.all()
    [<Taxonomy: echinops (homonym)>]
    >>> tree.commons.all()
    [<Taxonomy: black rat (common)>]

    # Getting all taxon list
    >>> tree.taxa.all()
    [<Taxonomy: black rat (common)>, <Taxonomy: echinops (homonym)>, <Taxonomy:
    mus (scientific name)>, <Taxonomy: nannomys (synonym)>, <Taxonomy: rattus
    (scientific name)>]

    # Getting ambiguous taxa (synonms, commons, homonyms...)
    >>> tree.ambiguous.all() [<Taxonomy: black rat (common)>, <Taxonomy: echinops (homonym)>, <Taxonomy: nannomys (synonym)>]
    """
    name = models.CharField( max_length = 80, null=True )
    delimiter = models.CharField( max_length = 5, default=' ' )
    source = models.TextField()
    rooted = models.NullBooleanField()
    description = models.TextField( null = True )
    created = models.DateTimeField()
    updated = models.DateTimeField()
    is_valid = models.BooleanField( default = False )
    _from_collection = models.BooleanField( default = False )
    column_error = models.IntegerField( null = True )
    #
    collection = models.ForeignKey( 'TreeCollection', related_name = 'trees', null = True)
    bad_taxa = models.ManyToManyField( BadTaxa, related_name = 'trees'  )
    taxon_ids = []

    def __unicode__( self ):
        return "%s" % ( self.name )

    @transaction.atomic
    def save( self, dont_generate = False, **kwargs ):
        """
        save the tree into the database. If dont_generate is True, the method
        will call generate_tree_info method.
        """
        regenerate = False
        if not self.id: # if instance is not in the database
            regenerate = True
        super( Tree, self ).save( **kwargs )
        if regenerate and not dont_generate:
            self.generate_tree_infos()

    def generate_tree_infos( self ):
        """
        This method will extract all taxa name from source  and wrap them
        into Taxonomy objects. Those taxonomy objects will be linked to the
        tree.
        """
        global get_taxonomy_toc, BADTAXA_TOC
        TAXONOMY_TOC = get_taxonomy_toc()
        if [i for i in ('(',')',',') if i in self.delimiter]:
            raise ValueError, '"%s" is a bad delimiter' % self.delimiter
        tree = self.source.lower()#.replace( self.delimiter, ' ' )
        self.newick_parser = NewickParser()
        try:
            self.newick_parser.parse_string( tree ) # XXX Verifier plantage
            self.is_valid = True
        except ParseException, err:
            self.column_error = err.column
        self.save( dont_generate = True )
        taxa_list = self.newick_parser.get_taxa()
        self.taxon_ids = {}
        if BADTAXA_TOC is None:
            BADTAXA_TOC = set([i[0] for i in BadTaxa.objects.all().values_list( 'name')])
        self.bad_taxon_ids = BADTAXA_TOC
        for taxon_name in taxa_list:
            if taxon_name.strip():
                user_taxon_name = taxon_name.strip()
                #VR taxon_name = self.strip_taxon_name( taxon_name.strip(), self.delimiter ).replace( "_", " " ).strip()
                taxon_name = self.strip_taxon_name( taxon_name.strip(), self.delimiter )
                taxo = TAXONOMY_TOC.get( taxon_name, '' )
                if self._from_collection:
                    if taxo:
                        self.taxon_ids[user_taxon_name] = taxo
                    else:
                        if user_taxon_name not in self.bad_taxon_ids:
                            BadTaxa.objects.create( name = user_taxon_name )
                            self.bad_taxon_ids.add( user_taxon_name )
                        self.taxon_ids[user_taxon_name] = ''
                else:
                    if not taxo:
                        t, created = BadTaxa.objects.get_or_create( name = user_taxon_name )
                        if user_taxon_name not in self.bad_taxon_ids:
                            print "creating %s in bad taxa" % user_taxon_name
                            self.bad_taxon_ids.add( t.id )
                        print "adding %s in bad taxa" % t.id
                        self.bad_taxa.add( t )
                    else:
                        taxon = Taxonomy.objects.get( id = taxo )
                        tto, created = TaxonomyTreeOccurence.objects.get_or_create( 
                          taxon = taxon, tree = self, user_taxon_name = user_taxon_name )
                        tto.nb_occurence += 1
                        tto.save()

    def __get_relation( self ):
        class Meta: pass
        model_rel =  type( 'RelTreeColTaxa%s' % self.collection.id, #name.capitalize(),
          (AbstractTreeColTaxa,), {'__module__': Taxonomy.__module__, 'Meta':Meta})
        return model_rel.objects.filter( tree = self )
    rel = property( __get_relation )

    def print_error( self ):
        """
        display where the error is in the source

        >>> tree.print_error()
        (,(mus, rattus)
         ^
        """
        print self.source
        print " "*(err.column-1) + "^"

    def get_taxa( self ):
        if not self._from_collection:
            return Taxonomy.objects.filter( taxonomy_occurences__tree = self )
        return Taxonomy.objects.extra( 
          tables = ['djangophylocore_reltreecoltaxa%s' % self.collection.id],
          where = ['djangophylocore_taxonomy.id = djangophylocore_reltreecoltaxa%s.taxon_id and djangophylocore_reltreecoltaxa%s.tree_id = %s' % ( 
            self.collection.id, self.collection.id, self.id )]
        ).distinct()
    taxa = property( get_taxa, None, None,
        """
        return a queryset of all taxa included ambiguous names and scientific names
        
        >>> tree.taxa.all()
        [<Taxonomy: black rat (common)>, <Taxonomy: echinops (homonym)>, <Taxonomy: mus (scientific name)>, <Taxonomy: nannomys (synonym)>, <Taxonomy: rattus (scientific name)>]
        """)

    def get_ambiguous( self ):
        """
        return a queryset of taxonomy objects wich are not scientific name
        """
        return self.taxa.exclude( type_name = 'scientific name' )
    ambiguous = property( get_ambiguous, None, None,
        """
        return a queryset of all ambiguous names (homonym, synonym and common names)
        
        >>> tree.ambiguous.all()
        [<Taxonomy: black rat (common)>, <Taxonomy: echinops (homonym)>, <Taxonomy: nannomys (synonym)>]
        """)

    def get_scientific_names( self ):
        return self.taxa.filter( type_name = 'scientific name' )
    scientifics = property( get_scientific_names, None, None,
        """
        return a queryset of all scientific names

        >>> tree.scientifics.all()
        [<Taxonomy: mus (scientific name)>, <Taxonomy: rattus (scientific name)>]
        """)
        
    def get_synonyms( self ):
        return self.taxa.filter( type_name = 'synonym' )
    synonyms = property( get_synonyms, None, None,
        """
        return a queryset of  all synonym names

        >>> tree.synonyms.all()
        [<Taxonomy: nannomys (synonym)>]
        """)

    def get_commons( self ):
        return self.taxa.filter( type_name = 'common' )
    commons = property( get_commons, None, None,
        """
        return a queryset of  all common names

        >>> tree.commons.all()
        [<Taxonomy: black rat (common)>]
        """)

    def get_homonyms( self ):
        return self.taxa.filter( type_name = 'homonym' )
    homonyms = property( get_homonyms, None, None,
        """
        return a queryset of all homonym names

        >>> tree.homonyms.all()
        [<Taxonomy: echinops (homonym)>]
        """)

#    def __generate_arborescence_old( self, tree=None ):
#        if tree is None:
#            # Init attributes
#            import networkx as NX
#            self.__networkx_tree = NX.DiGraph()
#            tree = self.source
#            self.__children = []
#            self.__last_child = ""
#            self.__rel_name = {}
#            self.__miss_spelled = {}
#            tree = self.newick_parser.tree
#        if tree:
#            if tree == self.source:
#                parent_name = Taxonomy.objects.get( name = 'root' )
#            else:
#                if not self.__rel_name.has_key( tree ):
#                    taxa_list = [Taxonomy.objects.filter( name = i )[0] for i in self.newick_parser.get_taxa( tree ) if self.is_valid_name( i )]
#                    parent_name = self.get_first_common_parent(taxa_list)
#                else:
#                    parent_name = self.__rel_name[tree]
#            children = tree
#            for child_name in children:
#                if type(child_name) is list: # child is a node
#                    taxa_list = [Taxonomy.objects.filter( name = i )[0] for i in self.newick_parser.get_taxa( child_name ) if self.is_valid_name( i )]
#                    child = self.get_first_common_parent(taxa_list)
#                    if child is None:
#                        child = parent_name
#                    if child not in self.__children:
#                        self.__children.append( parent_name )
#                    self.__rel_name[child_name] = child
#                    self.__children.append( child )
#                    self.__last_child = child
#                    self.__networkx_tree.add_edge( parent_name, child )
#                    self.__generate_arborescence( tree = child_name )
#                else: # child is a taxon
#                    child_name = child_name.replace( self.delimiter, ' ' )
#                    child = self.get_object_from_name( child_name )
#                    if child is None:
#                        child = parent_name
#                    self.__networkx_tree.add_edge( parent_name, child )

    def __generate_arborescence( self ):
        import networkx as NX
        tree = NX.DiGraph()
        parser = NewickParser()
        tree_list = parser.parse_string( self.source )
        current_list = parser.get_taxa()
        rel_list = self.collection.rel.filter( user_taxon_name__in = current_list )
        d_infos = {}
        for rel in rel_list:
            try:
                d_infos[rel.user_taxon_name] = rel.taxon
            except:
                d_infos[rel.user_taxon_name] = BadTaxa.objects.get( name = rel.user_taxon_name )
        self.__generate_arborescence_rec( tree, tree_list , d_infos, True )
        return tree

    def __generate_arborescence_rec( self, tree, current_list, d_infos, root=False ):
        class MockObj( object ):
            """
            object wich represent a parent of bad taxa
            """
            def __init__( self ):
                self.name = ""
                self.parents = []
        if type(current_list) is list:
            is_scientific = True
            son_list = []
            scientific_son = []
            for i in current_list[:]:
                son = self.__generate_arborescence_rec( tree, i, d_infos )
                son_list.append( son )
                if son[1]: #it's a scientific son
                    scientific_son.append( son )
                else:
                    is_scientific = False
            parent = self.get_first_common_parent( [i[2] for i in scientific_son] )
            if parent is None:
                parent = MockObj() 
            node = ( parent.name, is_scientific, parent, len(tree)+1 )
            if root:
                tree.add_edge( ("root", True, Taxonomy.objects.get(name = "root"), 0), node ) 
            for son in son_list:
                tree.add_edge( node, son )
        else:
            taxon = d_infos[current_list]
            if hasattr( taxon, "type_name" ):
                is_scientific = taxon.type_name == "scientific name"
            else:
                is_scientific = False
            node = ( current_list, is_scientific, taxon, len(tree)+1 )
            tree.add_node( node )
        return node

    def get_arborescence( self ):
        return  self.__generate_arborescence()
    arborescence = property( get_arborescence, None, None,
        """
        return a NetworkX oriented graph from the tree
        """)

    def __graph2nwk_rec( self, tree, current_node, internal_label ):
        str_current =""
        sons_of_current_node = tree.successors(current_node)
        if len(sons_of_current_node) > 1 and type( sons_of_current_node ) is list:
            str_current += "("
            for s in sons_of_current_node:
                str_current += self.__graph2nwk_rec(tree, s, internal_label) +","
            str_current = str_current[:-1] # on a une "," en trop
            str_current += ")"
            if internal_label:
                str_current += "["
                if not current_node[1]:
                    if current_node[2].name:
                        str_current +="?" + current_node[2].name + "?"
                    else:
                        str_current += "?"
                else:
                    str_current += current_node[2].name
                str_current +=  "]"
        elif len(sons_of_current_node) == 1 and type( sons_of_current_node ) is list:
            # intgernal node of degre 2
            return self.__graph2nwk_rec(tree, sons_of_current_node[0], internal_label)
        else: # it's a leave
            str_current = current_node[0]
            if not current_node[1]:
                str_current +="?"
        return str_current

    def get_tree_as_nwk( self, internal_label=False ):
        """
        return the newick format of the tree. 
        If internal_label is True, then display internal_label.
        """
        tree = self.get_arborescence()
        return self.__graph2nwk_rec( tree,
          ( "root", True, Taxonomy.objects.get( name = "root"), 0 ),
          internal_label )+";"

    def get_reference_tree_as_nwk( self, internal_label = True ):
        """
        return the reference arborescence in a newick string
        If internal_label is True, then display internal_label.
        """
        tree = self.get_reference_arborescence()
        return self.__graph2nwk_rec( tree,
          ( "root", True, Taxonomy.objects.get( name = "root"), 0 ),
          internal_label )+";"

    def get_reference_arborescence( self ):
        """
        return a netorkX oriented graph wich represente the arborescence of
        the reference taxonomy (itis, ncbi...) of the list of taxa
        """
        return self.get_reference_graph( self.scientifics.all() )

    def get_nb_taxa_from_parent( self, parent_name ):#, taxa_occurence = None ):
        """
        return the number of taxa wich have `parent_name` for parent
        """
        global get_taxonomy_toc
        TAXONOMY_TOC = get_taxonomy_toc()
        if self._from_collection:
            cursor = connection.cursor()
            #parent = Taxonomy.objects.get( name = parent_name )
            parent_id = TAXONOMY_TOC[parent_name]
           # cur = cursor.execute( "select tree_id, count(taxon_id) from djangophylocore_reltreecoltaxa%s where (taxon_id IN (select taxon_id from djangophylocore_parentsrelation where parent_id = %s) or taxon_id = %s) and tree_id = %s GROUP BY  tree_id ;" % ( self.collection.id, parent_id, parent_id, self.id ) ) 
            cur = cursor.execute( "select tree_id, count(djangophylocore_reltreecoltaxa%s.taxon_id) from djangophylocore_reltreecoltaxa%s as rel, djangophylocore_parentsrelation as par where djangophylocore_reltreecoltaxa%s.taxon_id = djangophylocore_parentsrelation.taxon_id and djangophylocore_parentsrelation.parent_id = %s and djangophylocore_reltreecoltaxa%s.tree_id = %s GROUP BY tree_id;" % ( 
                self.collection.id, self.collection.id, self.collection.id, parent_id, self.collection.id, self.id )
            ) 
            if settings.DATABASE_ENGINE == 'sqlite3':
                results = cur.fetchone()[1]
            else:
                results = cursor.fetchone()[1]
            return results
        return self.taxa.filter( parents_relation_taxa__parent__name = parent_name ).count()

    def eval_query( self, query, usertaxa_list=[] ):
        """
        test if a query match the tree. The query format is a python
        boolean expression with taxa name beetween braces :

        tree.eval_query( "{muridae} > 2 and {primates}" )

        will return true if tree have more than 2 taxa wich have muridae as parents
        and at least 1 taxon wich have a primate as parents.

        if a taxa_list is not null, the query can have another variable
        {usertaxa}. this variable represente all taxa passed in the list.

        tree.eval_query( "{muridae} => 4 and {usertaxa} > 2", ['rattus', 'mus', 'pan', 'boss'] )

        will return true if tree have at least 4 taxa wich are muridae and
        more than 2 taxa wich are in the usertaxa_list 
        """
        res = query.strip()
        for pattern in re.findall("{([^}]+)}", query):
            striped_pattern = pattern.strip().lower()
            if not self.is_valid_name( striped_pattern ):
                raise NameError, striped_pattern
            else:
                nb_occurence = self.get_nb_taxa_from_parent( striped_pattern )
            res = res.replace("{"+pattern+"}", str(nb_occurence) )
        if res:
            try:
                return eval( res )
            except:
                raise SyntaxError, "bad query 3 %s" % query
        raise SyntaxError, "bad query 4 %s" % query


##################################################
#               TreeCollection                   #
##################################################

class TreeCollection( models.Model, TaxonomyReference ):
    """
    This object represent a collection of phylogenetic trees.

    Two formats are supported : phylip and nexus. The format is automatiquely
    detected.

    This object herites from TaxonomyReference

    This object is a Django model. See the django documentation for more
    details.

    >>> simple_col = "(mus,nannomys,black rat,echinops,blabla);(mus, black rat);"
    >>> col = TreeCollection.objects.create( name = 'simple', source = simple_col )
    >>> col.format
    'phylip'

    # Deeling with taxa
    >>> col.taxa.all()
    [<Taxonomy: black rat (common)>, <Taxonomy: echinops (homonym)>, <Taxonomy: mus (scientific name)>, <Taxonomy: nannomys (synonym)>]
    >>> col.ambiguous.all()
    [<Taxonomy: black rat (common)>, <Taxonomy: echinops (homonym)>, <Taxonomy: nannomys (synonym)>]
    >>> col.scientifics.all()
    [<Taxonomy: mus (scientific name)>]
    >>> col.homonyms.all()
    [<Taxonomy: echinops (homonym)>]
    >>> col.synonyms.all()
    [<Taxonomy: nannomys (synonym)>]
    >>> col.commons.all()
    [<Taxonomy: black rat (common)>]
    >>> col.bad_taxa.all()
    [<BadTaxa: blabla (0)>]

    # working with tree
    >>> col.trees.all()
    [<Tree: 1>, <Tree: 2>]
    >>> col.trees.count()
    2L
    """
    name = models.CharField( max_length = 80, null= True )
    source = models.TextField( null = True ) # source
    delimiter = models.CharField( max_length = 5, default=' ' )
    description = models.TextField( null = True)
    format = models.CharField( max_length = 20, null = True )
    created = models.DateTimeField()
    updated = models.DateTimeField()

    def __unicode__( self ):
        return "%s (%s)" % ( self.name, self.format )

    @transaction.atomic
    def save( self, collection_changed = False, dont_regenerate = False,  **kwargs ):
        """
        save the collection into the database. If dont_regenerate is False, it
        will call the regenerate_from_source method
        """
        collection_string_changed = False
        self.source = self.source.lower()#.replace( self.delimiter, '' )
        if not collection_changed and not dont_regenerate:
            if self.id: # if instance is in the database
                if TreeCollection.objects.get( id = self.id).source != self.source:
                    collection_string_changed = True
            else:
                collection_string_changed = True
        super( TreeCollection, self ).save( **kwargs )
        if self.source and collection_string_changed and not dont_regenerate:
            self.regenerate_from_source()
            self.__get_relation().model.objects.filter( taxon = '' ).update( taxon = None )

    def regenerate_collection_string_from_trees( self ):
        """
        regenerate source from trees in the collection
        """
        self.source = self.get_collection_string()
        self.save( collection_changed = True )

    def __create_relation( self, name, dump ):
        cursor = connection.cursor()
        if settings.DATABASE_ENGINE == 'sqlite3':
            cursor.execute( """ CREATE TABLE "djangophylocore_reltreecoltaxa%s" (
                "id" integer NOT NULL PRIMARY KEY,
                "collection_id" integer NOT NULL REFERENCES "djangophylocore_treecollection" ("id"),
                "tree_id" integer NOT NULL REFERENCES "djangophylocore_tree" ("id"),
                "taxon_id" integer NULL REFERENCES "djangophylocore_taxonomy" ("id"),
                "user_taxon_name" varchar(200) NULL
            );""" % name )
        elif settings.DATABASE_ENGINE == 'mysql':
            cursor.execute( """ CREATE TABLE `djangophylocore_reltreecoltaxa%s` (
                `id` integer NOT NULL PRIMARY KEY,
                `collection_id` integer NOT NULL REFERENCES `djangophylocore_treecollection` (`id`),
                `tree_id` integer NOT NULL REFERENCES `djangophylocore_tree` (`id`),
                `taxon_id` integer NULL REFERENCES `djangophylocore_taxonomy` (`id`),
                `user_taxon_name` varchar(200) NULL
            );""" % name )
        codecs.open( '/tmp/rel_%s.dmp' % name, encoding='utf-8', mode='w').write( ''.join( dump ) )
        if settings.DATABASE_NAME == ':memory:':
            db_name = settings.TEST_DATABASE_NAME
        else:
            db_name = settings.DATABASE_NAME
        if settings.DATABASE_ENGINE == 'sqlite3':
            os.system( "sqlite3 -separator '|' %s '.import /tmp/rel_%s.dmp djangophylocore_reltreecoltaxa%s'" % (
              db_name, name, name ) )
        elif settings.DATABASE_ENGINE == 'mysql':
            cmd = """mysql -u %s -p%s %s -e "LOAD DATA LOCAL INFILE '/tmp/rel_%s.dmp' INTO TABLE djangophylocore_reltreecoltaxa%s FIELDS TERMINATED BY '|';" """ % ( settings.DATABASE_USER, settings.DATABASE_PASSWORD, db_name, name, name )
            os.system( cmd )
        else:
            raise RuntimeError, "%s engine not supported" % settings.DATABASE_ENGINE
        if settings.DATABASE_ENGINE == 'mysql':
            cursor = connection.cursor()
            cursor.execute( """CREATE INDEX djangophylocore_reltreecoltaxa%s_taxon_id ON djangophylocore_reltreecoltaxa%s (taxon_id);""" % (self.id, self.id ))
            cursor.execute( """CREATE INDEX djangophylocore_reltreecoltaxa%s_tree_id ON djangophylocore_reltreecoltaxa%s (tree_id);""" % (self.id, self.id ))
            cursor.close()
        os.system( 'rm /tmp/rel_%s.dmp' % name )

    def __get_relation( self ):
        class Meta: pass
        model_rel =  type( 'RelTreeColTaxa%s' % self.id, #name.capitalize(),
          (AbstractTreeColTaxa,), {'__module__': Taxonomy.__module__, 'Meta':Meta})
        return model_rel.objects.filter( collection = self )
    rel = property( __get_relation )

    def regenerate_from_source( self ):
        """
        This method will parse the source in order to extract all trees and
        taxa. Those objects are linked to the collection.
        """
        global BADTAXA_TOC
        BADTAXA_TOC = set([i[0] for i in BadTaxa.objects.all().values_list( 'name')])
        # TODO mettre un signal sur cette method quand
        # source change
        #if self.trees.count():
        #    Tree.objects.filter( collections = self ).delete()
        if settings.DEBUG and settings.DATABASE_ENGINE == 'sqlite3':
            if hasattr( settings, 'BOOST_SQLITE' ):
                if settings.BOOST_SQLITE:
                    cursor = connection.cursor()
                    cursor.execute('PRAGMA temp_store = MEMORY;')
                    cursor.execute('PRAGMA synchronous=OFF')
        index = 0
        dump = []
        nwk_collection = remove_nexus_comments( self.source )
        # Nexus collection
        if nwk_collection.lower().strip()[:6] == "#nexus":
            nex = Nexus( nwk_collection )
            self.format = 'nexus'
            self.save( dont_regenerate = True )
            for name, (tree, rooted) in nex.collection.iteritems():
                t = Tree( name = name, source = tree, rooted = rooted,
                  delimiter = self.delimiter, _from_collection = True,
                  collection = self )
                t.save()
                if not t.taxon_ids:
                    index += 1
                    dump.append( '%s|%s|%s||\n' % (index, self.id, t.id ) )
                else:
                    for user_taxon_name in t.taxon_ids:
                        index += 1
                        # index, collection_id, tree_id, taxon_id, user_taxon_name
                        dump.append(  '%s|%s|%s|%s|%s\n' % (index, self.id,
                          t.id, t.taxon_ids[user_taxon_name], user_taxon_name ) )
        # Phylip collection
        else:
            self.format = 'phylip'
            self.save( dont_regenerate = True )
            name = 0
            l_trees = nwk_collection.strip().split(";")
            for nwktree in l_trees:
                tree = nwktree.lower().strip()
                if tree:
                    name += 1
                    t = Tree( name = "tree_%s" % name, source = tree, rooted = False, 
                      delimiter = self.delimiter, _from_collection = True,
                      collection = self )
                    t.save()
                    if not t.taxon_ids:
                        index += 1
                        dump.append( '%s|%s|%s||\n' % (index, self.id, t.id ) )
                    else:
                        for user_taxon_name in t.taxon_ids:
                            index += 1
                            # index, collection_id, tree_id, taxon_id, user_taxon_name
                            dump.append(  '%s|%s|%s|%s|%s\n' % (index,
                              self.id, t.id, t.taxon_ids[user_taxon_name], user_taxon_name ) )
        self.__create_relation( str(self.id), dump )

    def get_collection_string( self ):
        """
        Generate from trees the collection_string into the source format.
        """
        result = []
        for (name, source, rooted) in self.trees.values_list( 'name', 'source', 'rooted' ):
            if rooted: rooted = '[&R]'
            else: rooted = ''
            if self.format == 'nexus':
                result.append( 'TREE %s = %s %s' % (name, rooted, source ) )
            else:
                result.append( source.strip() )
        if self.format == 'nexus':
            return "#NEXUS\n\nBEGIN TREES;\n\n"+";\n".join( result )+"\n\nEND;\n"
        else:
            return ";\n".join( result )+';'

    def get_taxa( self ):
        return Taxonomy.objects.extra( 
          tables = ['djangophylocore_reltreecoltaxa%s' % self.id], 
          where = ['djangophylocore_taxonomy.id = djangophylocore_reltreecoltaxa%s.taxon_id' % self.id]
        ).distinct()
    taxa = property( get_taxa, None, None,
        """
        return a queryset of all taxa from the collection (included ambiguous one)

        >>> col.taxa.all()
        [<Taxonomy: black rat (common)>, <Taxonomy: echinops (homonym)>, <Taxonomy: mus (scientific name)>, <Taxonomy: nannomys (synonym)>]
        """)

    def get_user_taxon_names( self ):
        """
        return the user taxon names list
        """
        return set([i[0] for i in self.rel.values_list( 'user_taxon_name' ) if i])

    def get_ambiguous( self ):
        return self.taxa.exclude( type_name = 'scientific name' )
    ambiguous = property( get_ambiguous, None, None,
        """
        return a queryset of non scientific name taxonomy objects
        
        >>> col.ambiguous.all()
        [<Taxonomy: black rat (common)>, <Taxonomy: echinops (homonym)>, <Taxonomy: nannomys (synonym)>]
        """)

    def get_bad_taxa( self ):
        bad_taxa = BadTaxa.objects.extra(
          tables = ['djangophylocore_reltreecoltaxa%s' % self.id],
          where = ['djangophylocore_badtaxa.name = djangophylocore_reltreecoltaxa%s.user_taxon_name' % self.id]
        ).distinct()
        return [taxa for taxa in bad_taxa if not self.is_valid_name(self.strip_taxon_name(taxa.name,self.delimiter))]
             
    bad_taxa = property( get_bad_taxa, None, None,
        """
        return a queryset of bad taxa

        >>> col.bad_taxa.all()
        [<BadTaxa: blabla (0)>]
        """)

    def get_scientific_names( self ):
        return self.taxa.filter( type_name = 'scientific name').distinct()
    scientifics = property( get_scientific_names, None, None,
        """
        return a queryset of all scientific names of the collection

        >>> col.scientifics.all()
        [<Taxonomy: mus (scientific name)>]
        """)

    def get_synonyms( self ):
        return self.taxa.filter( type_name = 'synonym').distinct()
    synonyms = property( get_synonyms, None, None,
        """
        return a queryset of all synonym names of the collection

        >>> col.synonyms.all()
        [<Taxonomy: nannomys (synonym)>]
        """)

    def get_homonyms( self ):
        return self.taxa.filter( type_name = 'homonym' ).distinct()
    homonyms = property( get_homonyms, None, None,
        """
        return a queryset of all homonym names of the collection

        >>> col.homonyms.all()
        [<Taxonomy: echinops (homonym)>]
        """)

    def get_common_names( self ):
        return self.taxa.filter( type_name = 'common' ).distinct()
    commons = property( get_common_names, None, None,
        """
        return a queryset of all common names of the collection

        >>> col.commons.all()
        [<Taxonomy: black rat (common)>]
        """)

    def get_bad_trees( self ):
        return self.trees.filter( is_valid = False )
    bad_trees = property( get_bad_trees, None, None,
        """
        return a queryset of all bad (misformed) trees present in the
        collection

        >>> col.bad_trees.all()
        [<Tree: 2>, <Tree: 4>]
        """)


    def _query( self, query, treebase ):
        """
        """
        global get_taxonomy_toc
        TAXONOMY_TOC = get_taxonomy_toc()
        res = query.strip()
        d_trees = {}
        cursor = connection.cursor()
        l_patterns = re.findall("{([^}]+)}", query)
        if not l_patterns:
            raise ValueError, "malformed request"
        for pattern in l_patterns:
            striped_pattern = pattern.strip().lower()
            if not striped_pattern == 'usertaxa' and not self.is_valid_name( striped_pattern ):
                raise NameError, striped_pattern+" not found in taxonomy"
        for pattern in l_patterns:
            striped_pattern = pattern.strip().lower()
            if not self.is_scientific_name( striped_pattern ) and striped_pattern != "usertaxa":
                related_scientific_names = Taxonomy.objects.get( name = striped_pattern ).scientifics
                if len(related_scientific_names) == 1:
                    striped_pattern = related_scientific_names[0].name
                else:
                    raise ValueError, striped_pattern+" is ambiguous : %s" % [str(i.name) for i in related_scientific_names]
            if 'usertaxa' == striped_pattern and treebase:
                cur = cursor.execute( "select tree_id, count(taxon_id) from djangophylocore_reltreecoltaxa1 where taxon_id IN (select taxon_id from djangophylocore_reltreecoltaxa%s ) GROUP BY tree_id;" % (self.id ) )
            else:
                baseId=self.id;
                if treebase:
                    baseId=1;
                parent_id = TAXONOMY_TOC[striped_pattern]
                
                if settings.DATABASE_ENGINE == 'sqlite3':
                    cur = cursor.execute( "select tree_id, count(taxon_id) from djangophylocore_reltreecoltaxa%s where taxon_id IN (select taxon_id from djangophylocore_parentsrelation where parent_id = %s) or taxon_id = %s GROUP BY  tree_id ;" % ( baseId, parent_id, parent_id ) )
                else:
                    cur = cursor.execute( "select tu.tree_id, count(tu.taxon_id) from ((select tb.tree_id, tb.taxon_id from djangophylocore_reltreecoltaxa%s as tb, djangophylocore_parentsrelation as par where (tb.taxon_id = par.taxon_id and par.parent_id =%s)  ) union (select tb.tree_id, tb.taxon_id from djangophylocore_reltreecoltaxa%s as tb where (tb.taxon_id = %s) )) as tu group by tu.tree_id ;"% ( baseId, parent_id, baseId, parent_id ) )
            if settings.DATABASE_ENGINE == 'sqlite3':
                result = cur.fetchall()
            else:
                result = cursor.fetchall()
            for (tree_id, nb_occurence) in result:
                if tree_id not in d_trees:
                    d_trees[tree_id] = {}
                d_trees[tree_id][pattern] =  nb_occurence
        l_trees_id = set()
        
        #VR il faut faire sur tout les arbres et non sur ceux present dans le rsultat (a cause des 0)
        # if a tree with none of the taxa can be a solution we have to tests alltrees not only those in d_trees
        # ex : {rodent}<4  ex2: {rodent}<4 and ({primate}>3 or {primate} <4)
        res0 = res;
        needAllTrees = False;
        for pattern in l_patterns:
            res0 = res0.replace( '{'+pattern+'}', '0' )
        if res0:
                try:
                    if eval( res0 ):
                       needAllTrees = True; 
                except:
                    raise SyntaxError, "bad query 0 %s" % query
        
        if needAllTrees: 
            baseId=self.id;
            if treebase:
                baseId=1;
                
            cur = cursor.execute( "select tree_id,tree_id from djangophylocore_reltreecoltaxa%s;" %baseId)
            if settings.DATABASE_ENGINE == 'sqlite3':
                result = cur.fetchall()
            else:
                result = cursor.fetchall()
            
            for (tree_id,tmp) in result:
                if tree_id not in d_trees:
                    d_trees[tree_id] = {}
            
        #VR
        for tree_id in d_trees:
            result = res
            for pattern in l_patterns:
                if pattern in d_trees[tree_id]:
                    result = result.replace( '{'+pattern+'}', str(d_trees[tree_id][pattern]) )
                else:
                    result = result.replace( '{'+pattern+'}', '0' )
            if result:
                try:
                    if eval( result ):
                       l_trees_id.add( tree_id ) 
                except:
                    raise SyntaxError, "bad query 2 %s %s" %( query, result)
      
        return Tree.objects.filter( id__in = l_trees_id )

    def query( self, query ):
        """
        query the collection in order to extract a targeted trees list.

        The query format is a python boolean expression with taxa name
        beetween braces :

        >>> col.query( "{muridae} > 2 and {primates}" )
        [<Tree: 1>, <Tree: 4>]

        will return a list of trees wich have more than 2 taxa wich have
        muridae as parents and at least 1 taxon wich have a primate as
        parents.

        if a taxa_list is not null, the query can have another variable
        {usertaxa}. this variable represente all taxa passed in the list.

        tree.eval_query( "{muridae} => 4 and {usertaxa} > 2", ['rattus', 'mus', 'pan', 'boss'] )

        will return true if tree have at least 4 taxa wich are muridae and
        more than 2 taxa wich are in the usertaxa_list 
        """
        return self._query( query, False )

    def query_treebase( self, query ):
        """
        query the treebase collection in order to extract a targeted trees list.

        The query format is a python boolean expression with taxa name
        beetween braces :

        >>> col.query( "{muridae} > 2 and {primates}" )
        [<Tree: 1>, <Tree: 4>]

        will return a list of treebase's trees wich have more than 2 taxa wich have
        muridae as parents and at least 1 taxon wich have a primate as
        parents.

        The special keyword {usertaxa} represent all the taxa present in the
        user collection.

        >>> tree.eval_query( "{muridae} => 4 and {usertaxa} > 2" )

        will return a treebase's trees list wich have at least 4 muridae and
        at least 2 taxa which are present in the user collection.
        """
        return self._query( query, True )

    def get_collection_from_query( self, query, treebase = False ):
        """
        return a new collection with all trees that match the query
        """
        format = self.format
        if treebase:
            format = "nexus"
        if format == 'phylip':
            return TreeCollection.objects.create( delimiter = self.delimiter, 
              source = '\n;'.join( [i.source for i in self._query( query, treebase )] ) )
        else:
            source = "#NEXUS\n\nBEGIN TREES;\n\n"
            trees_list = self._query( query, treebase )
            for tree in trees_list:
                source += "TREE %s = [&R] %s;\n" % ( tree.name, tree.source )
            source += "\nEND;\n"
            return TreeCollection.objects.create( delimiter = self.delimiter,
              source = source )
        
    def get_tree_size_distribution( self ):
        """ return stat of Tree Size Distribution """
        stat = {}
        cursor = connection.cursor()
        #cur = cursor.execute( " select tree_id, count(taxon_id) from djangophylocore_reltreecoltaxa%s GROUP BY tree_id;" % (self.id))
        #VR
        #cur = cursor.execute( " select tree_id, count(rel.user_taxon_name) from djangophylocore_reltreecoltaxa%s as rel, djangophylocore_taxonomy as taxonomy where taxonomy.id = rel.taxon_id GROUP BY tree_id;" % (self.id))
        cur = cursor.execute( " select tree_id, count(rel.user_taxon_name) from djangophylocore_reltreecoltaxa%s as rel GROUP BY tree_id;" % (self.id))
      
        if settings.DATABASE_ENGINE == 'sqlite3':
            result = cur.fetchall()
            cur.close()
        else:
            result = cursor.fetchall()
        cursor.close()
        for (tree_id, nbtaxa) in result:
            if not stat.has_key( nbtaxa ):
                stat[nbtaxa] = 0
            stat[nbtaxa] += 1
        if not stat.keys():
            raise ValueError, "your collection must have trees"
        nbmax = max( stat.keys() )
        ratio = int( nbmax*10.0/100) or 1 # 1 to prevent crash in xrange
        result_stat = {}
        for i in xrange( 0, nbmax+1, ratio):
            result_stat[i] = 0
        for i,j in stat.iteritems():
            for key in result_stat.keys():
                if key <= i < key+ratio:
                    result_stat[key] += j
        return result_stat

    def get_taxon_frequency_distribution( self ):
        """
        return stat of taxon frequency distribution
        """
        stat = {}
        cursor = connection.cursor()
        #VR cur = cursor.execute( "select rel.tree_id, taxonomy.name from djangophylocore_reltreecoltaxa%s as rel, djangophylocore_taxonomy as taxonomy where taxonomy.id = rel.taxon_id" % self.id )
        cur = cursor.execute( " select rel.tree_id, rel.user_taxon_name from djangophylocore_reltreecoltaxa%s as rel" % (self.id))
        if settings.DATABASE_ENGINE == 'sqlite3':
            results = cur.fetchall()
            cur.close()
        else:
            results = cursor.fetchall()
        cursor.close()
        d_results = {}
        for (tree, taxa_name) in results:
            if tree not in d_results:
                d_results[tree] = []
            d_results[tree].append( taxa_name )
        for tree in d_results:
            already_done = set()
            for taxon in d_results[tree]:
                if taxon not in already_done:
                    if not stat.has_key( taxon ):
                        stat[taxon] = 0
                    stat[taxon] += 1    
                    already_done.add( taxon )
        if not stat.values():
            raise ValueError, "your collection must have trees"
        nbmax = max( stat.values() )
        ratio = int( nbmax*10.0/100) or 1
        result_stat = {}
        for i in xrange( 0, nbmax+1, ratio):
            result_stat[i] = 0
        for i in stat.values():
            for key in result_stat.keys():
                if key <= i < key+ratio:
                    result_stat[key] += 1
        return result_stat

    def get_taxon_name_from_parents( self, parent_name ):
        """
        return a taxon names list wich have `parent_name` for parent
        """
        return [i.name for i in self.get_taxon_from_parents( parent_name )] 

    def get_taxon_from_parents( self, parent_name ):
        """
        return taxon in collection wich are for parent 'parent_name'
        """
        global get_taxonomy_toc, BADTAXA_TOC
        TAXONOMY_TOC = get_taxonomy_toc()
        assert parent_name in TAXONOMY_TOC, "%s does not exist in the current taxonomy" % parent_name
        parent_id = Taxonomy.objects.get( name = parent_name ).id
        #return Taxonomy.objects.extra( 
        #  tables = ["djangophylocore_reltreecoltaxa%s" % self.id, "djangophylocore_parentsrelation"],
        #  where = ["djangophylocore_taxonomy.id = djangophylocore_reltreecoltaxa%s.taxon_id  and djangophylocore_reltreecoltaxa%s.taxon_id = djangophylocore_parentsrelation.taxon_id and djangophylocore_parentsrelation.parent_id = %s" % (
        #    self.id, self.id, parent_id )]
        #).distinct()
        return Taxonomy.objects.extra( where = ["id IN (select taxon_id from djangophylocore_reltreecoltaxa%s where taxon_id IN (select taxon_id from djangophylocore_parentsrelation where parent_id = %s))" % (self.id, Taxonomy.objects.get( name = parent_name ).id)])

    def get_nb_trees( self, taxon ):
        """
        return the number of trees in collection wich contain taxon
        """
        return len(set([i[0] for i in self.rel.filter( taxon = taxon ).values_list( 'tree' )]))

    def get_filtered_collection_string( self, taxon_name_list ):
        """
        return a collections string wich have been striped of all taxa present
        in the taxon_name_list
        """
        if self.format == 'nexus':
            new_col = "#NEXUS\nBEGIN TREES;\n"
        else:
            new_col = ''
        parser = NewickParser()
        trees_list = self.trees.all()
        for tree in trees_list:
            if tree.is_valid:
                parser.parse_string( tree.source )
                filtered_tree = parser.filter( taxon_name_list )
                if filtered_tree:
                    filtered_tree = self._list2nwk( filtered_tree )
                    # Recreate nexus collection
                    if parser.get_taxa():
                        if self.format == 'nexus':
                            new_col += "Tree "+str(tree.name)+" = "+filtered_tree+";\n"
                        else:
                            new_col += filtered_tree+";\n"
        if self.format == 'nexus':
            new_col += "END;\n"
        return new_col

    def get_restricted_collection( self, taxon_name_list, keep = True ):
        """
        return a collection string wich contains only the taxon present in
        taxon_name_list
        """
        #print "XXXXXXXXXXXX"
        #for n in taxon_name_list:
        #        print n
        if keep:
            #remove_taxon_list = [i.user_taxon_name for i in self.rel.exclude( taxon__name__in = taxon_name_list )]
            #taxon_name_list = remove_taxon_list
            #VR
            remove_taxon_list = [i.user_taxon_name for i in self.rel]
            taxon_name_list2 =[i for i in remove_taxon_list if i not in taxon_name_list]
            taxon_name_list=taxon_name_list2
            #VR
            #print ">>>>>>>"
            #for n in taxon_name_list:
            #    print n
        #print "XXXXXXXXXXXX"
        new_nwk = self.get_filtered_collection_string( taxon_name_list )
        return TreeCollection.objects.create( delimiter = self.delimiter, source = new_nwk )

    def get_corrected_collection_string( self, correction ):
        """
        return a corrected source string where `correction` is a dictionary
        
            {'bad_name': 'good_name'}

        new_source_string = col.get_corrected_collection( {'echinops': 'echinops <plant>', 'ratis': 'rattus' } )
        """
        parser = NewickParser()
        trees_list = []
        for i in self.trees.all():
            if i.is_valid:
                parser.parse_string( i.source )
                corrected_tree = self._list2nwk( parser.correct_tree( correction ) )
#                corrected_tree = str( parser.correct_tree( correction ) )
#                corrected_tree = corrected_tree.replace( "[u'", "['").replace(", u'", ", '" )
#                corrected_tree = corrected_tree.replace( "['", "[" ).replace( "']", "]" )
#                corrected_tree = corrected_tree.replace( "',", "," ).replace( ", '", ", " )
#                corrected_tree = corrected_tree.replace("[", "(" ).replace("]", ")")
#                corrected_tree = corrected_tree.replace("('","(").replace("')",")")
#                corrected_tree = corrected_tree.replace("',", ",").replace(", '", ",").replace( ", ", ",")
            else:
                corrected_tree = i.source
            trees_list.append( (i.name, corrected_tree ) )
        if self.format == 'nexus':
            source = "#NEXUS\n\nBEGIN TREES;\n\n"
        else:
            source = ""
        for tree in trees_list:
            if self.format == 'nexus' :
                source += "TREE %s = %s;\n" % ( tree[0], tree[1] )
            else:
                source += "%s;\n" % tree[1]
        if self.format == 'nexus':
            source += "\nEND;\n"
        return source

    def get_corrected_collection( self, correction ):
        """
        return a corrected collection where `correction` is a dictionary
        
            {'bad_name': 'good_name'}

        newcol = col.get_corrected_collection( {'echinops': 'echinops <plant>', 'ratis': 'rattus' } )
        """
        new_nwk = self.get_corrected_collection_string( correction )
        return TreeCollection.objects.create( delimiter = self.delimiter, source = new_nwk )

    def get_autocorrected_collection( self ):
        """
        get a new collection wich have all single ambiguous taxa corrected
        """
        list_correction =  [(i.user_taxon_name, i.taxon.scientifics.get().name)\
          for i in self.rel.extra( where = ["taxon_id in (select id from djangophylocore_taxonomy where type_name != 'scientific name')"] )\
            if i.taxon.scientifics.count() == 1]
        #list_correction =  [(i.user_taxon_name, i.taxon.scientifics.get().name)\
        #  for i in self.rel.filter( taxon__in = self.taxa.exclude( type_name = 'scientific name' ) )\
        #    if i.taxon.scientifics.count() == 1] 
        return self.get_corrected_collection( list_correction ), list_correction

    def get_reference_arborescence( self ):
        """
        Take a taxon list, search in reference all parents names and
        return a networkx.DiGraph tree.
        """
        return self.get_reference_graph( self.taxa.filter( type_name = 'scientific name' ).iterator() )

    def get_statistics( self ):
        """
        return a dictionary with useful statistic informations.

        The dictionary structure take this form:
           
           {
              taxon_id: {
               'degree': degree_number,
               'scientific_taxon_list': set([]),
               'trees_list': set([]),
               'user_taxon_list': set([])
              },
            }

        i.e.:
        {
             349722: {'degree': 1,
                'scientific_taxon_list': set([u'carpomys']),
                'trees_list': set([6216, 6204, 6206]),
                'user_taxon_list': set([u'carpomys'])
             },
             359160: {'degree': 2,
                'scientific_taxon_list': set([u'oryza sativa']),
                'trees_list': set([6214]),
                'user_taxon_list': set([u'oryza sativa lhb'])
             },
            376913: {'degree': 2,
                'scientific_taxon_list': set([u'homo sapiens']),
                'trees_list': set([6214]),
                'user_taxon_list': set([u'homo sapiens cyg', u'homo sapiens mb',
                  u'homo sapiens ngb', u'homo sapiens hbt', u'homo sapiens hbb',
                  u'homo sapiens hba', u'homo sapiens hbg', u'homo sapiens hbd',
                  u'homo sapiens hbe'])
            },
        }
        """
        global TAXONOMY_TOC
	if not hasattr( self, "taxon_occurence" ):
            taxon_occurence = {}
            # initialisation of user taxon
            cursor = connection.cursor()
            cur = cursor.execute( "select rel.tree_id, rel.taxon_id, rel.user_taxon_name, taxonomy.name from djangophylocore_reltreecoltaxa%s as rel, djangophylocore_taxonomy as taxonomy where taxonomy.id = rel.taxon_id" % self.id )
            if settings.DATABASE_ENGINE == 'sqlite3':
                results = cur.fetchall()
                cur.close()
            else:
                results = cursor.fetchall()
            cursor.close()
            for ( tree_id, taxon_id, user_taxon_name, scientific_name ) in results:
                if taxon_id not in taxon_occurence:
                    taxon_occurence[taxon_id] = {"trees_list": set([]),
                      'user_taxon_list':set([]), "scientific_taxon_list":set([]), "degree":0}
                taxon_occurence[taxon_id]['trees_list'].add( tree_id )
                taxon_occurence[taxon_id]['user_taxon_list'].add( user_taxon_name )
                taxon_occurence[taxon_id]['scientific_taxon_list'].add( scientific_name )
            tree = self.get_reference_arborescence()
            if len( tree ):
                self.__compute_stats_arborescence( taxon_occurence, tree, Taxonomy.objects.get( name = 'root' ) )
            self.taxon_occurence = taxon_occurence
        return self.taxon_occurence

    def __compute_stats_arborescence( self, taxon_occurence, tree, node ):
        if node.id not in taxon_occurence:
            taxon_occurence[node.id] = {"trees_list": set([]),
              'user_taxon_list':set([]), "scientific_taxon_list":set([]), "degree":0}
        taxon_occurence[node.id]['degree'] = len( tree.successors( node ) ) + len( tree.predecessors( node ) )
        if tree.successors( node ): # if internal node
            for child in tree.successors( node ):
                self.__compute_stats_arborescence( taxon_occurence, tree, child )
                taxon_occurence[node.id]['trees_list'].update( taxon_occurence[child.id]['trees_list'] )
                taxon_occurence[node.id]['user_taxon_list'].update( taxon_occurence[child.id]['user_taxon_list'] )
                taxon_occurence[node.id]['scientific_taxon_list'].update( taxon_occurence[child.id]['scientific_taxon_list'] )

#    def _nxgraph2list( self, tree, node, d="" ):
#        if not d:
#            d = {}
#            for i in tree.edges():
#                if not d.has_key( i[0] ):
#                    d[i[0]] = []
#                d[i[0]].append( i[1] )
#            for i in tree.edges():
#                if not d.has_key( i[1] ):
#                    d[i[1]] = i[1]        
#        m = []
#        taxa_list = self.taxa.all()
#        for node in tree.successors( node ):
#            if isinstance( d[node], list ):
#                if node in taxa_list:
#                    m.append( node.name )
#                m.append( self._nxgraph2list( tree, node, d ) )
#            else:
#                m.append( d[node].name )
#        return m
#
    def _list2nwk( self, l ):
        NewickParser().remove_singleton( l )
        result = str( l )
        result = result.replace( "[u'", "['").replace(", u'", ", '" )
        result = result.replace( "['", "[" ).replace( "']", "]" )
        result = result.replace( "',", "," ).replace( ", '", ", " )
        result = result.replace("[", "(" ).replace("]", ")")
        result = result.replace("('","(").replace("')",")")
        result = result.replace("',", ",").replace(", '", ",").replace( ", ", ",")
        return result

#    def _nxgraph2nwk( self, tree, root ):
#        l = self._nxgraph2list( tree, root )
#        return self._list2nwk( l )[1:-1]+";" # removing extras parenthesis
#        
    def __graph2nwk_rec( self, tree, current_node, internal_label ):
        str_current =""
        sons_of_current_node = tree.successors(current_node)
        if len(sons_of_current_node) > 1:
            str_current += "("
            for s in sons_of_current_node:
                str_current += self.__graph2nwk_rec(tree, s, internal_label) +","
            str_current = str_current[:-1] # on a une "," en trop
            str_current += ")"
            if internal_label:
                str_current += "[" + current_node.name + "]"
            #if (currentNode.fiable == false)
            #    str_current +="?"
        elif len(sons_of_current_node) == 1:
            # intgernal node of degre 2
            return self.__graph2nwk_rec(tree, sons_of_current_node[0], internal_label)
        else: # it's a leave
            #if current_node.fiable == False:
            #    str_current +="?"
            str_current = current_node.name
        return str_current

    def get_reference_tree_as_nwk( self, internal_label = True ):
        """
        return the reference arborescence in a newick string
        """
        tree = self.get_reference_arborescence()
        if len(tree):
            return self.__graph2nwk_rec( tree, Taxonomy.objects.get( name = "root"), internal_label )+";"
        return ""

    def get_matrix( self ):
        """
        return the matrix wich described the presence of taxa in trees
        """
        matrix = {}
        trees_list = self.trees.all()
        stat = self.get_statistics()
        taxa_list = self.taxa.all()
        for taxa in taxa_list:
            matrix[taxa.id] = {}
            for tree in trees_list:
                tree_id = tree.id
                if tree_id in stat[taxa.id]['trees_list']:
                    matrix[taxa.id][tree_id] = 1
                else:
                    matrix[taxa.id][tree_id] = 0
        return matrix

    def get_matrix_csv( self ):
        """
        return the matrix in csv format
        """
        nb_trees = self.trees.count()
        nb_taxa = self.taxa.count()
        stat = self.get_statistics()
        matrix = self.get_matrix()
        # Sorting matrix
        sort_info = {'trees':{}, 'taxa':{}}
        for taxa, trees in matrix.iteritems():
            sort_info["taxa"][taxa] = len( [i for i in trees.values() if i] )
            for tree, presence in trees.iteritems():
                if presence:
                    if not tree in sort_info["trees"]:
                        sort_info["trees"][tree] = 0 
                    sort_info["trees"][tree] += 1 
        #
        taxa_list = [i[1] for i in sorted([(i,v) for v,i in sort_info['taxa'].items()])]
        tree_list = [i[1] for i in sorted([(i,v) for v,i in sort_info['trees'].items()])]
        d_tree_id_name = dict( Tree.objects.filter( id__in = tree_list).values_list( 'id','name' ) )
        csv = ["taxa,"+",".join([d_tree_id_name[tree] for tree in tree_list])]
        for taxa in taxa_list:
            tmp = matrix[taxa]
            line = [list(stat[taxa]['scientific_taxon_list'])[0]]
            for tree in tree_list:
                line.append( str(tmp[tree]) )
            csv.append( ",".join( line ) )
        return "\n".join( csv )




class AbstractTreeColTaxa( models.Model ):
    collection = models.ForeignKey( TreeCollection )#, related_name = 'rel' )
    tree = models.ForeignKey( Tree )#, related_name = 'rel' )
    taxon = models.ForeignKey( Taxonomy )#, related_name = 'rel', null = True )
    user_taxon_name = models.CharField( max_length = 200, null = True )
                
    class Meta:
        abstract = True
    def __unicode__( self ):
        return u"%s|%s|%s" % (self.collection.id, self.tree.id, self.user_taxon_name)

#############################################
#                Signals                    #
#############################################

def fill_created_updated_fields( sender, instance, signal, *args, **kwargs ):
    if not instance.id:
        instance.created = datetime.datetime.now()
    instance.updated = datetime.datetime.now()
signals.pre_save.connect(fill_created_updated_fields, sender=Tree)
signals.pre_save.connect(fill_created_updated_fields, sender=TreeCollection)

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

