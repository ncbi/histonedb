#
# Tree
#

def test_tree():
    return """
>>> from djangophylocore.models import *

# In this tree, there is only scientific names

>>> nwk_tree = " ( echinops <plantae>, (rattus, ( mus,(mus musculus))))"
>>> tree = Tree.objects.create( tree_string = nwk_tree, name = "bla")

>>> tree.is_valid
True
>>> tree.taxas.all()
[<Taxa: echinops <plantae>>, <Taxa: mus>, <Taxa: mus musculus>, <Taxa: rattus>]
>>> tree.bad_taxas.all()
[]
>>> tree.synonyms.all()
[]
>>> tree.homonyms.all()
[]
>>> tree.commons.all()
[]

# In this tree, there's bad taxas, common names, synonyms and homonyms

>>> nwk_tree = "( echinops, (rattus, ( mus, azerty, rat noir ), antilocapra anteflexa ))"
>>> tree = Tree.objects.create( tree_string = nwk_tree, name = "bla")

>>> tree.is_valid
True
>>> tree.taxas.all()
[<Taxa: mus>, <Taxa: rattus>]
>>> tree.bad_taxas.all()
[<BadTaxa: azerty (1)>]
>>> tree.synonyms.all()
[<SynonymName: antilocapra anteflexa>]
>>> tree.homonyms.all()
[<HomonymName: echinops>]
>>> tree.commons.all()
[<CommonName: rat noir (french)>]

# Getting all taxa list
>>> tree.taxonomy_objects.all()
[<Taxonomy: antilocapra anteflexa (synonym)>, <Taxonomy: echinops (homonym)>, <Taxonomy: mus (scientific name)>, <Taxonomy: rat noir (common)>, <Taxonomy: rattus (scientific name)>]

# Getting ambiguous taxa (synonms, commons, homonyms...)
>>> tree.ambiguous.all()
[<Taxonomy: antilocapra anteflexa (synonym)>, <Taxonomy: echinops (homonym)>, <Taxonomy: rat noir (common)>]


# getting arborescence

>>> tree.arborescence.edges()
[(<Taxa: root>, <Taxa: muridae>), (<Taxa: root>, <HomonymName: echinops>), (<Taxa: muridae>, <Taxa: mus>), (<Taxa: muridae>, <BadTaxa: azerty (1)>), (<Taxa: muridae>, <SynonymName: antilocapra anteflexa>), (<Taxa: muridae>, <CommonName: rat noir (french)>), (<Taxa: muridae>, <Taxa: rattus>)]

# we can know the number of occurence of each bad taxa

>>> nwk_tree = "(rattus, ( azerty, rat ), badname)"
>>> tree = Tree.objects.create( tree_string = nwk_tree, name = "bla")

>>> tree.bad_taxas.all()
[<BadTaxa: azerty (2)>, <BadTaxa: badname (1)>, <BadTaxa: rat (1)>]

# and getting arborescence even with bad taxa

>>> tree.arborescence.edges()
[(<Taxa: root>, <BadTaxa: azerty (2)>), (<Taxa: root>, <BadTaxa: rat (1)>), (<Taxa: root>, <BadTaxa: badname (1)>), (<Taxa: root>, <Taxa: rattus>)]

Going further
-------------

Passing bad newick string
~~~~~~~~~~~~~~~~~~~~~~~~~

If a bad newick string is passed, the creation of the tree won't crash.
This is useful to get statistics of the tree anyway.

>>> bad_nwk = "(mus,(,("
>>> bad_tree = Tree.objects.create( name = 'badtree', tree_string = bad_nwk )
>>> bad_tree.taxonomy_objects.all()
[<Taxonomy: mus (scientific name)>]

To know if a tree is valid or not, check the `is_valid` attribute

>>> bad_tree.is_valid
False

Working with delimiter
~~~~~~~~~~~~~~~~~~~~~~~

>>> nwk_tree = "(rattus_rattus, (mus, mus_musculus ))"
>>> tree2 = Tree.objects.create( tree_string = nwk_tree, name = "blu", delimiter = '_')
>>> tree2.taxonomy_objects.all()
[<Taxonomy: mus (scientific name)>, <Taxonomy: mus musculus (scientific name)>, <Taxonomy: rattus rattus (scientific name)>]

# delimiter can't contain  '(', ')' or ','

>>> nwk_tree = "(murinae, (mus, mus(,)musculus ))"
>>> try:
...     tree2 = Tree.objects.create( tree_string = nwk_tree, name = "blo", delimiter = '(,)')
... except ValueError, e:
...     e
ValueError('"(,)" is a bad delimiter',)

Making queries
~~~~~~~~~~~~~~
>>> nwk_tree = " ( echinops <plantae>, (rattus, ( mus,(mus musculus))))"
>>> tree = Tree.objects.create( tree_string = nwk_tree, name = "filtered")
>>> tree.get_nb_taxa_from_parent('muridae')
3L
>>> tree.get_nb_taxa_from_parent('plantae')
1L
>>> tree.eval_query('{plantae}==1')
True
>>> tree.eval_query('{plantae}>1')
False
>>> tree.eval_query('not {plantae}>1')
True
>>> tree.eval_query('{plantae}>1 and {murinae}>2')
Traceback (most recent call last):
 File "<console>", line 1, in <module>
    ...
    raise NameError, striped_pattern
NameError: murinae
>>> tree.eval_query('{plantae}==1 and {muridae}>2')
True
>>> tree.eval_query('{plantae}==1 and not {muridae}<2')
True
>>> tree.eval_query( '{plantae}==1 and {usertaxa} > 1', usertaxa_list=['rattus', 'mus'] )
True
>>> tree.eval_query( '{usertaxa} > 1', usertaxa_list=['mus'] )
False
>>> tree.eval_query( '{usertaxa} > 1', usertaxa_list=[] )
False

###########################
# TESTING PRIVATE METHODS #
###########################

>>> try:
...     tree._Tree__get_django_objects_from_nwk( '(rattus,bla)')
... except ValueError, e:
...     e
ValueError(u'bla not found in the database',)

>>> taxa_list = tree._Tree__get_django_objects_from_nwk( '(rattus,(mus,(echinops,badname)))')
>>> taxa_list
[<Taxa: rattus>, <Taxa: mus>, <HomonymName: echinops>, <BadTaxa: badname (1)>]
>>> tree._Tree__get_scientific_taxa( taxa_list )
[<Taxa: rattus>, <Taxa: mus>]



"""


