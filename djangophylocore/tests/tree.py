#
# Tree
#

def test_tree():
    return """
>>> TEST_LAUNCHED = True
>>> from djangophylocore.models import *

>>> TAXONOMY_TOC = get_taxonomy_toc(True)

# In this tree, there is only scientific names

>>> nwk_tree = " (echinops <plant>, (rattus, (mus,(mus musculus))))"
>>> tree = Tree.objects.create(source = nwk_tree, name = "1")

>>> tree.is_valid
True
>>> tree.taxa.all()
[<Taxonomy: echinops <plant> (scientific name)>, <Taxonomy: mus (scientific name)>, <Taxonomy: mus musculus (scientific name)>, <Taxonomy: rattus (scientific name)>]
>>> tree.bad_taxa.all()
[]
>>> tree.synonyms.all()
[]
>>> tree.homonyms.all()
[]
>>> tree.commons.all()
[]

# In this tree, there's bad taxa, common names, synonyms and homonyms

>>> nwk_tree = "(echinops, (rattus, (mus, azerty, black rat), nannomys))"
>>> tree = Tree.objects.create(source = nwk_tree, name = "2")

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
[<Taxonomy: black rat (common)>, <Taxonomy: echinops (homonym)>, <Taxonomy: mus (scientific name)>, <Taxonomy: nannomys (synonym)>, <Taxonomy: rattus (scientific name)>]

# Getting ambiguous taxa (synonms, commons, homonyms...)
>>> tree.ambiguous.all()
[<Taxonomy: black rat (common)>, <Taxonomy: echinops (homonym)>, <Taxonomy: nannomys (synonym)>]

# getting arborescence
>>> tree.arborescence.edges()
[(<Taxa: root>, <HomonymName: echinops>), (<Taxa: root>, <Taxa: murinae>), (<Taxa: murinae>, <Taxa: mus>), (<Taxa: murinae>, <BadTaxa: azerty (0)>), (<Taxa: murinae>, <Taxa: rattus>), (<Taxa: murinae>, <SynonymName: nannomys>), (<Taxa: murinae>, <CommonName: black rat (english)>)]

# we can know the number of occurence of each bad taxa

>>> nwk_tree = "(rattus, (azerty, black rat), badname)"
>>> tree = Tree.objects.create(source = nwk_tree, name = "3")

>>> tree.bad_taxa.all()
[<BadTaxa: azerty (0)>, <BadTaxa: badname (0)>]

# and getting arborescence even with bad taxa

>>> tree.arborescence.edges()
[(<Taxa: root>, <CommonName: black rat (english)>), (<Taxa: root>, <BadTaxa: azerty (0)>), (<Taxa: root>, <Taxa: rattus>), (<Taxa: root>, <BadTaxa: badname (0)>)]

Going further
-------------

Passing bad newick string
~~~~~~~~~~~~~~~~~~~~~~~~~

If a bad newick string is passed, the creation of the tree won't crash.
This is useful to get statistics of the tree anyway.

>>> bad_nwk = "(mus,(,("
>>> bad_tree = Tree.objects.create(name = 'badtree', source = bad_nwk)
>>> bad_tree.taxa.all()
[<Taxonomy: mus (scientific name)>]

To know if a tree is valid or not, check the `is_valid` attribute

>>> bad_tree.is_valid
False

Working with delimiter
~~~~~~~~~~~~~~~~~~~~~~~

>>> nwk_tree = "(rattus_rattus, (mus, mus_musculus))"
>>> tree2 = Tree.objects.create(source = nwk_tree, name = "4", delimiter = '_')
>>> tree2.taxa.all()
[<Taxonomy: mus (scientific name)>, <Taxonomy: mus musculus (scientific name)>, <Taxonomy: rattus rattus (scientific name)>]

# delimiter can't contain  '(', ')' or ','

>>> nwk_tree = "(murinae, (mus, mus(,)musculus))"
>>> try:
...     tree2 = Tree.objects.create(source = nwk_tree, name = "5", delimiter = '(,)')
... except ValueError, e:
...     e
ValueError('"(,)" is a bad delimiter',)

Making queries
~~~~~~~~~~~~~~
>>> nwk_tree = " (echinops <plant>, (rattus, (mus,(mus musculus))))"
>>> tree = Tree.objects.create(source = nwk_tree, name = "filtered")
>>> tree.get_nb_taxa_from_parent('murinae')
3L
>>> tree.get_nb_taxa_from_parent('cardueae')
1L
>>> tree.eval_query('{cardueae}==1')
True
>>> tree.eval_query('{cardueae}>1')
False
>>> tree.eval_query('not {cardueae}>1')
True
>>> tree.eval_query('{cardueae}>1 and {foobar}>2')
Traceback (most recent call last):
 File "<console>", line 1, in <module>
    ...
    raise NameError(striped_pattern)
NameError: foobar
>>> tree.eval_query('{cardueae}==1 and {murinae}>2')
True
>>> tree.eval_query('{cardueae}==1 and not {murinae}<2')
True

"""


