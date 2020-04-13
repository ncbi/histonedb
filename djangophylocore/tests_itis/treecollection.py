#
# TreeCollection
#

def test_tree_collection():
    return """
>>> from djangophylocore.models import *

Creation of a simple collection
-------------------------------

>>> simple_col = "(mus,antilocapra anteflexa,rat noir,echinops,blabla);(mus, rat noir);"
>>> col = TreeCollection.objects.create(name = 'simple', original_collection_string = simple_col)
>>> col.taxonomy_objects.all()
[<Taxonomy: antilocapra anteflexa (synonym)>, <Taxonomy: echinops (homonym)>, <Taxonomy: mus (scientific name)>, <Taxonomy: rat noir (common)>]
>>> col.ambiguous.all()
[<Taxonomy: antilocapra anteflexa (synonym)>, <Taxonomy: echinops (homonym)>, <Taxonomy: rat noir (common)>]
>>> col.taxas.all()
[<Taxa: mus>]
>>> col.homonyms.all()
[<HomonymName: echinops>]
>>> col.synonyms.all()
[<SynonymName: antilocapra anteflexa>]
>>> col.commons.all()
[<CommonName: rat noir (french)>]
>>> col.bad_taxas.all()
[<BadTaxa: blabla (1)>]

>>> col.trees.all()
[<Tree: 1>, <Tree: 2>]
>>> col.trees.count()
2L

Creation of collection from nexus (translated) format
-----------------------------------------------------

# if there is problem with encoding, try using icon :
# iconv -f iso-8859-15 -t utf-8 nexus_example > nexus_example_utf-8
# or the following lines
#>>> import codecs
#>>> nwk_col = codecs.open("./djangophylocore/tests/data/nexus_example", encoding='iso-8859-15').read()

>>> import djangophylocore, os
>>> nwk_col = open(os.path.join(djangophylocore.__path__[0],"tests","data","nexus_example_utf-8")).read()

After creating collection, you have to call
`generate_from_original_collection_string` juste after. This method will build
all trees from the collection string

>>> tree_col = TreeCollection.objects.create(name="test_nexus", original_collection_string=nwk_col)

The string format is detected

>>> tree_col.format
'nexus'

Dealing with collection
~~~~~~~~~~~~~~~~~~~~~~~

# Tree names are taken from nexus file

>>> tree_col.trees.all()
[<Tree: Fig._3>]

>>> tree_col.bad_taxas
[<BadTaxa: phytophthora_cajani (1)>, <BadTaxa: phytophthora_cinnamomi (1)>, <BadTaxa: phytophthora_drechsleri-like (1)>, <BadTaxa: phytophthora_melonis (1)>, <BadTaxa: phytophthora_pistaciae_a (1)>, <BadTaxa: phytophthora_pistaciae_b (1)>, <BadTaxa: phytophthora_sinensis (1)>, <BadTaxa: phytophthora_sojae_a (1)>, <BadTaxa: phytophthora_sojae_b (1)>, <BadTaxa: phytophthora_vignae (1)>]

>>> tree_col.bad_taxas.count()
10L

>>> tree_col.taxonomy_objects.all()
[<BadTaxa: phytophthora_melonis (1)>, <BadTaxa: phytophthora_drechsleri-like (1)>, <BadTaxa: phytophthora_sojae_a (1)>, <BadTaxa: phytophthora_sojae_b (1)>, <BadTaxa: phytophthora_sinensis (1)>, <BadTaxa: phytophthora_cinnamomi (1)>, <BadTaxa: phytophthora_pistaciae_a (1)>, <BadTaxa: phytophthora_pistaciae_b (1)>, <BadTaxa: phytophthora_vignae (1)>, <BadTaxa: phytophthora_cajani (1)>]


Export collection to string
~~~~~~~~~~~~~~~~~~~~~~~~~~~

>>> tree_col.get_collection_string()
u'#NEXUS\\n\\nBEGIN TREES;\\n\\nTREE Fig._3 = [&R] (((Phytophthora_sojae_a,Phytophthora_sojae_b),(((Phytophthora_vignae,Phytophthora_cajani),(Phytophthora_melonis,(Phytophthora_drechsleri-like,Phytophthora_sinensis))),(Phytophthora_pistaciae_a,Phytophthora_pistaciae_b))),Phytophthora_cinnamomi)\\n\\nEND;\\n'


Adding tree to the collection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

>>> tree = Tree.objects.create(name = 'newtree', tree_string = "(mus musculus, rattus rattus)")
>>> tree_col.trees.add(tree)
>>> tree_col.bad_taxas
[<BadTaxa: phytophthora_cajani (1)>, <BadTaxa: phytophthora_cinnamomi (1)>, <BadTaxa: phytophthora_drechsleri-like (1)>, <BadTaxa: phytophthora_melonis (1)>, <BadTaxa: phytophthora_pistaciae_a (1)>, <BadTaxa: phytophthora_pistaciae_b (1)>, <BadTaxa: phytophthora_sinensis (1)>, <BadTaxa: phytophthora_sojae_a (1)>, <BadTaxa: phytophthora_sojae_b (1)>, <BadTaxa: phytophthora_vignae (1)>]

>>> tree_col.taxonomy_objects.all()
[<Taxonomy: mus musculus (scientific name)>, <Taxonomy: rattus rattus (scientific name)>]

>>> tree_col.get_collection_string()
u'#NEXUS\\n\\nBEGIN TREES;\\n\\nTREE Fig._3 = [&R] (((Phytophthora_sojae_a,Phytophthora_sojae_b),(((Phytophthora_vignae,Phytophthora_cajani),(Phytophthora_melonis,(Phytophthora_drechsleri-like,Phytophthora_sinensis))),(Phytophthora_pistaciae_a,Phytophthora_pistaciae_b))),Phytophthora_cinnamomi);\\nTREE newtree =  (mus musculus, rattus rattus)\\n\\nEND;\\n'

Creation of collection from phylip format
-----------------------------------------

>>> nwk_col = "(echinops <plantae>,(rattus,(mus,(mus musculus)));\\n(rattus,(azerty,rat));\\n(echinops, (rattus, (mus, azerty, rat noir), antilocapra anteflexa));\\n(rattus, echinops, mus);"
>>> tree_col = TreeCollection.objects.create(name="test_phylip", original_collection_string=nwk_col)
>>> tree_col.format
'phylip'

# Tree names are generated

>>> tree_col.trees.all()
[<Tree: 1>, <Tree: 2>, <Tree: 3>, <Tree: 4>]

Dealing with collections
~~~~~~~~~~~~~~~~~~~~~~~~

>>> tree_col.taxas
[<Taxa: echinops <plantae>>, <Taxa: mus>, <Taxa: mus musculus>, <Taxa: rattus>]
>>> tree_col.taxas.count()
4L
>>> tree_col.bad_taxas
[<BadTaxa: azerty (4)>, <BadTaxa: rat (2)>]
>>> tree_col.homonyms
[<HomonymName: echinops>]
>>> tree_col.commons
[<CommonName: rat noir (french)>]
>>> tree_col.synonyms
[<SynonymName: antilocapra anteflexa>]
>>> tree_col.bad_taxas.count()
2L
>>> tree_col.get_collection_string()
u'(echinops <plantae>,(rattus,(mus,(mus musculus)));\\n(rattus,(azerty,rat));\\n(echinops,(rattus,(mus,azerty,rat noir),antilocapra anteflexa));\\n(rattus,echinops,mus);'

Adding trees in collection
~~~~~~~~~~~~~~~~~~~~~~~~~~

>>> tree
<Tree: newtree>
>>> tree_col.trees.add(tree)
>>> tree_col.trees.all()
[<Tree: newtree>, <Tree: 1>, <Tree: 2>, <Tree: 3>, <Tree: 4>]
>>> tree_col.get_collection_string()
u'(mus musculus, rattus rattus);\\n(echinops <plantae>,(rattus,(mus,(mus musculus)));\\n(rattus,(azerty,rat));\\n(echinops,(rattus,(mus,azerty,rat noir),antilocapra anteflexa));\\n(rattus,echinops,mus);'


Specify the delimiter
---------------------
>>> tree_col = TreeCollection.objects.create(name="test_delimiter",
...   original_collection_string="(rattus_rattus,(mus, mus_musculus))", delimiter = '_')
>>> tree_col.bad_taxas.all()
[]
>>> tree_col.taxonomy_objects.all()
[<Taxonomy: mus (scientific name)>, <Taxonomy: mus musculus (scientific name)>, <Taxonomy: rattus rattus (scientific name)>]


Getting bad trees from a collection
-----------------------------------

>>> bad_tree_nwk_1 = '(mus,(murinae,rattus)'
>>> bad_tree_1 = Tree.objects.create(name = 'bad_tree_1', tree_string = bad_tree_nwk_1)
>>> bad_tree_1.is_valid
False

>>> bad_tree_nwk_2 = '(mus,(murinae,rattus)0)'
>>> bad_tree_2 = Tree.objects.create(name = 'bad_tree_2', tree_string = bad_tree_nwk_2)
>>> bad_tree_2.is_valid
False

>>> bad_tree_nwk_3 = '(mus,('
>>> bad_tree_3 = Tree.objects.create(name = 'bad_tree_3', tree_string = bad_tree_nwk_3)
>>> bad_tree_3.is_valid
False

>>> good_tree_nwk_1 = '(mus,(rattus, mus musculus)) '
>>> good_tree_1 = Tree.objects.create(name = 'good_tree_1', tree_string = good_tree_nwk_1)
>>> good_tree_1.is_valid
True

>>> col = TreeCollection(name = 'with_bad_trees')
>>> col.save()
>>> col.trees.add(good_tree_1)
>>> col.trees.add(bad_tree_1)
>>> col.trees.add(bad_tree_2)
>>> col.trees.add(bad_tree_3)

>>> col.bad_trees.count()
3
>>> col.bad_trees.all()
[<Tree: bad_tree_1>, <Tree: bad_tree_2>, <Tree: bad_tree_3>]

Making queries
--------------

>>> simple_col = "(mus,echinops <plantae>,rattus);(mus, rattus);"
>>> col = TreeCollection.objects.create(name = 'query_col', original_collection_string = simple_col)
>>> col.query('{muridae} > 1')
[<Tree: 1>, <Tree: 2>]
>>> col.query('{muridae}< 1 or {plantae}>0')
[<Tree: 1>]
>>> col.query('{muridae}> 1 and not {plantae}')
[<Tree: 2>]
>>> col.query('{usertaxa} == 1', usertaxa_list = ['echinops <plantae>'])
[<Tree: 1>]
>>> col.query('{usertaxa} == 1', usertaxa_list = ['mus', 'rattus'])
[]
>>> col.query('{usertaxa} > 1', usertaxa_list = ['mus', 'rattus'])
[<Tree: 1>, <Tree: 2>]
>>> col.trees.all()
[<Tree: 1>, <Tree: 2>]
>>> col.id
6L
>>> filtered_col = col.get_collection_from_query('{plantae}')
>>> filtered_col.trees.all()
[<Tree: 1>]
>>> filtered_col.id
7L
>>> col.trees.all()
[<Tree: 1>, <Tree: 2>]


"""
