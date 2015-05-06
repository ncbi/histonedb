#
# TreeCollection
#

def test_tree_collection():
    return """
>>> from djangophylocore.models import *

Creation of a simple collection
-------------------------------

>>> TAXONOMY_TOC = get_taxonomy_toc( True )
>>> simple_col = "(mus,nannomys,black rat,echinops,blabla);(mus, black rat);"
>>> col = TreeCollection.objects.create( name = 'simple', source = simple_col )
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
#>>> nwk_col = codecs.open( "./djangophylocore/tests/data/nexus_example", encoding='iso-8859-15' ).read()

>>> import djangophylocore, os
>>> nwk_col = open( os.path.join( djangophylocore.__path__[0],"tests","data","nexus_example_utf-8" )).read()

After creating collection, you have to call
`generate_from_source` juste after. This method will build
all trees from the collection string

>>> tree_col = TreeCollection.objects.create( name="test_nexus", source=nwk_col, delimiter = '_' )

The string format is detected

>>> tree_col.format
'nexus'

Dealing with collection
~~~~~~~~~~~~~~~~~~~~~~~

# Tree names are taken from nexus file

>>> tree_col.trees.all()
[<Tree: fig._3>]

>>> tree_col.bad_taxa.all()
[<BadTaxa: phytophthora_cajani (0)>, <BadTaxa: phytophthora_cinnamomi (0)>, <BadTaxa: phytophthora_drechsleri-like (0)>, <BadTaxa: phytophthora_melonis (0)>, <BadTaxa: phytophthora_pistaciae_a (0)>, <BadTaxa: phytophthora_pistaciae_b (0)>, <BadTaxa: phytophthora_sinensis (0)>, <BadTaxa: phytophthora_sojae_a (0)>, <BadTaxa: phytophthora_sojae_b (0)>, <BadTaxa: phytophthora_vignae (0)>]

>>> tree_col.bad_taxa.count()
10L

>>> tree_col.taxa.all()
[]


Export collection to string
~~~~~~~~~~~~~~~~~~~~~~~~~~~

>>> tree_col.get_collection_string()
u'#NEXUS\\n\\nBEGIN TREES;\\n\\nTREE fig._3 =  (((phytophthora_sojae_a,phytophthora_sojae_b),(((phytophthora_vignae,phytophthora_cajani),(phytophthora_melonis,(phytophthora_drechsleri-like,phytophthora_sinensis))),(phytophthora_pistaciae_a,phytophthora_pistaciae_b))),phytophthora_cinnamomi)\\n\\nEND;\\n'


Creation of collection from phylip format
-----------------------------------------

>>> nwk_col = "(echinops <plant>,(rattus,( mus,(mus musculus))));\\n(rattus,(azerty,ratis));\\n(echinops, (rattus, ( mus, azerty, black rat ), nannomys ));\\n(rattus, echinops, mus);"
>>> tree_col = TreeCollection.objects.create( name="test_phylip", source=nwk_col )
>>> tree_col.format
'phylip'

# Tree names are generated

>>> tree_col.trees.all()
[<Tree: 1>, <Tree: 2>, <Tree: 3>, <Tree: 4>]

Dealing with collections
~~~~~~~~~~~~~~~~~~~~~~~~

>>> tree_col.scientifics.all()
[<Taxonomy: echinops <plant> (scientific name)>, <Taxonomy: mus (scientific name)>, <Taxonomy: mus musculus (scientific name)>, <Taxonomy: rattus (scientific name)>]
>>> tree_col.scientifics.count()
4L
>>> tree_col.bad_taxa
[<BadTaxa: azerty (0)>, <BadTaxa: ratis (0)>]
>>> tree_col.homonyms
[<Taxonomy: echinops (homonym)>]
>>> tree_col.commons
[<Taxonomy: black rat (common)>]
>>> tree_col.synonyms
[<Taxonomy: nannomys (synonym)>]
>>> tree_col.bad_taxa.count()
2L
>>> tree_col.get_collection_string()
u'(echinops <plant>,(rattus,( mus,(mus musculus))));\\n(rattus,(azerty,ratis));\\n(echinops, (rattus, ( mus, azerty, black rat ), nannomys ));\\n(rattus, echinops, mus);'

Specify the delimiter
---------------------
>>> tree_col = TreeCollection.objects.create( name="test_delimiter",
...   source="(rattus_rattus,(mus, mus_musculus))", delimiter = '_' )
>>> tree_col.bad_taxa.all()
[]
>>> tree_col.taxa.all()
[<Taxonomy: mus (scientific name)>, <Taxonomy: mus musculus (scientific name)>, <Taxonomy: rattus rattus (scientific name)>]


Getting bad trees from a collection
-----------------------------------

>>> good_tree_nwk_1 = '(mus,(rattus, mus musculus )) '
>>> bad_tree_nwk_1 = '(mus,( murinae,rattus)'
>>> good_tree_nwk_2 = '(mus,( murinae,rattus)0)'
>>> bad_tree_nwk_2 = '(mus,( '
>>> nwk_col = ';'.join( [good_tree_nwk_1, bad_tree_nwk_1, good_tree_nwk_2, bad_tree_nwk_2] )
>>> col = TreeCollection( name = 'with_bad_trees', source = nwk_col)
>>> col.save()

>>> col.bad_trees.count()
2
>>> col.bad_trees.all()
[<Tree: 2>, <Tree: 4>]

Making queries
--------------

>>> simple_col = "(mus,echinops <plant>,rattus);(mus, rattus);"
>>> col = TreeCollection.objects.create( name = 'query_col', source = simple_col )
>>> col.query( '{murinae} > 1' )
[<Tree: 1>, <Tree: 2>]
>>> col.query( '{murinae}< 1 or {cardueae}>0' )
[<Tree: 1>]
>>> col.query( '{murinae}> 1 and not {cardueae}' )
[<Tree: 2>]
>>> col.query_treebase( '{usertaxa} == 1' )
[<Tree: 1>]
>>> col.query_treebase( '{usertaxa} == 1' )
[]
>>> col.query_treebase( '{usertaxa} > 1' )
[<Tree: 1>, <Tree: 2>]
>>> col.trees.all()
[<Tree: 1>, <Tree: 2>]
>>> col.id
6L
>>> filtered_col = col.get_collection_from_query( '{cardueae}' )
>>> filtered_col.trees.all()
[<Tree: 1>]
>>> filtered_col.id
7L
>>> col.trees.all() # 2
[<Tree: 1>, <Tree: 2>]

Filter and restrict collection
------------------------------

>>> simple_col = "(mus musculus, rattus, glis);( glis, mus, rattus rattus);(mus, rattus);"
>>> col = TreeCollection.objects.create( source = simple_col )
>>> col.get_filtered_collection_string( ['mus', 'glis'] )
'(mus musculus,rattus);\\n(rattus rattus);\\n(rattus);\\n'
>>> col.get_filtered_collection_string( ['mus', 'rattus'] )
'(mus musculus,glis);\\n(glis,rattus rattus);\\n'
>>> '' is col.get_filtered_collection_string( ['mus musculus', 'rattus rattus', 'glis', 'mus', 'rattus'] )
True

>>> restricted_col = col.get_restricted_collection( ['mus', 'glis'] )
>>> restricted_col.taxa.all()
[<Taxonomy: mus (scientific name)>]
>>> restricted_col.bad_taxa.all()
[<BadTaxa: glis (0)>]
>>> restricted_col = col.get_restricted_collection( ['mus musculus', 'rattus'] )
>>> restricted_col.taxa.all()
[<Taxonomy: mus musculus (scientific name)>, <Taxonomy: rattus (scientific name)>]


Correct the collection
----------------------

>>> simple_col = "(mus musculus, (rattus, echinops));((echinops, mis france), rattis rattus);(mus, (rattus));"
>>> col = TreeCollection.objects.create( source = simple_col )
>>> col.taxa.all()
[<Taxonomy: echinops (homonym)>, <Taxonomy: mus (scientific name)>, <Taxonomy: mus musculus (scientific name)>, <Taxonomy: rattus (scientific name)>]
>>> col.bad_taxa.all()
[<BadTaxa: mis (0)>, <BadTaxa: rattis (0)>]
>>> col.homonyms.all()
[<Taxonomy: echinops (homonym)>]
>>> col.get_corrected_collection_string( {'mis': 'mus', 'rattis rattus':'rattus rattus', 'echinops': 'echinops <plant>'} )
'(mus musculus,(rattus,echinops <plant>));\\n((echinops <plant>,mus france),rattus rattus);\\n(mus,rattus);\\n'

>>> corrected_col = col.get_corrected_collection({'mis': 'mus', 'rattis rattus':'rattus rattus', 'echinops': 'echinops <plant>'})
>>> corrected_col.bad_taxa.all()
[]
>>> corrected_col.homonyms.all()
[]
>>> corrected_col.taxa.all()
[<Taxonomy: echinops <plant> (scientific name)>, <Taxonomy: mus (scientific name)>, <Taxonomy: mus musculus (scientific name)>, <Taxonomy: rattus (scientific name)>, <Taxonomy: rattus rattus (scientific name)>]

"""
