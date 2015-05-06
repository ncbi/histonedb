#
# TaxonomyReference
#

def test_taxonomy_reference():
    return """

All objects that inerits from TaxonomyReference will have some more methods.
Those methods provide informations from the taxonomy.

>>> from djangophylocore.models import *
>>> taxoref = TaxonomyReference()


Working with names
------------------

>>> taxoref.is_valid_name( 'rattus' )
True
>>> taxoref.is_valid_name( 'rattus other' )
False

>>> taxoref.is_scientific_name( 'rattus' )
True
>>> taxoref.is_scientific_name( 'antilocapra anteflexa' )
False

>>> taxoref.is_common( 'rat noir' )
True
>>> taxoref.is_common( 'rat' )
False

>>> taxoref.is_synonym( 'antilocapra anteflexa' )
True
>>> taxoref.is_synonym( 'rattus' )
False

>>> taxoref.is_homonym( 'echinops' )
True
>>> taxoref.is_homonym( 'echinops <plantae>' )
False

>>> taxoref.is_bad_taxa( 'rattus' )
False
>>> taxoref.is_bad_taxa( 'abadname' )
True
>>> abadname = BadTaxa.objects.create( name = 'abadname' )
>>> taxoref.is_bad_taxa( 'abadname' )
True

>>> taxoref.get_name_from_common( 'souris commune' )
[<Taxa: mus musculus>]
>>> taxoref.get_name_from_common( 'mus' )
[]

>>> taxoref.get_name_from_synonym( 'antilocapra anteflexa' )
[<Taxa: antilocapra americana>]
>>> taxoref.get_name_from_synonym( 'mus' )
[]

>>> taxoref.get_name_from_homonym( 'echinops' )
[<Taxa: echinops <animalia>>, <Taxa: echinops <plantae>>]
>>> taxoref.get_name_from_homonym( 'homo' )
[]

>>> taxoref.strip_taxa_name( "rattus" )
'rattus'
>>> taxoref.strip_taxa_name( "rattus_france" )
'rattus'
>>> taxoref.strip_taxa_name( "rattus france, delimiter=' '" )
'rattus'
>>> taxoref.strip_taxa_name( "mus_musculus" )
'mus musculus'
>>> taxoref.strip_taxa_name( "mus_musculus_france" )
'mus musculus'

>>> taxoref.strip_taxa_name( "rattus rattus france", delimiter=' ')
'rattus rattus'

'mus something foo' is not in our test database so we're falling down to 'mus'
>>> taxoref.strip_taxa_name( "mus something foo", delimiter=' ')
'mus'

>>> taxoref.correct( 'rattus' ) is None
True
>>> taxoref.correct( 'house mouse' )
[<Taxa: mus musculus>]
>>> taxoref.correct( 'echinops' )
[<Taxa: echinops <animalia>>, <Taxa: echinops <plantae>>]
>>> taxoref.correct( 'antilocapra anteflexa' )
[<Taxa: antilocapra americana>]
>>> taxoref.correct( 'taxa not in database' )
[0]

Getting django objects from name
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
>>> taxoref.get_object_from_name( 'mus' )
<Taxa: mus>
>>> taxoref.get_object_from_name( 'antilocapra anteflexa' )
<SynonymName: antilocapra anteflexa>
>>> taxoref.get_object_from_name( 'echinops' )
<HomonymName: echinops>
>>> taxoref.get_object_from_name( 'badname' )
<BadTaxa: badname (1)>

The name must be present into the database

>>> try:
...     taxoref.get_object_from_name( 'very badname' )
... except ValueError, e:
...     e
ValueError('very badname not found in the database',)

Working with taxa
-----------------

>>> root = Taxa.objects.get( name = 'root' )
>>> muridae = Taxa.objects.get( name = 'muridae' )
>>> mus = Taxa.objects.get( name = 'mus' )
>>> mus_musculus = Taxa.objects.get( name = 'mus musculus' )
>>> rattus = Taxa.objects.get( name = 'rattus' )
>>> echinops_plantae = Taxa.objects.get( name = 'echinops <plantae>' )
>>> mammalia = Taxa.objects.get( name = 'mammalia' )

Getting interval parents
~~~~~~~~~~~~~~~~~~~~~~~~

>>> try:
...    taxoref.get_interval_parents( mus, rattus )
... except AssertionError,e:
...    str(e)
'mus is not a parent of rattus'
>>> try:
...     taxoref.get_interval_parents( mus_musculus, muridae )
... except AssertionError,e:
...     str(e)
'mus musculus is not a parent of muridae'
>>> taxoref.get_interval_parents( muridae, mus_musculus )
[<Taxa: mus>]
>>> taxoref.get_interval_parents( mammalia, mus_musculus )
[<Taxa: mus>, <Taxa: muridae>, <Taxa: rodentia>]


Getting common parents
~~~~~~~~~~~~~~~~~~~~~~

>>> taxoref.get_common_parents( [root] )
[]
>>> taxoref.get_common_parents( [muridae] )
[<Taxa: rodentia>, <Taxa: mammalia>, <Taxa: chordata>, <Taxa: animalia>, <Taxa: root>]
>>> taxoref.get_common_parents( [muridae, mus] )
[<Taxa: rodentia>, <Taxa: mammalia>, <Taxa: chordata>, <Taxa: animalia>, <Taxa: root>]
>>> taxoref.get_common_parents( [mus, mus_musculus] )
[<Taxa: muridae>, <Taxa: rodentia>, <Taxa: mammalia>, <Taxa: chordata>, <Taxa: animalia>, <Taxa: root>]
>>> taxoref.get_common_parents( [mus_musculus, mus] )
[<Taxa: muridae>, <Taxa: rodentia>, <Taxa: mammalia>, <Taxa: chordata>, <Taxa: animalia>, <Taxa: root>]

Getting the first common parent
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

>>> taxoref.get_first_common_parent( [mus] )
<Taxa: muridae>
>>> taxoref.get_first_common_parent( [mus, rattus] )
<Taxa: muridae>
>>> taxoref.get_first_common_parent( [rattus, mus] )
<Taxa: muridae>
>>> taxoref.get_first_common_parent( [mus, rattus, muridae] )
<Taxa: rodentia>
>>> taxoref.get_first_common_parent( [mus, muridae] )
<Taxa: rodentia>

Root has no parents
>>> taxoref.get_first_common_parent( [root] ) is None
True
>>> taxoref.get_first_common_parent( [root, mus] ) is None
True

Working with networkx.DiGraph
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

>>> graph = taxoref.get_reference_arborescence( [mus, rattus, echinops_plantae] )
>>> graph.nodes()
[<Taxa: root>, <Taxa: animalia>, <Taxa: magnoliophyta>, <Taxa: asterales>, <Taxa: asteraceae>, <Taxa: magnoliopsida>, <Taxa: mammalia>, <Taxa: mus>, <Taxa: muridae>, <Taxa: chordata>, <Taxa: plantae>, <Taxa: rattus>, <Taxa: rodentia>, <Taxa: echinops <plantae>>]
>>> graph.edges()
[(<Taxa: root>, <Taxa: animalia>), (<Taxa: root>, <Taxa: plantae>), (<Taxa: animalia>, <Taxa: chordata>), (<Taxa: magnoliophyta>, <Taxa: magnoliopsida>), (<Taxa: asterales>, <Taxa: asteraceae>), (<Taxa: asteraceae>, <Taxa: echinops <plantae>>), (<Taxa: magnoliopsida>, <Taxa: asterales>), (<Taxa: mammalia>, <Taxa: rodentia>), (<Taxa: muridae>, <Taxa: mus>), (<Taxa: muridae>, <Taxa: rattus>), (<Taxa: chordata>, <Taxa: mammalia>), (<Taxa: plantae>, <Taxa: magnoliophyta>), (<Taxa: rodentia>, <Taxa: muridae>)]

>>> graph = taxoref.get_reference_arborescence( [mus, rattus, muridae] )
>>> graph.nodes()
[<Taxa: root>, <Taxa: animalia>, <Taxa: mammalia>, <Taxa: mus>, <Taxa: muridae>, <Taxa: chordata>, <Taxa: rattus>, <Taxa: rodentia>]
>>> graph.edges()
[(<Taxa: root>, <Taxa: animalia>), (<Taxa: animalia>, <Taxa: chordata>), (<Taxa: mammalia>, <Taxa: rodentia>), (<Taxa: muridae>, <Taxa: mus>), (<Taxa: muridae>, <Taxa: rattus>), (<Taxa: chordata>, <Taxa: mammalia>), (<Taxa: rodentia>, <Taxa: muridae>)]


"""
