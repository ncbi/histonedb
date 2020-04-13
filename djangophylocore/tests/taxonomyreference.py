#
# TaxonomyReference
#

def test_taxonomy_reference():
    return """

All objects that inerits from TaxonomyReference will have some more methods.
Those methods provide informations from the taxonomy.

>>> from djangophylocore.models import *
>>> TAXONOMY_TOC = get_taxonomy_toc(True)
>>> taxoref = TaxonomyReference()


Working with names
------------------

>>> taxoref.is_valid_name('rattus')
True
>>> taxoref.is_valid_name('rattus other')
False

>>> taxoref.is_scientific_name('rattus')
True
>>> taxoref.is_scientific_name('antilocapra anteflexa')
False

>>> taxoref.is_common('black rat')
True
>>> taxoref.is_common('rat noir')
False

>>> taxoref.is_synonym('nannomys')
True
>>> taxoref.is_synonym('rattus')
False

>>> taxoref.is_homonym('echinops')
True
>>> taxoref.is_homonym('echinops <plant>')
False

>>> taxoref.is_bad_taxon('rattus')
False
>>> taxoref.is_bad_taxon('abadname')
True
>>> abadname = BadTaxa.objects.create(name = 'abadname')
>>> taxoref.is_bad_taxon('abadname')
True

>>> taxoref.get_name_from_common('mouse')
[<Taxonomy: mus musculus (scientific name)>]
>>> taxoref.get_name_from_common('mus')
[]

>>> taxoref.get_name_from_synonym('nannomys')
[<Taxonomy: mus (scientific name)>]
>>> taxoref.get_name_from_synonym('mus')
[]

>>> taxoref.get_name_from_homonym('echinops')
[<Taxonomy: echinops <mammal> (scientific name)>, <Taxonomy: echinops <plant> (scientific name)>]
>>> taxoref.get_name_from_homonym('homo')
[]

>>> taxoref.strip_taxon_name("rattus")
'rattus'
>>> taxoref.strip_taxon_name("rattus_france")
'rattus'
>>> taxoref.strip_taxon_name("rattus france, delimiter=' '")
'rattus'
>>> taxoref.strip_taxon_name("mus_musculus")
'mus musculus'
>>> taxoref.strip_taxon_name("mus_musculus_france")
'mus musculus'

>>> taxoref.strip_taxon_name("rattus rattus france", delimiter=' ')
'rattus rattus'

'mus something foo' is not in our test database so we're falling down to 'mus'
>>> taxoref.strip_taxon_name("mus something foo", delimiter=' ')
'mus'

>>> taxoref.correct('rattus') is None
True
>>> taxoref.correct('house mouse')
[<Taxonomy: mus musculus (scientific name)>]
>>> taxoref.correct('echinops')
[<Taxonomy: echinops <mammal> (scientific name)>, <Taxonomy: echinops <plant> (scientific name)>]
>>> taxoref.correct('nannomys')
[<Taxonomy: mus (scientific name)>]
>>> taxoref.correct('taxon not in database')
[0]

Getting django objects from name
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
>>> taxoref.get_object_from_name('mus')
<Taxonomy: mus (scientific name)>
>>> taxoref.get_object_from_name('nannomys')
<Taxonomy: nannomys (synonym)>
>>> taxoref.get_object_from_name('echinops')
<Taxonomy: echinops (homonym)>
>>> taxoref.get_object_from_name('badname')
<BadTaxa: badname (0)>

The name must be present into the database

>>> try:
...     taxoref.get_object_from_name('very badname')
... except ValueError, e:
...     e
ValueError('very badname not found in the database',)

Working with taxa
-----------------

>>> root = Taxonomy.objects.get(name = 'root')
>>> muridae = Taxonomy.objects.get(name = 'muridae')
>>> mus = Taxonomy.objects.get(name = 'mus')
>>> mus_musculus = Taxonomy.objects.get(name = 'mus musculus')
>>> rattus = Taxonomy.objects.get(name = 'rattus')
>>> echinops_plantae = Taxonomy.objects.get(name = 'echinops <plant>')
>>> mammalia = Taxonomy.objects.get(name = 'mammalia')

Getting interval parents
~~~~~~~~~~~~~~~~~~~~~~~~

>>> try:
...    taxoref.get_interval_parents(mus, rattus)
... except AssertionError,e:
...    str(e)
'mus (scientific name) is not a parent of rattus (scientific name)'
>>> try:
...     taxoref.get_interval_parents(mus_musculus, muridae)
... except AssertionError,e:
...     str(e)
'mus musculus (scientific name) is not a parent of muridae (scientific name)'

>>> taxoref.get_interval_parents(muridae, mus_musculus)
[<Taxonomy: mus (scientific name)>, <Taxonomy: murinae (scientific name)>]

>>> taxoref.get_interval_parents(mammalia, mus_musculus)
[<Taxonomy: mus (scientific name)>, <Taxonomy: murinae (scientific name)>, <Taxonomy: muridae (scientific name)>, <Taxonomy: muroidea (scientific name)>, <Taxonomy: sciurognathi (scientific name)>, <Taxonomy: rodentia (scientific name)>, <Taxonomy: glires (scientific name)>, <Taxonomy: euarchontoglires (scientific name)>, <Taxonomy: eutheria (scientific name)>, <Taxonomy: theria (scientific name)>]


Getting common parents
~~~~~~~~~~~~~~~~~~~~~~

>>> taxoref.get_common_parents([root])
[]
>>> taxoref.get_common_parents([muridae])
[<Taxonomy: muroidea (scientific name)>, <Taxonomy: sciurognathi (scientific name)>, <Taxonomy: rodentia (scientific name)>, <Taxonomy: glires (scientific name)>, <Taxonomy: euarchontoglires (scientific name)>, <Taxonomy: eutheria (scientific name)>, <Taxonomy: theria (scientific name)>, <Taxonomy: mammalia (scientific name)>, <Taxonomy: amniota (scientific name)>, <Taxonomy: tetrapoda (scientific name)>, <Taxonomy: sarcopterygii (scientific name)>, <Taxonomy: euteleostomi (scientific name)>, <Taxonomy: teleostomi (scientific name)>, <Taxonomy: gnathostomata <vertebrate> (scientific name)>, <Taxonomy: vertebrata (scientific name)>, <Taxonomy: craniata <chordata> (scientific name)>, <Taxonomy: chordata (scientific name)>, <Taxonomy: deuterostomia (scientific name)>, <Taxonomy: coelomata (scientific name)>, <Taxonomy: bilateria (scientific name)>, <Taxonomy: eumetazoa (scientific name)>, <Taxonomy: metazoa (scientific name)>, <Taxonomy: fungi/metazoa group (scientific name)>, <Taxonomy: eukaryota (scientific name)>, <Taxonomy: cellular organisms (scientific name)>, <Taxonomy: root (scientific name)>]

>>> taxoref.get_common_parents([muridae, mus])
[<Taxonomy: muroidea (scientific name)>, <Taxonomy: sciurognathi (scientific name)>, <Taxonomy: rodentia (scientific name)>, <Taxonomy: glires (scientific name)>, <Taxonomy: euarchontoglires (scientific name)>, <Taxonomy: eutheria (scientific name)>, <Taxonomy: theria (scientific name)>, <Taxonomy: mammalia (scientific name)>, <Taxonomy: amniota (scientific name)>, <Taxonomy: tetrapoda (scientific name)>, <Taxonomy: sarcopterygii (scientific name)>, <Taxonomy: euteleostomi (scientific name)>, <Taxonomy: teleostomi (scientific name)>, <Taxonomy: gnathostomata <vertebrate> (scientific name)>, <Taxonomy: vertebrata (scientific name)>, <Taxonomy: craniata <chordata> (scientific name)>, <Taxonomy: chordata (scientific name)>, <Taxonomy: deuterostomia (scientific name)>, <Taxonomy: coelomata (scientific name)>, <Taxonomy: bilateria (scientific name)>, <Taxonomy: eumetazoa (scientific name)>, <Taxonomy: metazoa (scientific name)>, <Taxonomy: fungi/metazoa group (scientific name)>, <Taxonomy: eukaryota (scientific name)>, <Taxonomy: cellular organisms (scientific name)>, <Taxonomy: root (scientific name)>]

>>> taxoref.get_common_parents([mus, mus_musculus])
[<Taxonomy: murinae (scientific name)>, <Taxonomy: muridae (scientific name)>, <Taxonomy: muroidea (scientific name)>, <Taxonomy: sciurognathi (scientific name)>, <Taxonomy: rodentia (scientific name)>, <Taxonomy: glires (scientific name)>, <Taxonomy: euarchontoglires (scientific name)>, <Taxonomy: eutheria (scientific name)>, <Taxonomy: theria (scientific name)>, <Taxonomy: mammalia (scientific name)>, <Taxonomy: amniota (scientific name)>, <Taxonomy: tetrapoda (scientific name)>, <Taxonomy: sarcopterygii (scientific name)>, <Taxonomy: euteleostomi (scientific name)>, <Taxonomy: teleostomi (scientific name)>, <Taxonomy: gnathostomata <vertebrate> (scientific name)>, <Taxonomy: vertebrata (scientific name)>, <Taxonomy: craniata <chordata> (scientific name)>, <Taxonomy: chordata (scientific name)>, <Taxonomy: deuterostomia (scientific name)>, <Taxonomy: coelomata (scientific name)>, <Taxonomy: bilateria (scientific name)>, <Taxonomy: eumetazoa (scientific name)>, <Taxonomy: metazoa (scientific name)>, <Taxonomy: fungi/metazoa group (scientific name)>, <Taxonomy: eukaryota (scientific name)>, <Taxonomy: cellular organisms (scientific name)>, <Taxonomy: root (scientific name)>]

>>> taxoref.get_common_parents([mus_musculus, mus])
[<Taxonomy: murinae (scientific name)>, <Taxonomy: muridae (scientific name)>, <Taxonomy: muroidea (scientific name)>, <Taxonomy: sciurognathi (scientific name)>, <Taxonomy: rodentia (scientific name)>, <Taxonomy: glires (scientific name)>, <Taxonomy: euarchontoglires (scientific name)>, <Taxonomy: eutheria (scientific name)>, <Taxonomy: theria (scientific name)>, <Taxonomy: mammalia (scientific name)>, <Taxonomy: amniota (scientific name)>, <Taxonomy: tetrapoda (scientific name)>, <Taxonomy: sarcopterygii (scientific name)>, <Taxonomy: euteleostomi (scientific name)>, <Taxonomy: teleostomi (scientific name)>, <Taxonomy: gnathostomata <vertebrate> (scientific name)>, <Taxonomy: vertebrata (scientific name)>, <Taxonomy: craniata <chordata> (scientific name)>, <Taxonomy: chordata (scientific name)>, <Taxonomy: deuterostomia (scientific name)>, <Taxonomy: coelomata (scientific name)>, <Taxonomy: bilateria (scientific name)>, <Taxonomy: eumetazoa (scientific name)>, <Taxonomy: metazoa (scientific name)>, <Taxonomy: fungi/metazoa group (scientific name)>, <Taxonomy: eukaryota (scientific name)>, <Taxonomy: cellular organisms (scientific name)>, <Taxonomy: root (scientific name)>]

Getting the first common parent
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

>>> taxoref.get_first_common_parent([mus])
<Taxonomy: murinae (scientific name)>
>>> taxoref.get_first_common_parent([mus, rattus])
<Taxonomy: murinae (scientific name)>
>>> taxoref.get_first_common_parent([rattus, mus])
<Taxonomy: murinae (scientific name)>
>>> taxoref.get_first_common_parent([mus, rattus, muridae])
<Taxonomy: muroidea (scientific name)>
>>> taxoref.get_first_common_parent([mus, muridae])
<Taxonomy: muroidea (scientific name)>

Root has no parents
>>> taxoref.get_first_common_parent([root]) is None
True
>>> taxoref.get_first_common_parent([root, mus]) is None
True

Working with networkx.DiGraph
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

>>> graph = taxoref.get_reference_arborescence([mus, rattus, echinops_plantae])
>>> graph.nodes()
[<Taxonomy: root (scientific name)>, <Taxonomy: fungi/metazoa group (scientific name)>, <Taxonomy: eutheria (scientific name)>, <Taxonomy: rodentia (scientific name)>, <Taxonomy: tetrapoda (scientific name)>, <Taxonomy: amniota (scientific name)>, <Taxonomy: theria (scientific name)>, <Taxonomy: sciurognathi (scientific name)>, <Taxonomy: euteleostomi (scientific name)>, <Taxonomy: streptophytina (scientific name)>, <Taxonomy: muroidea (scientific name)>, <Taxonomy: chordata (scientific name)>, <Taxonomy: euarchontoglires (scientific name)>, <Taxonomy: glires (scientific name)>, <Taxonomy: coelomata (scientific name)>, <Taxonomy: streptophyta (scientific name)>, <Taxonomy: tracheophyta (scientific name)>, <Taxonomy: spermatophyta (scientific name)>, <Taxonomy: echinops <plant> (scientific name)>, <Taxonomy: euphyllophyta (scientific name)>, <Taxonomy: core eudicotyledons (scientific name)>, <Taxonomy: eumetazoa (scientific name)>, <Taxonomy: carduoideae (scientific name)>, <Taxonomy: bilateria (scientific name)>, <Taxonomy: vertebrata (scientific name)>, <Taxonomy: teleostomi (scientific name)>, <Taxonomy: murinae (scientific name)>, <Taxonomy: magnoliophyta (scientific name)>, <Taxonomy: eukaryota (scientific name)>, <Taxonomy: eudicotyledons (scientific name)>, <Taxonomy: viridiplantae (scientific name)>, <Taxonomy: cardueae (scientific name)>, <Taxonomy: metazoa (scientific name)>, <Taxonomy: muridae (scientific name)>, <Taxonomy: embryophyta (scientific name)>, <Taxonomy: sarcopterygii (scientific name)>, <Taxonomy: gnathostomata <vertebrate> (scientific name)>, <Taxonomy: mammalia (scientific name)>, <Taxonomy: deuterostomia (scientific name)>, <Taxonomy: mus (scientific name)>, <Taxonomy: campanulids (scientific name)>, <Taxonomy: cellular organisms (scientific name)>, <Taxonomy: asterales (scientific name)>, <Taxonomy: asteraceae (scientific name)>, <Taxonomy: craniata <chordata> (scientific name)>, <Taxonomy: asterids (scientific name)>, <Taxonomy: rattus (scientific name)>]

>>> graph.edges()
[(<Taxonomy: root (scientific name)>, <Taxonomy: cellular organisms (scientific name)>), (<Taxonomy: fungi/metazoa group (scientific name)>, <Taxonomy: metazoa (scientific name)>), (<Taxonomy: eutheria (scientific name)>, <Taxonomy: euarchontoglires (scientific name)>), (<Taxonomy: rodentia (scientific name)>, <Taxonomy: sciurognathi (scientific name)>), (<Taxonomy: tetrapoda (scientific name)>, <Taxonomy: amniota (scientific name)>), (<Taxonomy: amniota (scientific name)>, <Taxonomy: mammalia (scientific name)>), (<Taxonomy: theria (scientific name)>, <Taxonomy: eutheria (scientific name)>), (<Taxonomy: sciurognathi (scientific name)>, <Taxonomy: muroidea (scientific name)>), (<Taxonomy: euteleostomi (scientific name)>, <Taxonomy: sarcopterygii (scientific name)>), (<Taxonomy: streptophytina (scientific name)>, <Taxonomy: embryophyta (scientific name)>), (<Taxonomy: muroidea (scientific name)>, <Taxonomy: muridae (scientific name)>), (<Taxonomy: chordata (scientific name)>, <Taxonomy: craniata <chordata> (scientific name)>), (<Taxonomy: euarchontoglires (scientific name)>, <Taxonomy: glires (scientific name)>), (<Taxonomy: glires (scientific name)>, <Taxonomy: rodentia (scientific name)>), (<Taxonomy: coelomata (scientific name)>, <Taxonomy: deuterostomia (scientific name)>), (<Taxonomy: streptophyta (scientific name)>, <Taxonomy: streptophytina (scientific name)>), (<Taxonomy: tracheophyta (scientific name)>, <Taxonomy: euphyllophyta (scientific name)>), (<Taxonomy: spermatophyta (scientific name)>, <Taxonomy: magnoliophyta (scientific name)>), (<Taxonomy: euphyllophyta (scientific name)>, <Taxonomy: spermatophyta (scientific name)>), (<Taxonomy: core eudicotyledons (scientific name)>, <Taxonomy: asterids (scientific name)>), (<Taxonomy: eumetazoa (scientific name)>, <Taxonomy: bilateria (scientific name)>), (<Taxonomy: carduoideae (scientific name)>, <Taxonomy: cardueae (scientific name)>), (<Taxonomy: bilateria (scientific name)>, <Taxonomy: coelomata (scientific name)>), (<Taxonomy: vertebrata (scientific name)>, <Taxonomy: gnathostomata <vertebrate> (scientific name)>), (<Taxonomy: teleostomi (scientific name)>, <Taxonomy: euteleostomi (scientific name)>), (<Taxonomy: murinae (scientific name)>, <Taxonomy: mus (scientific name)>), (<Taxonomy: murinae (scientific name)>, <Taxonomy: rattus (scientific name)>), (<Taxonomy: magnoliophyta (scientific name)>, <Taxonomy: eudicotyledons (scientific name)>), (<Taxonomy: eukaryota (scientific name)>, <Taxonomy: fungi/metazoa group (scientific name)>), (<Taxonomy: eukaryota (scientific name)>, <Taxonomy: viridiplantae (scientific name)>), (<Taxonomy: eudicotyledons (scientific name)>, <Taxonomy: core eudicotyledons (scientific name)>), (<Taxonomy: viridiplantae (scientific name)>, <Taxonomy: streptophyta (scientific name)>), (<Taxonomy: cardueae (scientific name)>, <Taxonomy: echinops <plant> (scientific name)>), (<Taxonomy: metazoa (scientific name)>, <Taxonomy: eumetazoa (scientific name)>), (<Taxonomy: muridae (scientific name)>, <Taxonomy: murinae (scientific name)>), (<Taxonomy: embryophyta (scientific name)>, <Taxonomy: tracheophyta (scientific name)>), (<Taxonomy: sarcopterygii (scientific name)>, <Taxonomy: tetrapoda (scientific name)>), (<Taxonomy: gnathostomata <vertebrate> (scientific name)>, <Taxonomy: teleostomi (scientific name)>), (<Taxonomy: mammalia (scientific name)>, <Taxonomy: theria (scientific name)>), (<Taxonomy: deuterostomia (scientific name)>, <Taxonomy: chordata (scientific name)>), (<Taxonomy: campanulids (scientific name)>, <Taxonomy: asterales (scientific name)>), (<Taxonomy: cellular organisms (scientific name)>, <Taxonomy: eukaryota (scientific name)>), (<Taxonomy: asterales (scientific name)>, <Taxonomy: asteraceae (scientific name)>), (<Taxonomy: asteraceae (scientific name)>, <Taxonomy: carduoideae (scientific name)>), (<Taxonomy: craniata <chordata> (scientific name)>, <Taxonomy: vertebrata (scientific name)>), (<Taxonomy: asterids (scientific name)>, <Taxonomy: campanulids (scientific name)>)]

>>> graph = taxoref.get_reference_arborescence([mus, rattus, muridae])
>>> graph.nodes()
[<Taxonomy: root (scientific name)>, <Taxonomy: fungi/metazoa group (scientific name)>, <Taxonomy: eutheria (scientific name)>, <Taxonomy: rodentia (scientific name)>, <Taxonomy: tetrapoda (scientific name)>, <Taxonomy: amniota (scientific name)>, <Taxonomy: theria (scientific name)>, <Taxonomy: sciurognathi (scientific name)>, <Taxonomy: euteleostomi (scientific name)>, <Taxonomy: muroidea (scientific name)>, <Taxonomy: chordata (scientific name)>, <Taxonomy: euarchontoglires (scientific name)>, <Taxonomy: glires (scientific name)>, <Taxonomy: coelomata (scientific name)>, <Taxonomy: eumetazoa (scientific name)>, <Taxonomy: bilateria (scientific name)>, <Taxonomy: vertebrata (scientific name)>, <Taxonomy: teleostomi (scientific name)>, <Taxonomy: murinae (scientific name)>, <Taxonomy: eukaryota (scientific name)>, <Taxonomy: metazoa (scientific name)>, <Taxonomy: muridae (scientific name)>, <Taxonomy: sarcopterygii (scientific name)>, <Taxonomy: gnathostomata <vertebrate> (scientific name)>, <Taxonomy: mammalia (scientific name)>, <Taxonomy: deuterostomia (scientific name)>, <Taxonomy: mus (scientific name)>, <Taxonomy: cellular organisms (scientific name)>, <Taxonomy: craniata <chordata> (scientific name)>, <Taxonomy: rattus (scientific name)>]


>>> graph.edges()
[(<Taxonomy: root (scientific name)>, <Taxonomy: cellular organisms (scientific name)>), (<Taxonomy: fungi/metazoa group (scientific name)>, <Taxonomy: metazoa (scientific name)>), (<Taxonomy: eutheria (scientific name)>, <Taxonomy: euarchontoglires (scientific name)>), (<Taxonomy: rodentia (scientific name)>, <Taxonomy: sciurognathi (scientific name)>), (<Taxonomy: tetrapoda (scientific name)>, <Taxonomy: amniota (scientific name)>), (<Taxonomy: amniota (scientific name)>, <Taxonomy: mammalia (scientific name)>), (<Taxonomy: theria (scientific name)>, <Taxonomy: eutheria (scientific name)>), (<Taxonomy: sciurognathi (scientific name)>, <Taxonomy: muroidea (scientific name)>), (<Taxonomy: euteleostomi (scientific name)>, <Taxonomy: sarcopterygii (scientific name)>), (<Taxonomy: muroidea (scientific name)>, <Taxonomy: muridae (scientific name)>), (<Taxonomy: chordata (scientific name)>, <Taxonomy: craniata <chordata> (scientific name)>), (<Taxonomy: euarchontoglires (scientific name)>, <Taxonomy: glires (scientific name)>), (<Taxonomy: glires (scientific name)>, <Taxonomy: rodentia (scientific name)>), (<Taxonomy: coelomata (scientific name)>, <Taxonomy: deuterostomia (scientific name)>), (<Taxonomy: eumetazoa (scientific name)>, <Taxonomy: bilateria (scientific name)>), (<Taxonomy: bilateria (scientific name)>, <Taxonomy: coelomata (scientific name)>), (<Taxonomy: vertebrata (scientific name)>, <Taxonomy: gnathostomata <vertebrate> (scientific name)>), (<Taxonomy: teleostomi (scientific name)>, <Taxonomy: euteleostomi (scientific name)>), (<Taxonomy: murinae (scientific name)>, <Taxonomy: mus (scientific name)>), (<Taxonomy: murinae (scientific name)>, <Taxonomy: rattus (scientific name)>), (<Taxonomy: eukaryota (scientific name)>, <Taxonomy: fungi/metazoa group (scientific name)>), (<Taxonomy: metazoa (scientific name)>, <Taxonomy: eumetazoa (scientific name)>), (<Taxonomy: muridae (scientific name)>, <Taxonomy: murinae (scientific name)>), (<Taxonomy: sarcopterygii (scientific name)>, <Taxonomy: tetrapoda (scientific name)>), (<Taxonomy: gnathostomata <vertebrate> (scientific name)>, <Taxonomy: teleostomi (scientific name)>), (<Taxonomy: mammalia (scientific name)>, <Taxonomy: theria (scientific name)>), (<Taxonomy: deuterostomia (scientific name)>, <Taxonomy: chordata (scientific name)>), (<Taxonomy: cellular organisms (scientific name)>, <Taxonomy: eukaryota (scientific name)>), (<Taxonomy: craniata <chordata> (scientific name)>, <Taxonomy: vertebrata (scientific name)>)]

"""
