#
# Taxa
#

def test_taxa():
    return """
>>> from djangophylocore.models import *

>>> taxa = Taxa.objects.get(name = 'echinops <plantae>')
>>> taxa.homonyms.all()
[<HomonymName: echinops>]
>>> taxa.rank
<Rank: genus>

>>> taxa = Taxa.objects.get(name = 'antilocapra americana')
>>> taxa.synonyms.all()
[<SynonymName: antilocapra anteflexa>, <SynonymName: antilope americana>]

>>> taxa = Taxa.objects.get(name = 'mus musculus')
>>> taxa
<Taxa: mus musculus>
>>> taxa.homonyms.all()
[]
>>> taxa.synonyms.all()
[]
>>> taxa.commons.all()
[<CommonName: house mouse (english)>, <CommonName: souris commune (french)>]
>>> taxa.rank
<Rank: species>

Dealing with parents
~~~~~~~~~~~~~~~~~~~~

Get the first parent

>>> taxa.parent
<Taxa: mus>

Parent can be chained

>>> taxa.parent.parent
<Taxa: muridae>

Get all parents by order (closest to furthest). 'parents' attribute is a
simple liste containing Taxa

>>> taxa.parents
[<Taxa: mus>, <Taxa: muridae>, <Taxa: rodentia>, <Taxa: mammalia>, <Taxa: chordata>, <Taxa: animalia>, <Taxa: root>]
>>> list(reversed(taxa.parents))
[<Taxa: root>, <Taxa: animalia>, <Taxa: chordata>, <Taxa: mammalia>, <Taxa: rodentia>, <Taxa: muridae>, <Taxa: mus>]

'root' doesn't have parents

>>> root = Taxa.objects.get(name = 'root')
>>> root.parents
[]

# get info from sources
#>>> taxa.get_id_in_source('ncbi')
#u'2334'

Dealing with rank
-----------------

>>> species = Rank.objects.get(name = 'species')
>>> species.taxas.all()
[<Taxa: antilocapra americana>, <Taxa: candida tropicalis>, <Taxa: echinops ritro>, <Taxa: echinops telfairi>, <Taxa: ephemerellomyces aquilonius>, <Taxa: extatosoma popa>, <Taxa: monarcha axillaris>, <Taxa: mus musculus>, <Taxa: myxilla mariana>, <Taxa: rattus rattus>]

>>> genus = Rank.objects.get(name = 'genus')
>>> genus.taxas.all()
[<Taxa: antilocapra>, <Taxa: candida>, <Taxa: echinops <animalia>>, <Taxa: echinops <plantae>>, <Taxa: ephemerellomyces>, <Taxa: extatosoma>, <Taxa: monarcha <animalia>>, <Taxa: mus>, <Taxa: myxilla>, <Taxa: rattus>]

"""

