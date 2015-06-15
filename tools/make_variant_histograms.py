from browse.models import Histone
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from django.db.models import Min, Max, Count

brewer = list(sns.color_palette("Set2", 7))
print brewer

#Set2 Brewer, used in variant colors
brewer_colors = [
    "#66c2a5",
    "#fc8d62",
    "#8da0cb",
    "#e78ac3",
    "#a6d854",
    "#ffd92f",
    "#e5c494"
    ]

f, axes = plt.subplots(2,3, figsize=(8, 6), sharey=True)
axes[0,0].set_ylabel("Counts")

for i, histone in enumerate([u'H2A', u'H2B', u'H3', u'H4', u'H1']):
    try:
        core_histone = Histone.objects.get(id=histone)
    except:
        continue
    variants = core_histone.variants.annotate(num_sequences=Count('sequences')).order_by("id").all().values_list("id", "num_sequences")
    print variants
    variants = [(v[0], v[1], color) for v, color in zip(variants, brewer_colors)]
    print variants
    names, counts, colors = zip(*variants)
    print names
    print np.array(names)
    ax = axes[int(i>=3), i%3]
    sns.set_palette("Set2", 7)
    g = sns.barplot(names, counts, ci=None, hline=.1, order=names, ax=ax)
    #g.set_xticklabels(ax, rotation=30)
    ax.set_title(histone)


sns.despine(bottom=True)
plt.setp(f.axes, yticks=[])
plt.tight_layout(h_pad=3)
plt.show()




