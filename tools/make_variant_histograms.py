from browse.models import Histone
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from django.db.models import Min, Max, Count

#Set2 Brewer, used in variant colors
colors = [
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

for i, histone in enumerate(["H2A", "H2B", "H3", "H4", "H1"]):
	variants = Histone.objects.get(id=histone).variants.annotate(num_sequences=Count('sequences')).order_by("id").all().values_list("id", "num_sequences")
	variants = [(id, num, color) for (id, num), color in zip(sorted(variants, key=lambda v:v[0]), colors)]
	names, counts, colors = zip(*variants)
	print names
	sns.barplot(np.array(names), np.array(counts), color=np.array(colors), ci=None, hline=.1, ax=axes[int(i>=3), i%3])
	axes[int(i>3), i%3].set_title(histone)

sns.despine(bottom=True)
plt.setp(f.axes, yticks=[])
plt.tight_layout(h_pad=3)
plt.show()




