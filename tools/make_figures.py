from browse.models import *
from django.db.models import Q

def make_table_1():
	table = """\\documentclass[10pt,a4]{article}
% better looking tables with `\\toprule`,`\\midrule`,`\\bottomrule`:
\\usepackage{booktabs}
\\begin{document}
  \\begin{table}[h!]
    \\begin{center}
      \\begin{tabular}{lcccc}
        \\toprule
        \\textbf{Variant} & \\textbf{Num Sequences} & \\textbf{Num Features} & \\textbf{Taxonomic Span} \\\\
"""
	for hist_type in Histone.objects.all():
		table += "        {} & & & \n".format(hist_type.id)
		for variant in hist_type.variants.all().order_by("id"):
			sequence_counts = "\\begin{tabular}{@{}c@{}}"
                        sequence_counts += "Curated set: {} \\\\".format(variant.sequences.filter(reviewed=True).count())
                        sequence_counts += "Automically extracted set: {}".format(variant.sequences.count())
                        sequence_counts += "\\end{tabular}"
			feature_counts = "{}".format(Feature.objects.filter(Q(name="General{}".format(hist_type.id))|Q(name=variant.id)).count())
			table += "        {} & {} & {} & {} \\\\\n".format(variant.id, sequence_counts, feature_counts, variant.taxonomic_span)
	table += """        \\bottomrule
      \end{tabular}
    \end{center}
  \end{table}
\end{document}"""
	print table

make_table_1()