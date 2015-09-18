from browse.models import *
from django.db.models import Q

def make_table_1():
    table = """\\documentclass[a4paper,landscape]{article}
    \usepackage[landscape]{geometry}
% better looking tables with `\\toprule`,`\\midrule`,`\\bottomrule`:
\\usepackage{booktabs}
\\begin{document}
\\begin{table}[h!]
\\begin{center}
\\begin{tabular}{lccccc}
\\toprule
\\textbf{Variant} & \\textbf{Curated Sequences} & \\textbf{Auto Sequences} & \\textbf{Num Features} & \\textbf{Taxonomic Span} \\\\
"""
    for hist_type in Histone.objects.all():
        table += "\\toprule \n"
        for variant in hist_type.variants.all().order_by("id"):
            sequence_counts = "{} &".format(variant.sequences.filter(reviewed=True).count())
            sequence_counts += "{} ".format(variant.sequences.count())
            feature_counts = "{}".format(Feature.objects.filter(Q(id__contains="General{}".format(hist_type.id))|Q(id__contains=variant.id)).count())
            table += "        {} & {} & {} & {} \\\\\n".format(variant.id.replace('_',' '), sequence_counts, feature_counts, variant.taxonomic_span)
    table += """\\bottomrule
\end{tabular}
\end{center}
\end{table}
\end{document}
"""
    return table

with open("paper/data_table.tex","w") as f:
    f.write(make_table_1())

