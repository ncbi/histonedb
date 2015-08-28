def make_table_1():
	table = """\documentclass[10pt,a4]{article}
% better looking tables with `\toprule`,`\midrule`,`\bottomrule`:
\usepackage{booktabs}
\begin{table}[h!]
\begin{center}
\begin{tabular}{lcccc}
\toprule
\textbf{Variant} & \textbf{# of Sequences} & \textbf{# of Features} & \textbf{Taxonomic Span} \\
"""
	for hist_type in Histone.objects.all():
		table += "{} & & & ".format(hist_type.id)
		for variant in hist_type.variants.all().order_by("id"):
			sequence_counts = "\begin{tabular}{@\{{}}c@\{{}}}Curated set: {} \\ Automically extracted set: {}\end{tabular}".format(
				Sequence.objects.filter(reviewed=True).count(),
				Sequence.objects.all().count())
			feature_counts = "{}".format(Feature.objects.filter(Q(name="General{}".format(hist_type.id))|Q(name=variant.id)).count())
			table += "   {} & {} & {} & {} \\".format(variant.id, sequence_counts, feature_counts, variant.taxonomic_span)
	table += """\bottomrule
\end{tabular}
\end{center}
\end{table}
\end{document}"""
	print table