from browse.models import *
from django.db.models import Q

#This script is used to export tables from database for futher use in research

with open('seqs.txt','w') as f:
    f.write("gi,hist_type,hist_var,taxid,curated\n")
    for seq in Sequence.objects.all():
        f.write("%s,%s,%s,%s,%s\n"%(seq.gi,seq.variant.hist_type,seq.variant,seq.taxonomy_id,seq.reviewed))

with open('scores.txt','w') as f:
    f.write("gi,hist_var,score,used_for_classification\n")
    for s in Score.objects.all():
        f.write("%s,%s,%s,%s\n"%(s.sequence.gi,s.variant,s.score,s.used_for_classification))

