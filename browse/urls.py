from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^projects/HistoneDB2\.0/$', views.browse_types),
    url(r'^$', views.browse_types),
    url(r'^browse/$', views.browse_types),

    url(r'^type/([a-zA-Z0-9]+)/$', views.browse_variants),
    url(r'^type/([a-zA-Z0-9]+)/variant/([a-zA-Z0-9\.]+)/(\d*)$', views.browse_variant),

    url(r'^search/$', views.search),
    url(r'^analyze/$', views.analyze),
    url(r'^help/$', views.help),
    url(r'^basket/$', views.basket),

    #Parameters are stored as session variables a GET response
    url(r'^data/sequences/json$', views.get_sequence_table_data),
    url(r'^data/scores/json$', views.get_all_scores),
    url(r'^data/msa/json$', views.get_all_sequences),
    url(r'^data/features/gff$', views.get_sequence_features),
    url(r'^data/sequences\+features/json$', views.get_aln_and_features),
    url(r'^data/seed/([a-zA-Z0-9\.]+)$', views.get_seed_aln_and_features),
    url(r'^data/sunburst$', views.get_sunburst_json),
]