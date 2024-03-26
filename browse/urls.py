from django.urls import path
from . import views

app_name = "browse"

urlpatterns = [
    path('browse/', views.browse_types, name="browse_types"),

    path('type/<str:histone_type>/', views.browse_variants, name="browse_variants"),
    path('type/<str:histone_type>/variant/<str:variant>/', views.browse_variant, name="browse_variants_variant"),
    path('variant/<str:variant>/', views.browse_variant_clipped, name="browse_variant_clipped"),
    # url('^type/([a-zA-Z0-9]+)/variant/([a-zA-Z0-9\._]+)/(\d+)', views.browse_variant_with_highlighted_sequence),
    path('type/([a-zA-Z0-9]+)/variant/([a-zA-Z0-9\._]+)/([a-zA-Z0-9\._]+)/',
         views.browse_variant_with_highlighted_sequence),

    path('search/', views.search, name="search"),
    path('analyze/', views.analyze, name="analyze"),
    path('help/', views.help, name="help"),
    path('basket/', views.basket, name="basket"),
    path('human/', views.human, name="human"),

    # Parameters are stored as session variables a GET response
    path('data/sequences/json', views.get_sequence_table_data, name="get_sequence_table_data"),
    path('data/scores/json', views.get_all_scores, name="get_all_scores"),
    path('data/msa/json', views.get_all_sequences, name="get_all_sequences"),
    path('data/features/gff', views.get_sequence_features, name="get_sequence_features"),
    path('data/sequences\+features/json', views.get_aln_and_features, name="get_aln_and_features"),
    path('data/seed/<str:seed>', views.get_seed_aln_and_features, name="get_seed_aln_and_features"),
    path('data/sunburst', views.get_sunburst_json, name="get_sunburst_json"),
    path('', views.browse_types, name="index"),
]
