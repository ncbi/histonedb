from django.conf.urls import patterns, include, url
from django.contrib import admin
from django.conf.urls.static import static
from django.conf import settings


urlpatterns = patterns('browse.views',
    #url(r'^admin/', include(admin.site.urls)),
    url(r'^$', 'browse_types'),
    url(r'^browse/$', 'browse_types'),

    url(r'^type/([a-zA-Z0-9]+)/$', 'browse_variants'),
    url(r'^type/([a-zA-Z0-9]+)/variant/([a-zA-Z0-9\.]+)/$', 'browse_variant'),

    url(r'^search/$', 'search'),
    url(r'^analyze/$', 'analyze'),
    url(r'^help/$', 'help'),
    url(r'^basket/$', 'basket'),

    #Parameters are stored as session variables a GET response
    url(r'^data/sequences/json$', 'get_sequence_table_data'),
    url(r'^data/scores/json$', 'get_all_scores'),
    url(r'^data/msa/json$', 'get_all_sequences'),
    url(r'^data/features/gff$', 'get_sequence_features'),
    url(r'^data/sequences\+features/json$', 'get_aln_and_features'),
    url(r'^data/seed/([a-zA-Z0-9\.]+)$', 'get_seed_aln_and_features'),

    url(r'^data/(type)/json/([a-zA-Z0-9]+)/species/$', 'get_starburst_json'),
    url(r'^data/(variant)/json/([a-zA-Z0-9\.]+)/species/$', 'get_starburst_json'),
)

# urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
urlpatterns += [url(r'^admin/', include(admin.site.urls))]


if settings.DEBUG:
    import debug_toolbar
    urlpatterns += patterns('', url(r'^__debug__/', include(debug_toolbar.urls)),)
