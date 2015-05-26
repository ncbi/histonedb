from django.conf.urls import patterns, include, url
from django.contrib import admin
from django.conf.urls.static import static
from django.conf import settings


urlpatterns = patterns('browse.views',
    #url(r'^admin/', include(admin.site.urls)),
    url(r'^$', 'browse_types'),
    url(r'^browse/$', 'browse_types', name="browse_types"),

    url(r'^type/([a-zA-Z0-9]+)/$', 'browse_variants', name="browse_variants"),
    url(r'^variant/([a-zA-Z0-9\.]+)/$', 'browse_variant', name="browse_variant"),

    url(r'^search/$', 'search', name="search"),
    url(r'^upload/$', 'upload', name="upload"),

    #url(r'^data/type/phyloxml/([a-zA-Z0-9]+)/', 'get_phyloxml_of_type'),
    url(r'^data/(type)/json/([a-zA-Z0-9]+)/species/$', 'get_starburst_json'),
    url(r'^data/(type)/json/([a-zA-Z0-9]+)/all/$', 'get_sequence_table_data'),

    url(r'^data/(variant)/json/([a-zA-Z0-9\.]+)/species/$', 'get_starburst_json'),
    url(r'^data/(variant)/json/([a-zA-Z0-9\.]+)/all/$', 'get_sequence_table_data'),

    url(r'^treemap/$', 'treemap')


)

urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
