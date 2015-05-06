"""HistoneDB URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.8/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Add an import:  from blog import urls as blog_urls
    2. Add a URL to urlpatterns:  url(r'^blog/', include(blog_urls))
"""
from django.conf.urls import include, url, patterns
from django.contrib import admin
from django.conf.urls.static import static
from django.conf import settings

urlpatterns = patterns('browse.views',
    #url(r'^admin/', include(admin.site.urls)),
    url(r'^$', 'browse_types'),
    url(r'^browse/$', 'browse_types'),

    url(r'^type/([a-zA-Z0-9]+)/$', 'browse_variants'),
    url(r'^variant/([a-zA-Z0-9\.]+)/$', 'browse_variant'),

    url(r'^search/$', 'search'),
    url(r'^upload/$', 'upload'),

    url(r'^data/type/phyloxml/([a-zA-Z0-9]+)/', 'get_phyloxml_of_type'),
    url(r'^data/(type)/json/([a-zA-Z0-9]+)/species/$', 'get_starburst_json'),
    url(r'^data/(type)/json/([a-zA-Z0-9]+)/all/$', 'get_sequence_table_data'),

    url(r'^data/(variant)/json/([a-zA-Z0-9]+)/species/$', 'get_starburst_json'),
    url(r'^data/(variant)/json/([a-zA-Z0-9\.]+)/all/$', 'get_sequence_table_data'),

    url(r'search/([a-zA-Z0-9\.]+)/', 'search'),


)

urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
