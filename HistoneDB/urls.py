from django.conf.urls import patterns, include, url
from django.contrib import admin
from django.conf.urls.static import static
from django.conf import settings

urlpatterns = [
    url(r'', include('browse.urls')),
]

urlpatterns += [url(r'^admin/', include(admin.site.urls))]
