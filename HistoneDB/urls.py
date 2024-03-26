from django.urls import include, path
from django.conf import settings
from django.contrib import admin

urlpatterns = [
    path("", include('browse.urls', namespace="browse")),
    path(settings.ADMIN_URL, admin.site.urls),
]
