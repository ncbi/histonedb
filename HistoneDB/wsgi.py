"""
WSGI config for HistoneDB project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/1.8/howto/deployment/wsgi/
"""

import os, sys

# sys.path.append('/home/l_singh/.local/lib/python2.7/site-packages')
# output = ''
# output += 'sys.version = %s\n' % repr(sys.version)
# output += 'sys.prefix = %s\n' % repr(sys.prefix)
# output += 'sys.path = %s' % repr(sys.path)
# print output
#
# import django
# import json

from django.core.wsgi import get_wsgi_application

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "HistoneDB.settings")

GUNICORN = True if (os.getenv('GUNICORN', "0") == "1") else False

application = get_wsgi_application()
if(GUNICORN):
    from whitenoise.django import DjangoWhiteNoise
    application = DjangoWhiteNoise(application)
