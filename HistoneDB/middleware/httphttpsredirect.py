from django.http import HttpResponse, HttpResponsePermanentRedirect
import os
class redirect(object):
    # Check if client IP is allowed
    def process_request(self, request):
        if os.getenv('HTTP_X_FORWARDED_PROTO',0):
            return None
        else:
            newpath=request.build_absolute_uri().replace('http','https')
            response= HttpResponsePermanentRedirect(newpath)
            return response

#Add this line to middleware_classes in settings.py
#    'HistoneDB.middleware.httphttpsredirect.redirect',