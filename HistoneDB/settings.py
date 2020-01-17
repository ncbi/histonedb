"""
Django settings for HistoneDB project.

Generated by 'django-admin startproject' using Django 1.8.1.

For more information on this file, see
https://docs.djangoproject.com/en/1.8/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.8/ref/settings/
"""

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
import os
import dj_database_url

# import sys
# output = ''
# output += 'sys.version = %s\n' % repr(sys.version)
# output += 'sys.prefix = %s\n' % repr(sys.prefix)
# output += 'sys.path = %s' % repr(sys.path)
# print '------------------------------'
# print output
# import mod_wsgi

GUNICORN = True if (os.getenv('GUNICORN', "0") == "1") else False
if(GUNICORN):
    print("GUNICORN setup enabled")
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Set the MySQL dtabase information
NCBI_database_info = {}
def load_settings(path=os.path.join(BASE_DIR, "HistoneDB", "NCBI_database_info.txt")):
    with open(path) as NCBI_database_info_file:
        for line in NCBI_database_info_file:
            if line.startswith("#"): continue

            try:
                name, value = map(lambda s:s.strip(), line.strip().split("="))
            except:
                continue

            if name == "file" and value != path:
                load_settings(value)
            else:
                NCBI_database_info[name] = value
load_settings()

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.8/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = NCBI_database_info["SECRET_KEY"]
SESSION_COOKIE_HTTPONLY = True

# SECURITY WARNING: don't run with debug turned on in production!
#DEBUG = os.environ.get("NCBI_database_info_DEBUG", "True") == "True"
DEBUG = True
if not DEBUG:
    #X_FRAME_OPTIONS = "DENY"
    #CSRF_COOKIE_HTTPONLY = True
    #CSRF_COOKIE_SECURE = True
    SESSION_COOKIE_SECURE = True
    #SECURE_SSL_REDIRECT = True
    SECURE_BROWSER_XSS_FILTER = True
    SECURE_CONTENT_TYPE_NOSNIFF = True
    #SECURE_HSTS_SECONDS = 0

ALLOWED_HOSTS = ['*']


# Application definition

INSTALLED_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'browse',
    'djangophylocore',
    'django_extensions',
    # 'mod_wsgi.server',
    'analytics'
)

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.auth.middleware.SessionAuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'django.middleware.security.SecurityMiddleware',
)

ROOT_URLCONF = 'HistoneDB.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [os.path.join(BASE_DIR, "templates"),],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

TEMPLATE_DIRS = [os.path.join(BASE_DIR, "templates"),]

TEMPLATE_CONTEXT_PROCESSORS = (
    'django.core.context_processors.request',
)

WSGI_APPLICATION = 'HistoneDB.wsgi.application'

# Database
try:
    db_type=NCBI_database_info["db_type"]
except:
    db_type='mysql'
# https://docs.djangoproject.com/en/1.8/ref/settings/#databases
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.'+db_type, # Add 'postgresql_psycopg2', 'postgresql', 'mysql', 'sqlite3' or 'oracle'.
        'NAME': NCBI_database_info["name"],
        'USER': NCBI_database_info["user"],
        'PASSWORD': NCBI_database_info["password"],
        'HOST': NCBI_database_info["host"],
        'PORT':NCBI_database_info["port"],
        'CONN_MAX_AGE': 3600,
        'SSL_DISABLED': True
    }
}

DATABASE_ENGINE="mysql" #this is for djangophylocore

# Parse database configuration from $DATABASE_URL if available
if(dj_database_url.config()):
    DATABASES['default'] =  dj_database_url.config() # this is for cloud setup
# Internationalization
# https://docs.djangoproject.com/en/1.8/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_L10N = True

USE_TZ = True

GOOGLE_ANALYTICS_ID = NCBI_database_info["GOOGLE_ANALYTICS_ID"]


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.8/howto/static-files/

STATIC_URL = NCBI_database_info["STATIC_URL"]

STATIC_ROOT_AUX = os.path.join(BASE_DIR, "static")

if(not GUNICORN):
    STATICFILES_DIRS = [
        os.path.join(BASE_DIR, "static"),
     ]

if "FORCE_SCRIPT_NAME" in NCBI_database_info:
    FORCE_SCRIPT_NAME = NCBI_database_info["FORCE_SCRIPT_NAME"]

# Simplified static file serving.
# https://warehouse.python.org/project/whitenoise/
if(GUNICORN):
    STATICFILES_STORAGE = 'whitenoise.django.GzipManifestStaticFilesStorage'
    STATIC_ROOT = os.path.join(BASE_DIR, "static")


