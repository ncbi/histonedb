from django.conf import settings
from django.core.management.base import CommandError, BaseCommand
from django.db import connection
import logging
from optparse import make_option
from djangophylocore.models import *
        
class Command(BaseCommand):
    option_list = BaseCommand.option_list + (
        make_option('--noinput', action='store_false',
                    dest='interactive', default=True,
                    help='Tells Django to NOT prompt the user for input of any kind.'),
   )
    help = "Resets the database for this project."

    def handle(self, *args, **options):
        """
        Resets the database for this project.
    
        Note: Transaction wrappers are in reverse as a work around for
        autocommit, anybody know how to do this the right way?
        """
        if options.get('interactive'):
            confirm = input("""
You have requested a database clean up.
This will IRREVERSIBLY DESTROY
ALL TreeCollection in the database "%s".
Are you sure you want to do this?

Type 'yes' to continue, or 'no' to cancel: """ % (settings.DATABASE_NAME,))
        else:
            confirm = 'yes'
        if confirm != 'yes':
            print("Clean up cancelled.")
            return
        if TreeCollection.objects.all().count():
            cursor = connection.cursor()
            last_id = TreeCollection.objects.latest('id').id
            for i in range(1, last_id + 1):
                try:
                    cursor.execute(
                      'DROP TABLE djangophylocore_reltreecoltaxa%s;' % i)
                except Exception as e:
                    pass
            TreeCollection.objects.all().delete()
            cursor.execute("ALTER TABLE djangophylocore_treecollection AUTO_INCREMENT = 1;")
            BadTaxa.objects.all().delete()
            print("Clean up successful.")
        else:
            print("Nothing to clean.")
        os.system("python manage.py loadtreebase -v")
        print("TreeBase loaded")
