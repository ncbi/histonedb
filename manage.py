#!/usr/bin/env python
import os
import sys

if __name__ == "__main__":
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "HistoneDB.settings")

#    import sys
#    output = ''
#    output += 'sys.version = %s\n' % repr(sys.version)
#    output += 'sys.prefix = %s\n' % repr(sys.prefix)
#    output += 'sys.path = %s' % repr(sys.path)
#    print('------------------------------')
#    print(output)
#    import django

    from django.core.management import execute_from_command_line

    execute_from_command_line(sys.argv)
