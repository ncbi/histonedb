# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('browse', '0002_auto_20150527_1957'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='score',
            name='search_program',
        ),
        migrations.RemoveField(
            model_name='score',
            name='specificity',
        ),
        migrations.RemoveField(
            model_name='score',
            name='train_program',
        ),
        migrations.AddField(
            model_name='score',
            name='evalue',
            field=models.FloatField(default=1),
            preserve_default=False,
        ),
    ]
