# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.core.validators


class Migration(migrations.Migration):

    dependencies = [
        ('djangophylocore', '0001_initial'),
        ('browse', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='OldStyleVariant',
            fields=[
                ('name', models.CharField(max_length=255, serialize=False, primary_key=True)),
                ('gene', models.IntegerField(null=True, validators=[django.core.validators.MaxValueValidator(15), django.core.validators.MinValueValidator(1)])),
                ('splice', models.IntegerField(null=True, validators=[django.core.validators.MaxValueValidator(15), django.core.validators.MinValueValidator(1)])),
                ('taxonomy', models.ForeignKey(related_name='+', to='djangophylocore.Taxonomy')),
            ],
        ),
        migrations.CreateModel(
            name='Publication',
            fields=[
                ('id', models.IntegerField(serialize=False, primary_key=True)),
                ('cited', models.BooleanField()),
            ],
        ),
        migrations.RemoveField(
            model_name='score',
            name='evalue',
        ),
        migrations.AddField(
            model_name='score',
            name='above_threshold',
            field=models.BooleanField(default=False),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='variant',
            name='aucroc',
            field=models.IntegerField(null=True),
        ),
        migrations.AddField(
            model_name='variant',
            name='hmmthreshold',
            field=models.IntegerField(null=True),
        ),
        migrations.AddField(
            model_name='publication',
            name='variants',
            field=models.ManyToManyField(to='browse.Variant'),
        ),
        migrations.AddField(
            model_name='oldstylevariant',
            name='updated_variant',
            field=models.ForeignKey(related_name='old_names', to='browse.Variant'),
        ),
    ]
