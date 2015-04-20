# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Histone',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('core_type', models.CharField(max_length=25)),
                ('taxonomic_span', models.CharField(max_length=25)),
                ('description', models.CharField(max_length=255)),
            ],
        ),
        migrations.CreateModel(
            name='Sequence',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('gene', models.IntegerField()),
                ('splice', models.IntegerField()),
                ('GI', models.CharField(max_length=25)),
                ('score', models.FloatField()),
                ('evalue', models.FloatField()),
                ('header', models.CharField(max_length=255)),
                ('program', models.CharField(max_length=25)),
                ('sequence', models.TextField()),
                ('alphaN_start', models.IntegerField()),
                ('alphaN_end', models.IntegerField()),
                ('alpha1_start', models.IntegerField()),
                ('alpha1_end', models.IntegerField()),
                ('alpha1ext_start', models.IntegerField()),
                ('alpha1ext_end', models.IntegerField()),
                ('alpha2_start', models.IntegerField()),
                ('alpha2_end', models.IntegerField()),
                ('alpha3_start', models.IntegerField()),
                ('alpha3_end', models.IntegerField()),
                ('alpha3ext_start', models.IntegerField()),
                ('alpha3ext_end', models.IntegerField()),
                ('alphaC_start', models.IntegerField()),
                ('alphaC_end', models.IntegerField()),
                ('beta1_start', models.IntegerField()),
                ('beta1_end', models.IntegerField()),
                ('beta2_start', models.IntegerField()),
                ('beta2_end', models.IntegerField()),
                ('loopL1_start', models.IntegerField()),
                ('loopL1_end', models.IntegerField()),
                ('loopL2_start', models.IntegerField()),
                ('loopL2_end', models.IntegerField()),
                ('mgarg1_start', models.IntegerField()),
                ('mgarg1_end', models.IntegerField()),
                ('mgarg2_start', models.IntegerField()),
                ('mgarg2_end', models.IntegerField()),
                ('mgarg3_start', models.IntegerField()),
                ('mgarg3_end', models.IntegerField()),
                ('docking_domain_start', models.IntegerField()),
                ('docking_domain_end', models.IntegerField()),
                ('core', models.FloatField()),
                ('reviewed', models.BooleanField()),
            ],
        ),
        migrations.CreateModel(
            name='Taxon',
            fields=[
                ('id', models.IntegerField(serialize=False, primary_key=True)),
                ('species', models.CharField(max_length=255)),
                ('genus', models.CharField(max_length=255)),
                ('family', models.CharField(max_length=255)),
                ('phylum', models.CharField(max_length=255)),
                ('domain', models.CharField(max_length=255)),
                ('distance', models.FloatField()),
            ],
        ),
        migrations.CreateModel(
            name='Variant',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('taxonmic_span', models.CharField(max_length=25)),
                ('description', models.CharField(max_length=255)),
                ('core_type', models.ForeignKey(to='server.Histone')),
            ],
        ),
        migrations.AddField(
            model_name='sequence',
            name='species',
            field=models.ForeignKey(to='server.Taxon'),
        ),
        migrations.AddField(
            model_name='sequence',
            name='variant',
            field=models.ForeignKey(to='server.Variant'),
        ),
    ]
