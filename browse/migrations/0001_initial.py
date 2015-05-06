# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.core.validators


class Migration(migrations.Migration):

    dependencies = [
        ('djangophylocore', '__first__'),
    ]

    operations = [
        migrations.CreateModel(
            name='Histone',
            fields=[
                ('id', models.CharField(max_length=25, serialize=False, primary_key=True)),
                ('taxonomic_span', models.CharField(max_length=25)),
                ('description', models.CharField(max_length=255)),
            ],
        ),
        migrations.CreateModel(
            name='Score',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('score', models.FloatField()),
                ('evalue', models.FloatField()),
                ('train_program', models.CharField(max_length=25)),
                ('search_program', models.CharField(max_length=25)),
                ('specificity', models.IntegerField(null=True, validators=[django.core.validators.MaxValueValidator(100), django.core.validators.MinValueValidator(1)])),
                ('hmmStart', models.IntegerField()),
                ('hmmEnd', models.IntegerField()),
                ('seqStart', models.IntegerField()),
                ('seqEnd', models.IntegerField()),
            ],
        ),
        migrations.CreateModel(
            name='Sequence',
            fields=[
                ('id', models.CharField(max_length=25, serialize=False, primary_key=True)),
                ('gene', models.IntegerField(null=True, validators=[django.core.validators.MaxValueValidator(15), django.core.validators.MinValueValidator(1)])),
                ('splice', models.IntegerField(null=True, validators=[django.core.validators.MaxValueValidator(15), django.core.validators.MinValueValidator(1)])),
                ('header', models.CharField(max_length=255)),
                ('sequence', models.TextField()),
                ('reviewed', models.BooleanField()),
            ],
        ),
        migrations.CreateModel(
            name='Variant',
            fields=[
                ('id', models.CharField(max_length=25, serialize=False, primary_key=True)),
                ('taxonmic_span', models.CharField(max_length=25)),
                ('description', models.CharField(max_length=255)),
                ('core_type', models.ForeignKey(related_name='variants', to='browse.Histone')),
            ],
        ),
        migrations.CreateModel(
            name='Features',
            fields=[
                ('sequence', models.OneToOneField(related_name='features', primary_key=True, serialize=False, to='browse.Sequence')),
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
            ],
        ),
        migrations.CreateModel(
            name='Structure',
            fields=[
                ('sequence', models.OneToOneField(related_name='structures', primary_key=True, serialize=False, to='browse.Sequence')),
                ('pdb', models.CharField(max_length=25)),
                ('mmdb', models.CharField(max_length=25)),
                ('chain', models.CharField(max_length=25)),
            ],
        ),
        migrations.AddField(
            model_name='sequence',
            name='taxonomy',
            field=models.ForeignKey(to='djangophylocore.Taxonomy'),
        ),
        migrations.AddField(
            model_name='sequence',
            name='variant',
            field=models.ForeignKey(related_name='sequences', to='browse.Variant'),
        ),
        migrations.AddField(
            model_name='score',
            name='sequence',
            field=models.ForeignKey(related_name='scores', to='browse.Sequence'),
        ),
        migrations.AddField(
            model_name='score',
            name='variant',
            field=models.ForeignKey(related_name='+', to='browse.Variant'),
        ),
    ]
