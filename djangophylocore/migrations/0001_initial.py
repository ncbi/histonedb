# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import djangophylocore.models


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='BadTaxa',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=200)),
                ('nb_occurence', models.IntegerField(default=0)),
            ],
            options={
                'ordering': ['name'],
            },
        ),
        migrations.CreateModel(
            name='ParentsRelation',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('index', models.IntegerField()),
            ],
            options={
                'ordering': ['index'],
            },
        ),
        migrations.CreateModel(
            name='Rank',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=80)),
            ],
            options={
                'ordering': ['name'],
            },
        ),
        migrations.CreateModel(
            name='RelCommonTaxa',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('language', models.CharField(max_length=80, null=True)),
            ],
            options={
                'ordering': ['taxon', 'common'],
            },
        ),
        migrations.CreateModel(
            name='RelHomonymTaxa',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
            options={
                'ordering': ['taxon', 'homonym'],
            },
        ),
        migrations.CreateModel(
            name='RelSynonymTaxa',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
            options={
                'ordering': ['taxon', 'synonym'],
            },
        ),
        migrations.CreateModel(
            name='Taxonomy',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=200)),
                ('type_name', models.CharField(max_length=50)),
                ('_parents', models.ManyToManyField(related_name='children', through='djangophylocore.ParentsRelation', to='djangophylocore.Taxonomy')),
                ('parent', models.ForeignKey(related_name='direct_children', to='djangophylocore.Taxonomy', null=True)),
                ('rank', models.ForeignKey(related_name='taxa', to='djangophylocore.Rank', null=True)),
            ],
            options={
                'ordering': ['name'],
            },
        ),
        migrations.CreateModel(
            name='TaxonomyTreeOccurence',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('user_taxon_name', models.CharField(max_length=200, null=True)),
                ('nb_occurence', models.IntegerField(default=0)),
                ('taxon', models.ForeignKey(related_name='taxonomy_occurences', to='djangophylocore.Taxonomy')),
            ],
        ),
        migrations.CreateModel(
            name='Tree',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=80, null=True)),
                ('delimiter', models.CharField(default=b' ', max_length=5)),
                ('source', models.TextField()),
                ('rooted', models.NullBooleanField()),
                ('description', models.TextField(null=True)),
                ('created', models.DateTimeField()),
                ('updated', models.DateTimeField()),
                ('is_valid', models.BooleanField(default=False)),
                ('_from_collection', models.BooleanField(default=False)),
                ('column_error', models.IntegerField(null=True)),
                ('bad_taxa', models.ManyToManyField(related_name='trees', to='djangophylocore.BadTaxa')),
            ],
            bases=(models.Model, djangophylocore.models.TaxonomyReference),
        ),
        migrations.CreateModel(
            name='TreeCollection',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=80, null=True)),
                ('source', models.TextField(null=True)),
                ('delimiter', models.CharField(default=b' ', max_length=5)),
                ('description', models.TextField(null=True)),
                ('format', models.CharField(max_length=20, null=True)),
                ('created', models.DateTimeField()),
                ('updated', models.DateTimeField()),
            ],
            bases=(models.Model, djangophylocore.models.TaxonomyReference),
        ),
        migrations.AddField(
            model_name='tree',
            name='collection',
            field=models.ForeignKey(related_name='trees', to='djangophylocore.TreeCollection', null=True),
        ),
        migrations.AddField(
            model_name='taxonomytreeoccurence',
            name='tree',
            field=models.ForeignKey(related_name='taxonomy_occurences', to='djangophylocore.Tree'),
        ),
        migrations.AddField(
            model_name='relsynonymtaxa',
            name='synonym',
            field=models.ForeignKey(related_name='synonym_from_taxa', to='djangophylocore.Taxonomy'),
        ),
        migrations.AddField(
            model_name='relsynonymtaxa',
            name='taxon',
            field=models.ForeignKey(related_name='taxa_from_synonym', to='djangophylocore.Taxonomy'),
        ),
        migrations.AddField(
            model_name='relhomonymtaxa',
            name='homonym',
            field=models.ForeignKey(related_name='homonym_from_taxa', to='djangophylocore.Taxonomy'),
        ),
        migrations.AddField(
            model_name='relhomonymtaxa',
            name='taxon',
            field=models.ForeignKey(related_name='taxa_from_homonym', to='djangophylocore.Taxonomy'),
        ),
        migrations.AddField(
            model_name='relcommontaxa',
            name='common',
            field=models.ForeignKey(related_name='common_from_taxa', to='djangophylocore.Taxonomy'),
        ),
        migrations.AddField(
            model_name='relcommontaxa',
            name='taxon',
            field=models.ForeignKey(related_name='taxa_from_common', to='djangophylocore.Taxonomy'),
        ),
        migrations.AddField(
            model_name='parentsrelation',
            name='parent',
            field=models.ForeignKey(related_name='parents_relation_parents', to='djangophylocore.Taxonomy'),
        ),
        migrations.AddField(
            model_name='parentsrelation',
            name='taxon',
            field=models.ForeignKey(related_name='parents_relation_taxa', to='djangophylocore.Taxonomy'),
        ),
        migrations.AlterUniqueTogether(
            name='taxonomytreeoccurence',
            unique_together=set([('taxon', 'tree', 'user_taxon_name')]),
        ),
        migrations.AlterUniqueTogether(
            name='taxonomy',
            unique_together=set([('name', 'type_name')]),
        ),
        migrations.AlterUniqueTogether(
            name='relsynonymtaxa',
            unique_together=set([('synonym', 'taxon')]),
        ),
        migrations.AlterUniqueTogether(
            name='relhomonymtaxa',
            unique_together=set([('homonym', 'taxon')]),
        ),
        migrations.AlterUniqueTogether(
            name='relcommontaxa',
            unique_together=set([('common', 'taxon')]),
        ),
        migrations.AlterUniqueTogether(
            name='parentsrelation',
            unique_together=set([('taxon', 'parent')]),
        ),
    ]
