# Generated by Django 3.1.13 on 2021-12-21 13:48

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('cdl_rest_api', '0003_experiment_created'),
    ]

    operations = [
        migrations.CreateModel(
            name='Coincidences',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('c00', models.FloatField(blank=True, null=True)),
                ('c10', models.FloatField(blank=True, null=True)),
                ('c01', models.FloatField(blank=True, null=True)),
                ('c11', models.FloatField(blank=True, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='Countrates',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('d1', models.PositiveIntegerField(blank=True, null=True)),
                ('d2', models.PositiveIntegerField(blank=True, null=True)),
                ('d3', models.PositiveIntegerField(blank=True, null=True)),
                ('d4', models.PositiveIntegerField(blank=True, null=True)),
                ('d5', models.PositiveIntegerField(blank=True, null=True)),
                ('d6', models.PositiveIntegerField(blank=True, null=True)),
                ('d7', models.PositiveIntegerField(blank=True, null=True)),
                ('d8', models.PositiveIntegerField(blank=True, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='ExperimentData',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('countratePerDetector', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, to='cdl_rest_api.countrates')),
                ('encodedQubitMeasurements', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, to='cdl_rest_api.coincidences')),
            ],
        ),
    ]
