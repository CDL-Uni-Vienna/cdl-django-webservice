# Generated by Django 3.1.13 on 2021-12-21 09:53

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('cdl_rest_api', '0002_auto_20211218_1519'),
    ]

    operations = [
        migrations.AddField(
            model_name='experiment',
            name='created',
            field=models.DateTimeField(auto_now_add=True),
            preserve_default=False,
        ),
    ]