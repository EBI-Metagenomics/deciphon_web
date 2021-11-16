# Generated by Django 3.2.5 on 2021-11-15 10:25

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = []

    operations = [
        migrations.CreateModel(
            name="Job",
            fields=[
                ("id", models.AutoField(primary_key=True, serialize=False)),
                ("multi_hits", models.BooleanField(default=False)),
                ("hmmer3_compat", models.BooleanField(default=False)),
                (
                    "state",
                    models.CharField(
                        choices=[
                            ("pend", "Pending"),
                            ("run", "Running"),
                            ("done", "Done"),
                            ("fail", "Failed"),
                        ],
                        default="pend",
                        max_length=10,
                    ),
                ),
                ("error", models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                "db_table": "job",
                "abstract": False,
                "managed": False,
            },
        ),
        migrations.CreateModel(
            name="QuerySequence",
            fields=[
                ("id", models.AutoField(primary_key=True, serialize=False)),
                ("name", models.CharField(max_length=255)),
                ("data", models.TextField()),
            ],
            options={
                "db_table": "seq",
                "abstract": False,
                "managed": False,
            },
        ),
        migrations.CreateModel(
            name="Result",
            fields=[
                ("id", models.AutoField(primary_key=True, serialize=False)),
                ("match_id", models.IntegerField()),
                ("prof_name", models.CharField(max_length=255)),
                ("alphabet", models.CharField(db_column="abc_name", max_length=10)),
                ("loglik", models.FloatField()),
                ("null_loglik", models.FloatField()),
                ("model", models.CharField(max_length=10)),
                ("version", models.CharField(max_length=10)),
                ("match_data", models.TextField()),
            ],
            options={
                "db_table": "prod",
                "abstract": False,
                "managed": False,
            },
        ),
        migrations.CreateModel(
            name="TargetDb",
            fields=[
                ("id", models.AutoField(primary_key=True, serialize=False)),
                ("name", models.CharField(max_length=255)),
                ("filepath", models.CharField(max_length=255)),
            ],
            options={
                "db_table": "db",
                "abstract": False,
                "managed": False,
            },
        ),
    ]
