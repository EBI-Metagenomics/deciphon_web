from django.db import models


class DeciphonModel(models.Model):
    id = models.AutoField(primary_key=True)

    class Meta:
        abstract = True
        managed = False


class Alphabet(DeciphonModel):
    ALPHABET_TYPES = [
        ("dna", "dna"),
        ("rna", "rna"),
    ]

    name = models.CharField(max_length=15)
    size = models.IntegerField()
    sym_idx64 = models.CharField(max_length=255)
    symbols = models.CharField(max_length=31)
    creation = models.DateTimeField(auto_now=True)
    type = models.CharField(max_length=7, choices=ALPHABET_TYPES)
    any_symbol = models.CharField(max_length=1)

    class Meta(DeciphonModel.Meta):
        db_table = "abc"

    def __str__(self):
        return self.name


class Target(DeciphonModel):
    name = models.CharField(max_length=255)
    filepath = models.CharField(max_length=255)
    xxh3 = models.IntegerField()

    class Meta(DeciphonModel.Meta):
        db_table = "target"

    def __str__(self):
        return self.name or str(self.id)


class DeciphonUser(DeciphonModel):
    username = models.CharField(max_length=31)
    name = models.CharField(max_length=255)

    class Meta(DeciphonModel.Meta):
        db_table = "user"

    def __str__(self):
        return self.username


class Job(DeciphonModel):
    STATUSES = [
        ("pend", "Pending"),
        ("run", "Running"),
        ("done", "Done"),
        ("fail", "Failed"),
    ]
    sid = models.CharField(max_length=19)
    multiple_hits = models.BooleanField(default=False)
    hmmer3_compat = models.BooleanField(default=False)
    abc = models.ForeignKey(
        Alphabet, on_delete=models.DO_NOTHING, related_name="jobs", db_column="abc"
    )
    target = models.ForeignKey(
        Target, on_delete=models.DO_NOTHING, related_name="jobs", db_column="target"
    )
    status = models.CharField(max_length=10, choices=STATUSES)
    status_log = models.CharField(max_length=255, null=True)
    submission = models.DateTimeField(auto_created=True)
    exec_started = models.DateTimeField(null=True)
    exec_ended = models.DateTimeField(null=True)
    user = models.ForeignKey(
        DeciphonUser, on_delete=models.DO_NOTHING, related_name="jobs", db_column="user"
    )

    class Meta(DeciphonModel.Meta):
        db_table = "job"

    def __str__(self):
        return self.sid or str(self.id)


class Query(DeciphonModel):
    name = models.CharField(max_length=255)
    query = models.TextField()
    job = models.ForeignKey(
        Job, on_delete=models.DO_NOTHING, related_name="queries", db_column="job"
    )

    class Meta(DeciphonModel.Meta):
        db_table = "query"

    def __str__(self):
        return self.name or str(self.id)


class Result(DeciphonModel):
    job = models.ForeignKey(
        Job, on_delete=models.DO_NOTHING, related_name="results", db_column="job"
    )
    amino_faa = models.TextField()
    codon_fna = models.TextField()
    output_gff = models.TextField()

    class Meta(DeciphonModel.Meta):
        db_table = "result"
