import dataclasses
import time

from django.conf import settings
from django.db import models

from deciphon.result_utils import make_gff, make_fasta


class DeciphonModel(models.Model):
    objects = models.Manager()
    id = models.AutoField(primary_key=True)

    class Meta:
        abstract = True
        managed = getattr(settings, "UNIT_TESTING", False)


ALPHABETS = []


@dataclasses.dataclass
class Alphabet:
    name: str
    display_name: str
    symbols: str

    def __post_init__(self):
        ALPHABETS.append(self)


DNA = Alphabet(name="dna_iupac", display_name="DNA", symbols="actg")
RNA = Alphabet(name="rna_iupac", display_name="RNA", symbols="acug")


class TargetDb(DeciphonModel):
    name = models.CharField(max_length=255)
    filepath = models.CharField(max_length=255)

    class Meta(DeciphonModel.Meta):
        db_table = "db"

    def __str__(self):
        return self.name or str(self.id)


class TimestampField(models.IntegerField):
    def __init__(self, *args, **kwargs):
        self.auto_now = kwargs.pop("auto_now", False)
        super().__init__(*args, **kwargs)

    def pre_save(self, model_instance, add):
        if self.auto_now:
            value = time.time()
            setattr(model_instance, self.attname, value)
            return value
        else:
            return super().pre_save(model_instance, add)


class Job(DeciphonModel):
    PENDING = "pend"
    RUNNING = "run"
    DONE = "done"
    FAILED = "fail"
    STATES = [
        (PENDING, "Pending"),
        (RUNNING, "Running"),
        (DONE, "Done"),
        (FAILED, "Failed"),
    ]
    multi_hits = models.BooleanField(default=False)
    hmmer3_compat = models.BooleanField(default=False)
    target_db = models.ForeignKey(
        TargetDb, on_delete=models.DO_NOTHING, related_name="jobs", db_column="db_id"
    )
    state = models.CharField(max_length=10, choices=STATES, default="pend")
    error = models.CharField(max_length=255, blank=True, default="")
    submission = TimestampField(auto_now=True)
    exec_started = TimestampField(default=0)
    exec_ended = TimestampField(default=0)

    @property
    def gff(self):
        return make_gff(self)

    @property
    def amino_faa(self):
        return make_fasta(self, "amino")

    @property
    def codon_fna(self):
        return make_fasta(self, "codon")

    @property
    def frag_fasta(self):
        return make_fasta(self, "frag")

    class Meta(DeciphonModel.Meta):
        db_table = "job"

    def __str__(self):
        return str(self.id)


class QuerySequence(DeciphonModel):
    name = models.CharField(max_length=255)
    data = models.TextField()
    job = models.ForeignKey(
        Job, on_delete=models.DO_NOTHING, related_name="queries", db_column="job_id"
    )

    class Meta(DeciphonModel.Meta):
        db_table = "seq"

    def __str__(self):
        return self.name or str(self.id)


class Result(DeciphonModel):
    job = models.ForeignKey(
        Job, on_delete=models.DO_NOTHING, related_name="results", db_column="job_id"
    )
    seq = models.ForeignKey(
        QuerySequence,
        on_delete=models.DO_NOTHING,
        related_name="results",
        db_column="seq_id",
    )
    match_id = models.IntegerField()
    prof_name = models.CharField(max_length=255)
    alphabet = models.CharField(max_length=10, db_column="abc_name")
    loglik = models.FloatField()
    null_loglik = models.FloatField()
    model = models.CharField(max_length=10)
    version = models.CharField(max_length=10)
    match_data = models.TextField()

    class Meta(DeciphonModel.Meta):
        db_table = "prod"
