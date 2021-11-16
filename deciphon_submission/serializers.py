from django.db import transaction
from django.http import Http404
from django.urls import reverse
from rest_framework import serializers

from deciphon.models import Job, QuerySequence, TargetDb
from deciphon.serializers import DeciphonJobSerializer
from deciphon_submission.models import SubmittedJob


class SubmittedJobSerializer(serializers.ModelSerializer):
    job = DeciphonJobSerializer()

    result_urls = serializers.SerializerMethodField()

    def get_result_urls(self, obj):
        return {
            "gff": reverse("result", args=(obj.id, "gff")),
            "amino_faa": reverse("result", args=(obj.id, "faa")),
            "codon_fna": reverse("result", args=(obj.id, "fna")),
        }

    class Meta:
        model = SubmittedJob
        fields = ["id", "job", "result_urls"]
        read_only_fields = ["id", "result_urls"]

    def create(self, validated_data):
        job_data = validated_data.pop("job")
        queries_data = job_data.pop("queries")
        target_db_data = job_data.pop("target_db")
        try:
            target_db = TargetDb.objects.get(**target_db_data)
        except (TargetDb.DoesNotExist, TargetDb.MultipleObjectsReturned):
            raise Http404("Target database not specificied correctly")
        with transaction.atomic():
            job = Job.objects.create(target_db=target_db, **job_data)
            for query_data in queries_data:
                QuerySequence.objects.create(job=job, **query_data)
            submitted_job = SubmittedJob.objects.create(job_id=job.id, **validated_data)
        return submitted_job
