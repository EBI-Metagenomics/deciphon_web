import uuid

from django.db import models

from deciphon.models import Job


class SubmittedJob(models.Model):
    """
    Lookup model for user-facing Job UUIDs to internal deciphon integer IDs.
    This is **not** a foreign key because of cross-database foreign-key constraint complexity.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    job_id = models.IntegerField()

    @property
    def job(self):
        return Job.objects.filter(id=self.job_id).first()
