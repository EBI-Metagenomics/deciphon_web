from django.contrib import admin
from django.urls import reverse
from django.utils.html import format_html

from deciphon_submission.models import SubmittedJob


@admin.register(SubmittedJob)
class SubmittedJobAdmin(admin.ModelAdmin):
    def deciphon_job_exists(self, obj):
        return obj.job is not None

    def deciphon_job_link(self, obj):
        if self.deciphon_job_exists(obj):
            job = obj.job
            url = reverse(
                "admin:%s_%s_change" % (job._meta.app_label, job._meta.model_name),
                args=[job.id],
            )
            return format_html("<a href='{url}'>{url}</a>", url=url)
        else:
            return "⚠️"

    list_display = ("id", "deciphon_job_exists", "deciphon_job_link", "created")
    list_filter = (
        "id",
        "created",
    )
