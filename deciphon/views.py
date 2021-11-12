from django.http import Http404, HttpResponse
from django.shortcuts import get_object_or_404
from django.views.generic import TemplateView, DetailView
from django.views.generic.detail import BaseDetailView

from deciphon.models import Job, SubmittedJob


class IndexView(TemplateView):
    template_name = "index.html"


class ResultView(DetailView):
    template_name = "result.html"

    model = Job

    def get_object(self, queryset=None):
        return get_object_or_404(self.model, sid=self.kwargs["job_id"])


class ResultDownloadView(BaseDetailView):
    def get_object(self, queryset=None):
        # submission_job = get_object_or_404(SubmittedJob, sid=self.kwargs["sid"])
        # return get_object_or_404(Job, id=submission_job.job_id)
        return get_object_or_404(Job, id=int(self.kwargs["job_sid"]))

    def get(self, request, *args, **kwargs):
        job = self.get_object()
        # result = job.results.first()
        if not job.results.exists():
            raise Http404("No results available")

        file_format = self.kwargs["filetype"]
        if file_format not in ["faa", "fna", "gff"]:
            raise Http404("File format not accepted")

        # if file_format == "faa":
        #     response = HttpResponse(
        #         result.amino_faa,
        #         headers={
        #             "Content-Type": "application/plain",
        #             "Content-Disposition": f'attachment; filename="{job.sid}-{result.id}.faa"',
        #         },
        #     )
        #     return response
        #
        # if file_format == "fna":
        #     response = HttpResponse(
        #         result.codon_fna,
        #         headers={
        #             "Content-Type": "application/plain",
        #             "Content-Disposition": f'attachment; filename="{job.sid}-{result.id}.fna"',
        #         },
        #     )
        #     return response

        if file_format == "gff":
            response = HttpResponse(
                job.gff,
                headers={
                    "Content-Type": "application/plain",
                    "Content-Disposition": f'attachment; filename="{job.id}.gff"',
                },
            )
            return response
