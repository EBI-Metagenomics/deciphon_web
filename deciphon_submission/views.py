from django.http import Http404, HttpResponse
from django.shortcuts import get_object_or_404
from django.views.generic import TemplateView, DetailView
from django.views.generic.detail import BaseDetailView
from rest_framework import viewsets, mixins

from deciphon.models import TargetDb
from deciphon.serializers import TargetDbSerializer
from deciphon_submission.models import SubmittedJob
from deciphon_submission.serializers import SubmittedJobSerializer


class IndexView(TemplateView):
    template_name = "index.html"


class ResultView(DetailView):
    template_name = "result.html"

    model = SubmittedJob
    context_object_name = "submitted_job"


class ResultDownloadView(BaseDetailView):
    def get(self, request, *args, **kwargs):
        submitted_job = get_object_or_404(SubmittedJob, id=self.kwargs["pk"])
        deciphon_job = submitted_job.job
        if not deciphon_job:
            raise Http404()
        if not deciphon_job.results.exists():
            raise Http404("No results available")

        file_format = self.kwargs["filetype"]
        if file_format not in ["faa", "fna", "gff"]:
            raise Http404("File format not accepted")

        if file_format == "faa":
            response = HttpResponse(
                deciphon_job.amino_faa,
                headers={
                    "Content-Type": "application/plain",
                    "Content-Disposition": f'attachment; filename="{submitted_job.id}.faa"',
                },
            )
            return response

        if file_format == "fna":
            response = HttpResponse(
                deciphon_job.codon_fna,
                headers={
                    "Content-Type": "application/plain",
                    "Content-Disposition": f'attachment; filename="{submitted_job.id}.fna"',
                },
            )
            return response

        if file_format == "gff":
            response = HttpResponse(
                deciphon_job.gff,
                headers={
                    "Content-Type": "application/plain",
                    "Content-Disposition": f'attachment; filename="{submitted_job.id}.gff"',
                },
            )
            return response


class RestJobView(
    viewsets.GenericViewSet, mixins.CreateModelMixin, mixins.RetrieveModelMixin
):
    queryset = SubmittedJob.objects.all()
    serializer_class = SubmittedJobSerializer
    permission_classes = []


class RestTargetDbView(viewsets.ReadOnlyModelViewSet):
    queryset = TargetDb.objects.all()
    serializer_class = TargetDbSerializer
    permission_classes = []
