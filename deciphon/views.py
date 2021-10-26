from django.http import Http404, HttpResponse
from django.shortcuts import get_object_or_404
from django.views.generic import TemplateView, DetailView
from django.views.generic.detail import BaseDetailView

from deciphon.models import Job


class IndexView(TemplateView):
    template_name = "index.html"


class ResultView(DetailView):
    template_name = "result.html"

    model = Job

    def get_object(self, queryset=None):
        return get_object_or_404(self.model, sid=self.kwargs['job_sid'])


class ResultDownloadView(BaseDetailView):
    def get_object(self, queryset=None):
        return get_object_or_404(Job, sid=self.kwargs['job_sid'])

    def get(self, request, *args, **kwargs):
        job = self.get_object()
        result = job.results.first()
        if not result:
            raise Http404('No results available')

        file_format = self.kwargs['filetype']
        if file_format not in ['faa', 'fna', 'gff']:
            raise Http404('File format not accepted')

        if file_format == 'faa':
            response = HttpResponse(result.amino_faa, headers={
                'Content-Type': 'application/plain',
                'Content-Disposition': f'attachment; filename="{job.sid}-{result.id}.faa"',
            })
            return response

        if file_format == 'fna':
            response = HttpResponse(result.codon_fna, headers={
                'Content-Type': 'application/plain',
                'Content-Disposition': f'attachment; filename="{job.sid}-{result.id}.fna"',
            })
            return response

        if file_format == 'gff':
            response = HttpResponse(result.output_gff, headers={
                'Content-Type': 'application/plain',
                'Content-Disposition': f'attachment; filename="{job.sid}-{result.id}.gff"',
            })
            return response
