from io import StringIO
from typing import List, Optional

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from django.db import transaction
from django.db.models import QuerySet
from django.http import Http404, HttpResponseRedirect
from django.urls import reverse
from django_unicorn.components import UnicornView

from deciphon.models import Alphabet, Query, Job, Target, DeciphonUser
from deciphon.utils import create_memorable_job_name, alphabet_of_seqrecord


class SubmitJobView(UnicornView):
    template_name = "unicorn/submit-job.html"

    alphabet_options: QuerySet = None
    target_options: QuerySet = None

    alphabet_selected: Optional[str] = None
    target_selected: str = None

    queryText: str = None
    _seqs: List[SeqRecord] = []

    def mount(self):
        self.alphabet_options = Alphabet.objects.all()
        self.target_options = Target.objects.all()
        default_target = self.target_options.first()
        default_alphabet = self.alphabet_options.first()
        if default_target:
            self.target_selected = str(default_target.id)
        if default_alphabet:
            self.alphabet_selected = str(default_alphabet.id)

    def can_submit(self):
        return (
            len(self._seqs) > 0
            and self.alphabet_selected is not None
            and self.target_selected is not None
        )

    def detect_sequence(self, queryText):
        if not queryText:
            return

        try:
            self.queryText = queryText.replace("⏎", "\n")
            self._seqs = list(SeqIO.parse(StringIO(self.queryText), "fasta"))
            alphabet = alphabet_of_seqrecord(self._seqs[0])
            assert alphabet is not None
        except (AssertionError, IndexError):
            self.alphabet_selected = None
        else:
            self.alphabet_selected = str(alphabet.id)

    def submit(self):
        if not self.can_submit():
            return

        user = DeciphonUser.sentinels.first()
        if not user:
            raise Http404()

        job_name = create_memorable_job_name()

        with transaction.atomic():
            job = Job.objects.create(
                target_id=self.target_selected,
                abc_id=self.alphabet_selected,
                user=user,
                sid=job_name,
            )
            for sequence in self._seqs:
                best_query_name = (
                    sequence.description
                    if sequence.name == "Generated" and sequence.description is not None
                    else sequence.name
                )
                job.queries.create(name=best_query_name, seq=sequence.seq)

        return HttpResponseRedirect(reverse("result", args=(job_name,)))
