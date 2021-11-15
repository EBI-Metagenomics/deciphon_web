from django_unicorn.components import UnicornView, PollUpdate

from deciphon.models import Job


class ResultView(UnicornView):
    template_name = "unicorn/result.html"

    is_polling: bool = True

    def check_state(self):
        self.submitted_job.job.refresh_from_db()
        if self.submitted_job.job.state in [Job.FAILED, Job.DONE]:
            return PollUpdate(disable=True, method="check_state")
