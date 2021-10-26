import pytest
from Bio.SeqRecord import SeqRecord
from django.test import TestCase

from deciphon.models import DeciphonUser, Job, Target, Alphabet
from deciphon.utils import create_memorable_job_name, alphabet_of_seqrecord


class DeciphonTestCase(TestCase):
    pass


@pytest.mark.django_db
class TestSentinelUser(DeciphonTestCase):
    def setUp(self) -> None:
        super().setUp()
        DeciphonUser.objects.create(username="sentinel", name="Sentinel User")
        DeciphonUser.objects.create(username="spock", name="Ensign Spock")

    def test_sentinel_user(self):
        """There exists exactly one sentinel user"""
        sentinels = DeciphonUser.sentinels
        self.assertEqual(sentinels.count(), 1)
        self.assertEqual(sentinels.first().username, "sentinel")


@pytest.mark.django_db
class TestJobNaming(DeciphonTestCase):
    def setUp(self):
        super().setUp()
        self.user = DeciphonUser.objects.create(
            username="sentinel", name="Sentinel User"
        )
        self.target = Target.objects.create(
            name="trekkersdb", filepath="/to/boldly/go", xxh3=111111
        )
        self.alphabet = Alphabet.objects.create(
            name="dna",
            size=4,
            sym_idx64="yyyy",
            symbols="acgt",
            creation=1,
            type="dna",
            any_symbol="x",
        )

    def test_job_naming(self):
        # Can create a manual job name
        job = Job.objects.create(
            target=self.target, abc=self.alphabet, user=self.user, sid="enterprise"
        )

        # Can create auto named job
        name = create_memorable_job_name()
        job = Job.objects.create(
            target=self.target, abc=self.alphabet, user=self.user, sid=name
        )
        self.assertIn("-", job.sid)


@pytest.mark.django_db
class TestAlphabetDetection(DeciphonTestCase):
    def setUp(self) -> None:
        super().setUp()
        self.dna = Alphabet.objects.create(
            name="dna",
            size=4,
            sym_idx64="yyyy",
            symbols="acgt",
            creation=1,
            type="dna",
            any_symbol="x",
        )
        self.rna = Alphabet.objects.create(
            name="rna",
            size=4,
            sym_idx64="yyyy",
            symbols="acgu",
            creation=1,
            type="rna",
            any_symbol="x",
        )

    def test_alphabet_detection(self):
        record = SeqRecord("actgactgactg")
        alphabet = alphabet_of_seqrecord(record)
        self.assertEqual(alphabet, self.dna)

        record = SeqRecord("acguacgu")
        alphabet = alphabet_of_seqrecord(record)
        self.assertEqual(alphabet, self.rna)

        record = SeqRecord("abcdefg")
        alphabet = alphabet_of_seqrecord(record)
        self.assertIsNone(alphabet)
