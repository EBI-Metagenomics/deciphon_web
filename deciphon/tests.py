import os

import pytest
from Bio.SeqRecord import SeqRecord
from django.conf import settings
from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from django.test import TestCase
from selenium.webdriver import Keys, ActionChains
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.webdriver import WebDriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions
from selenium.webdriver.support.wait import WebDriverWait

from deciphon.models import Job, TargetDb, Result, DNA, RNA, QuerySequence
from deciphon.test_result_fixtures import MATCH1, MATCH2, GFF, FAA, FNA
from deciphon.utils import alphabet_of_seqrecord
from deciphon_submission.models import SubmittedJob


class DeciphonTestCase(TestCase):
    pass


@pytest.mark.django_db
class TestAlphabetDetection(DeciphonTestCase):
    def test_alphabet_detection(self):
        record = SeqRecord("actgactgactg")
        alphabet = alphabet_of_seqrecord(record)
        self.assertEqual(alphabet, DNA)

        record = SeqRecord("acguacgu")
        alphabet = alphabet_of_seqrecord(record)
        self.assertEqual(alphabet, RNA)

        record = SeqRecord("abcdefg")
        alphabet = alphabet_of_seqrecord(record)
        self.assertIsNone(alphabet)





@pytest.mark.django_db
class TestResultsFiles(DeciphonTestCase):
    def setUp(self) -> None:
        self.target_db = TargetDb.objects.create(name='pdb', filepath='/stairway/to/heaven')
        self.job = Job.objects.create(target_db=self.target_db, state=Job.DONE)
        self.seq = QuerySequence.objects.create(job=self.job, name='ZEP', data='actgatcg')
        Result.objects.create(
            job=self.job, seq=self.seq, alphabet=DNA.name, **MATCH1
        )
        Result.objects.create(
            job=self.job, seq=self.seq, alphabet=DNA.name, **MATCH2
        )

    def test_gff(self):
        self.assertEqual(self.job.gff.read(), GFF)

    def test_amino_faa(self):
        self.assertEqual(self.job.amino_faa.read(), FAA)

    def test_codon_dna(self):
        self.assertEqual(self.job.codon_fna.read(), FNA)


@pytest.mark.django_db
class InterfaceTests(StaticLiveServerTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        options = Options()
        options.headless = True
        prefs = {
            "download.default_directory": str(
                os.path.join(settings.BASE_DIR, "downloads")
            )
        }
        options.add_experimental_option("prefs", prefs)
        cls.selenium = WebDriver(options=options)
        cls.selenium.implicitly_wait(10)

    @classmethod
    def tearDownClass(cls):
        cls.selenium.quit()
        super().tearDownClass()

    def setUp(self) -> None:
        super().setUp()
        self.target = TargetDb.objects.create(
            name="trekkersdb", filepath="/to/boldly/go",
        )

    def test_query(self):
        self.selenium.get(self.live_server_url)
        wait = WebDriverWait(self.selenium, 10)

        # Dismiss GDPR banner
        self.selenium.find_element(By.ID, "data-protection-agree").click()

        query_input_component = self.selenium.find_element(By.ID, "queryText")
        query_input = query_input_component.find_element(
            By.XPATH, '//div[@contenteditable="true"]'
        )
        submit_button = self.selenium.find_element(By.ID, "submit")

        # Submit button disabled at first
        self.assertFalse(submit_button.is_enabled())

        # Enter a DNA sequence
        query_input.click()
        query_input.send_keys("actg actg hello")

        # Test the sequence cleanup
        self.selenium.find_element(By.ID, "check-query").click()
        self.assertIn("Generated Header", query_input.text)
        self.assertIn("actgactg", query_input.text)
        self.assertNotIn("hello", query_input.text)

        # Auto detected alphabet
        selected_alphabet = self.selenium.find_element(By.ID, "alphabet-select")
        self.assertEqual(str(DNA.name), selected_alphabet.get_attribute("value"))

        target_radio = self.selenium.find_element(By.ID, f"target_{self.target.id}")
        self.assertTrue(target_radio.is_selected())

        # Submit button enabled
        self.assertTrue(submit_button.is_enabled())

        # Enter an RNA
        query_input.click()
        query_input.clear()
        query_input.send_keys("> RNAQUERY")
        query_input.send_keys(Keys.RETURN)
        query_input.send_keys("acgu acgu rna now")
        self.selenium.find_element(By.ID, "check-query").click()
        self.assertIn("> RNAQUERY", query_input.text)
        self.assertIn("acguacgu", query_input.text)
        self.assertNotIn("now", query_input.text)

        # alphabet should change
        self.assertEqual(str(RNA.name), selected_alphabet.get_attribute("value"))
        self.assertTrue(submit_button.is_enabled())
        ActionChains(self.selenium).move_to_element(submit_button).perform()
        submit_button.click()

        # Should be forwarded to result page
        wait.until(expected_conditions.url_contains("/result/"))

        submitted_job = SubmittedJob.objects.order_by("-created").first()
        self.assertIsNotNone(submitted_job)

        self.assertIn(str(submitted_job.id), self.selenium.current_url)

        self.assertIn(
            "Job is pending", self.selenium.find_element(By.TAG_NAME, "body").text
        )

        job = submitted_job.job
        job.state = Job.RUNNING
        job.save()

        wait.until(
            expected_conditions.text_to_be_present_in_element(
                (By.TAG_NAME, "body"), "Job is running"
            )
        )

        Result.objects.create(
            job=job, seq=job.queries.first(), alphabet=DNA.name, **MATCH1
        )
        Result.objects.create(
            job=job, seq=job.queries.first(), alphabet=DNA.name, **MATCH2
        )
        job.state = Job.DONE
        job.save()
        wait.until(
            expected_conditions.text_to_be_present_in_element(
                (By.TAG_NAME, "body"), "Job complete"
            )
        )

        # Files should be downloaded
        for link in ["Download GFF", "Download FNA", "Download FAA"]:
            link = self.selenium.find_element(By.LINK_TEXT, link)
            link.click()

        #
        # for file_format in ["gff", "fna", "faa"]:
        #     self.assertTrue(
        #         os.path.isfile(
        #             f"{settings.BASE_DIR}/downloads/{str(submitted_job.id)}.{file_format}"
        #         )
        #     )
