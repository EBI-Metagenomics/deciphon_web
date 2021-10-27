import pytest
from Bio.SeqRecord import SeqRecord
from django.test import TestCase, LiveServerTestCase
from selenium.webdriver import Keys, ActionChains
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.webdriver import WebDriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions
from selenium.webdriver.support.wait import WebDriverWait

from deciphon.models import DeciphonUser, Job, Target, Alphabet, Result
from deciphon.test_result_fixtures import AMINO_FAA, CODON_FNA, GFF_OUTPUT
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


@pytest.mark.django_db
class InterfaceTests(LiveServerTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        options = Options()
        options.headless = True
        cls.selenium = WebDriver(options=options)
        cls.selenium.implicitly_wait(10)

    @classmethod
    def tearDownClass(cls):
        cls.selenium.quit()
        super().tearDownClass()

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
        self.user = DeciphonUser.objects.create(
            username="sentinel", name="Sentinel User"
        )
        self.target = Target.objects.create(
            name="trekkersdb", filepath="/to/boldly/go", xxh3=111111
        )

    def test_query(self):
        self.selenium.get(self.live_server_url)
        wait = WebDriverWait(self.selenium, 10)

        # Dismiss GDPR banner
        self.selenium.find_element(By.ID, 'data-protection-agree').click()

        query_input_component = self.selenium.find_element(By.ID, 'queryText')
        query_input = query_input_component.find_element(By.XPATH, '//div[@contenteditable="true"]')
        submit_button = self.selenium.find_element(By.ID, 'submit')

        # Submit button disabled at first
        self.assertFalse(submit_button.is_enabled())

        # Enter a DNA sequence
        query_input.click()
        query_input.send_keys('actg actg hello')

        # Test the sequence cleanup
        self.selenium.find_element(By.ID, 'check-query').click()
        self.assertIn("Generated Header", query_input.text)
        self.assertIn("actgactg", query_input.text)
        self.assertNotIn('hello', query_input.text)

        # Auto detected alphabet
        selected_alphabet = self.selenium.find_element(By.ID, 'alphabet-select')
        self.assertEqual(str(self.dna.id), selected_alphabet.get_attribute('value'))

        target_radio = self.selenium.find_element(By.ID, f'target_{self.target.id}')
        wait.until(expected_conditions.element_to_be_selected(target_radio))

        # Submit button enabled
        self.assertTrue(submit_button.is_enabled())

        # Enter an RNA
        query_input.click()
        query_input.clear()
        query_input.send_keys('> RNAQUERY')
        query_input.send_keys(Keys.RETURN)
        query_input.send_keys('acgu acgu rna now')
        self.selenium.find_element(By.ID, 'check-query').click()
        self.assertIn("> RNAQUERY", query_input.text)
        self.assertIn("acguacgu", query_input.text)
        self.assertNotIn('now', query_input.text)

        # alphabet should change
        self.assertEqual(str(self.rna.id), selected_alphabet.get_attribute('value'))
        self.assertTrue(submit_button.is_enabled())
        ActionChains(self.selenium).move_to_element(submit_button).perform()
        submit_button.click()

        # Should be forwarded to result page
        wait.until(expected_conditions.url_contains('/result/'))

        job = Job.objects.order_by('-id').first()
        self.assertIsNotNone(job)

        self.assertIn(job.sid, self.selenium.current_url)

        self.assertIn('Job is pending', self.selenium.find_element(By.TAG_NAME, 'body').text)

        job.status = Job.RUNNING
        job.save()

        wait.until(expected_conditions.text_to_be_present_in_element((By.TAG_NAME, 'body'), 'Job is running'))

        Result.objects.create(
            job=job,
            amino_faa=AMINO_FAA,
            codon_fna=CODON_FNA,
            output_gff=GFF_OUTPUT
        )
        job.status = Job.DONE
        job.save()
        wait.until(expected_conditions.text_to_be_present_in_element((By.TAG_NAME, 'body'), 'Job complete'))

        for link in ['Download GFF', 'Download FNA', 'Download FAA']:
            link = self.selenium.find_element(By.LINK_TEXT, link)
            link.click()