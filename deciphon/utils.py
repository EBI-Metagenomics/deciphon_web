from typing import Optional

import randomname
import uuid

from Bio.SeqRecord import SeqRecord

from deciphon.models import Job, Alphabet


def create_memorable_job_name():
    random_sid = randomname.get_name()
    allowed_length = Job._meta.get_field("sid").max_length
    if len(random_sid) > allowed_length:
        random_sid = str(uuid.uuid4())[:allowed_length]
    return random_sid


def alphabet_of_seqrecord(record: SeqRecord) -> Optional[Alphabet]:
    """
    Find the first Alphabet that matches a seqrecord's sequence.
    :param record: A BioPython Seqrecord containing the sequence for which an alphabet is needed
    :return: An Alphabet object or None if no match found.
    """
    chars_user = set(record.seq)
    for alphabet in Alphabet.objects.all():
        if chars_user.issubset(set(alphabet.symbols)):
            return alphabet
