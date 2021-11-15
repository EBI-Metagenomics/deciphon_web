from typing import Optional


from Bio.SeqRecord import SeqRecord

from deciphon.models import Job, ALPHABETS, Alphabet


def alphabet_of_seqrecord(record: SeqRecord) -> Optional[Alphabet]:
    """
    Find the first Alphabet that matches a seqrecord's sequence.
    :param record: A BioPython Seqrecord containing the sequence for which an alphabet is needed
    :return: An Alphabet object or None if no match found.
    """
    chars_user = set(record.seq)
    for alphabet in ALPHABETS:
        if chars_user.issubset(set(alphabet.symbols)):
            return alphabet
