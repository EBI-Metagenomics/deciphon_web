import io
from typing import TYPE_CHECKING

from Bio import SeqIO

if TYPE_CHECKING:
    from models import Job

from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

EPSILON = "0.01"


def _make_hits(rec, match_data, hit_id_offset, qualifiers):
    hit_start = 0
    hit_end = 0
    offset = 0
    frags = []
    codons = []
    aminos = []
    hit_start_found = False
    hit_end_found = False

    for frag_match in match_data.split(";"):
        frag, state, codon, amino = frag_match.split(",")
        frags.append(frag)
        codons.append(codon)
        aminos.append(amino)

        if not hit_start_found and (state.startswith("M") or state.startswith("I")):
            hit_start = offset
            hit_start_found = True

        if hit_start_found and not (state.startswith("M") or state.startswith("I")):
            hit_end = offset + len(frag)
            hit_end_found = True

        if hit_end_found:
            hit = SeqFeature(
                FeatureLocation(hit_start, hit_end, strand=None),
                type="CDS",
                qualifiers=dict(qualifiers, ID=hit_id_offset),
            )
            rec.features.append(hit)
            hit_start_found = False
            hit_end_found = False
            hit_id_offset += 1

        offset += len(frag)

    return hit_id_offset


def make_gff(job: "Job"):
    gff_records = []

    hit_id_offset = 1
    for sequence in job.queries.all():
        seq = Seq(sequence.data)
        rec = SeqRecord(seq, sequence.name)

        for result in sequence.results.all():
            lrt = -2 * (result.null_loglik - result.loglik)

            qualifiers = {
                "source": f"deciphon:{result.version}",
                "score": f"{lrt:.17g}",
                "Target_alph": result.alphabet,
                "Profile_acc": result.prof_name,
                "Epsilon": EPSILON,
            }

            hit_id_offset = _make_hits(
                rec, result.match_data, hit_id_offset, qualifiers
            )
            gff_records.append(rec)

        gff_io = io.StringIO()
        GFF.write(gff_records, gff_io, False)
        gff_io.seek(0)
        return gff_io


def make_fasta(job: "Job", fasta_type: str):
    assert fasta_type in ["amino", "frag", "codon"]

    records = []

    for sequence in job.queries.all():
        for result in sequence.results.all():
            start_found = False
            end_found = False

            data = result.match_data
            elements = []

            for match in data.split(";"):
                frag, state, codon, amino = match.split(",")
                if fasta_type == "frag":
                    elements.append(frag)
                if fasta_type == "amino":
                    elements.append(amino)
                if fasta_type == "codon":
                    elements.append(codon)

                if end_found:
                    continue

                if not start_found and state.startswith("M") or state.startswith("I"):
                    start_found = True

                if start_found and not (state.startswith("M") or state.startswith("I")):
                    end_found = True

            records.append(
                SeqRecord(
                    Seq("".join(elements)),
                    id=str(result.match_id),
                )
            )
    fasta_io = io.StringIO()
    SeqIO.write(records, fasta_io, "fasta")
    fasta_io.seek(0)
    return fasta_io
