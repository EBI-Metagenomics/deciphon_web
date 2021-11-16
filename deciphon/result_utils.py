import io
from typing import TYPE_CHECKING

from Bio import SeqIO

if TYPE_CHECKING:
    from models import Job

from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def make_gff(job: "Job"):
    gff_records = []

    for sequence in job.queries.all():
        seq = Seq(sequence.data)
        rec = SeqRecord(seq, sequence.name)

        for result in sequence.results.all():
            frags = []
            codons = []
            aminos = []
            start = 0
            start_found = False
            end = 0
            end_found = False
            offset = 0

            data = result.match_data

            for match in data.split(";"):
                frag, state, codon, amino = match.split(",")
                frags.append(frag)
                codons.append(codon)
                aminos.append(amino)

                if end_found:
                    continue

                if not start_found and state.startswith("M") or state.startswith("I"):
                    start = offset
                    start_found = True

                if start_found and not (state.startswith("M") or state.startswith("I")):
                    start = offset
                    end_found = True

            lrt = -2 * (result.null_loglik - result.loglik)

            gff_item = SeqFeature(
                FeatureLocation(start, end, strand=None),
                type="CDS",
                qualifiers={
                    "source": "deciphon:0.0.1",
                    "score": f"{lrt:.17g}",
                    "ID": result.match_id,
                    "Target_alph": result.alphabet,
                    "Profile_acc": result.prof_name,
                    "Epsilon": "0.01",
                },
            )
            rec.features.append(gff_item)
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
