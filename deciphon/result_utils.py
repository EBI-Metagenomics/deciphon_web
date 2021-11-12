import io
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from models import Job, QuerySequence

from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def make_gff(job: 'Job'):
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

            #
            #
            # frag_file.write_item(f"{match_id}", "".join(frags))
            # codon_file.write_item(f"{match_id}", "".join(codons))
            # amino_file.write_item(f"{match_id}", "".join(aminos))

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
                    "Epsilon": "0.01"
                }
            )
            rec.features.append(gff_item)
        gff_records.append(rec)
    gff_io = io.StringIO()
    GFF.write(gff_records, gff_io, True)
    gff_io.seek(0)
    return gff_io
        # for field in item.__dataclass_fields__:
        #     value = getattr(item, field)
        #     out_file.write(str(value))
        #     if field == "attributes":
        #         out_file.write("\n")
        #     else:
        #         out_file.write("\t")

#
# with write_fasta("fragment.fna") as frag_file:
#     with write_fasta("codon.fna") as codon_file:
#         with write_fasta("amino.faa") as amino_file:
#             with open("output.gff", "w") as out_file:
#                 write_friendly(frag_file, codon_file, amino_file, out_file)
#
