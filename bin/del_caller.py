import sys
import pysam
import pandas as pd
from collections import defaultdict

REF_SEQ = ("cactagaagctttattgcggtagtttatcacagttaaattgctaacgcagtcagtgggcctcgcggccaagctaggcaatccggtactgttggtaaagccaccatggtcttcacactcgaagatttcgttggggactggcgacagacagc")
REF_LEN = len(REF_SEQ)

##BQ is basecall quality, MQ is mapping quality
BAM_FILE  = sys.argv[1] 
CONTIG    = "ref"
MIN_BQ    = 20
MIN_MQ    = 30
COV_LEN_INC_PROP = 0.99

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

def complement(base):
    return base.translate(COMPLEMENT)

#normalise substitutions to a pyrimide ref
def pyrimidine_snv(ref, alt):
    ref, alt = ref.upper(), alt.upper()
    if ref in ("C", "T"):
        return ref, alt, False
    return complement(ref), complement(alt), True

##Collect reads 
def collect_reads(bam_path, contig=CONTIG, start=0, end=REF_LEN,
                  min_bq=MIN_BQ, min_mq=MIN_MQ):
    base_calls = {}
    deletions = {}
    insertions = {}
    read_spans = {}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(contig, start, end):

            #apply basic filters for processing.
            if (read.is_unmapped or read.is_secondary or read.is_supplementary
                    or read.is_duplicate or read.is_qcfail
                    or read.mapping_quality < min_mq):
                continue

            key   = (read.query_name, read.is_read1)
            quals = read.query_qualities
            seq   = read.query_sequence

            #call snvs before indels
            calls = {}
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if rpos is None or not (start <= rpos < end):
                    continue
                if quals is not None and quals[qpos] < min_bq:
                    continue
                base = seq[qpos].upper()
                if base in "ACGT":
                    calls[rpos] = base
            base_calls[key] = calls
            read_spans[key] = (read.reference_start, read.reference_end)

            #process CIGAR indels
            ref_pos = read.reference_start
            qpos = 0
            read_dels = []
            read_ins = []

            for op, length in (read.cigartuples or []):
                #del
                if op == 2:
                    ds, de = ref_pos, ref_pos + length
                    if de > start and ds < end:
                        cs, ce = max(ds, start), min(de, end)
                        read_dels.append({
                            "del_start": ds,
                            "del_end": de,
                            "del_size": length,
                            "deleted_seq": REF_SEQ[cs:ce].upper(),
                        })
                    ref_pos += length
                #ins
                elif op == 1:
                    if start <= ref_pos < end:
                        read_ins.append({
                            "ins_anchor": ref_pos,
                            "ins_size": length,
                            "ins_seq": seq[qpos:qpos + length].upper(),
                        })
                    qpos += length
                #soft clip
                elif op == 4:
                    qpos += length
                #target type
                elif op in (0, 7, 8):
                    ref_pos += length
                    qpos += length
                #Ns
                elif op == 3:
                    ref_pos += length

            deletions[key] = read_dels
            insertions[key] = read_ins

    return base_calls, deletions, insertions, read_spans

##produce a consensus from pairs based on annotations in collect_reads
def call_molecules(bam_path, contig=CONTIG, start=0, end=REF_LEN, min_bq=MIN_BQ, min_mq=MIN_MQ):
    base_calls, deletions, insertions, read_spans = collect_reads(bam_path, contig, start, end, min_bq, min_mq)
    names = set(name for name, _ in base_calls)
    molecules = []
    snv_rows  = []

    for name in sorted(names):
        r1_calls = base_calls.get((name, True), {})
        r2_calls = base_calls.get((name, False), {})
        r1_dels = deletions.get((name, True), [])
        r2_dels = deletions.get((name, False), [])
        r1_ins = insertions.get((name, True), [])
        r2_ins = insertions.get((name, False), [])

        #extent of coverage
        r1_span = read_spans.get((name, True), (None, None))
        r2_span = read_spans.get((name, False), (None, None))
        r1_start, r1_end = r1_span
        r2_start, r2_end = r2_span

        r1_present = r1_start is not None
        r2_present = r2_start is not None

        #aymmetric consensus
        coverage_asymmetric = False
        stunted_read = None
        r1_stunted = False
        r2_stunted = False

        if r1_present and r2_present:
            r1_len = r1_end - r1_start
            r2_len = r2_end - r2_start
            longer = max(r1_len, r2_len)
            if longer > 0:
                ratio = min(r1_len, r2_len) / longer
                if ratio < COV_LEN_INC_PROP:
                    coverage_asymmetric = True
                    if r1_len < r2_len:
                        stunted_read = "R1"
                        r1_stunted = True
                    else:
                        stunted_read = "R2"
                        r2_stunted = True
        elif r1_present and not r2_present:
            coverage_asymmetric = True
            stunted_read = "R2"
            r2_stunted = True
        elif r2_present and not r1_present:
            coverage_asymmetric = True
            stunted_read = "R1"
            r1_stunted = True

        #Merge deletions
        def del_key(d):
            return (d["del_start"], d["del_end"])

        r1_del_keys = {del_key(d) for d in r1_dels}
        r2_del_keys = {del_key(d) for d in r2_dels}

        del_lookup = {}
        for d in r1_dels + r2_dels:
            k = del_key(d)
            if k not in del_lookup:
                del_lookup[k] = dict(d)
            del_lookup[k]["pair_confirmed"] = (k in r1_del_keys and k in r2_del_keys)
            del_lookup[k]["confirmation_note"] = "both_reads"

        all_deletions = sorted(del_lookup.values(), key=lambda d: d["del_size"], reverse=True)

        focal = all_deletions[0] if all_deletions else None
        secondary = all_deletions[1:]

        deleted_positions = set()
        for d in all_deletions:
            deleted_positions.update(range(d["del_start"], d["del_end"]))

        #merge insertions
        def ins_key(i):
            return (i["ins_anchor"], i["ins_seq"])

        r1_ins_keys = {ins_key(i) for i in r1_ins}
        r2_ins_keys = {ins_key(i) for i in r2_ins}
        ins_lookup = {}

        for i in r1_ins + r2_ins:
            k = ins_key(i)
            if k not in ins_lookup:
                ins_lookup[k] = dict(i)
            ins_lookup[k]["pair_confirmed"] = (k in r1_ins_keys and k in r2_ins_keys)

        all_insertions = sorted(ins_lookup.values(), key=lambda i: i["ins_size"], reverse=True)

        ##classify if reads reach consensus
        all_positions = set(r1_calls) | set(r2_calls)
        concordant = {}
        r1_only = {}
        r2_only = {}

        for pos in all_positions:
            b1, b2 = r1_calls.get(pos), r2_calls.get(pos)
            if b1 and b2:
                if b1 == b2:
                    concordant[pos] = b1
            elif b1:
                r1_only[pos] = b1
            else:
                r2_only[pos] = b2

        ##call snvs
        snvs = []

        def evidence_label(pos, which):
            if r2_stunted and which == "r1" and r2_end is not None and pos >= r2_end:
                return "r1_only"
            if r1_stunted and which == "r2" and r1_end is not None and pos >= r1_end:
                return "r2_only"
            return f"{which}_only"

        def make_snv(pos, base, evidence):
            ref_base = REF_SEQ[pos].upper()
            if base == ref_base:
                return None
            norm_ref, norm_alt, comped = pyrimidine_snv(ref_base, base)
            return {
                "template": name,
                "pos": pos,
                "pos_1based": pos + 1,
                "ref_base": ref_base,
                "alt_base": base,
                "norm_ref": norm_ref,
                "norm_alt": norm_alt,
                "mutation_type": f"{norm_ref}>{norm_alt}",
                "complemented": comped,
                "concordant": evidence == "concordant",
                "evidence": evidence,
            }

        for source, label_fn in ((concordant, lambda p: "concordant"), (r1_only, lambda p: evidence_label(p, "r1")), (r2_only, lambda p: evidence_label(p, "r2"))):
            for pos, base in source.items():
                if pos in deleted_positions:
                    continue
                snv = make_snv(pos, base, label_fn(pos))
                if snv:
                    snvs.append(snv)
                    snv_rows.append(snv)

##classify failures
        #dels confirmed by both reads
        def strictly_confirmed(d):
            return (d.get("pair_confirmed", False)
                and d.get("confirmation_note", "") == "both_reads")

        fail_reasons = []

        #abberant consensus coverage
        if not r1_present or not r2_present or coverage_asymmetric:
            fail_reasons.append("stunted_or_single_reads")

        #abberant consensus support
        if ((focal is not None and not strictly_confirmed(focal))
            or any(not strictly_confirmed(d) for d in secondary)
            or any(not i.get("pair_confirmed", False)
                   for i in all_insertions)
            or any(s["evidence"] != "concordant" for s in snvs)
        ): fail_reasons.append("variant_mismatch")
        fail = "; ".join(fail_reasons) if fail_reasons else None

#read record output
        mol = {
            "template": name,
            "fdel_start": focal["del_start"] if focal else None,
            "fdel_end": focal["del_end"] if focal else None,
            "fdel_size": focal["del_size"] if focal else 0,
            "fdel_seq": focal["deleted_seq"] if focal else None,
            "secondary_deletions": secondary,
            "n_secondary_dels": len(secondary),
            "ins": all_insertions,
            "n_ins": len(all_insertions),
            "snvs": snvs,
            "n_snvs": len(snvs),
            "mutation_class": mutation_class(focal, snvs, all_insertions),
            "r1_present": r1_present,
            "r2_present": r2_present,
            "coverage_asymmetric": coverage_asymmetric,
            "stunted_read": stunted_read,
            "r1_span": r1_span,
            "r2_span": r2_span,
            "fail": fail,        
        }
        molecules.append(mol)

    return molecules, snv_rows

##classify reads by variants
def mutation_class(focal, snvs, insertions):
    has_del = focal is not None
    has_snv = len(snvs) > 0
    has_ins = len(insertions) > 0

    if has_del and has_snv:
        base = "deletion+SNV"
    elif has_del:
        base = "deletion_only"
    elif has_snv:
        base = "SNV_only"
    else:
        base = "reference"

    if has_ins:
        base += "+ins"
    return base

##output consensus calls and snv df
def to_dataframes(molecules, snv_rows):
    flat = []
    for m in molecules:

        sec_summary = "; ".join(
            f"{d['del_start']}-{d['del_end']}({d['del_size']}bp)"
            for d in m["secondary_deletions"]
        ) or None

        ins_summary = "; ".join(
            f"{i['ins_anchor']}+{i['ins_seq']}({i['ins_size']}bp)"
            for i in m["ins"]
        ) or None

        snv_summary = "; ".join(
            f"{s['pos_1based']}{s['mutation_type']}[{s['evidence']}]"
            for s in m["snvs"]
        ) or None

        flat.append({
            "template": m["template"],
            "mut_class": m["mutation_class"],
            "fail": m["fail"],
            "fdel_start": m["fdel_start"],
            "fdel_end": m["fdel_end"],
            "fdel_size": m["fdel_size"],
            "fdel_seq": m["fdel_seq"],
            "n_secondary_dels": m["n_secondary_dels"],
            "secondary_dels": sec_summary,
            "n_ins": m["n_ins"],
            "ins": ins_summary,
            "n_snvs": m["n_snvs"],
            "snvs": snv_summary,
            "r1cov": str(m["r1_span"]),
            "r2cov": str(m["r2_span"]),
        })

    molecule_df = pd.DataFrame(flat)
    molecule_df["fdel_start"] = molecule_df["fdel_start"].astype("Int64")
    molecule_df["fdel_end"]   = molecule_df["fdel_end"].astype("Int64")

    snv_df = pd.DataFrame(snv_rows) if snv_rows else pd.DataFrame()
    return molecule_df, snv_df

if __name__ == "__main__":

    molecules, snv_rows = call_molecules(BAM_FILE)
    mol_df, snv_df = to_dataframes(molecules, snv_rows)
    ##apply filter for output
    mol_df = mol_df[
        (mol_df["mut_class"] == "deletion_only")
        & (mol_df["fail"].isna())
        & (mol_df["fdel_size"] > 1)
    ]
    mol_df.to_csv(sys.argv[2], index=False)



