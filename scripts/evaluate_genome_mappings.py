import multiprocessing as mp
from collections import defaultdict
from pathlib import Path
from pafpy import PafFile
from Bio import SeqIO
from tqdm import tqdm
import numpy as np
import intervaltree as itree

PROCESSES=80
MIN_COMPLEXITY = 0.8


def index_fasta(fasta_file):
    return {record.id: record.seq.upper() for record in (SeqIO.parse(fasta_file, "fasta"))}

def get_seq_complexity(seq):
    k=19
    kmers = set(str(seq[i:i+k]) for i in range(len(seq) - k + 1))
    return len(kmers) / (len(seq) - k + 1)

if __name__ == "__main__":
    chimp_index = index_fasta("/home/Users/blk6/Data/assemblies/chimpanzee/Clint_PTRv2.fna")
    macaque_index =  index_fasta("/home/Users/blk6/Data/assemblies/macaque/GCF_003339765.1_Mmul_10_genomic.fna")

    name_to_index = {
        "chimpanzee": chimp_index,
        "macaque": macaque_index
    }

    ref = "human"

    filtered = []
    lc = []
    overall = []

    for query in ["chimpanzee", "macaque"]:
        print(f"{ref}--{query}")
        print("-"*20)
        latex_str = f"{query} "
        for (ss, pi) in [(200, 95), (250, 90), (300, 85)]:
            latex_str = f"{query} & {pi}% \n\t"
            for tool in ["mashmap2", "mashmap3"]:
                paf_filepath=f"/home/Users/blk6/Contribute/wfmash/output-pafs/{ref}-{query}/{ref}-{query}.L5000.ss{ss}.p{pi}.{tool}.aligned.paf"
                if Path(paf_filepath).exists():
                    # MashMap outputs mutliple mappings when they have the same predicted ANI. 
                    # To avoid oversampling these regions, we only take the first mapping for each segment
                    with open(paf_filepath) as paf_handle:
                        sampled_records = sorted(
                            list(PafFile(paf_handle)),
                            key=lambda x: (x.qname, x.qstart)
                        )
                    sampled_records = [
                        sampled_records[i] for i in range(len(sampled_records)) 
                        if i==0 or sampled_records[i].qstart != sampled_records[i-1].qstart
                    ]
                else:
                    continue
                if len(sampled_records) < 1:
                    print(f"No records for {tool:<10}  pi={pi : <4}  ss={ss : <4}  ")
                    print()
                    continue


                lc_errs = []
                lc_lens = []
                filtered_errs = []
                filtered_lens = []
                contig_to_itree = defaultdict(itree.IntervalTree)
                lc_contig_to_itree = defaultdict(itree.IntervalTree)
                with mp.Pool(PROCESSES) as pool:
                    complexities = list(tqdm(pool.imap(
                        get_seq_complexity, 
                        [name_to_index[query][r.qname][r.qstart:r.qend] for r in sampled_records]
                    )))
                for idx, complexity in enumerate(complexities):
                    r = sampled_records[idx]
                    err = r.tags["md"].value - r.tags["gi"].value
                    if complexity >= MIN_COMPLEXITY:
                        filtered_errs.append(err)
                        filtered_lens.append(r.qend - r.qstart)
                        contig_to_itree[r.qname].add(
                            itree.Interval(r.qstart, r.qend)
                        )
                    else:
                        lc_errs.append(err)
                        lc_lens.append(r.qend - r.qstart)
                        lc_contig_to_itree[r.qname].add(
                            itree.Interval(r.qstart, r.qend)
                        )

                filtered_coverage = 0
                for contig, tree in contig_to_itree.items():
                    tree.merge_overlaps()
                    for interval in tree:
                        filtered_coverage += interval.end - interval.begin

                lc_coverage = 0
                for contig, tree in lc_contig_to_itree.items():
                    tree.merge_overlaps()
                    for interval in tree:
                        lc_coverage += interval.end - interval.begin

                lc_latex_str = latex_str
                overall_latex_str = latex_str

                latex_str += f"{filtered_coverage / 1e9 : .2f} & "
                lc_latex_str += f"{lc_coverage / 1e9 : .2f} & "
                overall_latex_str += f"{(lc_coverage + filtered_coverage) / 1e9 : .2f} & "

                latex_str += f"{np.average(filtered_errs, weights=filtered_lens):.2f} & "
                lc_latex_str += f"{np.average(lc_errs, weights=lc_lens):.2f} & "
                overall_latex_str += f"{np.average(lc_errs + filtered_errs, weights=(lc_lens + filtered_lens)) :.2f} & "

                latex_str += f"{np.average([abs(x) for x in filtered_errs], weights=filtered_lens):.2f} & \n\t"
                lc_latex_str += f"{np.average([abs(x) for x in lc_errs], weights=lc_lens):.2f} & \n\t"
                overall_latex_str += f"{np.average([abs(x) for x in (lc_errs + filtered_errs)]) :.2f} & \n\t"

                print(f"{tool:<10}  pi={pi : <4}  ss={ss : <4}  ")
                print(f"Low-complexity: {lc_coverage/1e9 :.2f} ME: {np.average(lc_errs, weights=lc_lens) : 0.2f}\t MAE: {np.average([abs(x) for x in lc_errs], weights=lc_lens) : 0.2f}")
                print(f"Filtered: {filtered_coverage / 1e9:.2f} ME: {np.average(filtered_errs, weights=filtered_lens) : 0.2f}\t MAE: {np.average([abs(x) for x in filtered_errs], weights=filtered_lens) : 0.2f}")
                print(f"Overall: {(lc_coverage + filtered_coverage) / 1e9:.2f} ME: {np.average(lc_errs + filtered_errs) : 0.2f}\t MAE: {np.average([abs(x) for x in (lc_errs + filtered_errs)]) : 0.2f}")
                print()

            filtered.append(latex_str)
            lc.append(lc_latex_str)
            overall.append(overall_latex_str)

    print("Low complexity:")
    for line in lc:
        print(line)

    print()
    print("Filtered:")
    for line in filtered:
        print(line)

    print()
    print("Overall")
    for line in overall:
        print(line)

