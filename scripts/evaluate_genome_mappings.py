import multiprocessing as mp
from collections import defaultdict, namedtuple
from datetime import timedelta
from pathlib import Path
from pafpy import PafFile
from Bio import SeqIO
from tqdm import tqdm
import numpy as np
import intervaltree as itree

PROCESSES=80
MIN_COMPLEXITY = 0.8

import re

PerformanceProfile = namedtuple("PerformanceProfile", ["cputime", "walltime", "memory"])
BasicAlignment =  namedtuple("BasicAlignment", ["length", "ANI_est", "ANI_wfmash"])


def file_to_profile(input_f):
    kbInGb = 1048576
    user, sys, walltime, cputime, memory = 0, 0, 0, 0, 0
    if not Path(input_f).exists():
        # print(input_f)
        return None
    with open(input_f) as profile_in:
        for line in profile_in:
            m = re.search(r'User time \(seconds\): (?P<time>\d+\D\d+)', line)
            if m:
                user = float(m.groupdict()["time"])
                continue
            m = re.search(r'System time \(seconds\): (?P<time>\d+\D\d+)', line)
            if m:
                sys = float(m.groupdict()["time"])
                continue
            m = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (?P<time>\d+:\d+:\d+)', line)
            if m:
                h, mm, ss = m.groupdict()["time"].split(":")
                walltime = timedelta(seconds=float(ss), minutes=float(mm), hours=float(h))
                continue
            m = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (?P<time>\d+:\d+.\d+)', line)
            if m:
                mm, ss = m.groupdict()["time"].split(":")
                walltime = timedelta(seconds=float(ss), minutes=float(mm))
                continue
            m = re.search(r'Maximum resident set size \(kbytes\): (?P<mem>\d+)', line)
            if m:
                memory = int(m.groupdict()["mem"])
                continue
            
            m = re.search(r'Exit status: (?P<status>\d+)', line)
            if m:
                if m.groupdict()["status"] != "0":
                    return None
                continue
        cputime = timedelta(seconds=user+sys)
    return PerformanceProfile(np.round(cputime.total_seconds() / 60, 2), np.round(walltime.total_seconds() / 60, 2), np.round(memory / kbInGb, 2)) if (cputime and walltime and memory) else None


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

    filtered = defaultdict(list)
    lc = defaultdict(list)
    overall = defaultdict(list)
    
    # header = ["ref", "query", "ss", "pi", "tool", "CPU", "RAM"]

    for query in ["chimpanzee", "macaque"]:
    # for query in ["macaque"]:
        print(f"{ref}--{query}")
        print("-"*20)
        latex_str = f"{query} "
        # for (ss, pi) in [(100, 95), (100, 90), (100, 85)]:
        for pi in [95, 90, 85]:
            for ss in [50, 100, 150, 200]:
                latex_str = f"{query} & {pi}% \n\t"
                for tool in ["mashmap2", "mashmap3"]:
                # for tool in ["mashmap3"]:
                    paf_filepath=f"/home/Users/blk6/Contribute/wfmash/output-pafs-nomerge/{ref}-{query}/{ref}-{query}.L5000.ss{ss}.p{pi}.{tool}.aligned.paf"
                    if Path(paf_filepath).exists():
                        # MashMap outputs mutliple mappings when they have the same predicted ANI. 
                        # To avoid oversampling these regions, we only take the first mapping for each segment
                        with open(paf_filepath) as paf_handle, open("/tmp/xyz.paf", 'w') as filtered_paf_h:
                            filtered_paf_h.writelines(filter(lambda line: (not ("INVALID" in line)), paf_handle))
                        with open("/tmp/xyz.paf") as paf_handle:
                            sampled_records = sorted(
                                list(PafFile(paf_handle)),
                                key=lambda x: (x.qname, x.qstart)
                            )
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

                    filtered_coverage = 0
                    lc_coverage = 0
                    for idx, complexity in enumerate(complexities):
                        r = sampled_records[idx]
                    # for r in sampled_records:
                        err = r.tags["md"].value - r.tags["gi"].value
                        if complexity >= MIN_COMPLEXITY:
                            filtered_errs.append(err)
                            filtered_lens.append(r.qend - r.qstart)
                            # contig_to_itree[r.qname].add(
                                # itree.Interval(r.qstart, r.qend)
                            # )
                            filtered_coverage += r.qend - r.qstart
                        else:
                            lc_errs.append(err)
                            lc_lens.append(r.qend - r.qstart)
                            # lc_contig_to_itree[r.qname].add(
                                # itree.Interval(r.qstart, r.qend)
                            # )
                            lc_coverage += r.qend - r.qstart

                    # for contig, tree in contig_to_itree.items():
                        # tree.merge_overlaps()
                        # for interval in tree:
                            # filtered_coverage += interval.end - interval.begin

                    # for contig, tree in lc_contig_to_itree.items():
                        # tree.merge_overlaps()
                        # for interval in tree:
                            # lc_coverage += interval.end - interval.begin

                    profile = file_to_profile(f"/home/Users/blk6/Contribute/wfmash/output-pafs-nomerge/{ref}-{query}/{ref}-{query}.L5000.ss{ss}.p{pi}.{tool}.approx.err")

                    lc_latex_str = latex_str
                    overall_latex_str = latex_str

                    latex_str += f"{filtered_coverage / 1e9 : .2f} & {profile.cputime} & {profile.memory} & "
                    lc_latex_str += f"{lc_coverage / 1e9 : .2f} & "
                    overall_latex_str += f"{(lc_coverage + filtered_coverage) / 1e9 : .2f} & "

                    latex_str += f"{np.average(filtered_errs):.2f} & "
                    lc_latex_str += f"{np.average(lc_errs):.2f} & "
                    overall_latex_str += f"{np.average(lc_errs + filtered_errs) :.2f} & "

                    latex_str += f"{np.average([abs(x) for x in filtered_errs]):.2f} & \n\t"
                    lc_latex_str += f"{np.average([abs(x) for x in lc_errs]):.2f} & \n\t"
                    overall_latex_str += f"{np.average([abs(x) for x in (lc_errs + filtered_errs)]) :.2f} & \n\t"

                    print(f"{tool:<10}  pi={pi : <4}  ss={ss : <4}  ")
                    print(profile)
                    print(f"Low-complexity: {lc_coverage/1e6 :.2f} ME: {np.average(lc_errs) : 0.2f}\t MAE: {np.average([abs(x) for x in lc_errs]) : 0.2f}")
                    print(f"Filtered: \t{filtered_coverage / 1e9:.2f}\tME: {np.average(filtered_errs) : 0.2f}\t MAE: {np.average([abs(x) for x in filtered_errs]) : 0.2f}")
                    print(f"Overall: {(lc_coverage + filtered_coverage) / 1e9:.2f} ME: {np.average(lc_errs + filtered_errs) : 0.2f}\t MAE: {np.average([abs(x) for x in (lc_errs + filtered_errs)]) : 0.2f}")
                    print()

                filtered[ss].append(latex_str)
                lc[ss].append(lc_latex_str)
                overall[ss].append(overall_latex_str)

    for ss in filtered:
        print(ss)
        print("Low complexity:")
        for line in lc[ss]:
            print(line)

        print()
        print("Filtered:")
        for line in filtered[ss]:
            print(line)

        print()
        print("Overall")
        for line in overall[ss]:
            print(line)

