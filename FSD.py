import subprocess
import time
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
import pysam

RANGE_START = 100
RANGE_END = 220


def filter_bam_to_bed(
    chrom, bam_file, output_bed, range_start=RANGE_START, range_end=RANGE_END
):
    """Filter BAM to BED."""
    ranges = []
    unpaired = {}
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        for read in samfile.fetch(chrom):
            if (
                read.is_paired
                and read.mapping_quality >= 30
                and read.mate_is_mapped
                and not read.is_duplicate
                and not read.is_secondary
            ):
                qname = read.query_name
                unpaired_read = unpaired.get(qname, None)
                if unpaired_read is None:
                    unpaired[qname] = [read.reference_start, read.reference_end]
                else:
                    paired_reads = unpaired_read + [
                        read.reference_start,
                        read.reference_end,
                    ]
                    read_max = max(paired_reads)
                    read_min = min(paired_reads)
                    read_width = read_max - read_min
                    if range_start <= read_width <= range_end:
                        ranges.append([chrom, read_min, read_max])
                    del unpaired[qname]
    bed_cols = ["chrom", "chromStart", "chromEnd"]
    df_ranges = pd.DataFrame(ranges, columns=bed_cols).sort_values("chromStart")
    df_ranges.to_csv(output_bed, sep="\t", index=False, header=False)


def filter_and_overlap_bed(chrom, bed_file, output_dir):
    """Calculate GC; filter and count overlaps."""
    gc_bed_file = f"{output_dir}/{chrom}.gc.bed"
    with open(gc_bed_file, "w", encoding="ascii") as out_bed:
        subprocess.run(
            ["bedtools", "nuc", "-fi", "/path/to/genome.fa", "-bed", bed_file],
            stdout=out_bed,
            check=True,
        )

    filtered_bed_file = f"{output_dir}/{chrom}.gc.filtered.bed"
    with open(filtered_bed_file, "w", encoding="ascii") as out_bed:
        intersected = subprocess.Popen(
            [
                "bedtools",
                "intersect",
                "-a",
                gc_bed_file,
                "-b",
                f"/path/to/{chrom}.filters.bed",
                "-wa",
                "-v",
            ],
            stdout=subprocess.PIPE,
        )
        subprocess.run(
            ["cut", "-f1,2,3,8,15"],
            stdin=intersected.stdout,
            stdout=out_bed,
            check=True,
        )

    overlap_file = f"{output_dir}/{chrom}.overlap.tsv"
    with open(overlap_file, "w", encoding="ascii") as out_tsv:
        overlaped = subprocess.Popen(
            [
                "bedtools",
                "intersect",
                "-a",
                f"/path/to/{chrom}.AB.bed",
                "-b",
                filtered_bed_file,
                "-wa",
                "-wb",
                "-loj",
            ],
            stdout=subprocess.PIPE,
        )
        subprocess.run(
            ["cut", "-f1,4,10,11"], stdin=overlaped.stdout, stdout=out_tsv, check=True
        )


def frag_stats(chrom):
    """Process statistics on fragments."""
    bed_file = f"{OUTPUT_DIR}/{chrom}.bed"
    filter_bam_to_bed(chrom, BAM_FILE, bed_file)
    filter_and_overlap_bed(chrom, bed_file, OUTPUT_DIR)


def frag_arm(sample_id, output_dir):
    AB_num = {
        1: 2188,
        2: 2330,
        3: 1942,
        4: 1861,
        5: 1746,
        6: 1663,
        7: 1513,
        8: 1407,
        9: 1080,
        10: 1271,
        11: 1296,
        12: 1296,
        13: 949,
        14: 870,
        15: 778,
        16: 750,
        17: 763,
        18: 733,
        19: 553,
        20: 587,
        21: 328,
        22: 332,
    }

    overlap_all = f"{output_dir}/overlap.tsv"
    Path(overlap_all).unlink(missing_ok=True)
    overlap_file_list = sorted(
        list(output_dir.glob("*.overlap.tsv")), key=lambda x: int(x.name.split(".")[0])
    )
    for ovlp in overlap_file_list:
        with open(ovlp, "r", encoding="ascii") as file_in:
            with open(overlap_all, "a", encoding="ascii") as file_out:
                file_out.write(file_in.read())

    # count overlaps & calc mean(gc) by AB (and width)
    overlap = pd.read_table(
        overlap_all,
        header=None,
        names=["chrom", "AB", "gc", "width"],
        dtype={"chrom": "Int64", "AB": "Int64", "gc": "Float64", "width": "Int64"},
        na_values={"gc": ".", "width": "-1"},
    )
    overlap = overlap.dropna()

    zero_overlap_AB = pd.DataFrame(
        [(chr, i, 0) for chr, n_ab in AB_num.items() for i in range(1, n_ab + 1)],
        columns=["chrom", "AB", "nul"],
    )
    zero_overlap_AB = zero_overlap_AB.set_index(["chrom", "AB"])

    counts = overlap.groupby(["chrom", "AB"]).count()
    counts = counts.drop(columns="width")

    counts = pd.merge(
        counts, zero_overlap_AB, how="outer", left_index=True, right_index=True
    )
    counts = counts.fillna(0)
    counts["gc"].to_csv(f"{output_dir}/counts.csv", index=False, header=False)

    zero_overlap_AB_width = pd.DataFrame(
        [
            (chr, i, j, 0)
            for chr, n_ab in AB_num.items()
            for i in range(1, n_ab + 1)
            for j in range(RANGE_START, RANGE_END + 1)
        ],
        columns=["chrom", "AB", "width", "nul"],
    )
    zero_overlap_AB_width = zero_overlap_AB_width.set_index(["chrom", "AB", "width"])

    counts_by_width = overlap.groupby(["chrom", "AB", "width"]).count()
    counts_by_width = pd.merge(
        counts_by_width,
        zero_overlap_AB_width,
        how="outer",
        left_index=True,
        right_index=True,
    )
    counts_by_width = counts_by_width.fillna(0)
    counts_by_width["gc"] = counts_by_width["gc"].astype(int)
    counts_by_width_wide = counts_by_width.reset_index().pivot(
        index=["chrom", "AB"], columns="width", values="gc"
    )
    counts_by_width_wide.to_csv(
        f"{output_dir}/counts_by_width.csv", index=True, header=False
    )

    zero_bin_gc = pd.DataFrame(
        [(chr, i, -1) for chr, n_ab in AB_num.items() for i in range(1, n_ab + 1)],
        columns=["chrom", "AB", "nul"],
    )
    zero_bin_gc = zero_bin_gc.set_index(["chrom", "AB"])

    bin_gc = overlap.groupby(["chrom", "AB"]).mean().drop(columns="width")
    bin_gc = pd.merge(
        bin_gc, zero_bin_gc, how="outer", left_index=True, right_index=True
    )
    bin_gc["gc"].to_csv(
        f"{output_dir}/bingc.csv", index=False, header=False, na_rep="NA"
    )
    dic = {}
    dic["sample_id"] = sample_id
    dic["output_dir"] = output_dir
    dic["RANGE_START"] = RANGE_START
    dic["RANGE_END"] = RANGE_END
    subprocess.check_call(
        "Rscript fragment_analysis.r --sample {sample_id} --output {output_dir} --range_start {RANGE_START} --range_end {RANGE_END} --step 5 --binsize 5 --ab AB.rds".format(
            **dic
        ),
        shell=True,
    )
    print("Analysis complete.")


if __name__ == "__main__":
    import sys

    start_time = time.time()

    BAM_FILE = Path(sys.argv[1])
    BAM_ID = sys.argv[2]
    OUTPUT_DIR = Path(sys.argv[3])
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    N_THREADS = int(sys.argv[4])

    arg_list = [str(i) for i in range(1, 23)]
    with Pool(processes=N_THREADS) as p:
        p.map(frag_stats, arg_list)

    frag_arm(BAM_ID, OUTPUT_DIR)

    print("run time:", time.time() - start_time)
