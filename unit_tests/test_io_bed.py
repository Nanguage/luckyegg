import os
from os.path import join
import random

from luckyegg.io.bed import *


def create_sample(name:str,
        bed_type:BEDLike,
        lines:int=10,
        header:bool=True,
        sample_dir:str="/tmp") -> str:
    path = join(sample_dir, name)
    with open(path, 'w') as f:
        if header:
            # write 3 example header
            f.write("#header\n")
            f.write("track name=pairedReads description=\"Clone Paired Reads\" useScore=1\n")
            f.write("browser position chr7:127471196-127495720")
        for _ in range(lines):
            line = example_bed_line(bed_type)
            f.write(line+"\n")
    return path


def example_bed_line(bed_type:BEDLike) -> str:
    def fmt_line(fields):
        return "\t".join([str(i) for i in fields])
    CHROMS = ['chr'+str(i) for i in range(1, 23)]
    MAX_LENGTH = 100000
    INTVAL_LEN = 100
    MAX_VAL = 100
    chrom = random.choices(CHROMS)[0]
    start = random.randint(1, MAX_LENGTH - INTVAL_LEN)
    end = random.randint(INTVAL_LEN, MAX_LENGTH)
    # bed graph
    if bed_type == BedGraph:
        val = random.randint(0, MAX_VAL)
        return fmt_line([chrom, start, end, val])
    name = "."
    strand = random.choice(['.', '+', '-'])
    score = random.randint(0, MAX_VAL)
    # bed6
    if bed_type == Bed6:
        return fmt_line([chrom, start, end, name, score, strand])
    thick_start = random.randint(1, MAX_LENGTH - INTVAL_LEN)
    thick_end = random.randint(INTVAL_LEN, MAX_LENGTH)
    item_rgb = ",".join([str(random.randint(0, 255)) for _ in range(3)])
    # bed9
    if bed_type == Bed9:
        return fmt_line([chrom, start, end, name, score, strand, thick_start, thick_end, item_rgb])
    blockc = 2
    block_sizes = "567,488,"
    block_starts = "0,3512"
    # bed12
    if bed_type == Bed12:
        return fmt_line([chrom, start, end, name, score, strand, thick_start, thick_end, item_rgb, blockc, block_sizes, block_starts])


def test_infer_bed_type():
    bedg = create_sample('example_bedgraph', BedGraph)
    assert infer_bed_type(bedg) == BedGraph
    os.remove(bedg)
    bed6 = create_sample('example_bed6', Bed6)
    assert infer_bed_type(bed6) == Bed6
    os.remove(bed6)
    bed9 = create_sample('example_bed9', Bed9)
    assert infer_bed_type(bed9) == Bed9
    os.remove(bed9)
    bed12 = create_sample('example_bed12', Bed12)
    assert infer_bed_type(bed12) == Bed12
    os.remove(bed12)


def test_infer_header_rows():
    bed = create_sample('example_bed', Bed6, header=True)
    assert infer_header_rows(bed) == 3
    bed = create_sample('example_bed', Bed6, header=False)
    assert infer_header_rows(bed) == 0
    os.remove(bed)


def test_read_bed():

    def test_(bed_type):
        bed_path = create_sample('example_bed', bed_type)
        header_rows = infer_header_rows(bed_path)
        with open(bed_path) as f:
            header_rows = [f.readline() for _ in range(header_rows)]
            for bed, line in zip(read_bed(bed_path), f):
                assert isinstance(bed, bed_type)
                assert len(bed) == len(bed_type.fields)
                assert str(bed) == line.strip()
                items = line.strip().split()
                genome_range_str = f"{items[0]}:{items[1]}-{items[2]}"
                assert str(bed.genome_range) == genome_range_str
        os.remove(bed_path)

    test_(BedGraph)
    test_(Bed6)
    test_(Bed9)
    test_(Bed12)

