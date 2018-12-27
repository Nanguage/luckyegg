import os
from os.path import join
import pytest

from luckyegg.genome import *


def test_genome_range():
    gr = genome_range("chr1:1000-2000")
    assert (gr.chrom, gr.start, gr.end) == ("chr1", 1000, 2000)
    gr = genome_range("chr1:1000")
    assert (gr.chrom, gr.start, gr.end) == ("chr1", 1000, None)
    gr = genome_range("chr1", 1000)
    assert (gr.chrom, gr.start, gr.end) == ("chr1", 1000, None)
    gr = genome_range("chr1", 1000, 2000)
    assert (gr.chrom, gr.start, gr.end) == ("chr1", 1000, 2000)

    with pytest.raises(AssertionError) as excinfo:
        gr = genome_range(11, 1000, 2000)
    assert "chrom expect" in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        gr = genome_range("chr1", None, 1)
    assert "if start is None, end must also None." in str(excinfo.value)

    with pytest.raises(AssertionError) as excinfo:
        gr = genome_range("chr1", "111")
    assert "start expect instance of int" in str(excinfo.value)

    with pytest.raises(AssertionError) as excinfo:
        gr = genome_range("chr1", "111", 111)
    assert "start expect instance of int" in str(excinfo.value)

    with pytest.raises(AssertionError) as excinfo:
        gr = genome_range("chr1", 111, "111")
    assert "end expect instance of int" in str(excinfo.value)


def test_GenomeRange_in():
    gr1 = GenomeRange("chr1", 1000, 2000)
    gr2 = GenomeRange("chr1", 0, 3000)
    gr3 = GenomeRange("chr2", 0, 3000)
    assert gr1 in gr2
    assert gr2 not in gr1
    assert gr1 not in gr3
    assert "chr1" in gr1
    assert 1000 in gr1
    assert 2000 in gr1
    assert 0 not in gr1
    assert genome_range("chr1:1000") in gr1
    assert genome_range("chr1:2000") not in gr1
    assert genome_range("chr1:1000") in genome_range("chr1:1000")
    assert genome_range("chr1:1000") not in genome_range("chr1:2000")
    assert genome_range("chr1") not in genome_range("chr1:2000")
    assert genome_range("chr1") not in gr1
    assert gr1 in genome_range("chr1")


def test_GenomeRange_length():
    gr1 = genome_range("chr1:1000-2000")
    print(len(gr1))
    assert len(gr1) == 1000
    gr2 = genome_range("chr1:1000")
    assert len(gr2) == 1
    gr3 = genome_range("chr1")
    with pytest.raises(ValueError):
        len(gr3)


def test_GenomeRange_convert():
    gr1 = GenomeRange("chr1", 0, 1001)
    gbr1 = gr1.to_bin(binsize=1000)
    assert isinstance(gbr1, GenomeBinRange)
    assert (gbr1.start, gbr1.end) == (0, 2)
    gr2 = gbr1.to_bp(binsize=1000)
    assert not isinstance(gr2, GenomeBinRange)
    assert (gr2.start, gr2.end) == (0, 2000)
    gr3 = genome_range("chr1:1000")
    gbr3 = gr3.to_bin(binsize=1000)
    assert isinstance(gbr3, GenomeBinRange)
    assert (gbr3.start, gbr3.end) == (1, None)
    gr4 = genome_range("chr1")
    gbr4 = gr4.to_bin(binsize=1000)
    assert isinstance(gbr4, GenomeBinRange)
    assert (gbr4.start, gbr4.end) == (None, None)
    assert genome_range("chr1:1500").to_bin(binsize=1000).start == 1
    assert genome_range("chr1:1000").to_bin(binsize=1000).start == 1
    assert genome_range("chr1:999").to_bin(binsize=1000).start == 0
    assert genome_range("chr1:1000-1002").to_bin(binsize=1000).end == 2
    assert genome_range("chr1:100-800").to_bin(binsize=1000).end == 1
    assert genome_range("chr1:100-1000").to_bin(binsize=1000).end == 1


def gen_example_chromsizes_file(tmp_dir="/tmp"):
    file_ = join(tmp_dir, "examl")
    with open(file_, 'w') as f:
        f.write("chr1\t10000\n")
        f.write("chr2\t20000\n")
        f.write("chr3\t15000\n")
        f.write("chr4\t2000\n")
    return file_


def test_ChromSizes_from_file():
    file_ = gen_example_chromsizes_file()
    chrsizes = ChromSizes.from_file(file_)
    assert chrsizes.unit == "bp"
    assert "chr1" in chrsizes.sizes
    assert "chr1" in chrsizes
    assert isinstance(chrsizes["chr1"], int)
    with pytest.raises(TypeError):
        1 in chrsizes
    os.remove(file_)


def test_ChromSizes_contains():
    chrsizes = ChromSizes({
        "chr1": 10000,
        "chr2": 200000,
    })
    assert "chr1" in chrsizes
    assert "chr3" not in chrsizes
    assert genome_range("chr1") in chrsizes
    assert genome_range("chr3") not in chrsizes
    assert genome_range("chr1:1000") in chrsizes
    assert genome_range("chr1:3000000") not in chrsizes
    assert genome_range("chr3:1000-2000") not in chrsizes
    assert genome_range("chr1:100-200") in chrsizes
    

def test_ChromSizes_convert():
    chrsizes = ChromSizes({
        "chr1": 15000,
        "chr2": 20000,
    })
    chrsizes_bin = chrsizes.to_bin(binsize=10000)
    assert chrsizes_bin.unit == "bin"
    assert chrsizes_bin["chr1"] == 2
    assert chrsizes_bin["chr2"] == 2

