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