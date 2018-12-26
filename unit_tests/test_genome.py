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

