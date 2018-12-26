import re
from collections import namedtuple
import doctest


def change_chromname(chrom:str) -> str:
    if chrom.startswith("chr"):
        return chrom.replace("chr", "")
    else:
        return "chr" + chrom


GenomeRange_ = namedtuple("GenomeRange", ["chrom", "start", "end"])

class GenomeRange(GenomeRange_):
    """
    Objects for represent a genome range, a genome position or a chromosome.

    Attributes
    ----------
    chrom : str
        chromosome name
    start : {int, None}
        Range start position, 0 based.
    end : {int, None}
        Range end position, 1 based.
    """

    def __str__(self) -> str:
        if (self.start is not None) and (self.end is not None):
            return f"{self.chrom}:{self.start}-{self.end}"
        elif (self.start is not None) and (self.end is None):
            return f"{self.chrom}:{self.start}"
        else:
            return self.chrom

    def __repr__(self) -> str:
        return f"GenomeRange({self.chrom}, {self.start}, {self.end})"

    def change_chromname(self) -> 'GenomeRange':
        chrom_ = change_chromname(self.chrom)
        return GenomeRange(chrom_, self.start, self.end)

    @property
    def length(self) -> int:
        return self.end - self.start

    def __contains__(self, another:'GenomeRange') -> bool:
        if another.chrom != self.chrom:
            return False
        if another.start < self.start:
            return False
        if another.end > self.end:
            return False
        return True

    @staticmethod
    def from_str(region:str) -> 'GenomeRange':
        if '-' in region:
            chr_, s, e = re.split("[:-]", region)[:3]
            s, e = int(s), int(e)
            grange = GenomeRange(chr_, s, e)
        else:
            chr_ = region
            grange = GenomeRange(chr_, None, None)
        return grange

    @staticmethod
    def check(grange:'GenomeRange') -> bool:
        """
        Check a GenomeRange object is valid or not.
        """
        chrom, start, end = grange
        msg = f"GenomeRange object: {repr(grange)} is not valid. "
        if not isinstance(chrom, str):
            raise ValueError(msg + f"chrom expect instance of str, get {type(chrom)}")

        if start is None:
            if end is not None:
                raise ValueError(msg + "if start is None, end must also None.")
        else:  # start is not None
            if end is None:
                assert isinstance(start, int), msg + f"start expect instance of int, get {type(start)}"
            else:
                assert isinstance(start, int), msg + f"start expect instance of int, get {type(start)}"
                assert isinstance(end, int), msg + f"end expect instance of int, get {type(end)}"

                if start > end:
                    raise ValueError(msg + "end must >= start.")

        return True


class GenomeBinRange(GenomeRange):
    pass


def genome_range(*args) -> GenomeRange:
    """
    A convenient function for construct the GenomeRange object.

    >>> genome_range("chr1:1000-2000")
    GenomeRange(chr1, 1000, 2000)
    >>> genome_range("chr1")
    GenomeRange(chr1, None, None)
    >>> genome_range("chr1", 1000, 2000)
    GenomeRange(chr1, 1000, 2000)
    >>> genome_range("chr1", 1000)
    GenomeRange(chr1, 1000, None)
    """
    if len(args) == 1:
        region_str = args[0]
        result = GenomeRange.from_str(region_str)
    else:
        if len(args) == 2:
            chrom, start = args
            end = None
        else:
            chrom, start, end = args[:3]
        result = GenomeRange(chrom, start, end)
    GenomeRange.check(result)
    return result


if __name__ == "__main__":
    doctest.testmod()
