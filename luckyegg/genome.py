import re
from collections import namedtuple
from typing import Union, Any
import doctest


def change_chromname(chrom:str) -> str:
    if chrom.startswith("chr"):
        return chrom.replace("chr", "")
    else:
        return "chr" + chrom


GenomeRange_ = namedtuple("GenomeRange", ["chrom", "start", "end"])

class GenomeRange(GenomeRange_):
    """
    Object for represent a genome range, a genome position or a chromosome.

    Attributes
    ----------
    chrom : str
        chromosome name.
    start : {int, None}
        Range start position, 0 based.
    end : {int, None}
        Range end position, 1 based.
    range_type : {'chromosome', 'range', 'point'}
        Type of the GenomeRange.
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

    @property
    def range_type(self) -> str:
        """
        Genome range type:

        >>> GenomeRange("chr1", None, None).range_type
        'chromosome'
        >>> GenomeRange("chr1", 1000, None).range_type
        'point'
        >>> GenomeRange("chr1", 0, 1000).range_type
        'range'
        """
        if self.start is None and self.end is None:
            return "chromosome"
        elif self.start is not None and self.end is None:
            return "point"
        else:
            return "range"

    def change_chromname(self) -> 'GenomeRange':
        """
        Change chromosome name style.

        >>> GenomeRange("chr1", 1000, 2000).change_chromname()
        GenomeRange(1, 1000, 2000)
        >>> GenomeRange("1", 1000, 2000).change_chromname()
        GenomeRange(chr1, 1000, 2000)
        """
        chrom_ = change_chromname(self.chrom)
        return GenomeRange(chrom_, self.start, self.end)

    def __len__(self) -> int:
        if self.end is not None:
            return abs(self.end - self.start)
        else:
            if self.start is not None:
                return 1
            else:
                raise ValueError("If start and end both None, can not compute length")

    def __contains__(self, another:Union['GenomeRange', Any]) -> bool:
        if isinstance(another, GenomeRange):
            if another.chrom != self.chrom:
                return False

            if self.range_type == "range" and another.range_type == "range":
                if another.start < self.start:
                    return False
                if another.end > self.end:
                    return False
                return True
            elif self.range_type == "range" and another.range_type == "point":
                return self.start <= another.start < self.end
            elif self.range_type == "chromosome" and (another.range_type == "point" or another.range_type == "range"):
                return True
            elif self.range_type == "point" and another.range_type == "point":
                return self == another
            else:
                return False

        else:
            return another in tuple(self)

    @staticmethod
    def from_str(region:str) -> 'GenomeRange':
        if '-' in region:
            chr_, s_, e_ = re.split("[:-]", region)[:3]
            s, e = int(s_), int(e_)
            grange = GenomeRange(chr_, s, e)
        elif ':' in region:
            chr_, s_ = region.split(":")
            s = int(s_)
            grange = GenomeRange(chr_, s, None)
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
        assert isinstance(chrom, str), msg + f"chrom expect instance of str, get {type(chrom)}"

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

    def to_bin(self, binsize:int) -> 'GenomeBinRange':
        """
        Convert to GenomeBinRange, unit in 'bin'.

        >>> GenomeRange("chr1", 0, 1001).to_bin(binsize=1000)
        GenomeBinRange(chr1, 0, 1)
        """
        if self.start is None:
            return GenomeBinRange(self.chrom, self.start, self.end)
        elif self.end is None:
            return GenomeBinRange(self.chrom, self.start//binsize, self.end)
        else:
            return GenomeBinRange(self.chrom, self.start//binsize, self.end//binsize)


class GenomeBinRange(GenomeRange):
    """
    Similar to GenomeRange, but the unit is the number of 'bin'.
    """

    def __repr__(self) -> str:
        return f"GenomeBinRange({self.chrom}, {self.start}, {self.end})"

    def to_bp(self, binsize:int) -> GenomeRange:
        """
        Convert to GenomeRange, unit in 'bp'

        >>> GenomeBinRange("chr1", 0, 1).to_bp(1000)
        GenomeRange(chr1, 0, 1000)
        """
        if self.start is None:
            return GenomeRange(self.chrom, self.start, self.end)
        elif self.end is None:
            return GenomeRange(self.chrom, self.start*binsize, self.end)
        else:
            return GenomeRange(self.chrom, self.start*binsize, self.end*binsize)


def genome_range(*args) -> GenomeRange:
    """
    A convenient function for construct a valid GenomeRange object.

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


class ChromSizes(object):
    """
    Object for represent the Chromosomes's length of a Genome.

    Attributes
    ----------
    sizes : dict
        A dict store Chromosome's name to it's size.
    unit : {'bp', 'bin'}
        The unit of chromosome size.
    """
    def __init__(self, chromsizes:dict, unit:str='bp') -> None:
        self.sizes = chromsizes
        self.unit = unit

    def to_bin(self, binsize) -> 'ChromSizes':
        """
        Return a ChromSize object unit in 'bin'
        """
        if self.unit != 'bin':
            sizes = {}
            for key, value in self.sizes.items():
                sizes[key] = value // binsize
            return ChromSizes(sizes, "bin")
        else:
            return self

    def contain_range(self, grange:GenomeRange) -> bool:
        """
        Judge a GenomeRange within the genome or not.
        """
        if self.unit == "bin" and not isinstance(grange, GenomeBinRange):
            return False
        if self.unit != "bin" and isinstance(grange, GenomeBinRange):
            return False
        if grange.chrom not in self.sizes:
            return False
        chr_len = self.sizes[grange.chrom]
        chrom_range = GenomeRange(grange.chrom, 0, chr_len)
        return grange in chrom_range

    def __contains__(self, another:Union[GenomeRange, str]) -> bool:
        if isinstance(another, GenomeRange):
            if another.range_type == "chromosome":
                return another.chrom in self
            else:
                return self.contain_range(another)
        elif isinstance(another, str):
            return another in self.sizes
        else:
            raise TypeError("ChromSizes can only contains GenomeRange or str object.")

    def __getitem__(self, key:str):
        return self.sizes[key]

    @staticmethod
    def from_file(path:str) -> 'ChromSizes':
        """
        From chromosome length file.
        """
        sizes = {}
        with open(path) as f:
            for line in f:
                items = line.strip().split()
                sizes[items[0]] = int(items[1])
        return ChromSizes(sizes, unit="bp")


if __name__ == "__main__":
    doctest.testmod()
