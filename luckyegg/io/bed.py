from collections import namedtuple
from typing import Iterable, Union, Type, List, NewType, NamedTuple, TypeVar

from luckyegg.genome import GenomeRange

class BED6_(NamedTuple):
    chrom: str
    start: int
    end: int
    name: str
    score: Union[int, float, str]
    strand: str

class BED9_(NamedTuple):
    chrom: str
    start: int
    end: int
    name: str
    score: Union[int, float, str]
    strand: str
    thickStart: str
    thickEnd: str
    itemRGB: str

class BED12_(NamedTuple):
    chrom: str
    start: int
    end: int
    name: str
    score: Union[int, float, str]
    strand: str
    thickStart: str
    thickEnd: str
    itemRGB: str
    blockCount: str
    blockSizes: str
    blockStarts: str

class BEDGraph_(NamedTuple):
    chrom: str
    start: int
    end: int
    value: Union[int, float, str]


class BEDLike(object):
    @classmethod
    def from_line(cls, line):
        """
        compose BED record from a line.
        """
        items = line.split()
        items[1] = int(items[1])  # cast start and end to int
        items[2] = int(items[2])
        return cls(*items)

    @property
    def genome_range(self):
        return GenomeRange(self.chrom, self.start, self.end)

    def __str__(self):
        return "\t".join([str(i) for i in self])


class Bed6(BEDLike, BED6_):
    pass

class Bed9(BEDLike, BED9_):
    pass

class Bed12(BEDLike, BED12_):
    pass

class BedGraph(BEDLike, BEDGraph_):
    pass


def read_bed(path: str) -> Iterable[BEDLike]:
    bed_type = infer_bed_type(path)
    header_rows = infer_header_rows(path)
    with open(path) as f:
        [f.readline() for _ in range(header_rows)]
        for line in f:
            line = line.strip()
            yield bed_type.from_line(line)


def infer_bed_type(path:str) -> Type[BEDLike]:
    with open(path) as f:
        while True:
            line = f.readline()
            if not is_header(line):
                break
        else:
            raise IOError(f"Bed-like file {path} don't have enough content.")
    items = line.strip().split()
    if len(items) == 4:
        return BedGraph
    elif len(items) == 6:
        return Bed6
    elif len(items) == 9:
        return Bed9
    else:
        return Bed12


def infer_header_rows(path:str) -> int:
    with open(path) as f:
        i = 0
        while True:
            line = f.readline()
            if not is_header(line):
                break
            i += 1
    return i


def is_header(line) -> bool:
    if line.startswith("#") or\
       line.startswith("header") or\
       line.startswith("track") or\
       line.startswith("browser"):
        return True
    else:
        return False
