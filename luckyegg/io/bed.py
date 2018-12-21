from collections import namedtuple
from typing import Iterable, Union, Type, List

from luckyegg.genome import GenomeRange

BED_BASIC_FIELDS = ['chrom', 'start', 'end']
BED6_FIELDS = BED_BASIC_FIELDS + ['name', 'score', 'strand']
BED9_FIELDS = BED6_FIELDS + ['thickStart', 'thichEnd', 'itemRGB']
BED12_FIELDS = BED9_FIELDS + ['blockCount', 'blockSizes', 'blockStarts']
BED_GRAPH_FIELDS = BED_BASIC_FIELDS + ['value']

def bed_like(name:str, fields:List[str]) -> Type:
    base_class = namedtuple(name+"_", fields)

    def from_line(cls, line):
        """
        compose BED record from a line.
        """
        items = line.split()
        items[1] = int(items[1])  # cast start and end to int
        items[2] = int(items[2])
        return cls(*items)

    class_ = type(name, (base_class,), {
        'fields': fields,
        'genome_range': property(
            lambda self: GenomeRange(self.chrom, self.start, self.end)),
        '__str__': lambda self: "\t".join([str(i) for i in self]),
        'from_line': classmethod(from_line),
    })
    return class_


Bed6 = bed_like("Bed6", BED6_FIELDS)
Bed9 = bed_like("Bed9", BED9_FIELDS)
Bed12 = bed_like("Bed12", BED12_FIELDS)
BedGraph = bed_like("BedGraph", BED_GRAPH_FIELDS)

BEDLike = Union[Bed6, Bed9, Bed12, BedGraph]


def read_bed(path: str) -> Iterable[BEDLike]:
    bed_type = infer_bed_type(path)
    header_rows = infer_header_rows(path)
    with open(path) as f:
        header_rows = [f.readline() for _ in range(header_rows)]
        for line in f:
            line = line.strip()
            items = line.split()
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
