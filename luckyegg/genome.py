import re
from collections import namedtuple


def change_chromname(chrom:str) -> str:
    if chrom.startswith("chr"):
        return chrom.replace("chr", "")
    else:
        return "chr" + chrom


GenomeRange_ = namedtuple("GenomeRange", ["chrom", "start", "end"])

class GenomeRange(GenomeRange_):
    def __str__(self) -> str:
        if (self.start is not None) and (self.end is not None):
            return "{}:{}-{}".format(self.chrom, self.start, self.end)
        else:
            return self.chrom

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


def region_str2genome_range(region:str) -> GenomeRange:
    if '-' in region:
        chr_, s, e = re.split("[:-]", region)[:3]
        grange = GenomeRange(chr_, s, e)
    else:
        chr_ = region
        grange = GenomeRange(chr_, None, None)
    return grange
