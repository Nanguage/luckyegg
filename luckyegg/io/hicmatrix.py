from luckyegg.genome import GenomeRange

from cooler.api import Cooler


class MatrixSelector(object):
    """
    Selector for fetch the matrix from cool file.

    Parameters
    ----------
    cool : `cooler.api.Cooler`
        cool object.
    balance : bool
        balance matrix or not.
    """
    def __init__(self, cool:Cooler, balance:bool=True) -> None:
        self.cool = cool
        self.balance = balance

    @property
    def chromsizes(self):
        return self.cool.chromsizes.to_dict()

    @property
    def binsize(self):
        return self.cool.binsize

