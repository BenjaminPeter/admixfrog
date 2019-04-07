from abc import ABC, abstractmethod
from collections import defaultdict, Counter
import itertools
from .utils import AUTOSOMES
from math import ceil
from pybedtools import BedTool


class Blocks(ABC):
    """Abstract Base class to define coordinates for sampling and resampling

    The :class:`Blocks`-class is used to handle

        - which sites are to be read (from bed-file)
        - how sites are grouped into blocks for processing of data. 

    This can be used
    for e.g. calculating a statistic in a sliding window,
    or for implementing a (blocked) jackknive or
    bootstrat resampling.

    Its key method is to provide an iterator over all target coordinates,
    using :meth:`blocks` (or, equivalently :meth:`__iter__`).

    .. testsetup:: *

       from pgdirect.blocks import *

    .. testcode:: 

       blocks = SNPBlocks("data/sample.bed", block_size=2, max_snps=7)
       for block, (chrom, pos) in blocks:
           print("chr%s:%s is in block %s" %(chrom, pos, block))

    .. testoutput:: 

       chr1:567136 is in block 0
       chr1:752565 is in block 0
       chr1:768447 is in block 0
       chr1:776545 is in block 1
       chr1:777121 is in block 1
       chr11:300363 is in block 2
       chr11:306790 is in block 2


    Blocks are integer-indexed and indices are disjoint by chromosome, for easy
    parallel handling.

    For actual use, use one of the subclasses

        - :class:`NoBlocks` : No blocks (assign zero to all SNP)
        - :class:`SNPBlocks` : Block by target coordinates. If bed defines a
            SNP-file, this creates blocks of a certain size
        - :class:`BPBlocks` : Block by genomic coordinate (bp), irrespective of
            whether sites are in the target position
        - :class:`MorganBlocks` : Block by genomic coordinate (cM), irrespective of
            whether sites are in the target position
        - :class:`SingleBlocks` : Each site is its own block
        - :class:`ChromBlocks` : Each chromosome is its own block
    
    Args:
        bed (:obj:`pybedtools.bed`): a bed handle that contains all 
            positions to be processed
        chroms (iterable of :obj:`str`, optional) : The chromosome to be handled.
            Defaults to human autosomes ('1' ... '22')
        starts (:obj:`dict` of int, optional) : 
            The start positions for each chromosome. Defaults to 0
        ends  (:obj:`dict` of int, optional) : 
            The end positions for each chromosome. Defaults to 10^12
        snp_alleles (bool) : 
            if snp alleles are stored in columns 4/5, they can be iterated over.
            for longer fragments this assumes that entry 4, 5 are strings with the ref/alt sequence
            (default True)
        max_snps (int) : Maximum number of SNP to be read
            
    """

    def __init__(
        self,
        file,
        chroms=AUTOSOMES,
        starts=defaultdict(lambda: 0),
        ends=defaultdict(lambda: int(1e12)),
        snp_alleles=True,
        max_snps=None,
    ):
        self.file = file
        self.starts = starts
        self.ends = ends
        self.has_snp_alleles = snp_alleles
        self.max_snps = max_snps

    @abstractmethod
    def block_range(self):
        """generates the block-id for each consecutive SNP

        Yields:
            block_id (int): id of the block the next snp is contained in
        Return("return block id for each SNP")
        """

    def __iter__(self):
        return self.blocks()

    @abstractmethod
    def blocks(self):
        pass


class BedBlocks(Blocks):
    """
    blocks by bed file
    """

    def __init__(self, file, *args, **kwargs):
        super.__init__(file, *args, **kwargs)
        if type(file) is not BedTool:
            self.file = BedTool(file)


    def blocks(self):
        """Joint iterator over all sites and the corresponding block id.
        
        Yields
        -------
        tuple (block, (chrom, pos, ref, alt)) of genomic coordinates and block assigned
        to it
        """
        snp_positions = itertools.chain.from_iterable(self._get_snps_interval())
        block_ids = self.block_range()

        if self.max_snps is None:
            return zip(block_ids, snp_positions)
        else:
            return itertools.islice(zip(block_ids, snp_positions), self.max_snps)


    def _get_snps_interval(self):
        """iterate over positions of a single interval in bed files

        Yields: 
            (chrom, pos) : The chromosome and position of the next site
        """
        for b in self.bed:
            s, e = self.starts[b.chrom], self.ends[b.chrom]
            if s <= b.start <= e or s <= b.end <= e:
                cur_s, cur_e = max(s, b.start), min(e, b.end)
                if self.has_snp_alleles:
                    ref, alt = b.name, b.score
                    pos = zip(itertools.cycle([b.chrom]), range(cur_s, cur_e), ref, alt)
                    yield pos
                else:
                    pos = zip(itertools.cycle([b.chrom]), range(cur_s, cur_e))
                    yield pos
        return

class SNPBlocks(BedBlocks):
    """Genome-wide blocking
    
    .. testsetup:: *

       from pgdirect.blocks import *

    .. testcode:: 

       blocks = SNPBlocks("data/sample.bed", block_size=2, max_snps=7)
       for block, (chrom, pos) in blocks:
           print("chr%s:%s is in block %s" %(chrom, pos, block))

    .. testoutput:: 

       chr1:567136 is in block 0
       chr1:752565 is in block 0
       chr1:768447 is in block 0
       chr1:776545 is in block 1
       chr1:777121 is in block 1
       chr11:300363 is in block 2
       chr11:306790 is in block 2
    
    
    """

    def __init__(self, *args, block_size, **kwargs):
        super().__init__(*args, **kwargs)
        self.block_size = block_size

    def _count_chrom(self):
        """_count_chrom

        counts the number of callable bases per chromosome
        
        Parameters
        ----------
        bed : pybedtools.bed
            bed file to calculate blocks on
        start : int
            absolute start of interval
        end : type
            absolute end of interval
        
        Yields
        -------
        chromosome id for each site
        """

        for b in self.bed:
            if (
                self.starts[b.chrom] <= b.start <= self.ends[b.chrom]
                or self.starts[b.chrom] <= b.end <= self.ends[b.chrom]
            ):
                s = max(self.starts[b.chrom], b.start)
                e = min(self.ends[b.chrom], b.end)
                for i in range(s, e):
                    yield b.chrom

    def _setup(self):
        """_setup
        Counts number of SNP per chromosome, to make evenly-sized windows
        """

        self.snp_per_chrom = Counter(self._count_chrom())
        self.blocks_per_chrom = dict(
            (k, max(1, v // self.block_size)) for (k, v) in self.snp_per_chrom.items()
        )
        self.sizes = dict(
            (k, ceil(self.snp_per_chrom[k] / self.blocks_per_chrom[k]))
            for k in self.snp_per_chrom
        )
        self.block_start = dict()
        b = 0
        for k, v in self.blocks_per_chrom.items():
            self.block_start[k] = b
            b += v

    def block_range(self):
        """block_range

        Yields
        ======
        given current parametrization, returns the block id of the i-th SNP
        i.e. if block-size is 5, then return is 0,0,0,0,0,1,1,1,1,1,2, ...
        """
        if not hasattr(self, "snp_per_chrom"):
            self._setup()
        for chrom in self.chroms:
            for i in range(self.snp_per_chrom[chrom]):
                yield self.block_start[chrom] + (i // self.sizes[chrom])


class BPBlocks(BedBlocks):
    """Genome-wide blocking by base-pairs

    See :class:`Blocks` for details. Note that the specified
    `block_size` is approximate, in that it is slightly increased to avoid
    having "rounded-off"-blocks at the end of the chromosome

    .. testsetup:: *

       from pybedtools import BedTool
       from pgdirect.blocks import *

    .. testcode:: 

       blocks = BPBlocks(BedTool("data/sample.bed"), 
           block_size=100000, max_snps=7)
       for block, (chrom, pos) in blocks:
           print("chr%s:%s is in block %s" %(chrom, pos, block))

    .. testoutput:: 

       chr1:567136 is in block 0
       chr1:752565 is in block 1
       chr1:768447 is in block 1
       chr1:776545 is in block 1
       chr1:777121 is in block 1
       chr11:300363 is in block 2
       chr11:306790 is in block 2

    Args:
        block_size (int) : the number of bases each block is to contain.

    """

    def __init__(self, *args, block_size, **kwargs):
        super().__init__(*args, **kwargs)
        self.block_size = block_size

    def _chrom_range(self):
        """count_chrom

        gets min and max position per chrom
        
        Parameters
        ----------
        bed : pybedtools.bed
            bed file to calculate blocks on
        start : int
            absolute start of interval
        end : type
            absolute end of interval
        
        Yields
        -------
        chromosome id for each site
        """

        min_ = defaultdict(lambda: int(1e12))
        max_ = defaultdict(lambda: 0)

        for b in self.bed:
            if b.start < min_[b.chrom]:
                min_[b.chrom] = b.start
            if b.end > max_[b.chrom]:
                max_[b.chrom] = b.end

        range_ = dict((i, max_[i] - min_[i]) for i in min_)
        return min_, max_, range_

    def _setup(self):
        """_setup
        Counts number of SNP per chromosome, to make evenly-sized windows
        """
        self.min, self.max, self.range = self._chrom_range()
        self.blocks_per_chrom = dict(
            (k, max(1, v // self.block_size)) for (k, v) in self.range.items()
        )

        self.sizes = dict(
            (k, ceil(self.range[k] / self.blocks_per_chrom[k])) for k in self.range
        )
        self.block_start = dict()
        b = 0
        for k, v in self.blocks_per_chrom.items():
            self.block_start[k] = b
            b += v

    def block_range(self):
        """block_range

        Yields
        ======
        given current parametrization, returns the block id of the i-th SNP
        i.e. if block-size is 5, then return is 0,0,0,0,0,1,1,1,1,1,2, ...
        """
        if not hasattr(self, "blocks_per_chrom"):
            self._setup()
        for b in self.bed:
            for pos in range(b.start, b.end):
                offset = (pos - self.min[b.chrom]) // self.sizes[b.chrom]
                yield int(self.block_start[b.chrom] + offset)


class NoBlocks(BedBlocks):
    """All SNP are assigned to one block

    See :class:`Blocks` for details

    .. testsetup:: *

       from pybedtools import BedTool
       from pgdirect.blocks import *

    .. testcode:: 

       blocks = NoBlocks(BedTool("data/sample.bed"), max_snps=7)
       for block, (chrom, pos) in blocks:
           print("chr%s:%s is in block %s" %(chrom, pos, block))

    .. testoutput:: 

       chr1:567136 is in block 0
       chr1:752565 is in block 0
       chr1:768447 is in block 0
       chr1:776545 is in block 0
       chr1:777121 is in block 0
       chr11:300363 is in block 0
       chr11:306790 is in block 0

    """

    def block_range(self):
        while True:
            yield 0


class SingleBlocks(BedBlocks):
    """Each SNP is its own block

    .. testsetup:: *

       from pybedtools import BedTool
       from pgdirect.blocks import *

    .. testcode:: 

       blocks = SingleBlocks(BedTool("data/sample.bed"), max_snps=7)
       for block, (chrom, pos) in blocks:
           print("chr%s:%s is in block %s" %(chrom, pos, block))

    .. testoutput:: 

       chr1:567136 is in block 0
       chr1:752565 is in block 1
       chr1:768447 is in block 2
       chr1:776545 is in block 3
       chr1:777121 is in block 4
       chr11:300363 is in block 5
       chr11:306790 is in block 6

    See :class:`Blocks` for details
    """

    def block_range(self):
        block = 0
        while True:
            yield block
            block += 1


class ChromBlocks(BedBlocks):
    """each chromosome is its own block

    See :class:`Blocks` for details


    .. testsetup:: *

       from pgdirect.blocks import *

    .. testcode:: 

       blocks = ChromBlocks("data/sample.bed", max_snps=7)
       for block, (chrom, pos) in blocks:
           print("chr%s:%s is in block %s" %(chrom, pos, block))

    .. testoutput:: 

       chr1:567136 is in block 1
       chr1:752565 is in block 1
       chr1:768447 is in block 1
       chr1:776545 is in block 1
       chr1:777121 is in block 1
       chr11:300363 is in block 11
       chr11:306790 is in block 11
        
    """

    def block_range(self):
        for b in self.bed:
            for _ in range(b.start, b.end):
                yield b.chrom


class MorganBlocks(BedBlocks):
    """Not yet implemented"""

    def __init__(self):
        raise NotImplemented


__all__ = [
    "Blocks",
    "BPBlocks",
    "SNPBlocks",
    "NoBlocks",
    "SingleBlocks",
    "ChromBlocks",
    "MorganBlocks",
]
