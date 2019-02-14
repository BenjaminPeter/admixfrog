from enum import IntEnum
from collections import namedtuple

AUTOSOMES = [str(i) for i in range(1, 23)]

rev_complement = {"A": "T", "C": "G", "G": "C", "T": "A"}

unambiguous = {
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "N": "N",
}

DEAMINATION = IntEnum("DEAMINATION", ("NONE", "THREE_BASES", "FIRST_BASE", "DOUBLE"))

SAMPLING = IntEnum(
    "SAMPLING",
    ("RANDOM", "RANDOM_RG", "RANDOM_RGDEAM", "RANDOM_DEAM", "DIPLOID", "HAPLOID"),
)

NA = "N"

Read = namedtuple("Read", ("RG", "deam", "base", "len", "bq", "is_reverse", "pos_in_read", "deam_full"))
Read.__new__.__defaults__ = ("", (-1, -1), NA, -1, 0, False, -1, [])
NA_READ = Read("", (-1, -1), NA, -1, 0, False, -1, [])

__all__ = ["NA", "NA_READ", "SAMPLING", "DEAMINATION", "Read"]
