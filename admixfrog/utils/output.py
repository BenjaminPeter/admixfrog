import logging
import yaml
from numba import njit
from collections import Counter
from scipy.stats import binom
import pandas as pd
import numpy as np
import itertools

from .utils import posterior_table, posterior_table_slug


# yaml writer with indent
class IndentDumper(yaml.Dumper):
    def increase_indent(self, flow=False, indentless=False):
        return super(IndentDumper, self).increase_indent(flow, False)


def write_pars_table(pars, outname=None):
    try:
        P = pars.__dict__.copy()
        for k, v in P.items():
            try:
                P[k] = v.tolist()
            except:
                pass

    except AttributeError:
        P = pars
    # s = yaml.safe_dump(P, default_flow_style=False, indent=4)
    s = yaml.dump(P, Dumper=IndentDumper, default_flow_style=False, indent=4)

    if outname is not None:
        with open(outname, "wt") as f:
            f.write(s)

    return s
