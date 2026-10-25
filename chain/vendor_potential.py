"""The plateau of the chain.  SumPotential is a thin re-export of the hierarchical plateau class of cadena_cierre
(components/cadena_cierre/common.py, Potential), so that every stage uses literally the same code; PolePotential is the order-4
pole of U22 (chain/pole_potential.py).  Potential is the one chosen by WHALEFALL_PLATEAU (ESPEC_U26; default 'pole')."""
from ucommon import PLATEAU
from common import Potential as SumPotential  # components/cadena_cierre/common.py (chain/ucommon.py puts it on sys.path)
from pole_potential import PolePotential

Potential = PolePotential if PLATEAU == "pole" else SumPotential
