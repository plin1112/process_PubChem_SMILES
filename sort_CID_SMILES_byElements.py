import pandas as pd
from collections import defaultdict
import re

# modified from pysmiles:
# https://github.com/pckroon/pysmiles
import enum

ISOTOPE_PATTERN = r'(?P<isotope>[\d]+)?'
ELEMENT_PATTERN = r'(?P<element>b|c|n|o|s|p|\*|[A-Z][a-z]{0,2})'
STEREO_PATTERN = r'(?P<stereo>@|@@|@TH[1-2]|@AL[1-2]|@SP[1-3]|@OH[\d]{1,2}|'\
                  r'@TB[\d]{1,2})?'
HCOUNT_PATTERN = r'(?P<hcount>H[\d]?)?'
CHARGE_PATTERN = r'(?P<charge>(-|\+)(\++|-+|[\d]{1,2})?)?'
CLASS_PATTERN = r'(?::(?P<class>[\d]+))?'
ATOM_PATTERN = re.compile(r'^\[' + ISOTOPE_PATTERN + ELEMENT_PATTERN +
                          STEREO_PATTERN + HCOUNT_PATTERN + CHARGE_PATTERN +
                          CLASS_PATTERN + r'\]$')

VALENCES = {"B": (3,), "C": (4,), "N": (3, 5), "O": (2,), "P": (3, 5),
            "S": (2, 4, 6), "F": (1,), "Cl": (1,), "Br": (1,), "I": (1,)}

AROMATIC_ATOMS = "B C N O P S Se As *".split()

class TokenType(enum.Enum):
    """Possible SMILES token types"""
    ATOM = 1
    BOND_TYPE = 2
    BRANCH_START = 3
    BRANCH_END = 4
    RING_NUM = 5
    EZSTEREO = 6


def _tokenize(smiles):
    """
    Iterates over a SMILES string, yielding tokens.

    Parameters
    ----------
    smiles : iterable
        The SMILES string to iterate over

    Yields
    ------
    tuple(TokenType, str)
        A tuple describing the type of token and the associated data
    """
    organic_subset = 'B C N O P S F Cl Br I * b c n o s p'.split()
    smiles = iter(smiles)
    token = ''
    peek = None
    while True:
        char = peek if peek else next(smiles, '')
        peek = None
        if not char:
            break
        if char == '[':
            token = char
            for char in smiles:
                token += char
                if char == ']':
                    break
            yield TokenType.ATOM, token
        elif char in organic_subset:
            peek = next(smiles, '')
            if char + peek in organic_subset:
                yield TokenType.ATOM, char + peek
                peek = None
            else:
                yield TokenType.ATOM, char
        elif char in '-=#$:.':
            yield TokenType.BOND_TYPE, char
        elif char == '(':
            yield TokenType.BRANCH_START, '('
        elif char == ')':
            yield TokenType.BRANCH_END, ')'
        elif char == '%':
            # If smiles is too short this will raise a ValueError, which is
            # (slightly) prettier than a StopIteration.
            yield TokenType.RING_NUM, int(next(smiles, '') + next(smiles, ''))
        elif char in '/\\':
            yield TokenType.EZSTEREO, char
        elif char.isdigit():
            yield TokenType.RING_NUM, int(char)

def parse_charge(charge_str):
    """
    Parses a SMILES charge specification.

    Parameters
    ----------
    charge_str : str
        The charge specification to parse.

    Returns
    -------
    int
        The charge.
    """
    if not charge_str:
        return 0
    signs = {'-': -1, '+': 1}
    sign = signs[charge_str[0]]
    if len(charge_str) > 1 and charge_str[1].isdigit():
        charge = sign * int(charge_str[1:])
    else:
        charge = sign * charge_str.count(charge_str[0])
    return charge

def parse_hcount(hcount_str):
    """
    Parses a SMILES hydrogen count specifications.

    Parameters
    ----------
    hcount_str : str
        The hydrogen count specification to parse.

    Returns
    -------
    int
        The number of hydrogens specified.
    """
    if not hcount_str:
        return 0
    if hcount_str == 'H':
        return 1
    return int(hcount_str[1:])

def parse_atom(atom):
    """
    Parses a SMILES atom token, and returns a dict with the information.

    Note
    ----
    Can not deal with stereochemical information yet. This gets discarded.

    Parameters
    ----------
    atom : str
        The atom string to interpret. Looks something like one of the
        following: "C", "c", "[13CH3-1:2]"

    Returns
    -------
    dict
        A dictionary containing at least 'element', 'aromatic', and 'charge'. If
        present, will also contain 'hcount', 'isotope', and 'class'.
    """
    defaults = {'charge': 0, 'hcount': 0, 'aromatic': False}
    if not atom.startswith('[') and not atom.endswith(']'):
        if atom != '*':
            # Don't specify hcount to signal we don't actually know anything
            # about it
            return {'element': atom.capitalize(), 'charge': 0,
                    'aromatic': atom.islower()}
        else:
            return defaults.copy()
    match = ATOM_PATTERN.match(atom)
    if match is None:
        raise ValueError('The atom {} is malformatted'.format(atom))
    out = defaults.copy()
    out.update({k: v for k, v in match.groupdict().items() if v is not None})

    if out.get('element', 'X').islower():
        out['aromatic'] = True

    parse_helpers = {
        'isotope': int,
        'element': str.capitalize,
        'stereo': lambda x: x,
        'hcount': parse_hcount,
        'charge': parse_charge,
        'class': int,
        'aromatic': lambda x: x,
    }

    for attr, val_str in out.items():
        out[attr] = parse_helpers[attr](val_str)

    if out['element'] == '*':
        del out['element']

    if out.get('element') == 'H' and out.get('hcount', 0):
        pass
        # raise ValueError("A hydrogen atom can't have hydrogens")

    if 'stereo' in out:
        pass
        # LOGGER.warning('Atom "%s" contains stereochemical information that will be discarded.', atom)

    return out


def read_smiles(smiles, explicit_hydrogen=False, zero_order_bonds=False,
                reinterpret_aromatic=False):
    """
    Parses a SMILES string.

    Parameters
    ----------
    smiles : iterable
        The SMILES string to parse. Should conform to the OpenSMILES
        specification.
    explicit_hydrogen : bool
        Whether hydrogens should be explicit nodes in the outout graph, or be
        implicit in 'hcount' attributes.
    reinterprit_aromatic : bool
        Whether aromaticity should be determined from the created molecule,
        instead of taken from the SMILES string.

    Returns
    -------
    type of systems: (complex or monomer)
    total chemical formula (heavy atoms only)
    charged or not
    has isotopes or not
    """
    total_formula = defaultdict(int)
    is_monomer = True 
    is_charged = False
    has_isotope = False
    bond_to_order = {'-': 1, '=': 2, '#': 3, '$': 4, ':': 1.5, '.': 0}
    for tokentype, token in _tokenize(smiles):
        if tokentype == TokenType.ATOM:
            i_atom = parse_atom(token)
            if 'element' in i_atom.keys():
                total_formula[i_atom['element']] += 1
            if i_atom['charge'] != 0:
                is_charged = True
            if 'isotope' in i_atom.keys():
                has_isotope = True
        elif tokentype == TokenType.BRANCH_START:
            pass 
        elif tokentype == TokenType.BRANCH_END:
            pass 
        elif tokentype == TokenType.BOND_TYPE:
            if bond_to_order[token] == 0:
                is_monomer = False
        elif tokentype == TokenType.RING_NUM:
            pass 
        elif tokentype == TokenType.EZSTEREO:
            pass
            # LOGGER.warning('E/Z stereochemical information, which is specified by "%s", will be discarded', token)

    return total_formula, is_monomer, is_charged, has_isotope

BUFFER_SIZE_LINES = 100000  # Maximum number of lines to buffer in memory

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, help='Input CID-SMILES filename')
args = parser.parse_args()

r = open(args.input, "r")
f_isotope = args.input + '-isotopes'
f_HO_nC_mon = args.input + '-HO-noCharge-monomer'
f_HO_nC_com = args.input + '-HO-noCharge-complex'
f_HO_wC_mon = args.input + '-HO-wCharge-monomer'
f_HO_wC_com = args.input + '-HO-wCharge-complex'
f_HCNO_nC_mon = args.input + '-HCNO-noCharge-monomer'
f_HCNO_nC_com = args.input + '-HCNO-noCharge-complex'
f_HCNO_wC_mon = args.input + '-HCNO-wCharge-monomer'
f_HCNO_wC_com = args.input + '-HCNO-wCharge-complex'
f_HCNOFPSCl_nC_mon = args.input + '-HCNOFPSCl-noCharge-monomer'
f_HCNOFPSCl_nC_com = args.input + '-HCNOFPSCl-noCharge-complex'
f_HCNOFPSCl_wC_mon = args.input + '-HCNOFPSCl-wCharge-monomer'
f_HCNOFPSCl_wC_com = args.input + '-HCNOFPSCl-wCharge-complex'
f_remaining_nC_mon = args.input + '-remaining-noCharge-monomer'
f_remaining_nC_com = args.input + '-remaining-noCharge-complex'
f_remaining_wC_mon = args.input + '-remaining-wCharge-monomer'
f_remaining_wC_com = args.input + '-remaining-wCharge-complex'
w_isotope = open(f_isotope, "w")
w_HO_nC_mon = open(f_HO_nC_mon, "w")
w_HO_nC_com = open(f_HO_nC_com, "w")
w_HO_wC_mon = open(f_HO_wC_mon, "w")
w_HO_wC_com = open(f_HO_wC_com, "w")
w_HCNO_nC_mon = open(f_HCNO_nC_mon, "w")
w_HCNO_nC_com = open(f_HCNO_nC_com, "w")
w_HCNO_wC_mon = open(f_HCNO_wC_mon, "w")
w_HCNO_wC_com = open(f_HCNO_wC_com, "w")
w_HCNOFPSCl_nC_mon = open(f_HCNOFPSCl_nC_mon, "w")
w_HCNOFPSCl_nC_com = open(f_HCNOFPSCl_nC_com, "w")
w_HCNOFPSCl_wC_mon = open(f_HCNOFPSCl_wC_mon, "w")
w_HCNOFPSCl_wC_com = open(f_HCNOFPSCl_wC_com, "w")
w_remaining_nC_mon = open(f_remaining_nC_mon, "w")
w_remaining_nC_com = open(f_remaining_nC_com, "w")
w_remaining_wC_mon = open(f_remaining_wC_mon, "w")
w_remaining_wC_com = open(f_remaining_wC_com, "w")

buf_isotope = ""
buf_HO_nC_mon = ""
buf_HO_nC_com = ""
buf_HO_wC_mon = ""
buf_HO_wC_com = ""
buf_HCNO_nC_mon = ""
buf_HCNO_nC_com = ""
buf_HCNO_wC_mon = ""
buf_HCNO_wC_com = ""
buf_HCNOFPSCl_nC_mon = ""
buf_HCNOFPSCl_nC_com = ""
buf_HCNOFPSCl_wC_mon = ""
buf_HCNOFPSCl_wC_com = ""
buf_remaining_nC_mon = ""
buf_remaining_nC_com = ""
buf_remaining_wC_mon = ""
buf_remaining_wC_com = ""

bufLines = 0
for lineIn in r:
    total_formula, is_monomer, is_charged, has_isotope = read_smiles(lineIn.split()[1])
    # print(molinfo.inchi_str)
    # print(molinfo.charge)
    # print(molinfo.total_formula)
    if has_isotope:
        buf_isotope += lineIn
    elif set(total_formula.keys()).issubset({"H", "O"}):
        if is_charged:
            if is_monomer:
                buf_HO_wC_mon += lineIn
            else:
                buf_HO_wC_com += lineIn
        else:
            if is_monomer:
                buf_HO_nC_mon += lineIn
            else:
                buf_HO_nC_com += lineIn
    elif set(total_formula.keys()).issubset({"H", "C", "N", "O"}):
        if is_charged:
            if is_monomer:
                buf_HCNO_wC_mon += lineIn
            else:
                buf_HCNO_wC_com += lineIn
        else:
            if is_monomer:
                buf_HCNO_nC_mon += lineIn
            else:
                buf_HCNO_nC_com += lineIn
    elif set(total_formula.keys()).issubset(
            {"H", "C", "N", "O", "F", "P", "S", "Cl"}
        ):
        if is_charged:
            if is_monomer:
                buf_HCNOFPSCl_wC_mon += lineIn
            else:
                buf_HCNOFPSCl_wC_com += lineIn
        else:
            if is_monomer:
                buf_HCNOFPSCl_nC_mon += lineIn
            else:
                buf_HCNOFPSCl_nC_com += lineIn
    else:
        if is_charged:
            if is_monomer:
                buf_remaining_wC_mon += lineIn
            else:
                buf_remaining_wC_com += lineIn
        else:
            if is_monomer:
                buf_remaining_nC_mon += lineIn
            else:
                buf_remaining_nC_com += lineIn

    bufLines += 1
    if bufLines >= BUFFER_SIZE_LINES:
        # Flush buffer to disk
        if buf_isotope:
            w_isotope.write(buf_isotope)
        if buf_HO_nC_mon:
            w_HO_nC_mon.write(buf_HO_nC_mon)
        if buf_HO_nC_com:
            w_HO_nC_com.write(buf_HO_nC_com)
        if buf_HO_wC_mon:
            w_HO_wC_mon.write(buf_HO_wC_mon)
        if buf_HO_wC_com:
            w_HO_wC_com.write(buf_HO_wC_com)
        if buf_HCNO_nC_mon:
            w_HCNO_nC_mon.write(buf_HCNO_nC_mon)
        if buf_HCNO_nC_com:
            w_HCNO_nC_com.write(buf_HCNO_nC_com)
        if buf_HCNO_wC_mon:
            w_HCNO_wC_mon.write(buf_HCNO_wC_mon)
        if buf_HCNO_wC_com:
            w_HCNO_wC_com.write(buf_HCNO_wC_com)
        if buf_HCNOFPSCl_nC_mon:
            w_HCNOFPSCl_nC_mon.write(buf_HCNOFPSCl_nC_mon)
        if buf_HCNOFPSCl_nC_com:
            w_HCNOFPSCl_nC_com.write(buf_HCNOFPSCl_nC_com)
        if buf_HCNOFPSCl_wC_mon:
            w_HCNOFPSCl_wC_mon.write(buf_HCNOFPSCl_wC_mon)
        if buf_HCNOFPSCl_wC_com:
            w_HCNOFPSCl_wC_com.write(buf_HCNOFPSCl_wC_com)
        if buf_remaining_nC_mon:
            w_remaining_nC_mon.write(buf_remaining_nC_mon)
        if buf_remaining_nC_com:
            w_remaining_nC_com.write(buf_remaining_nC_com)
        if buf_remaining_wC_mon:
            w_remaining_wC_mon.write(buf_remaining_wC_mon)
        if buf_remaining_wC_com:
            w_remaining_wC_com.write(buf_remaining_wC_com)

        buf_isotope = ""
        buf_HO_nC_mon = ""
        buf_HO_nC_com = ""
        buf_HO_wC_mon = ""
        buf_HO_wC_com = ""
        buf_HCNO_nC_mon = ""
        buf_HCNO_nC_com = ""
        buf_HCNO_wC_mon = ""
        buf_HCNO_wC_com = ""
        buf_HCNOFPSCl_nC_mon = ""
        buf_HCNOFPSCl_nC_com = ""
        buf_HCNOFPSCl_wC_mon = ""
        buf_HCNOFPSCl_wC_com = ""
        buf_remaining_nC_mon = ""
        buf_remaining_nC_com = ""
        buf_remaining_wC_mon = ""
        buf_remaining_wC_com = ""
        bufLines = 0


# Flush remaining buffer to disk
if buf_isotope:
    w_isotope.write(buf_isotope)
if buf_HO_nC_mon:
    w_HO_nC_mon.write(buf_HO_nC_mon)
if buf_HO_nC_com:
    w_HO_nC_com.write(buf_HO_nC_com)
if buf_HO_wC_mon:
    w_HO_wC_mon.write(buf_HO_wC_mon)
if buf_HO_wC_com:
    w_HO_wC_com.write(buf_HO_wC_com)
if buf_HCNO_nC_mon:
    w_HCNO_nC_mon.write(buf_HCNO_nC_mon)
if buf_HCNO_nC_com:
    w_HCNO_nC_com.write(buf_HCNO_nC_com)
if buf_HCNO_wC_mon:
    w_HCNO_wC_mon.write(buf_HCNO_wC_mon)
if buf_HCNO_wC_com:
    w_HCNO_wC_com.write(buf_HCNO_wC_com)
if buf_HCNOFPSCl_nC_mon:
    w_HCNOFPSCl_nC_mon.write(buf_HCNOFPSCl_nC_mon)
if buf_HCNOFPSCl_nC_com:
    w_HCNOFPSCl_nC_com.write(buf_HCNOFPSCl_nC_com)
if buf_HCNOFPSCl_wC_mon:
    w_HCNOFPSCl_wC_mon.write(buf_HCNOFPSCl_wC_mon)
if buf_HCNOFPSCl_wC_com:
    w_HCNOFPSCl_wC_com.write(buf_HCNOFPSCl_wC_com)
if buf_remaining_nC_mon:
    w_remaining_nC_mon.write(buf_remaining_nC_mon)
if buf_remaining_nC_com:
    w_remaining_nC_com.write(buf_remaining_nC_com)
if buf_remaining_wC_mon:
    w_remaining_wC_mon.write(buf_remaining_wC_mon)
if buf_remaining_wC_com:
    w_remaining_wC_com.write(buf_remaining_wC_com)

w_isotope.close()
w_HO_nC_mon.close()
w_HO_nC_com.close()
w_HO_wC_mon.close()
w_HO_wC_com.close()
w_HCNO_nC_mon.close()
w_HCNO_nC_com.close()
w_HCNO_wC_mon.close()
w_HCNO_wC_com.close()
w_HCNOFPSCl_nC_mon.close()
w_HCNOFPSCl_nC_com.close()
w_HCNOFPSCl_wC_mon.close()
w_HCNOFPSCl_wC_com.close()
w_remaining_nC_mon.close()
w_remaining_nC_com.close()
w_remaining_wC_mon.close()
w_remaining_wC_com.close()

r.close()

