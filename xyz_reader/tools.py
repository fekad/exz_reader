import re
from pathlib import Path
from typing import Union
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
import numpy as np

from ase.utils import basestring
from ase.spacegroup.spacegroup import Spacegroup

UNPROCESSED_KEYS = ['uid']


def key_val_str_to_dict(string, sep=None):
    """
    Parse an xyz properties string in a key=value and return a dict with
    various values parsed to native types.

    Accepts brackets or quotes to delimit values. Parses integers, floats
    booleans and arrays thereof. Arrays with 9 values are converted to 3x3
    arrays with Fortran ordering.

    If sep is None, string will split on whitespace, otherwise will split
    key value pairs with the given separator.

    """
    # store the closing delimiters to match opening ones
    delimiters = {
        "'": "'",
        '"': '"',
        '(': ')',
        '{': '}',
        '[': ']',
    }

    # Make pairs and process afterwards
    kv_pairs = [
        [[]]]  # List of characters for each entry, add a new list for new value
    delimiter_stack = []  # push and pop closing delimiters
    escaped = False  # add escaped sequences verbatim

    # parse character-by-character unless someone can do nested brackets
    # and escape sequences in a regex
    for char in string.strip():
        if escaped:  # bypass everything if escaped
            kv_pairs[-1][-1].extend(['\\', char])
            escaped = False
        elif delimiter_stack:  # inside brackets
            if char == delimiter_stack[-1]:  # find matching delimiter
                delimiter_stack.pop()
            elif char in delimiters:
                delimiter_stack.append(delimiters[char])  # nested brackets
            elif char == '\\':
                escaped = True  # so escaped quotes can be ignored
            else:
                kv_pairs[-1][-1].append(char)  # inside quotes, add verbatim
        elif char == '\\':
            escaped = True
        elif char in delimiters:
            delimiter_stack.append(delimiters[char])  # brackets or quotes
        elif (sep is None and char.isspace()) or char == sep:
            if kv_pairs == [[[]]]:  # empty, beginning of string
                continue
            elif kv_pairs[-1][-1] == []:
                continue
            else:
                kv_pairs.append([[]])
        elif char == '=':
            if kv_pairs[-1] == [[]]:
                del kv_pairs[-1]
            kv_pairs[-1].append([])  # value
        else:
            kv_pairs[-1][-1].append(char)

    kv_dict = {}

    for kv_pair in kv_pairs:
        if len(kv_pair) == 0:  # empty line
            continue
        elif len(kv_pair) == 1:  # default to True
            key, value = ''.join(kv_pair[0]), 'T'
        else:  # Smush anything else with kv-splitter '=' between them
            key, value = ''.join(kv_pair[0]), '='.join(
                ''.join(x) for x in kv_pair[1:])

        if key.lower() not in UNPROCESSED_KEYS:
            # Try to convert to (arrays of) floats, ints
            split_value = re.findall(r'[^\s,]+', value)
            try:
                try:
                    numvalue = np.array(split_value, dtype=int).tolist()
                except (ValueError, OverflowError):
                    # don't catch errors here so it falls through to bool
                    numvalue = np.array(split_value, dtype=float).tolist()
                if len(numvalue) == 1:
                    numvalue = numvalue[0]  # Only one number
                # elif len(numvalue) == 9:
                #     # special case: 3x3 matrix, fortran ordering
                #     numvalue = np.array(numvalue).reshape((3, 3), order='F')
                value = numvalue
            except (ValueError, OverflowError):
                pass  # value is unchanged

            # Parse boolean values: 'T' -> True, 'F' -> False,
            #                       'T T F' -> [True, True, False]
            if isinstance(value, basestring):
                str_to_bool = {'T': True, 'F': False}

                try:
                    boolvalue = [str_to_bool[vpart] for vpart in
                                 re.findall(r'[^\s,]+', value)]
                    if len(boolvalue) == 1:
                        value = boolvalue[0]
                    else:
                        value = boolvalue
                except KeyError:
                    pass  # value is unchanged

        kv_dict[key] = value

    return kv_dict


def key_val_dict_to_str(d, sep=' '):
    """
    Convert atoms.info dictionary to extended XYZ string representation
    """
    if len(d) == 0:
        return ''
    s = ''
    type_val_map = {(bool, True): 'T',
                    (bool, False): 'F',
                    (np.bool_, True): 'T',
                    (np.bool_, False): 'F'}

    s = ''
    for key in d.keys():
        val = d[key]
        if isinstance(val, dict):
            continue
        if hasattr(val, '__iter__'):
            val = np.array(val)
            val = ' '.join(str(type_val_map.get((type(x), x), x))
                           for x in val.reshape(val.size, order='F'))
            val.replace('[', '')
            val.replace(']', '')
        elif isinstance(val, Spacegroup):
            val = val.symbol
        else:
            val = type_val_map.get((type(val), val), val)

        if val is None:
            s = s + '%s%s' % (key, sep)
        elif isinstance(val, basestring) and ' ' in val:
            s = s + '%s="%s"%s' % (key, val, sep)
        else:
            s = s + '%s=%s%s' % (key, str(val), sep)

    return s.strip()


def new_convert(text: str):
    """Convert string to python object by guessing its type"""

    # array-like object
    elements = text.split()
    if len(elements) > 1:
        return [new_convert(el) for el in elements]

    try:
        return int(text)
    except ValueError:
        pass

    try:
        return float(text)
    except ValueError:
        pass

    if text == 'T' or text == 'True':
        return True
    elif text == 'F' or text == 'False':
        return False
    else:
        # return as it is
        return text

# https://regex101.com/r/52W6cm/4
# \s*(?P<key>[\w]+[\w-]*)
# (?:
#   (?:\s*=\s*
#     (?:
#       (?P<false>[Ff]|false)|
#       (?P<true>[Tt]|true)|
#       [\[\{\(\"]+(?P<bool_list>\s*(?:(?:[tf]|true|false)\s*[,\s]\s*(?:[tf]|true|false)?)+)[\]\}\)\"]+|
#       [\[\{\(\"]+(?P<number_list>\s*(?:(?:[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*[,\s]\s*(?:[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)?)+)[\]\}\)\"]+|
#       \"(?P<quoted>(?:[^\\"]|\\\")*)\"|
#       (?P<nested1>\[(?:\[[^\]]*\]|[^\[\]])*\])|
#       (?P<nested2>\{(?:\[[^\]]*\]|[^\[\]])*\})|
#       (?P<single_number>[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)|
#       (?P<single_string>[^\s]+)
#     )\s+
#   )
#   |
#   (?P<key_only>\s+)
# )
#
# KEY = r'[\w]+[\w-]*'
# NUMBER = r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'
# TRUE = r'[t]|true'
# FALSE = r'[f]|false'
#
# re_key = r'(?P<key>[A-Za-z_]+[A-Za-z0-9_-]*)'
# re_values = r'(?:' \
#             r'(?P<single_value>[^\s"]+)' \
#             r'|' \
#             r'["\{\}](?P<quoted_value>[^"\{\}]+)["\{\}]' \
#             r')\s*'
