import re

INT = r'[-+]?[0-9]+(?:[eE][-+]?[0-9]+)?'
FLOAT = r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'
KEY = r'[\w]+[\w-]*'
# WORD = r'[^\s]+'
WORD = r'\w+'
NUMBER = r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'
TRUE = r'[tT]|true|True|TRUE'
FALSE = r'[fF]|false|False|FALSE'
DELIMITER= r'\s*[,\s]\s*'
re_key = r'(?P<key>[A-Za-z_]+[A-Za-z0-9_-]*)'
re_values = r'(?:' \
            r'(?P<single_value>[^\s"]+)' \
            r'|' \
            r'["\{\}](?P<quoted_value>[^"\{\}]+)["\{\}]' \
            r')\s*'

complex_xyz_string = (
    ' '  # start with a separator
    'str=astring '
    'quot="quoted value" '
    u'quote_special="a_to_Z_$%%^&*\xfc\u2615" '
    r'escaped_quote="esc\"aped" '
    'true_value '
    'false_value = F '
    'integer=22 '
    'floating=1.1 '
    'int_array={1 2 3} '
    'float_array="3.3 4.4" '
    'a3x3_array="1 4 7 2 5 8 3 6 9" '  # fortran ordering
    'Lattice="  4.3  0.0 0.0 0.0  3.3 0.0 0.0 0.0  7.0 " '  # spaces in array
    'scientific_float=1.2e7 '
    'scientific_float_2=5e-6 '
    'scientific_float_array="1.2 2.2e3 4e1 3.3e-1 2e-2" '
    'not_array="1.2 3.4 text" '
    'nested_brackets=[[1,2],[3,4]] '  # gets flattented if not 3x3
    'bool_array={T F T F} '
    'bool_array_2=" T, F, T " '  # leading spaces
    'not_bool_array=[T F S] '
    # read and write
    u'\xfcnicode_key=val\xfce '
    u'unquoted_special_value=a_to_Z_$%%^&*\xfc\u2615 '
    '2body=33.3 '
    'hyphen-ated '
    # parse only
    'many_other_quotes=({[4 8 12]}) '
    'comma_separated="7, 4, -1" '
    'bool_array_commas=[T, T, F, T] '
    'Properties=species:S:1:pos:R:3 '
    'multiple_separators       '
    'double_equals=abc=xyz '
    'trailing'
)
# \s*(?P<key>[\w]+[\w-]*)
# (?:
#   (?:\s*=\s*
#     (?:
#       (?P<false>[f]|false)|
#       (?P<true>[t]|true)|
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

regex = (f'''
\s*(?P<key>{KEY})\s*
(?:
  (?:\s*=\s*
    (?:
      (?P<single_true>{TRUE})|
      (?P<single_false>{FALSE})|
      (?P<single_int>{INT})|
      (?P<single_float>{FLOAT})|
      \[(?P<bool_list>\s*(?:(?:{TRUE}|{FALSE}){DELIMITER}(?:{TRUE}|{FALSE})?)+)\]|
      (?P<array>\[(?:\[[^\]]*\]|[^\[\]])*\])|
      # "(?P<quoted_sting>(?:[^\\\"]|\\\")*)"|
      (?P<single_word>{WORD})
    )\s+
  )
  |
  (?P<key_only>\s+)
)
''')

if __name__ == '__main__':
    print(regex)
    print(complex_xyz_string)
    converter = re.compile(''.join(regex.split()))
    # r = converter.findall(complex_xyz_string)
    for match in converter.finditer(complex_xyz_string):
        print([(k, v) for k, v in match.groupdict().items() if v is not None])
