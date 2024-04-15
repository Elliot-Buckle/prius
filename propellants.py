cea_cards = {
    "ABS":"""
fuel C3.85H4.85N0.43 C 3.85 H 4.85 N 0.43      wt%=100
h,cal=15060     t(k)=293.15
""",

"PLA":"""
fuel C3H4O2 C 3 H 4 O 2      wt%=100
h,cal=-72180     t(k)=293.15
""",

"PMMA":"""
fuel C5H8O2 C 5 H 8 O 2      wt%=100
h,cal=-148700     t(k)=293.15
""",

"PC":"""
fuel C16H14O3 C 16 H 14 O 3      wt%=100
h,cal=-24620     t(k)=293.15
""",

"HDPE":"""
fuel C2H4 C 2 H 4     wt%=100
h,cal=-6118.546845   t(k)=293.15
"""
}

regression_coefficients = {
    "HTPB/GOX":{"a":0.000907, "n":0.67},
    "PMMA/GOX":{"a":0.000359, "n":0.615},
    "HDPE/GOX":{"a":0.000416, "n":0.498},
    "PLA/GOX":{"a":0.00001, "n":0.674},
    "ABS/GOX":{"a":0.000073, "n":0.351},
    "PC/GOX":{"a":0.000012, "n":0.723},
    "HDPE/N2O":{"a":0.000531, "n":0.331},
    "PMMA/N2O":{"a":0.000614, "n":0.335},
    "HTPB/N2O":{"a":0.000927, "n":0.347},
    "ABS/N2O":{"a":0.000103, "n":0.215},             
}

densities = {
    "HTPB":930,
    "PMMA":1180,
    "HDPE":961,
    "ABS":975,
    "PC":1200,
    "PLA":1240
}