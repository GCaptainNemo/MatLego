from ase.build import cut
import os
from pymatgen.ext.matproj import MPRester
import re
from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifWriter
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
import numpy as np



debuglogger = print
episilon = 0.000001

def is_number(s):
    # 判断一个字符串是不是数
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False


class crys_data:
    """ To store the data of each atoms(name, Femi_level, atom radii, atom color). """
    # 金属性，是金属为1，不是为0
    metal_or_not = [np.nan,
0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
1, 1, 1, 0, 0, 0, 0, 0, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 0, 0, 0, 0, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 0, 0, 0, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 0, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1]


    # （单位埃）
    vander_wals_radii = \
        [
        np.nan, 1.2, 1.43, 2.12,      # He 为推断值
        1.98, 1.91, 1.77, 1.66, 1.5,
        1.46, 1.58, 2.5, 2.51, 2.25,  # Ne 为推断值
        2.19, 1.9, 1.89, 1.82, 1.83,
        2.73, 2.62, 2.58, 2.46, 2.42,
        2.45, 2.45, 2.44, 2.4, 2.4,
        2.38, 2.39, 2.32, 2.29, 1.88,
        1.82, 1.86, 2.25, 3.21, 2.84,
        2.75, 2.52, 2.56, 2.45, 2.44,
        2.46, 2.44, 2.15, 2.53, 2.49,
        2.43, 2.42, 2.47, 1.99, 2.04,
        2.06, 3.48, 3.03, 2.98, 2.88,
        2.92, 2.95, np.nan, 2.9, 2.87,
        2.83, 2.79, 2.87, 2.81, 2.83,
        2.79, 2.8, 2.74, 2.63, 2.53,
        2.57, 2.49, 2.48, 2.41, 2.29,
        2.32, 2.45, 2.47, 2.6, 2.54,
        np.nan, np.nan, np.nan, np.nan, np.nan,
        2.8, 2.93, 2.88, 2.71, 2.82,
        2.81, 2.83, 3.05, 3.4, 3.05,
        2.7, np.nan, np.nan, np.nan, np.nan,
        np.nan, np.nan, np.nan, np.nan, np.nan,
        np.nan, np.nan, np.nan, np.nan, np.nan,
        np.nan, np.nan, np.nan, np.nan, np.nan
        ]

    # 元素的费米能级(单位：*10^-19 J)
    atom_Femi_level = [np.nan,  np.nan,  np.nan,   7.5,   22.6,   np.nan,
                       np.nan,  np.nan,   np.nan,    np.nan,    np.nan,
                       5.1,     11.4,     9.0,   np.nan,    np.nan,
                       np.nan,  np.nan,   np.nan,   3.3,     6.9,
                       10.9,    13.7,     16,    11.1,      17.4,
                       17.7,    18.4,     18.8,  11.3,      15.1,
                       8.0,     11.6,     np.nan,    np.nan,    np.nan,
                       np.nan,  2.9,      6.3,   9.0,       11.3,
                       8.5,     9.4,      np.nan,    10.2,    10.2,
                       9.7,     8.8,      12,    6.6,       10.3,
                       12.5,    np.nan,   np.nan,    np.nan,   2.5,
                       5.8,     8.29,     8.77,      8.74,   8.8,
                       np.nan,  9.0,      7.02,      9.0,    9.19,
                       9.28,    9.36,     9.46,      9.58,   7.77,
                       9.70,    11.6,     13.4,      14.71,  15.45,
                       15.96,   15.83,    9.56,      8.8,    np.nan,
                       6.2,     9.5,      11.3,      13.1,   np.nan,
                       np.nan,  np.nan,   5.2,       np.nan, 9.00,
                       10.9,    12.2,     8.6,       np.nan, 8.7,
                       np.nan,  np.nan,   np.nan,    np.nan, np.nan,
                       np.nan,  np.nan,   np.nan,    np.nan, np.nan,
                       np.nan,  np.nan,   np.nan,    np.nan, np.nan,
                       np.nan,  np.nan,   np.nan,    np.nan, np.nan,
                       np.nan,  np.nan,   np.nan]

    # 原子半径，单位(埃)
    atom_radii_lst = [np.nan, 0.53, 0.31, 1.67, 1.12, 0.87, 0.67, 0.56, 0.48, 0.42, 0.38,
                      1.90, 1.45, 1.18, 1.11, 0.98, 0.88, 0.79, 0.71, 2.43, 1.94, 1.84,
                      1.76, 1.71, 1.66, 1.61, 1.56, 1.52, 1.49, 1.45, 1.42, 1.36, 1.25,
                      1.14, 1.03, 0.94, 0.88, 2.65, 2.19, 2.12, 2.06, 1.98, 1.90, 1.83,
                      1.78, 1.73, 1.69, 1.65, 1.61, 1.56, 1.45, 1.33, 1.23, 1.15, 1.08,
                      2.98, 2.53, np.nan, np.nan, 2.47, 2.06, 2.05, 2.38, 2.31, 2.33,
                      2.25, 2.28, np.nan, 2.26, 2.22, 2.22, 2.17, 2.08, 2.00, 1.93, 1.88,
                      1.85, 1.80, 1.77, 1.74, 1.71, 1.56, 1.54, 1.43, 1.35, np.nan, 1.20,
                      np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                      np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                      np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                      np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

    # 原子的颜色
    atom_color = [(0, 0, 0, 0), (0.18039, 0.5451, 0.34118, 1.0), (0.12549, 0.69804, 0.66667, 1.0),
                  (0, 1, 0.49804), (0, 1, 0, 1.0), (0.46667, 0.53333, 0.6, 1.0),
                  (0.51765, 0.43922, 1, 1.0), (1, 0.87059, 0.67843, 1.0), (0.39216, 0.58431, 0.92941, 1.0),
                  (0.7451, 0.7451, 0.7451, 1.0),  (0.94118, 1, 0.94118, 1.0), (0.41961, 0.55686, 0.13725, 1.0), (0.93333, 0.9098, 0.66667, 1.0),
                  (0.72941, 0.33333, 0.82745, 1.0), (0.78039, 0.08235, 0.52157, 1.0), (0.81569, 0.12549, 0.56471, 1.0),
                  (1, 0, 1, 1.0), (0.93333, 0.5098, 0.93333, 1.0), (0.86667, 0.62745, 0.86667, 1.0),
                  (0.96078, 0.87059, 0.70196, 1.0), (0.93333, 0.9098, 0.80392, 1.0), (0.80392, 0.78431, 0.69412, 1.0),
                  (0.93333, 0.93333, 0.87843, 1.0), (0.80392, 0.80392, 0.75686, 1.0), (0.96078, 0.96078, 0.86275, 1.0),
                  (1, 1, 0.87843, 1.0), (0.2549, 0.41176, 0.88235, 1.0), (0.11765, 0.56471, 1, 1.0),
                  (0.52941, 0.80784, 0.921571, 1.0), (1, 0.84314, 0, 1.0), (0.93333, 0.86667, 0.5098, 1.0),
                  (0, 0.98039, 0.60392, 1.0), (0.67843, 1, 0.18431, 1.0), (0.19608, 0.80392, 0.19608, 1.0),
                  (0.60392, 0.80392, 0.19608, 1.0), (0.6902, 0.87843, 0.90196, 1.0), (0, 0.80784, 0.81961, 1.0),
                  (0.25098, 0.87843, 0.81569, 1.0), (0.87843, 1, 1, 1.0), (0.37255, 0.61961, 0.62745, 1.0),
                  (0.49804, 1, 0.83137, 1.0), (0.33333, 0.41961, 0.18431, 1.0), (0.94118, 1, 1, 1.0),
                  (0.90196, 0.90196, 0.98039, 1.0), (1, 0.89412, 0.88235, 1.0), (0.97255, 0.97255, 1, 1.0),
                  (0.86275, 0.86275, 0.86275, 1.0), (0.5451, 0.53725, 0.43922, 1.0), (0.51373, 0.5451, 0.51373, 1.0),
                  (1, 0.94118, 0.96078, 1.0), (0.93333, 0.87843, 0.89804, 1.0), (0.85882, 0.43922, 0.57647, 1.0),
                  (0.6902, 0.18824, 0.37647, 1.0), (1, 0.54902, 0, 1.0), (1, 0.49804, 0.31373, 1.0),
                  (0.94118, 0.50196, 0.50196, 1.0), (1, 1, 0.94118, 1.0), (1, 0.96078, 0.93333, 1.0),
                  (0.41176, 0.41176, 0.41176, 1.0),  (0, 0, 0.50196, 1.0), (0.98039, 0.94118, 0.90196, 1.0),
                  (1, 0.92157, 0.80392, 1.0), (0.62745, 0.32157, 0.17647, 1.0), (0.87059, 0.72157, 0.52941, 1.0),
                  (1, 0.89412, 0.76863, 1.0),  (0.13333, 0.5451, 0.13333, 1.0), (0.8549, 0.64706, 0.12549, 1.0),
                  (0.72157, 0.52549, 0.04314, 1.0), (0.73725, 0.56078, 0.56078, 1.0), (1, 0.71373, 0.75686, 1.0),
                  (0.80392, 0.71765, 0.61961, 1.0), (0.5451, 0.4902, 0.41961, 1.0), (0.80392, 0.70196, 0.5451, 1.0),
                  (0.5451, 0.47451, 0.36863, 1.0), (1, 0.98039, 0.80392, 1.0), (0.93333, 0.91373, 0.74902, 1.0),
                  (0.80392, 0.78824, 0.64706, 1.0),  (1, 0.38824, 0.27843, 1.0), (1, 0.27059, 0, 1.0),
                  (1, 0.41176, 0.70588, 1.0), (0.84706, 0.74902, 0.84706, 1.0), (1, 0.98039, 0.98039, 1.0),
                  (0.93333, 0.91373, 0.91373, 1.0), (0.80392, 0.78824, 0.78824, 1.0), (1, 0.96078, 0.93333, 1.0),
                  (0.80392,	0.77255, 0.74902, 1.0), (0.69804, 0.13333, 0.13333, 1.0), (0.64706, 0.16471, 0.16471, 1.0),
                  (0.91373, 0.58824, 0.47843, 1.0), (0.98039, 0.50196, 0.44706, 1.0), (1, 0.62745, 0.47843, 1.0),
                  (1, 0.64706, 0, 1.0), (0.5451, 0.5451, 0.51373, 1.0), (0.94118, 1, 0.94118, 1.0),
                  (0.87843, 0.93333, 0.87843, 1.0), (0.75686, 0.80392, 0.75686, 1.0), (0.6, 0.19608, 0.8, 1.0),
                  (0.58039, 0, 0.82745, 1.0), (0.54118, 0.16863, 0.88627, 1.0), (0.62745, 0.12549, 0.94118, 1.0),
                  (0.57647, 0.43922, 0.85882, 1.0), (0.95686, 0.64314, 0.37647, 1.0), (0.82353, 0.70588, 0.54902, 1.0),
                  (1, 0.87059, 0.67843, 1.0), (0.93333, 0.81176, 0.63137, 1.0), (0.5451, 0.52549, 0.5098, 1.0),
                  (1, 0.93725, 0.85882, 1.0), (0.93333, 0.87451, 0.8, 1.0), (0.80392, 0.75294, 0.6902, 1.0),
                  (0.5451, 0.51373, 0.47059, 1.0), (1, 0.89412, 0.76863, 1.0), (0.93333, 0.8359, 0.71765, 1.0),
                  (0.80392, 0.75686, 0.77255, 1.0), (0.5451, 0.51373, 0.52549, 1.0), (1, 0.89412, 0.88235, 1.0),
                  (0.82353, 0.41176, 0.11765, 1.0),  (1, 0.8549, 0.72549, 1.0), (0.80392, 0.68627, 0.58431, 1.0),
                 (0.93333, 0.83529, 0.82353, 1.0)]

    atom_formula = ('', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
                        'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V',
                        'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As',
                        'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc',
                        'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',
                        'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu',
                        'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
                        'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl',
                        'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa',
                        'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md',
                        'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg',
                        'Uub', 'Uut', 'Uuq', 'Uup', 'Uuh', 'Uus', 'Uuo')

    conduct2DMaterials = [
        "AgF2", "C", "CeSiI", "Co(OH)2", "CoI2", "CoO2", "CrSe2",
        "CuAgTe2", "CuTe", "DySBr", "DySI", "ErHCl", "ErSeI",
        "EuOBr", "EuOI", "FeBr2", "FeI2", "FeO2", "FeOCl",
        "FeS", "FeSe", "FeTe", "Hf3Te2", "HfSiTe", "HfTe2",
        "La2PBr2", "La2PI2", "La2TeI2", "LaCl", "LaGeI",
        "LaI2", "LuHCl", "MoS2", "MoTe2", "Na2PdH2",
        "Nb2CS2", "NbF4", "NbS2", "NbSe2",
        "NbSe2", "NdOBr", "OsOCl2", "PrOBr", "PrOI",
        "ReSe2", "RuOCl2", "Sc2NCl2", "ScCl", "ScHCl",
        "SiTe2", "SmOBr", "SmSI", "Ta2CS2",
        "TaS2", "TaSe2", "TbBr",
        "TbCl", "Ti2PTe2", "TiBr2", "TiCl2", "TiNI",
        "TiOBr", "TiOCl", "TiSe2", "TiTe2", "TmI2",
        "TmOI", "VOBr", "VOCl", "VS2", "VSe2",
        "VTe2", "W2N3", "WTe2", "YbOBr", "YbOCl",
        "YCBr", "YCI", "YGaI", "Zr2PTe2", "ZrBr",
        "ZrCl", "ZrSiTe", "ZrTe2", "ZrTiSe4", "ZrTiTe4"]

    conductor_Elemental_composition = \
        [
        ["Ag", "F"], ["C"], ["Ce", "Si", "I"], ["Co", "O", "H"], ["Co", "I"], ["Co", "O"], ["Cr", "Se"],
        ["Cu", "Ag", "Te"], ["Cu", "Te"], ["Dy", "S", "Br"], ["Dy", "S", "I"], ["Er", "H", "Cl"], ["Er", "Se", "I"],
        ["Eu", "O", "Br"], ["Eu", "O", "I"], ["Fe", "Br"], ["Fe", "I"], ["Fe", "O"], ["Fe", "O", "Cl"],
        ["Fe", "S"], ["Fe", "Se"], ["Fe", "Te"], ["Hf", "Te"], ["Hf", "Si", "Te"], ["Hf", "Te"],
        ["La", "P", "Br"], ["La", "P", "I"], ["La", "Te", "I"], ["La", "Cl"], ["La", "Ge", "I"],
        ["La", "I"], ["Lu", "H", "Cl"], ["Mo", "S"], ["Mo", "Te"], ["Na", "Pd", "H"],
        ["Nb", "C", "S"], ["Nb", "F"], ["Nb", "S"], ["Nb", "Se"], ["Nd", "O", "Br"], ["Os", "O", "Cl"],
        ["Pr", "O", "Br"], ["Pr", "O", "I"],
        ["Re", "Se"], ["Ru", "O", "Cl"], ["Sc", "N", "Cl"], ["Sc", "Cl"], ["Sc", "H", "Cl"],
        ["Si", "Te"], ["Sm", "O", "Br"], ["Sm", "S", "I"], ["Ta", "C", "S"], ["Ta", "S"],
        ["Ta", "Se"], ["Tb", "Br"],
        ["Tb", "Cl"], ["Ti", "P", "Te"], ["Ti", "Br"], ["Ti", "Cl"], ["Ti", "N", "I"],
        ["Ti", "O", "Br"], ["Ti", "O", "Cl"], ["Ti", "Se"], ["Ti", "Te"], ["Tm", "I"],
        ["Tm", "O", "I"], ["V", "O", "Br"], ["V", "O", "Cl"], ["V", "S"], ["V", "Se"],
        ["V", "Te"], ["W", "N"], ["W", "Te"], ["Yb", "O", "Br"], ["Yb", "O", "Cl"],
        ["Y", "C", "Br"], ["Y", "C", "I"], ["Y", "Ga", "I"], ["Zr", "P", "Te"], ["Zr", "Br"],
        ["Zr", "Cl"], ["Zr", "Si", "Te"], ["Zr", "Te"], ["Zr", "Ti", "Se"], ["Zr", "Ti", "Te"]
    ]



    semiconductor2DMaterial = \
        [
        ["TiS2", 0.1], ["CoBr2", 0.2], ["CoCl2", 0.2], ["Cu2Te", 0.2],
        ["CuCl2", 0.2], ["ErSCl", 0.2], ["CdOCl", 0.3], ["NiI2", 0.3],
        ["CrSBr", 0.4], ["La2GeI2", 0.4], ["ZrI2", 0.4], ["CrOBr", 0.5],
        ["HoSI", 0.5], ["Mn(OH)2", 0.5], ["ZrSe2", 0.5], ["Bi", 0.6],
        ["CrOCl", 0.6], ["HfSe2", 0.6], ["KAgSe", 0.6], ["LaBr", 0.6],
        ["LaBr2", 0.6], ["TiNBr", 0.6], ["TiNCl", 0.6], ["As", 0.7],
        ["CrI2", 0.7], ["Sb2Te2Se", 0.7], ["Sb2Te3", 0.7], ["Sb2TeSe2", 0.7],
        ["SnTe", 0.7], ["YbI2", 0.7], ["ZrNI", 0.7], ["CrBr2", 0.8],
        ["GaGeTe", 0.8], ["NiBr2", 0.8], ["Sb2TeSe2", 0.8], ["SnSe2", 0.8],
        ["Bi2Se3", 0.9], ["Bi2Te2Se", 0.9], ["CLu2Cl2", 0.9], ["FeCl2", 0.9],
        ["LiAlTe2", 0.9], ["P", 0.9], ["PdCl2", 0.9], ["SbTeI", 0.9],
        ["Sc2CCl2", 0.9], ["VOBr2", 0.9], ["VOCl2", 0.9], ["Bi2Te2S", 1.0],
        ["Bi2Te3", 1.0], ["HfNI", 1.0], ["In2Se3", 1.0], ["LiAuI4", 1.0],
        ["Tl2O", 1.0], ["ZrCl2", 1.0], ["Bi2TeSe2", 1.1], ["GeSe", 1.1],
        ["MoTe2", 1.1], ["NiCl2", 1.1], ["RhTeCl", 1.1], ["WTe2", 1.1],
        ["AuSe", 1.2], ["Bi2Te2S", 1.2], ["PdS2", 1.2], ["Sb", 1.2],
        ["VI2", 1.2], ["ZrS2", 1.2], ["AgBr", 1.3], ["HfS2", 1.3],
        ["NiO2", 1.3], ["PtSe2", 1.3], ["VBr2", 1.3], ["ZrNI", 1.3],
        ["GaTe", 1.3], ["MnI2", 1.4], ["Tl2S", 1.4], ["VCl2", 1.4],
        ["BiOI", 1.5], ["CuBr", 1.5], ["InSe", 1.5], ["KTlO", 1.5],
        ["MoSe2", 1.5], ["SnS2", 1.5], ["AgNO2", 1.6], ["BiTeI", 1.6],
        ["GeS", 1.6], ["MoS2", 1.6], ["PbTe", 1.6], ["WSe2", 1.6],
        ["InSe", 1.7], ["PtO2", 1.7], ["ZrNBr", 1.7], ["BiTeCl", 1.8],
        ["GaSe", 1.8], ["HgI2", 1.8], ["PtS2", 1.8], ["WS2", 1.8],
        ["ZrNBr", 1.8], ["MnBr2", 1.9], ["P", 1.9], ["ZrNCl", 1.9],
        ["AuBr", 2.0], ["AuI", 2.0], ["CuI", 2.0], ["GeI2", 2.0],
        ["HfNBr", 2.0], ["ZnI2", 2.0], ["CuI", 2.1], ["GeI2", 2.1],
        ["MnCl2", 2.1], ["SnO", 2.1], ["AgI", 2.2], ["GaS", 2.2],
        ["HfNBr", 2.2], ["SiH", 2.2], ["Cd(OH)2", 2.3], ["GaS", 2.3],
        ["GaTeCl", 2.3], ["InOBr", 2.3], ["PbIF", 2.3], ["BiOBr", 2.4],
        ["CdI2", 2.4], ["HfNCl", 2.4], ["PbO", 2.4], ["PbF4", 2.5],
        ["PbI2", 2.5], ["BiOCl", 2.7], ["KTlCl4", 2.8], ["NaOH", 2.8],
        ["ZrNBr", 2.8], ["AgO4Cl", 2.9]
    ]

    semiconductor_Element_composition =\
        [['Ti', 'S'], ['Co', 'Br'], ['Co', 'Cl'], ['Cu', 'Te'], ['Cu', 'Cl'], ['Er', 'S', 'Cl'], ['Cd', 'O', 'Cl'],
     ['Ni', 'I'], ['Cr', 'S', 'Br'], ['La', 'Ge', 'I'], ['Zr', 'I'], ['Cr', 'O', 'Br'], ['Ho', 'S', 'I'],
     ['Mn', 'O', 'H'], ['Zr', 'Se'], ['Bi'], ['Cr', 'O', 'Cl'], ['Hf', 'Se'], ['K', 'Ag', 'Se'], ['La', 'Br'],
     ['La', 'Br'], ['Ti', 'N', 'Br'], ['Ti', 'N', 'Cl'], ['As'], ['Cr', 'I'], ['Sb', 'Te', 'Se'], ['Sb', 'Te'],
     ['Sb', 'Te', 'Se'], ['Sn', 'Te'], ['Yb', 'I'], ['Zr', 'N', 'I'], ['Cr', 'Br'],
     ['Ga', 'Ge', 'Te'], ['Ni', 'Br'], ['Sb', 'Te', 'Se'], ['Sn', 'Se'], ['Bi', 'Se'], ['Bi', 'Te', 'Se'],
     ['C', 'Lu', 'Cl'], ['Fe', 'Cl'], ['Li', 'Al', 'Te'], ['P'], ['Pd', 'Cl'], ['Sb', 'Te', 'I'], ['Sc', 'C', 'Cl'],
     ['V', 'O', 'Br'], ['V', 'O', 'Cl'], ['Bi', 'Te', 'S'], ['Bi', 'Te'], ['Hf', 'N', 'I'], ['In', 'Se'],
     ['Li', 'Au', 'I'], ['Tl', 'O'], ['Zr', 'Cl'], ['Bi', 'Te', 'Se'], ['Ge', 'Se'], ['Mo', 'Te'],
     ['Ni', 'Cl'], ['Rh', 'Te', 'Cl'], ['W', 'Te'], ['Au', 'Se'],
     ['Bi', 'Te', 'S'], ['Pd', 'S'], ['Sb'], ['V', 'I'], ['Zr', 'S'], ['Ag', 'Br'], ['Hf', 'S'],
     ['Ni', 'O'], ['Pt', 'Se'], ['V', 'Br'], ['Zr', 'N', 'I'], ['Ga', 'Te'], ['Mn', 'I'], ['Tl', 'S'],
     ['V', 'Cl'], ['Bi', 'O', 'I'], ['Cu', 'Br'], ['In', 'Se'], ['K', 'Tl', 'O'], ['Mo', 'Se'], ['Sn', 'S'],
     ['Ag', 'N', 'O'], ['Bi', 'Te', 'I'], ['Ge', 'S'], ['Mo', 'S'], ['Pb', 'Te'], ['W', 'Se'], ['In', 'Se'],
     ['Pt', 'O'], ['Zr', 'N', 'Br'], ['Bi', 'Te', 'Cl'], ['Ga', 'Se'], ['Hg', 'I'], ['Pt', 'S'], ['W', 'S'],
     ['Zr', 'N', 'Br'], ['Mn', 'Br'], ['P'], ['Zr', 'N', 'Cl'], ['Au', 'Br'],
     ['Au', 'I'], ['Cu', 'I'], ['Ge', 'I'], ['Hf', 'N', 'Br'], ['Zn', 'I'], ['Cu', 'I'],
     ['Ge', 'I'], ['Mn', 'Cl'], ['Sn', 'O'], ['Ag', 'I'], ['Ga', 'S'], ['Hf', 'N', 'Br'],
     ['Si', 'H'], ['Cd', 'O', 'H'], ['Ga', 'S'], ['Ga', 'Te', 'Cl'], ['In', 'O', 'Br'], ['Pb', 'I', 'F'],
     ['Bi', 'O', 'Br'], ['Cd', 'I'], ['Hf', 'N', 'Cl'], ['Pb', 'O'], ['Pb', 'F'], ['Pb', 'I'], ['Bi', 'O', 'Cl'],
     ['K', 'Tl', 'Cl'], ['Na', 'O', 'H'], ['Zr', 'N', 'Br'], ['Ag', 'O', 'Cl']]

    # > 3.0eV
    Insulator2DMaterials = \
    [
        ["PbBrF", 3.0], ["SnO", 3.0], ["CdBr2", 3.2], ["ScOBr", 3.2],
        ["Mg(OH)2", 3.3], ["OLuI", 3.3], ["ZrNCl", 3.3], ["BaHI", 3.4],
        ["LaOI", 3.4], ["ZnBr2", 3.4], ["PbClF", 3.5], ["MgI2", 3.6],
        ["TlF", 3.6], ["Ca(OH)2", 3.7], ["CaHI", 3.7], ["CaI2", 3.8],
        ["SrI2", 3.8], ["CdCl2", 3.9], ["LiOH", 3.9], ["SnF4", 3.9],
        ["SrHI", 3.9], ["LaOBr", 4.0], ["CaHBr", 4.2], ["ZnCl2", 4.3],
        ["OLuBr", 4.4], ["SrHBr", 4.4], ["YOCl", 4.4], ["ZnCl2", 4.5],
        ["RbCl", 4.6], ["BN", 4.7], ["MgBr2", 4.8], ["NaCN", 4.8],
        ["SrBrF", 5.3], ["AlOCl", 5.8], ["MgCl2", 6.0], ["LiBH4", 6.4]
    ]

    Insulator_Element_composition = \
        [
            ["Pb", "Br", "F"], ["Sn", "O"], ["Cd", "Br"], ["Sc", "O", "Br"],
            ["Mg", "O", "H"], ["O", "Lu", "I"], ["Zr", "N", "Cl"], ["Ba", "H", "I"],
            ["La", "O", "I"], ["Zn", "Br"], ["Pb", "Cl", "F"], ["Mg", "I"],
            ["Tl", "F"], ["Ca", "O", "H"], ["Ca", "H", "I"], ["Ca", "I"],
            ["Sr", "I"], ["Cd", "Cl"], ["Li", "O", "H"], ["Sn", "F"],
            ["Sr", "H", "I"], ["La", "O", "Br"], ["Ca", "H", "Br"], ["Zn", "Cl"],
            ["O", "Lu", "Br"], ["Sr", "H", "Br"], ["Y", "O", "Cl"], ["Zn", "Cl"],
            ["Rb", "Cl"], ["B", "N"], ["Mg", "Br"], ["Na", "C", "N"],
            ["Sr", "Br", "F"], ["Al", "O", "Cl"], ["Mg", "Cl"], ["Li", "B", "H"]
        ]

class MathOperation:
    """ Different MathOperation"""
    @staticmethod
    def found_periodic_function(function):
        for i in range(1, len(function)):  # i是可能周期数
            for j in range(0, len(function)):
                if function[j] == function[j % i]:
                    if j != len(function) - 1:
                        function[j] = function[j % i]
                        continue
                    else:
                        print("该函数周期为:", i)
                        return [function[x] for x in range(i)]
                else:
                    break
        print("该函数周期为：", len(function))
        return function

    @staticmethod
    def get_theat(vector1, vector2, rad=True):
        c = vector1 @ vector2
        amod = np.linalg.norm(vector1)
        bmod = np.linalg.norm(vector2)
        if rad == True:
            return(np.arccos(c / amod / bmod))
        else:
            return np.arccos(c / amod / bmod) / np.pi * 180

class FileOperations:
    """ File Operation."""
    def __init__(self, file_dir):
        self.file_dir = file_dir

    def file(self):                                    # 返回文件名，后缀列表
        s = []
        for root, dirs, files in os.walk(self.file_dir):
            for i in files:
                s.append(i)
        return s

    def file_dirction(self):             # 提取文件目录,返回file_dir文件夹中所有文件地址列表
        s = []
        for root, dirs, files in os.walk(self.file_dir):
            for i in files:
                l = (root + "\\" + i)
                s.append(l)
        return s

    def read_files(self):             # 读取文件夹中每一文件，从而得到Atoms对象列表
        dir_lst = self.file_dirction()
        lst = []
        for i in dir_lst:
            try:
                lst.append(read(i, index=None, format=None, parallel=True))
            except Exception as e:
                print(e)
        return lst


class CrystalOperationsNum1:
    """ Unary operation of crys."""
    def __init__(self, crys):
        self.crys = crys

    def add_cell(self):                  # 沿着z轴扩胞
        return self.crys.repeat((1, 1, 2))

    def gamaToObtuseangle(self):            # 把晶胞betta变成钝角方便比较切应变, 注：不改变self.crys
        cell_par_lst = self.crys.get_cell_lengths_and_angles()
        if cell_par_lst[5] >= 90:
            return self.crys
        else:
            vector = self.crys.get_cell()
            b = cell_par_lst[1]
            gama = cell_par_lst[5]
            mubiaob = [b * np.cos((180 - gama) / 180 * np.pi), b * np.sin((180 - gama) / 180 * np.pi), 0]
            A = vector.T
            r = np.linalg.solve(A, mubiaob)  # 求解线性方程组，直角坐标系下----用晶胞坐标系表示
            real = cut(self.crys, a=[1, 0, 0], b=r.T, c=[0, 0, 1], origo=(0, 0, 0))
            return real

    def gamaToAcuteangle(self):
        cell_par_lst = self.crys.get_cell_lengths_and_angles()
        if cell_par_lst[5] <= 90:
            return self.crys
        else:
            vector = self.crys.get_cell()
            b = cell_par_lst[1]
            gama = cell_par_lst[5]
            mubiaob = [b * np.cos((180 - gama) / 180 * np.pi), b * np.sin((180 - gama) / 180 * np.pi), 0]
            A = vector.T
            r = np.linalg.solve(A, mubiaob)  # 求解线性方程组，直角坐标系下----用晶胞坐标系表示
            real = cut(self.crys, a=[1, 0, 0], b=r.T, c=[0, 0, 1], origo=(0, 0, 0))

            return real

    # def SixToTri(self):
    #     crys_par_lst = self.crys.get_cell_lengths_and_angles()
    #     vector = self.crys.get_cell()
    #     A = np.array(vector).T
    #     b = crys_par_lst[1]
    #     gama = crys_par_lst[5]
    #     mubiaob = np.array([0, 2 * b * np.sin(gama / 180 * np.pi), 0])
    #     r = np.linalg.solve(A, mubiaob)  # 求解线性方程组，直角坐标系下----用晶胞坐标系表示
    #     layer_view = cut(self.crys, a=[1, 0, 0], b=r.T, c=[0, 0, 1], origo=(0, 0, 0))
    #     return layer_view

    def to_cif(self):     # 输出crys.cif文件
        f = AseAtomsAdaptor.get_structure(self.crys)
        print(AseAtomsAdaptor.get_structure(self.crys))
        f1 = str(f) + '\n'
        j_pv_lst = re.findall('abc\s\s\s:(.*?)\n', f1)[0]   # abc   :  19.257300  19.569178  21.133988
        j1_pv_lst = j_pv_lst.split(' ')                   # abc   :   6.419100   6.523059   7.044663
        while '' in j1_pv_lst:
            j1_pv_lst.remove('')
        par_lst_matrix = self.crys.get_cell()
        # 物质种类（比如：Lu2 Al4）
        y1 = re.findall('Full\sFormula\s(.*?)\n', f1)[0]
        material = re.findall('Full\sFormula\s(.*?)\n', f1)[0].lstrip('(').rstrip(')')
        self.mat_name = material
        elements = material.split(' ')
        print(elements)              # (Re4 S8)\(Re108 S216)
        zmb_lst = [chr(i) for i in range(97,123)] + [chr(i) for i in range(65,91)]
        szb_lst = [str(i) for i in range(0,10)]
        element_lst = []
        number_lst =[]
        for element in elements:       # 这个循环之后，得到原子与对应的原子数两个列表
            letter_lst = list(element)
            symbol_lst = []
            num_lst = []
            element_lst1 =[]
            number_lst1 =[]
            for i in range(len(letter_lst)):
                if letter_lst[i] in szb_lst:
                    num_lst.append(letter_lst[i])
                if letter_lst[i] in zmb_lst:
                    symbol_lst.append(letter_lst[i])
                element1 = ''.join(symbol_lst)
                element_lst1.append(element1)
                number1 = ''.join(num_lst)
                number_lst1.append(number1)
            ys = 'a'     # 元素
            gs = '0'     # 个数
            for i in element_lst1:
                if len(i) >= len(ys):
                    ys = i
            element_lst.append(i)
            for i in number_lst1:
                if len(i) >= len(gs):
                    gs = i
            number_lst.append(i)
        print(element_lst)
        print(number_lst)
        par_lst_species = []                   # 用于Cifwrite参数(species)的
        for i in range(len(element_lst)):
            num = int(number_lst[i])
            for j in range(num):
                par_lst_species.append(element_lst[i])
        print(par_lst_species)
        # 每个原子的坐标
        ord_lst = []  # 最终Cifwriter所需要的coords参数
        ord_lst2 = []  # 储存的形式为
        for element in element_lst:
            ord_lst1 = re.findall(element + '\s\s\s\s(.*?)\n', f1)
            for ord in ord_lst1:
                ord1 = ord.split(' ')
                while '' in ord1:
                    ord1.remove('')
                ord_lst2.append(ord1)
        for ord in ord_lst2:
            ord1 = []
            for string in ord:
                num = float(string)
                ord1.append(num)
                if len(ord1) == 3:
                    ord2 = ord1
            ord_lst.append(ord2)
        par_lst_coords = ord_lst
        # 构建Structure类
        structure = Structure(par_lst_matrix, par_lst_species, par_lst_coords)
        slab = CifWriter(structure,
                         write_magmoms=True)  # struct (Structure) – structure to write; symprec (float) – If not none, finds the symmetry of the structure and writes the cif with symmetry information. Passes symprec to the SpacegroupAnalyzer; write_magmoms (bool) – If True, will write magCIF file. Incompatible with symprec
        slab.write_file(r"C:\Users\wang1\Desktop\2D\2dmaterialsstack_supercell\{}.cif".format(self.mat_name))

class IdToCif:             # 用Id在pymatgen上下载Cif文件
    """ Download cif from www.materialproject.com with ID."""
    def __init__(self, mat_id):
        self.id = mat_id

    def to_cif(self):
        m = MPRester("iEe4KC8U9QBczglh")
        structure = m.query(self.id, ['initial_structure'])  # 改变不同的晶体ID录入不同的数字
        f = structure[0]['initial_structure']
        f1 = str(f) + '\n'
        j_pv_lst = re.findall('abc\s\s\s:\s\s\s(.*?)\n', f1)[0]
        j1_pv_lst = j_pv_lst.split('   ')
        a = float(j1_pv_lst[0])
        b = float(j1_pv_lst[1])
        c = float(j1_pv_lst[2])
        par_lst_matrix = [[a, 0, 0],
                          [0, b, 0],
                          [0, 0, c]]
        print(par_lst_matrix)
        # 物质种类（比如：Lu2 Al4）
        y1 = re.findall('Full\sFormula\s(.*?)\n', f1)[0]
        # y1_lst = y1.lstrip('(').rstrip(')')
        material = re.findall('Full\sFormula\s(.*?)\n', f1)[0].lstrip('(').rstrip(')')
        self.mat_name = material
        elements = material.split(' ')
        element_lst = []
        num_lst = []
        for element in elements:       # 这个循环之后，得到原子与对应的原子数两个列表
            letter_lst = list(element)
            num_lst.append(letter_lst[-1])
            letter_lst.pop(-1)
            element1 = ''.join(letter_lst)
            element_lst.append(element1)
        par_lst_species = []                   # 用于Cifwrite参数(species)的
        for i in range(len(num_lst)):
            num = int(num_lst[i])
            for j in range(num):
                par_lst_species.append(element_lst[i])
        print(par_lst_species)
        # 每个原子的坐标
        ord_lst = []   # 最终Cifwriter所需要的coords参数
        ord_lst2 = []  # 储存的形式为
        for element in element_lst:
            ord_lst1 = re.findall(element+'\s\s\s\s(.*?)\n',f1)
            for ord in ord_lst1:
                ord1 = ord.split(' ')
                while '' in ord1:
                    ord1.remove('')
                ord_lst2.append(ord1)
        for ord in ord_lst2:
            ord1 = []
            for string in ord:
                num = float(string)
                ord1.append(num)
                if len(ord1) == 3:
                    ord2 = ord1
            ord_lst.append(ord2)
        par_lst_coords = ord_lst
        print(par_lst_coords)
        # 构建Structure类
        structure = Structure(par_lst_matrix, par_lst_species, par_lst_coords)
        slab = CifWriter(structure, write_magmoms=True)   # struct (Structure) – structure to write; symprec (float) – If not none, finds the symmetry of the structure and writes the cif with symmetry information. Passes symprec to the SpacegroupAnalyzer; write_magmoms (bool) – If True, will write magCIF file. Incompatible with symprec
        slab.write_file(r'C:\Users\wang1\Desktop\2D\download vasp\{}.cif'.format(self.mat_name))
