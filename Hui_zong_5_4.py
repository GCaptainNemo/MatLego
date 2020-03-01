"""
The GUI of matlego 1.5.
"""
from PyQt5 import QtCore, QtGui, QtWidgets
from OpenGL.GLU import gluProject, gluUnProject
from OpenGL.GL import *
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifWriter
from pymatgen.ext.matproj import MPRester
from pymatgen.io.ase import AseAtomsAdaptor
from ase.build import cut, sort
from ase.spacegroup import get_spacegroup
from ase.io import read
from ase.atom import Atom
from re import findall
import os
import numpy as np
from fractions import Fraction
import qdarkstyle

from sys import argv, exit
import sqlite3
from operator import itemgetter
from copy import deepcopy
from ase.visualize import view


class Hetero_junction():
    Name = ...
    Rotate_angle = ...
    Optimal_orientation = ...
    Atomsobject = ...
    Normal_strain_limit = ...
    Shear_strain_limit = ...
    Area_error_limit = ...
    def __init__(self):
        ...

class MyaxisItem(pg.opengl.GLAxisItem):
    def __init__(self, size=None, antialias=True, glOptions='translucent', cell_coordinate=False, fixed=False):
        super(MyaxisItem, self).__init__(size, antialias, glOptions)
        self.cell_coordinate = cell_coordinate
        self.fixed = fixed

    def paint(self):
        self.setupGLState()
        if self.antialias:
            glEnable(GL_LINE_SMOOTH)
            glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
        if not self.cell_coordinate:
            if not self.fixed:
                glBegin(GL_LINES)
                x, y, z = self.size()
                glColor4f(0, 0, 1, .6)  # z is blue
                glVertex3f(0, 0, 0)
                glVertex3f(0, 0, z)
                glColor4f(0, 1, 0, .6)  # y is green
                glVertex3f(0, 0, 0)
                glVertex3f(0, y, 0)
                glColor4f(1, 0, 0, .6)  # x is red
                glVertex3f(0, 0, 0)
                glVertex3f(x, 0, 0)
                glEnd()
        else:
            if not self.fixed:
                glBegin(GL_LINES)
                glViewport(0, 0, 50, 50)
                glDisable(GL_DEPTH_TEST)
                x = self.cell_coordinate[0]
                y = self.cell_coordinate[1]
                z = self.cell_coordinate[2]
                glColor4f(0, 0, 1, .6)  # z is blue
                glVertex3f(0, 0, 0)
                glVertex3f(z[0], z[1], z[2])
                glColor4f(0, 1, 0, .6)  # y is green
                glVertex3f(0, 0, 0)
                glVertex3f(y[0], y[1], y[2])
                glColor4f(1, 0, 0, .6)  # x is red
                glVertex3f(0, 0, 0)
                glVertex3f(x[0], x[1], x[2])
                glEnd()


class ClassificationCrystalOperationsNum1:
    def __init__(self, crys, outputdirs):
        self.crys = crys
        self.outputdirs = outputdirs
        if self.cut_exf_layer():
            self.to_cif(self.crys)

    def add_cell(self):                   # 沿c方向扩胞两倍
        return self.crys.repeat((1, 1, 2))

    def cut_exf_layer(self):
        after_add_cell_self_ATO = self.add_cell()
        pos_lst = after_add_cell_self_ATO.get_positions().tolist()
        all_dis = after_add_cell_self_ATO.get_all_distances(vector=False)
        order = len(pos_lst)
        Atomic_number_lst = after_add_cell_self_ATO.get_atomic_numbers()
        vander_wals_matrix = np.diag([crys_data.vander_wals_radii[Atomic_number_lst[i]] for i in range(order)])
        vander_wals_matrix = all_dis + np.ones((order, order)) * 1.3 - \
                             np.transpose(np.ones((order, order)) @ vander_wals_matrix) - np.ones(
            (order, order)) @ vander_wals_matrix
        dis_or_not_matrix = (vander_wals_matrix > 0)
        gouzaolist = [pos_piece for pos_piece in pos_lst]
        for k in range(len(gouzaolist)):
            gouzaolist[k].append(k)
        gouzaolist.sort(key=itemgetter(2))  # 根据z轴由小到大排序
        height = 0
        exfoliat_height = 0
        index_lst = [gouzaolist[0][3]]
        for i in range(len(gouzaolist) - 1):
            if not dis_or_not_matrix[gouzaolist[i][3]][gouzaolist[i + 1][3]]:  # valence bond
                height += (gouzaolist[i + 1][2] - gouzaolist[i][2])
                index_lst.append(gouzaolist[i + 1][3])
            elif (gouzaolist[i + 1][2] - gouzaolist[i][2]) / \
                    all_dis[gouzaolist[i][3]][gouzaolist[i + 1][3]] < .5:
                height += (gouzaolist[i + 1][2] - gouzaolist[i][2])
                index_lst.append(gouzaolist[i + 1][3])
            else:
                exfoliat_height = gouzaolist[i + 1][2] - gouzaolist[i][2]
                break
        if exfoliat_height:
            return True

    def to_cif(self, Atomsobject):     # 输出crys.cif文件
        f = AseAtomsAdaptor.get_structure(Atomsobject)
        f1 = str(f) + '\n'
        j_pv_lst = findall('abc\s\s\s:(.*?)\n', f1)[0]   # abc   :  19.257300  19.569178  21.133988
        j1_pv_lst = j_pv_lst.split(' ')                   # abc   :   6.419100   6.523059   7.044663
        while '' in j1_pv_lst:
            j1_pv_lst.remove('')
        par_lst_matrix = self.crys.get_cell()
        # 物质种类（比如：Lu2 Al4）
        material = findall('Full\sFormula\s(.*?)\n', f1)[0].lstrip('(').rstrip(')')
        self.mat_name = material
        elements = material.split(' ')
        zmb_lst = [chr(i) for i in range(97, 123)] + [chr(i) for i in range(65,91)]
        szb_lst = [str(i) for i in range(0, 10)]
        element_lst = []
        number_lst =[]
        for element in elements:       # 这个循环之后，得到原子与对应的原子数两个列表
            letter_lst = list(element)
            symbol_lst = []
            num_lst = []
            element_lst1 =[]
            number_lst1 =[]
            for letter in letter_lst:
                if letter in szb_lst:
                    num_lst.append(letter)
                if letter in zmb_lst:
                    symbol_lst.append(letter)
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
        par_lst_species = []                   # 用于Cifwrite参数(species)的
        for i, element in enumerate(element_lst):
            num = int(number_lst[i])
            for j in range(num):
                par_lst_species.append(element)
        ord_lst = []  # 最终Cifwriter所需要的coords参数
        ord_lst2 = []  # 储存的形式为
        for element in element_lst:
            ord_lst1 = findall(element + '\s\s\s\s(.*?)\n', f1)
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
        slab.write_file(self.outputdirs + "\{}.cif".format(self.mat_name))

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

    electronegativity = [np.nan, 2.20, 3.89, .98, 1.57, 2.04, 2.55, 3.04,
                         3.44, 3.98, 3.67, .93, 1.31, 1.61, 1.90, 2.19,
                         2.58, 3.16, 3.3, .82, 1.00, 1.36, 1.54, 1.63,
                         1.66, 1.55, 1.83, 1.88, 1.91, 1.90, 1.65, 1.81,
                         2.01, 2.18, 2.55, 2.96, 3.00, .82, .95, 1.22,
                         1.33, 1.6, 2.16, 1.9, 2.2, 2.28, 2.20, 1.93, 1.69,
                         1.78, 1.96, 2.05, 2.1, 2.66, 2.67, .79, .89, 1.1,
                         1.12, 1.13, 1.14, 1.13, 1.17, 1.2, 1.2, 1.1, 1.22,
                         1.23, 1.24, 1.25, 1.1, 1.27, 1.3, 1.5, 2.36, 1.9, 2.2,
                         2.20, 2.28, 2.54, 2.00, 1.62, 2.33, 2.02, 2.0, 2.2, 2.2,
                         .7, .9, 1.1, 1.3, 1.5, 1.38, 1.36, 1.28, 1.13, 1.28,
                         1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.291]



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
    atom_radii_lst = [np.nan, .53, .31, 1.67, 1.12, .87, .67, .56, .48, .42, .38,
                      1.90, 1.45, 1.18, 1.11, .98, .88, .79, .71, 2.43, 1.94, 1.84,
                      1.76, 1.71, 1.66, 1.61, 1.56, 1.52, 1.49, 1.45, 1.42, 1.36, 1.25,
                      1.14, 1.03, .94, .88, 2.65, 2.19, 2.12, 2.06, 1.98, 1.90, 1.83,
                      1.78, 1.73, 1.69, 1.65, 1.61, 1.56, 1.45, 1.33, 1.23, 1.15, 1.08,
                      2.98, 2.53, np.nan, np.nan, 2.47, 2.06, 2.05, 2.38, 2.31, 2.33,
                      2.25, 2.28, np.nan, 2.26, 2.22, 2.22, 2.17, 2.08, 2.00, 1.93, 1.88,
                      1.85, 1.80, 1.77, 1.74, 1.71, 1.56, 1.54, 1.43, 1.35, np.nan, 1.20,
                      np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                      np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                      np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                      np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

    # 原子的颜色(jmol color)
    atom_color = [(0, 0, 0, 0), (1.0, 1.0, 1.0, .5), (.8509803921568627, 1.0, 1.0, .5),
                  (.8, .5019607843137255, 1.0, .5), (.7607843137254902, 1.0, 0, .5),
                  (1.0, .7098039215686275, .7098039215686275, .5), (.5647058823529412, .5647058823529412, .5647058823529412, .5),
                  (.18823529411764706, .3137254901960784, .9725490196078431, .5), (1.0, .050980392156862744, .050980392156862744, .5),
                  (.5647058823529412, .8784313725490196, .3137254901960784, .5), (.7019607843137254, .8901960784313725, .9607843137254902, .5),
                  (.6705882352941176, .3607843137254902, .9490196078431372, .5), (.5411764705882353, 1.0, 0, .5), (.7490196078431373, .6509803921568628, .6509803921568628, .5), (.9411764705882353, .7843137254901961, .6274509803921569, .5),
                  (1.0, .5019607843137255, 0, .5), (1.0, 1.0, .18823529411764706, .5), (.12156862745098039, .9411764705882353, .12156862745098039, .5), (.5019607843137255, .8196078431372549, .8901960784313725, .5), (.5607843137254902, .25098039215686274, .8313725490196079, .5),
                  (.23921568627450981, 1.0, 0, .5), (.9019607843137255, .9019607843137255, .9019607843137255, .5), (.7490196078431373, .7607843137254902, .7803921568627451, .5), (.6509803921568628, .6509803921568628, .6705882352941176, .5), (.5411764705882353, .6, .7803921568627451, .5),
                  (.611764705882353, .47843137254901963, .7803921568627451, .5), (.8784313725490196, .4, .2, .5), (.9411764705882353, .5647058823529412, .6274509803921569, .5), (.3137254901960784, .8156862745098039, .3137254901960784, .5), (.7843137254901961, .5019607843137255, .2, .5),
                  (.49019607843137253, .5019607843137255, .6901960784313725, .5), (.7607843137254902, .5607843137254902, .5607843137254902, .5), (.4, .5607843137254902, .5607843137254902, .5), (.7411764705882353, .5019607843137255, .8901960784313725, .5), (1.0, .6313725490196078, 0, .5),
                  (.6509803921568628, .1607843137254902, .1607843137254902, .5), (.3607843137254902, .7215686274509804, .8196078431372549, .5), (.4392156862745098, .1803921568627451, .6901960784313725, .5), (0, 1.0, 0, .5), (.5803921568627451, 1.0, 1.0, .5), (.5803921568627451, .8784313725490196, .8784313725490196, .5),
                  (.45098039215686275, .7607843137254902, .788235294117647, .5), (.32941176470588235, .7098039215686275, .7098039215686275, .5), (.23137254901960785, .6196078431372549, .6196078431372549, .5), (.1411764705882353, .5607843137254902, .5607843137254902, .5), (.0392156862745098, .49019607843137253, .5490196078431373, .5),
                  (0, .4117647058823529, .5215686274509804, .5), (.7529411764705882, .7529411764705882, .7529411764705882, .5), (1.0, .8509803921568627, .5607843137254902, .5), (.6509803921568628, .4588235294117647, .45098039215686275, .5), (.4, .5019607843137255, .5019607843137255, .5), (.6196078431372549, .38823529411764707, .7098039215686275, .5),
                  (.8313725490196079, .47843137254901963, 0, .5), (.5803921568627451, 0, .5803921568627451, .5), (.25882352941176473, .6196078431372549, .6901960784313725, .5), (.3411764705882353, .09019607843137255, .5607843137254902, .5), (0, .788235294117647, 0, .5), (.4392156862745098, .8313725490196079, 1.0, .5), (1.0, 1.0, .7803921568627451, .5),
                  (.8509803921568627, 1.0, .7803921568627451, .5), (.7803921568627451, 1.0, .7803921568627451, .5), (.6392156862745098, 1.0, .7803921568627451, .5), (.5607843137254902, 1.0, .7803921568627451, .5), (.3803921568627451, 1.0, .7803921568627451, .5), (.27058823529411763, 1.0, .7803921568627451, .5), (.18823529411764706, 1.0, .7803921568627451, .5),
                  (.12156862745098039, 1.0, .7803921568627451, .5), (0, 1.0, .611764705882353, .5), (0, .9019607843137255, .4588235294117647, .5), (0, .8313725490196079, .3215686274509804, .5), (0, .7490196078431373, .2196078431372549, .5), (0, .6705882352941176, .1411764705882353, .5), (.30196078431372547, .7607843137254902, 1.0, .5),
                  (.30196078431372547, .6509803921568628, 1.0, .5), (.12941176470588237, .5803921568627451, .8392156862745098, .5), (.14901960784313725, .49019607843137253, .6705882352941176, .5), (.14901960784313725, .4, .5882352941176471, .5), (.09019607843137255, .32941176470588235, .5294117647058824, .5), (.8156862745098039, .8156862745098039, .8784313725490196, .5),
                  (1.0, .8196078431372549, .13725490196078433, .5), (.7215686274509804, .7215686274509804, .8156862745098039, .5), (.6509803921568628, .32941176470588235, .30196078431372547, .5), (.3411764705882353, .34901960784313724, .3803921568627451, .5), (.6196078431372549, .30980392156862746, .7098039215686275, .5), (.6705882352941176, .3607843137254902, 0, .5),
                  (.4588235294117647, .30980392156862746, .27058823529411763, .5), (.25882352941176473, .5098039215686274, .5882352941176471, .5), (.25882352941176473, 0, .4, .5), (0, .49019607843137253, 0, .5), (.4392156862745098, .6705882352941176, .9803921568627451, .5), (0, .7294117647058823, 1.0, .5), (0, .6313725490196078, 1.0, .5), (0, .5607843137254902, 1.0, .5),
                  (0, .5019607843137255, 1.0, .5), (0, .4196078431372549, 1.0, .5), (.32941176470588235, .3607843137254902, .9490196078431372, .5), (.47058823529411764, .3607843137254902, .8901960784313725, .5), (.5411764705882353, .30980392156862746, .8901960784313725, .5), (.6313725490196078, .21176470588235294, .8313725490196079, .5), (.7019607843137254, .12156862745098039, .8313725490196079, .5),
                  (.7019607843137254, .12156862745098039, .7294117647058823, .5), (.7019607843137254, .050980392156862744, .6509803921568628, .5), (.7411764705882353, .050980392156862744, .5294117647058824, .5), (.7803921568627451, 0, .4, .5), (.8, 0, .34901960784313724, .5), (.8196078431372549, 0, .30980392156862746, .5), (.8509803921568627, 0, .27058823529411763, .5),
                  (.8784313725490196, 0, .2196078431372549, .5), (.9019607843137255, 0, .1803921568627451, .5), (.9215686274509803, 0, .14901960784313725, .5), (.80392, .75686, .77255, 1.0), (.5451, .51373, .52549, 1.0), (1, .89412, .88235, 1.0),
                  (.82353, .41176, .11765, 1.0),  (1, .8549, .72549, 1.0), (.80392, .68627, .58431, 1.0),  (1, .41176, .70588, 1.0), (.84706, .74902, .84706, 1.0), (1, .98039, .98039, 1.0),
                  (.93333, .91373, .91373, 1.0), (.80392, .78824, .78824, 1.0)]


    Element_formula = ['', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
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
                       'Uub', 'Uut', 'Uuq', 'Uup', 'Uuh', 'Uus', 'Uuo']

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
        ["TiS2", .1], ["CoBr2", .2], ["CoCl2", .2], ["Cu2Te", .2],
        ["CuCl2", .2], ["ErSCl", .2], ["CdOCl", .3], ["NiI2", .3],
        ["CrSBr", .4], ["La2GeI2", .4], ["ZrI2", .4], ["CrOBr", .5],
        ["HoSI", .5], ["Mn(OH)2", .5], ["ZrSe2", .5], ["Bi", .6],
        ["CrOCl", .6], ["HfSe2", .6], ["KAgSe", .6], ["LaBr", .6],
        ["LaBr2", .6], ["TiNBr", .6], ["TiNCl", .6], ["As", .7],
        ["CrI2", .7], ["Sb2Te2Se", .7], ["Sb2Te3", .7], ["Sb2TeSe2", .7],
        ["SnTe", .7], ["YbI2", .7], ["ZrNI", .7], ["CrBr2", .8],
        ["GaGeTe", .8], ["NiBr2", .8], ["Sb2TeSe2", .8], ["SnSe2", .8],
        ["Bi2Se3", .9], ["Bi2Te2Se", .9], ["CLu2Cl2", .9], ["FeCl2", .9],
        ["LiAlTe2", .9], ["P", .9], ["PdCl2", .9], ["SbTeI", .9],
        ["Sc2CCl2", .9], ["VOBr2", .9], ["VOCl2", .9], ["Bi2Te2S", 1.0],
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

    def to_cif(self):     # 输出crys.cif文件
        f = AseAtomsAdaptor.get_structure(self.crys)
        f1 = str(f) + '\n'
        j_pv_lst = findall('abc\s\s\s:(.*?)\n', f1)[0]   # abc   :  19.257300  19.569178  21.133988
        j1_pv_lst = j_pv_lst.split(' ')                   # abc   :   6.419100   6.523059   7.044663
        while '' in j1_pv_lst:
            j1_pv_lst.remove('')
        par_lst_matrix = self.crys.get_cell()
        # 物质种类（比如：Lu2 Al4）
        material = findall('Full\sFormula\s(.*?)\n', f1)[0].lstrip('(').rstrip(')')
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
            for letter in letter_lst:
                if letter in szb_lst:
                    num_lst.append(letter)
                if letter in zmb_lst:
                    symbol_lst.append(letter)
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
        par_lst_species = []                   # 用于Cifwrite参数(species)的
        for i, element in enumerate(element_lst):
            num = int(number_lst[i])
            for j in range(num):
                par_lst_species.append(element)
        # 每个原子的坐标
        ord_lst = []  # 最终Cifwriter所需要的coords参数
        ord_lst2 = []  # 储存的形式为
        for element in element_lst:
            ord_lst1 = findall(element + '\s\s\s\s(.*?)\n', f1)
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
        j_pv_lst = findall('abc\s\s\s:\s\s\s(.*?)\n', f1)[0]
        j1_pv_lst = j_pv_lst.split('   ')
        a = float(j1_pv_lst[0])
        b = float(j1_pv_lst[1])
        c = float(j1_pv_lst[2])
        par_lst_matrix = [[a, 0, 0],
                          [0, b, 0],
                          [0, 0, c]]
        # 物质种类（比如：Lu2 Al4）
        findall('Full\sFormula\s(.*?)\n', f1)[0]
        material = findall('Full\sFormula\s(.*?)\n', f1)[0].lstrip('(').rstrip(')')
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
        # 每个原子的坐标
        ord_lst = []   # 最终Cifwriter所需要的coords参数
        ord_lst2 = []  # 储存的形式为
        for element in element_lst:
            ord_lst1 = findall(element+'\s\s\s\s(.*?)\n',f1)
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
        slab = CifWriter(structure, write_magmoms=True)   # struct (Structure) – structure to write; symprec (float) – If not none, finds the symmetry of the structure and writes the cif with symmetry information. Passes symprec to the SpacegroupAnalyzer; write_magmoms (bool) – If True, will write magCIF file. Incompatible with symprec
        slab.write_file(r'C:\Users\wang1\Desktop\2D\download vasp\{}.cif'.format(self.mat_name))


class Database:
    info_data_base = ...
    def __init__(self, db_dir):
        '''
        Database类
        :param db_dir: 数据库路径
        '''
        print("dir = ", db_dir)
        self.path = db_dir
        self.conn = sqlite3.connect(self.path)
        self.cnn = self.conn.cursor()

    def set_up_properties(self, content):
        '''
        建立储存某种材料所有组合器件的物理性质的表(名字：Formula_Device)
        :param db_dir: 数据路径
        :param content: 材料的化学式
        :return: None
        '''
        cnn = self.conn.cursor()
        if content.isdigit():
            cnn.execute("SELECT Formula FROM test WHERE ID={}".format(int(content)))
            content = cnn.fetchone()[0]
        tbl_name = content + "_Device"
        cnn.execute('''CREATE TABLE {0}(
                       Splice_ID    INTEGER PRIMARY KEY, 
                       Splice                      TEXT,  
                       Device                      TEXT,
                       Hetero_junction             TEXT,
                       Optimal_Match               TEXT,
                       Layer                       TEXT,
                       "Binding_energy(eV)"        REAL,
                       "Schottky_barrier(eV)"      REAL,
                       Image_dir                   TEXT)'''.format(tbl_name))
        self.conn.commit()

    def insert_into_properties(self, content, splice_content):
        '''
        一对多实现（一种材料所有堆垛可能）
        :param content: 某种材料的ID或化学式
        :param splice_content: 用于堆垛的材料的 化学式或ID
        :return: None
        '''
        cnn = self.conn.cursor()
        if content.isdigit():
            cnn.execute("SELECT Formula FROM test WHERE ID={}".format(int(content)))
            r_formula = cnn.fetchone()[0]
        else:
            r_formula = content
        if splice_content.isdigit():
            cnn.execute("SELECT Formula FROM test WHERE ID={}".format(int(splice_content)))
            r_splice_formula = cnn.fetchone()[0]
            splice_ID = splice_content
        else:
            cnn.execute("SELECT Formula FROM test WHERE Formula={}".format(splice_content))
            splice_ID = cnn.fetchone()[0]
            r_splice_formula = splice_content
        tbl_name = r_formula + "_Device"
        print(tbl_name)
        splice_formula = "'" + r_splice_formula + "'"
        print('''INSERT INTO %s (Splice_ID, Splice) 
                 VALUES (%d, %s)''' % (tbl_name, int(splice_ID), splice_formula))
        cnn.execute('''INSERT INTO %s (Splice_ID, Splice) 
                       VALUES (%d, %s)''' % (tbl_name, int(splice_ID), splice_formula))
        self.conn.commit()



    def exhibit_all_info(self, fenkuai=False):
        '''
        展示数据库test表格中所有的信息
        :return: None
        '''

        self.cnn.execute("SELECT * FROM test")
        self.conn.commit()
        self.info_data_base = self.cnn.fetchall()
        if fenkuai == False:
            self.info_data_base.sort(key=itemgetter(1))     # 以ID排序，不是TrueID
        else:
            self.info_data_base.sort(key=itemgetter(1))     # 以ID排序，不是TrueID
            self.fenkuai_repeat()

    def fenkuai_repeat(self):
        lst = [self.info_data_base[0]]
        zong_lst = []
        for i in range(len(self.info_data_base) - 1):
            if self.info_data_base[i][1] == self.info_data_base[i + 1][1]:
                lst.append(self.info_data_base[i + 1])
            else:
                zong_lst.append(lst)
                lst = [self.info_data_base[i + 1]]
        zong_lst.append(lst)
        self.info_data_base = deepcopy(zong_lst)

    def exhibit_all_info2(self, content):
        '''
        显示某种材料所有堆垛后材料的性质（显示一对多所有情况）
        :param content: 材料（一） 的化学式或ID
        :return: None
        '''
        if content.isdigit():
            self.cnn.execute("SELECT Formula FROM test WHERE ID={0}".format(int(content)))
            formula = self.cnn.fetchone()[0]
        else:
            formula = content
        tbl_name = formula + "_Device"
        self.cnn.execute("SELECT * FROM {}".format(tbl_name))
        self.conn.commit()
        for x in self.cnn.fetchall():
            print(x)


    def insert_data(self, content, property, data, ID=True):
        '''
        向test表格中插入或修改数据
        :param content: 需要修改性质的材料的化学式或ID
        :param property: 相应物理性质
        :param data: 数据
        :return: None
        '''

        if ID == True:
            self.cnn.execute("UPDATE test SET {0} = {1} WHERE True_ID={2}".format(property, data, int(content)))
        self.conn.commit()
        print("edit", property, "to", data)

    def insert_data2(self, content, splice_content, property, data):
        '''
        在material_Device表中填入数据
        :param content: 材料的化学式或ID
        :param splice_content: 用于堆垛的材料的化学式或ID
        :param property: 需要更改那一项的名称
        :param data: 更改后的数据
        :return:
        '''
        if content.isdigit():
            self.cnn.execute("SELECT Formula FROM test WHERE ID={0}".format(int(content)))
            formula = self.cnn.fetchone()[0]
        else:
            formula = content
        tbl_name = formula + "_Device"
        if splice_content.isdigit():
            self.cnn.execute("UPDATE {0} SET {1} = {2} WHERE Splice_ID={3}".format(tbl_name, property, data, int(splice_content)))
        else:
            formula = "'" + splice_content + "'"
            self.cnn.execute("UPDATE {0} SET {1} = {2} WHERE Splice={3}".format(tbl_name, property, data, formula))
        self.conn.commit()

    def insert_row_into_test(self, True_ID, ID, formula):
        '''
        插入新的一行在test里
        :param ID: 插入材料的ID
        :param formula: 插入材料的化学式
        :return: None
        '''
        formula = "'" + formula + "'"
        self.cnn.execute('''INSERT INTO test (True_ID, ID, Formula)  
                            VALUES ({0}, {1}, {2})'''.format(True_ID, ID, formula))
        self.conn.commit()

    def delete_row(self, content, ID=True):
        '''
        用于删除某种材料在数据库中的记录
        :param content: 材料的化学式或者ID
        :return: None
        '''
        if ID==True:
            self.exhibit_all_info()
            print(self.info_data_base)
            content = int(content)
            for info_database_piece in self.info_data_base:
                if info_database_piece[0] == content:
                    formula = "'" + info_database_piece[2] + "'"
                    break
            print("formula = ", formula)
            self.cnn.execute("SELECT Formula from test WHERE True_ID={0}".format(int(content)))
            print("formula = ", formula)
            # self.cnn.execute("DELETE FROM test WHERE ID={0}".format(int(content)))
            self.cnn.execute("DELETE FROM test WHERE True_ID={0}".format(int(content)))

        else:
            formula = "'" + content + "'"
            self.cnn.execute("DROP TABLE {}".format(formula))
            self.cnn.execute("DELETE FROM test WHERE Formula={0}".format(formula))
        self.conn.commit()

    def set_up_properties(self, content):
        '''
        建立储存某种材料所有组合器件的物理性质的表(名字：Formula_Device)
        :param db_dir: 数据路径
        :param content: 材料的化学式
        :return: None
        '''
        cnn = self.conn.cursor()
        if content.isdigit():
            cnn.execute("SELECT Formula FROM test WHERE ID={}".format(int(content)))
            content = cnn.fetchone()[0]
        tbl_name = content + "_Device"
        cnn.execute('''CREATE TABLE {0}(
                          Splice_ID    INTEGER PRIMARY KEY, 
                          Splice                      TEXT,  
                          Device                      TEXT,
                          Hetero_junction             TEXT,
                          Optimal_Match               TEXT,
                          Layer                       TEXT,
                          "Binding_energy(eV)"        REAL,
                          "Schottky_barrier(eV)"      REAL,
                          Image_dir                   TEXT)'''.format(tbl_name))
        self.conn.commit()


    def close_db(self):
        '''
        关闭数据库
        :return: None
        '''
        self.conn.close()

class Thread(QtCore.QThread):
    def __init__(self):
        super(Thread, self).__init__()

class Myviewwidget(gl.GLViewWidget):      # 重写glviewwidget
    """ To rewrite the gl.GLViewwidget to achieve the drawrectangle function of self.view_widget. """
    num_mouse_track = 0
    signal_release = QtCore.pyqtSignal(tuple, tuple)
    set_font_size = ...
    def __init__(self, parent=None):
        super(Myviewwidget, self).__init__(parent)
        self.current_pos = None
        self.init_pos = None
        self.timer = QtCore.QTimer(self)
        self.timer.start(100)
        self.timer.timeout.connect(self.timeout_slot)
        global num
        num = 0


    def clear(self):
        """
        Remove all items from the scene.
        """
        for item in self.items:
            item._setView(None)
        self.items = []
        self.update()

    def setText_init(self, text1, text2):
        try:
            self.text1 = text1
            self.text2 = text2
            self.init = True
            self.update()
        except Exception as e:
            print(e)

    def setAtomstext(self, text_lst, pos_lst):
        self.text_lst = text_lst
        self.pos_lst = pos_lst
        self.init = False
        self.update()

    def setcoordinate_text(self, translate_pos, string, cell_array=None):
        self.translate_pos = translate_pos
        self.init = string
        self.cell_array = cell_array
        self.update()

    def paintGL(self, *args, **kwds):
        super().paintGL(*args, **kwds)
        if self.init is True:
            font1 = QtGui.QFont()
            font1.setPointSize(30)
            font1.setBold(True)
            font1.setWeight(75)
            font2 = QtGui.QFont()
            font2.setPointSize(20)
            font2.setBold(True)
            font2.setWeight(35)
            width = self.width() / 2
            height = self.height() / 2
            self.renderText(width - 100, height - 50, self.text1, font1)
            self.renderText(width - 170, height + 50, self.text2, font2)

        elif self.init is False:
            try:
                font2 = QtGui.QFont()
                font2.setBold(True)
                font2.setWeight(35)
                if self.set_font_size is None:
                    font2.setPointSize(20)
                else:
                    font2.setPointSize(self.set_font_size)
                for i, pos_piece in enumerate(self.pos_lst):
                    self.renderText(pos_piece[0], pos_piece[1], pos_piece[2],
                                    self.text_lst[i], font2)
            except Exception as e:
                print(e)
        elif self.init == "Cartestian":
            try:
                width = self.opts['distance'] * np.tan(self.opts['fov'] / 2 / 180 * np.pi)
                pos_lst = [list((np.array([0, 0, 0]) + self.translate_pos) / self.view_widget_orijin_objectwidth *
                                width),
                           list((np.array([2.5, 0, 0]) + self.translate_pos) / self.view_widget_orijin_objectwidth *
                                width),
                           list((np.array([0, 2.5, 0]) + self.translate_pos) / self.view_widget_orijin_objectwidth *
                                width),
                           list((np.array([0, 0, 2.5]) + self.translate_pos) / self.view_widget_orijin_objectwidth *
                                width)]
                screentuple_lst = []
                for pos_piece in pos_lst:
                    screentuple = list(gluProject(pos_piece[0], pos_piece[1], pos_piece[2]))
                    screentuple_lst.append(screentuple)
                width = screentuple_lst[0][0]
                height = screentuple_lst[0][1]
                worldtuple_lst = []
                for i, screentuple in enumerate(screentuple_lst):
                    screentuple[0] = screentuple[0] - width + 100
                    screentuple[1] = screentuple[1] - height + 100
                    worldtuple = gluUnProject(screentuple[0], screentuple[1], screentuple[2])
                    worldtuple_lst.append(list(worldtuple))
                orijin = worldtuple_lst[0]
                glBegin(GL_LINES)
                # glLineWidth(2.0)
                glColor4f(0, 0, 1, .6)  # z is blue
                glVertex3f(orijin[0], orijin[1], orijin[2])
                glVertex3f(worldtuple_lst[3][0], worldtuple_lst[3][1], worldtuple_lst[3][2])
                glColor4f(0, 1, 0, .6)  # y is green
                glVertex3f(orijin[0], orijin[1], orijin[2])
                glVertex3f(worldtuple_lst[2][0], worldtuple_lst[2][1], worldtuple_lst[2][2])
                glColor4f(1, 0, 0, .6)  # x is red
                glVertex3f(orijin[0], orijin[1], orijin[2])
                glVertex3f(worldtuple_lst[1][0], worldtuple_lst[1][1], worldtuple_lst[1][2])
                glEnd()
            except Exception as e:
                print(e)
        elif self.init == "Cell":
            try:
                width = self.opts['distance'] * np.tan(self.opts['fov'] / 2 / 180 * np.pi)
                pos_lst = [list((self.translate_pos) / self.view_widget_orijin_objectwidth *
                                width),
                           list((self.cell_array[0] / np.linalg.norm(self.cell_array[0]) * 2
                                 + self.translate_pos + np.array([.5, 0, 0])) / self.view_widget_orijin_objectwidth *
                                width),
                           list((self.cell_array[1] / np.linalg.norm(self.cell_array[1]) * 2
                                 + self.translate_pos + np.array([0, .5, 0])) / self.view_widget_orijin_objectwidth *
                                width),
                           list((self.cell_array[2] / np.linalg.norm(self.cell_array[2]) * 2
                                 + self.translate_pos + np.array([0, 0, .5])) / self.view_widget_orijin_objectwidth *
                                width)]
                screentuple_lst = []
                for pos_piece in pos_lst:
                    screentuple = list(gluProject(pos_piece[0], pos_piece[1], pos_piece[2]))
                    screentuple_lst.append(screentuple)
                width = screentuple_lst[0][0]
                height = screentuple_lst[0][1]
                worldtuple_lst = []
                for i, screentuple in enumerate(screentuple_lst):
                    screentuple[0] = screentuple[0] - width + 100
                    screentuple[1] = screentuple[1] - height + 100
                    worldtuple = gluUnProject(screentuple[0], screentuple[1], screentuple[2])
                    worldtuple_lst.append(list(worldtuple))
                orijin = worldtuple_lst[0]
                glBegin(GL_LINES)
                # glLineWidth(2.0)
                glColor4f(0, 0, 1, .6)  # z is blue
                glVertex3f(orijin[0], orijin[1], orijin[2])
                glVertex3f(worldtuple_lst[3][0], worldtuple_lst[3][1], worldtuple_lst[3][2])
                glColor4f(0, 1, 0, .6)  # y is green
                glVertex3f(orijin[0], orijin[1], orijin[2])
                glVertex3f(worldtuple_lst[2][0], worldtuple_lst[2][1], worldtuple_lst[2][2])
                glColor4f(1, 0, 0, .6)  # x is red
                glVertex3f(orijin[0], orijin[1], orijin[2])
                glVertex3f(worldtuple_lst[1][0], worldtuple_lst[1][1], worldtuple_lst[1][2])
                glEnd()
            except Exception as e:
                print(e)

    def mousePressEvent(self, event):
        try:
            if self.num_mouse_track == 1:   # drag rectangle
                self.init_pos = None
                self.current_pos = None
                if event.button() == 1:
                    self.init_pos = (event.pos().x(), event.pos().y())
                elif event.button() != 1:
                    self.init_pos = None
            elif self.num_mouse_track == 0:         # Normal
                return super().mousePressEvent(event)
            elif self.num_mouse_track == 2:
                if event.button() == 1:  # 鼠标左键按下
                    self.init_pos = (event.pos().x(), event.pos().y())
                    self.init_viewport = list(self.getViewport())
                elif event.button() != 1:  # 不是鼠标左键按下
                    self.init_pos = None
        except Exception as e:
            print(e)

    def mouseReleaseEvent(self, event):
        try:
            if self.num_mouse_track == 1:       # drag rectangle
                self.signal_release.emit(self.init_pos, self.current_pos)
                self.current_pos = deepcopy(self.init_pos)
                self.update()
                self.num_mouse_track = 0
                try:
                    self.repaint()
                except Exception as e:
                    print(e)
                self.num_mouse_track = 1
            elif self.num_mouse_track == 0:      # normal choose
                return super().mouseReleaseEvent(event)
            elif self.num_mouse_track == 2:        # translate viewport
                self.init_pos = None
                self.current_pos = None
                self.init_viewport = None
        except Exception as e:
            print(e)

    def mouseMoveEvent(self, event):
        try:
            if self.num_mouse_track == 1:       # drag rectangle
                self.current_pos = (event.pos().x(), event.pos().y())
                self.update()   # 调用paintevent
                self.num_mouse_track = 0
                self.repaint()
                self.num_mouse_track = 1
            elif self.num_mouse_track == 0:              # normal
                return super().mouseMoveEvent(event)
            elif self.num_mouse_track == 2:              # translate
                self.current_pos = (event.pos().x(), event.pos().y())
                vector = [self.current_pos[0] - self.init_pos[0],
                          self.init_pos[1] - self.current_pos[1]]
                self.set_new_viewport(vector)
        except Exception as e:
            print(e)

    def set_new_viewport(self, vector):
        try:
            width = self.width()
            height = self.height()
            viewport = [0, 0, 0, 0]
            viewport[0] = int(self.init_viewport[0] + vector[0])
            viewport[1] = int(self.init_viewport[1] + vector[1])
            viewport[2] = 2 * width - viewport[0]
            viewport[3] = 2 * height - viewport[1]
            self.opts['viewport'] = tuple(viewport)
            self.update()
        except Exception as e:
            print(e)

    def paintEvent(self, event):
        try:
            if self.num_mouse_track == 1:            # drag rectangle
                self.q = QtGui.QPainter(self)
                self.q.begin(self)
                self.drawRect(self.q)
                self.q.end()
            else:
                return super().paintEvent(event)
        except Exception as e:
            print(e)

    def drawRect(self, qp):
        try:
            if self.current_pos and self.init_pos:
                qp.setPen(QtGui.QColor(255, 255, 255))
                recwidth = self.current_pos[0] - self.init_pos[0]
                recheight = self.current_pos[1] - self.init_pos[1]
                self.num_mouse_track = 0
                self.repaint()
                qp.drawRect(self.init_pos[0], self.init_pos[1], recwidth, recheight)
                self.num_mouse_track = 1
        except Exception as e:
            print(e)

    def timeout_slot(self):
        """ To adjust the size of screen. """
        try:
            global num
            num += 1
            if num < 40:
                try:
                    height = self.height()
                    width = self.width()
                    self.opts['viewport'] = (-width, -height, 3 * width, 3 * height)
                except Exception as e:
                    print(e)
                self.view_widget_orijin_width = self.width()
                self.view_widget_orijin_height = self.height()
                self.view_widget_orijin_objectwidth = self.opts['distance'] * np.tan(self.opts['fov'] / 2 / 180 * np.pi)
        except:
            pass
        try:
            if self.num_mouse_track == 1:
                self.num_mouse_track = 0
                self.repaint()
                self.num_mouse_track = 1
        except:
            pass
        try:
            viewport = list(self.getViewport())
            width = self.width()
            height = self.height()
            viewport[0] = int(viewport[0] / self.view_widget_orijin_width * width)
            self.view_widget_orijin_width = width
            viewport[1] = int(viewport[1] / self.view_widget_orijin_height * height)
            self.view_widget_orijin_height = height
            viewport[2] = 2 * width - viewport[0]
            viewport[3] = 2 * height - viewport[1]
            self.opts['viewport'] = tuple(viewport)
        except Exception as e:
            print(e)


class Ui_MainWindow():  # 主窗口
    """ The layout of the mainwindow. """
    def setupUi(self, MainWindow):
        MainWindow.show()
        MainWindow.resize(2000, 1500)
        font = QtGui.QFont()
        font.setPointSize(9)
        MainWindow.setFont(font)
        # toolbar
        self.toolbarBox2 = QtGui.QToolBar()
        self.toolbarBox1 = QtGui.QToolBar()
        self.addToolBar(self.toolbarBox1)
        self.addToolBar(self.toolbarBox2)

        self.sonwindow = QtWidgets.QWidget()
        MainWindow.setCentralWidget(self.sonwindow)
        layout = QtWidgets.QHBoxLayout(self.sonwindow)
        splitter1 = QtWidgets.QSplitter(QtCore.Qt.Horizontal)    # 水平
        splitter2 = QtWidgets.QSplitter(QtCore.Qt.Vertical)     # 垂直
        # parametertree
        self.project_tree = QtWidgets.QTreeWidget()
        self.project_tree.setAnimated(True)
        self.project_tree.setDragDropMode(QtGui.QAbstractItemView.InternalMove)
        self.project_tree.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.project_tree.setColumnCount(4)
        self.project_tree.setHeaderLabels(['Files', 'Cell Formula', 'Heterostructure', 'Conductivity'])
        self.project_tree.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.project_tree.setStyleSheet("""
        QTreeView {
            background-color: rgb(35, 38, 41);
            alternate-background-color: rgb(200, 200, 200);
            color: rgb(238, 238, 238);
        }
        QLabel {
            color: rgb(238, 238, 238);
        }
        QTreeView::item:has-children {
            background-color: '#212627';
            color: rgb(233, 185, 110);
        }
        QTreeView::item:selected {
            background-color: rgb(92, 53, 102);
        }""")
        self.tab_tree = QtWidgets.QTabWidget()
        self.tab_tree.setFont(QtGui.QFont('Times New Roman', 10))

        self.tab_tree_project = QtWidgets.QWidget()
        self.tab_tree.addTab(self.tab_tree_project, "Project")
        tab_tree_widgetlayout = QtWidgets.QVBoxLayout(self.tab_tree_project)
        tab_tree_Filterlayout = QtWidgets.QHBoxLayout()
        # 给FileTab添加一个过滤器
        filterlabel1 = QtWidgets.QLabel("Filter from")
        self.filterLineEdit1 = QtWidgets.QLineEdit()
        self.filterLineEdit1.setMaximumWidth(200)
        self.filterLineEdit2 = QtWidgets.QLineEdit()
        self.filterLineEdit2.setMaximumWidth(200)
        filterlabel2 = QtWidgets.QLabel("to")
        self.filterbutton = QtWidgets.QPushButton('Ok')
        tab_tree_Filterlayout.addWidget(filterlabel1)
        tab_tree_Filterlayout.addWidget(self.filterLineEdit1)
        tab_tree_Filterlayout.addWidget(filterlabel2)
        tab_tree_Filterlayout.addWidget(self.filterLineEdit2)
        tab_tree_Filterlayout.addWidget(self.filterbutton)
        tab_tree_widgetlayout.addLayout(tab_tree_Filterlayout)
        tab_tree_widgetlayout.addWidget(self.project_tree)
        # 给projectTab设置一个搜索框
        searchlayout = QtWidgets.QHBoxLayout()
        self.searchLineEdit = QtWidgets.QLineEdit()
        self.searchLineEdit.setPlaceholderText("Search...")
        searchlayout.addWidget(self.searchLineEdit)
        tab_tree_widgetlayout.addLayout(searchlayout)
        # 左边控件 - object
        self.tab_tree_Object = QtWidgets.QWidget()
        self.tab_tree.addTab(self.tab_tree_Object, "Object")
        Verticallayout = QtWidgets.QVBoxLayout(self.tab_tree_Object)
        self.calculate_distance_button = QtWidgets.QPushButton('Distance')
        self.calculate_degree_button = QtWidgets.QPushButton('Degree')
        self.calculate_vector_button = QtWidgets.QPushButton('Vector')
        Horizonlayout = QtWidgets.QHBoxLayout()
        Horizonlayout.addWidget(self.calculate_distance_button)
        Horizonlayout.addWidget(self.calculate_vector_button)
        Horizonlayout.addWidget(self.calculate_degree_button)
        Verticallayout.addLayout(Horizonlayout)
        Splitter_object = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        self.object_tablewidget = QtWidgets.QTableWidget()
        self.object_tablewidget.setStyleSheet("color: rgb(205, 205, 205);background-color: rgb(35, 38, 41);")
        self.object_tablewidget.setColumnCount(3)
        self.object_tablewidget.setHorizontalHeaderLabels(["Atom", 'Color', 'Atom radii(Å)'])
        self.object_tablewidget.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.object_tablewidget.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.object_tablewidget.setColumnWidth(1, 70)
        Splitter_object.addWidget(self.object_tablewidget)
        self.objectTextwidget = QtWidgets.QTextEdit()
        self.objectTextwidget.setStyleSheet("color: rgb(205, 205, 205);background-color: rgb(35, 38, 41);")
        Splitter_object.addWidget(self.objectTextwidget)
        Verticallayout.addWidget(Splitter_object)
        self.ballfilling_button = QtWidgets.QPushButton('Space filling')
        self.ball_stick_button = QtWidgets.QPushButton('Ball-and-stick')
        self.stick_button = QtWidgets.QPushButton('Stick')
        Horizonlayout1 = QtWidgets.QHBoxLayout()
        Horizonlayout1.addWidget(self.ballfilling_button)
        Horizonlayout1.addWidget(self.ball_stick_button)
        Horizonlayout1.addWidget(self.stick_button)

        Verticallayout.addLayout(Horizonlayout1)

        splitter1.addWidget(self.tab_tree)
        ## 右边的控件
        self.view_widget = Myviewwidget()

        self.view_widget.set_font_size = 20
        self.view_widget.setMinimumSize(100, 100)
        self.view_widget.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        self.view_widget.setBackgroundColor((35, 38, 41, 0))
        self.widget_view_widget = QtWidgets.QWidget()
        view_widget_vertical_layout = QtWidgets.QVBoxLayout(self.widget_view_widget)
        view_widget_vertical_layout.addWidget(self.view_widget)
        self.splitterview_widget = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.splitterview_widget.addWidget(self.widget_view_widget)
        self.vertical_widget = QtWidgets.QWidget()
        database_vertical_layout = QtWidgets.QVBoxLayout(self.vertical_widget)
        self.datatable = QtWidgets.QTableWidget()
        self.datatable.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.datatable.setColumnCount(9)
        self.datatable.setHorizontalHeaderLabels(['ID', 'Formula', 'a(Å)', 'b(Å)', 'c(Å)',
                                                  'α(°)', 'β(°)', 'γ(°)', 'Volume'])
        self.show_image_table = QtWidgets.QTableWidget()
        self.show_image_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        viewwidget_tab = QtWidgets.QTabWidget()
        tab_view_widget_init = QtWidgets.QWidget()
        viewwidget_tab.setFont(QtGui.QFont('Times New Roman', 12))
        viewwidget_tab.addTab(tab_view_widget_init, "MatLego")
        tab_view_widgetlayout = QtWidgets.QVBoxLayout(tab_view_widget_init)
        # database 里的搜索框
        searchdatabase_layout = QtWidgets.QHBoxLayout()
        self.searchdatabase_LineEdit = QtWidgets.QLineEdit()
        self.searchdatabase_LineEdit.setPlaceholderText("Search...")
        searchdatabase_layout.addWidget(self.searchdatabase_LineEdit)
        # databaseradiobutton
        Hlayout0 = QtWidgets.QHBoxLayout()
        self.btn1 = QtWidgets.QRadioButton("Crystal")
        self.btn1.setFont(QtGui.QFont('Times New Roman', 10))
        self.btn1.setChecked(True)
        Hlayout0.addWidget(self.btn1)
        self.btn2 = QtWidgets.QRadioButton("Hetero-junction")
        self.btn2.setFont(QtGui.QFont('Times New Roman', 10))
        self.btn2.setChecked(False)
        Hlayout0.addWidget(self.btn2)
        self.btn3 = QtWidgets.QRadioButton("Device")
        self.btn3.setFont(QtGui.QFont('Times New Roman', 10))
        self.btn3.setChecked(False)
        Hlayout0.addWidget(self.btn3)
        database_vertical_layout.addLayout(Hlayout0)
        database_vertical_layout.addLayout(searchdatabase_layout)
        splitter_data_image = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        splitter_data_image.addWidget(self.datatable)
        splitter_data_image.addWidget(self.show_image_table)
        database_vertical_layout.addWidget(splitter_data_image)
        self.splitterview_widget.addWidget(self.vertical_widget)
        tab_view_widgetlayout.addWidget(self.splitterview_widget)
        splitter2.addWidget(viewwidget_tab)
        # 右下
        self.tab = QtWidgets.QTabWidget()
        self.Crystal_tab = QtWidgets.QWidget()
        self.Atom_tab = QtWidgets.QWidget()
        self.Hetero_device_tab = QtWidgets.QWidget()
        self.tab.addTab(self.Crystal_tab, "Crystal information")
        self.tab.addTab(self.Atom_tab, "Atom information")
        self.tab.addTab(self.Hetero_device_tab, "Hetero-junction / Device information")
        self.tab.setFont(QtGui.QFont('Times New Roman', 10))

        # tab1的布局
        self.Crystal_textwidget = QtWidgets.QTextEdit()
        self.Crystal_textwidget.setReadOnly(True)
        self.Crystal_textwidget.resize(24, 20)
        self.Crystal_textwidget.setStyleSheet("color: rgb(205, 205, 205);background-color: rgb(35, 38, 41);")
        tab1layout = QtWidgets.QVBoxLayout(self.Crystal_tab)
        tab1layout.addWidget(self.Crystal_textwidget)
        # tab2的布局
        self.Atom_textwidget = QtWidgets.QTextEdit()
        self.Atom_textwidget.setReadOnly(True)
        self.Atom_textwidget.resize(24, 20)
        self.Atom_textwidget.setStyleSheet("color: rgb(205, 205, 205);background-color: rgb(35, 38, 41);")
        tab2layout = QtWidgets.QVBoxLayout(self.Atom_tab)
        tab2layout.addWidget(self.Atom_textwidget)
        # tab3的布局
        self.Hetero_tablewidget = QtWidgets.QTableWidget()
        self.Hetero_tablewidget.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.Hetero_tablewidget.setStyleSheet("color: rgb(255, 255, 255);background-color: rgb(35, 38, 41);")
        self.Hetero_tablewidget.setColumnCount(5)
        self.Hetero_tablewidget.setHorizontalHeaderLabels(['Hetero-junction', 'Optimal Match', '#Layer',
                                                           'Binding  energy (eV/MoS2)', 'Schottky barrier (eV)'])
        self.Hetero_tablewidget.resizeColumnsToContents()
        self.Device_tablewidget = QtWidgets.QTableWidget()
        self.Device_tablewidget.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.Device_tablewidget.setStyleSheet("color: rgb(255, 255, 255);background-color: rgb(35, 38, 41);")
        self.Device_tablewidget.setColumnCount(5)
        self.Device_tablewidget.setHorizontalHeaderLabels(['Device', 'Optimal Match', '#Layer',
                                                           'Schottky barrier (eV)', 'I-V curve (Theory)'])
        self.Device_tablewidget.resizeColumnsToContents()
        tab3layout = QtWidgets.QVBoxLayout(self.Hetero_device_tab)
        tab3splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)    # 水平
        tab3splitter.addWidget(self.Hetero_tablewidget)
        tab3splitter.addWidget(self.Device_tablewidget)
        tab3layout.addWidget(tab3splitter)
        splitter2.addWidget(self.tab)
        splitter2.setStretchFactor(0, 5)
        splitter2.setStretchFactor(1, 2)       # 实现控件可以拉拽
        splitter1.addWidget(splitter2)
        splitter1.setStretchFactor(0, 2)
        splitter1.setStretchFactor(1, 8)
        # 以上均为主界面布局
        self.view_widget.setWindowTitle('view_cell')
        self.view_widget.setCameraPosition(distance=15, azimuth=-90)
        layout.addWidget(splitter1)
        self.menuBar = QtWidgets.QMenuBar(MainWindow)
        self.document_menu = QtWidgets.QMenu(self.menuBar)
        self.function_menu = QtWidgets.QMenu(self.menuBar)
        self.traversal_menu = QtWidgets.QMenu(self.menuBar)

        self.view_menu = QtWidgets.QMenu(self.menuBar)
        self.database_menu = QtWidgets.QMenu(self.menuBar)
        self.Setting_menu = QtWidgets.QMenu(self.menuBar)
        self.help_menu = QtWidgets.QMenu(self.menuBar)
        MainWindow.setMenuBar(self.menuBar)
        self.mainToolBar = QtWidgets.QToolBar(MainWindow)
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.mainToolBar)
        # New
        self.actionNew = QtWidgets.QAction(MainWindow)
        # self.actionNewToolBar = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\New.png"), "New", MainWindow)
        self.actionNewToolBar = QtWidgets.QAction(QtGui.QIcon("New.png"), "New", MainWindow)

        self.toolbarBox1.addAction(self.actionNewToolBar)
        # Open
        self.actionOpen = QtWidgets.QAction(MainWindow)
        # self.actionOpenToolBar = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\open1.png"), "Open", MainWindow)
        self.actionOpenToolBar = QtWidgets.QAction(QtGui.QIcon("open1.png"), "Open", MainWindow)
        self.toolbarBox1.addAction(self.actionOpenToolBar)
        self.actionExit_toolbar = QtWidgets.QAction(QtGui.QIcon("Exit4.png"), "Exit", MainWindow)
        self.toolbarBox1.addAction(self.actionExit_toolbar)
        # Save
        self.actionSave = QtWidgets.QAction(MainWindow)
        # exit
        self.actionExit = QtWidgets.QAction(MainWindow)
        # Cut
        self.actionCut = QtWidgets.QAction(MainWindow)
        self.actionTraversal_Cut = QtWidgets.QAction(MainWindow)
        # Cut ToolBar
        # self.actionCutToolBar = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\Cut7.png"), "Exfoliate", MainWindow)
        self.actionCutToolBar = QtWidgets.QAction(QtGui.QIcon("Cut7.png"), "Exfoliate", MainWindow)
        self.toolbarBox1.addAction(self.actionCutToolBar)
        # Judge3layer
        self.actionGuide = QtWidgets.QAction(MainWindow)
        self.action_setting_atom_par = QtWidgets.QAction(MainWindow)
        self.action_setting_displayed_par = QtWidgets.QAction(MainWindow)
        # Stack
        self.actionStack1 = QtWidgets.QAction(MainWindow)
        self.actionStack2 = QtWidgets.QAction(MainWindow)
        # self.actionStackToolbar = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\Stack2.png"), "Stack", MainWindow)
        self.actionStack2Toolbar = QtWidgets.QAction(QtGui.QIcon("Stack2.png"), "Stack(2)", MainWindow)
        self.toolbarBox1.addAction(self.actionStack2Toolbar)
        # Rotate
        self.actionRotate = QtWidgets.QAction(MainWindow)
        # drag_rectangle
        self.actiondrag_rectangle = QtWidgets.QAction(QtGui.QIcon("dragrectangle0.png"), "Drag rectangle", MainWindow)
        self.actionNormalchoose = QtWidgets.QAction(QtGui.QIcon("drag6.png"), "Drag", MainWindow)
        self.actiontranslate_choose = QtWidgets.QAction(QtGui.QIcon("drag6.png"), "Translate", MainWindow)
        self.actionSetlayer = QtWidgets.QAction(QtGui.QIcon("setlayer.png"), "Set layer", MainWindow)
        self.actionMovelayer = QtWidgets.QAction(QtGui.QIcon("movelayer.png"), "Move layer", MainWindow)
        self.actionleft_screen = QtWidgets.QAction(QtGui.QIcon("left.png"), "Left screen", MainWindow)
        self.actionright_screen = QtWidgets.QAction(QtGui.QIcon("right.png"), "Right screen", MainWindow)
        self.actionup_screen = QtWidgets.QAction(QtGui.QIcon("up.png"), "Up screen", MainWindow)
        self.actiondown_screen = QtWidgets.QAction(QtGui.QIcon("down.png"), "Down screen", MainWindow)
        self.actionmiddle_screen = QtWidgets.QAction(QtGui.QIcon("middle.png"), "Middle screen", MainWindow)

        self.toolbarBox1.addAction(self.actiondrag_rectangle)
        self.toolbarBox1.addAction(self.actionNormalchoose)
        self.toolbarBox1.addAction(self.actionSetlayer)
        self.toolbarBox1.addAction(self.actionMovelayer)
        self.toolbarBox2.addAction(self.actionleft_screen)
        self.toolbarBox2.addAction(self.actionright_screen)
        self.toolbarBox2.addAction(self.actionup_screen)
        self.toolbarBox2.addAction(self.actiondown_screen)
        self.toolbarBox2.addAction(self.actionmiddle_screen)
        self.toolbarBox2.addAction(self.actiontranslate_choose)

        self.actionTraversal_make_devices = QtWidgets.QAction(MainWindow)
        self.actionTraversal_rotate = QtWidgets.QAction(MainWindow)
        self.actionTraversal_Stack1 = QtWidgets.QAction(MainWindow)
        self.actionTraversal_Stack2 = QtWidgets.QAction(MainWindow)
        self.actiontwist_little_angle = QtWidgets.QAction(MainWindow)
        self.actionMulti_layer_opitmization = QtWidgets.QAction(MainWindow)
        self.actionClassification = QtWidgets.QAction(MainWindow)
        self.actionImport = QtWidgets.QAction(MainWindow)
        self.actionExport = QtWidgets.QAction(MainWindow)
        self.action_viewfrom_a = QtWidgets.QAction(MainWindow)
        self.action_viewfrom_b = QtWidgets.QAction(MainWindow)
        self.action_viewfrom_c = QtWidgets.QAction(MainWindow)
        # self.actionViewCell = QtWidgets.QAction(QtGui.QIcon("dot.png"), "View cell", MainWindow)
        self.actionViewCell = QtWidgets.QAction("View cell", MainWindow)
        self.actionDaytimemode = QtWidgets.QAction("Day time mode", MainWindow)
        self.actionNightmode = QtWidgets.QAction("Night mode", MainWindow)
        self.actionNightmode.setObjectName('Night')
        self.actionDaytimemode.setObjectName('Day')
        self.actionDaytimemode.setEnabled(False)

        self.actionNo_coordinate_system = QtWidgets.QAction(MainWindow)
        self.actionViewCoordinate_System = QtWidgets.QAction(MainWindow)
        self.actionViewCoordinate_cell = QtWidgets.QAction(MainWindow)
        self.action3D_coordinate = QtWidgets.QAction(MainWindow)
        self.actionFixed_coordinate = QtWidgets.QAction(MainWindow)
        self.actionshow_None = QtWidgets.QAction(MainWindow)
        self.actionshow_atom_index = QtWidgets.QAction(MainWindow)
        self.actionshow_atom_element = QtWidgets.QAction(MainWindow)
        # self.actionshow_coordinate_system = QtWidgets.QAction(MainWindow)

        self.actionViewGrid = QtWidgets.QAction(MainWindow)
        # view tool window
        self.actionViewDatabase = QtWidgets.QAction(MainWindow)
        self.actionViewProject = QtWidgets.QAction(MainWindow)
        self.actionViewText = QtWidgets.QAction(MainWindow)
        self.actionPerspective = QtWidgets.QAction(MainWindow)
        self.actionOrthogonal = QtWidgets.QAction(MainWindow)
        # document menu
        self.document_menu.addAction(self.actionNew)
        self.document_menu.addAction(self.actionOpen)
        self.document_menu.addSeparator()
        self.document_menu.addAction(self.actionSave)
        self.document_menu.addSeparator()
        self.document_menu.addAction(self.actionImport)
        self.document_menu.addAction(self.actionExport)
        self.document_menu.addSeparator()
        self.document_menu.addAction(self.actionExit)
        # function menu
        self.function_menu.addAction(self.actionCut)
        self.function_menu.addSeparator()
        self.function_menu.addAction(self.actionStack2)

        self.function_menu.addAction(self.actionStack1)
        self.function_menu.addAction(self.actionRotate)
        self.function_menu.addSeparator()
        self.function_menu.addAction(self.actiontwist_little_angle)
        self.function_menu.addAction(self.actionMulti_layer_opitmization)
        self.function_menu.addSeparator()
        self.function_menu.addAction(self.actionClassification)
        # Traversal-menu
        self.traversal_menu.addAction(self.actionTraversal_Cut)
        self.traversal_menu.addSeparator()
        self.traversal_menu.addAction(self.actionTraversal_Stack2)
        self.traversal_menu.addAction(self.actionTraversal_Stack1)
        self.traversal_menu.addAction(self.actionTraversal_rotate)
        self.traversal_menu.addSeparator()
        self.traversal_menu.addAction(self.actionTraversal_make_devices)
        # view menu
        self.toolwindow_menu = self.view_menu.addMenu("Tool Window")
        self.toolwindow_menu.addAction(self.actionViewDatabase)
        self.toolwindow_menu.addAction(self.actionViewProject)
        self.toolwindow_menu.addAction(self.actionViewText)
        #
        self.viewport_menu = self.view_menu.addMenu("Viewport")
        self.viewport_menu.addAction(self.actionPerspective)
        self.viewport_menu.addAction(self.actionOrthogonal)
        self.view_from = self.view_menu.addMenu("view from")
        self.view_from.addAction(self.action_viewfrom_a)
        self.view_from.addAction(self.action_viewfrom_b)
        self.view_from.addAction(self.action_viewfrom_c)
        self.view_menu.addSeparator()
        self.view_menu.addAction(self.actionViewCell)
        self.coordinate_system_menu = self.view_menu.addMenu("Coordinate System")
        self.coordinate_system_menu.addAction(self.actionNo_coordinate_system)
        self.coordinate_system_menu.addSeparator()
        self.coordinate_system_menu.addAction(self.actionViewCoordinate_System)
        self.coordinate_system_menu.addAction(self.actionViewCoordinate_cell)
        self.coordinate_system_menu.addSeparator()
        self.coordinate_system_menu.addAction(self.actionFixed_coordinate)
        self.coordinate_system_menu.addAction(self.action3D_coordinate)
        # self.action3D_coordinate.setEnabled(False)
        # self.actionFixed_coordinate.setEnabled(False)
        self.show_label_menu = self.view_menu.addMenu("Show label")
        self.show_label_menu.addAction(self.actionshow_None)
        self.show_label_menu.addSeparator()

        self.show_label_menu.addAction(self.actionshow_atom_index)
        self.show_label_menu.addAction(self.actionshow_atom_element)
        # self.show_label_menu.addAction(self.actionshow_coordinate_system)

        self.view_menu.addAction(self.actionViewGrid)
        self.view_menu.addSeparator()
        self.style_menu = self.view_menu.addMenu('Style')
        self.style_menu.addAction(self.actionDaytimemode)
        self.style_menu.addAction(self.actionNightmode)
        # help menu
        self.help_menu.addAction(self.actionGuide)
        # Setting menu
        self.Setting_menu.addAction(self.action_setting_atom_par)
        self.Setting_menu.addAction(self.action_setting_displayed_par)
        # database_menu
        self.action_db_connect = QtWidgets.QAction(MainWindow)
        self.action_add_data =  QtWidgets.QAction(MainWindow)
        self.action_cif_to_db = QtWidgets.QAction(MainWindow)
        self.database_menu.addAction(self.action_db_connect)
        self.database_menu.addAction(self.action_add_data)
        self.database_menu.addAction(self.action_cif_to_db)
        #
        self.menuBar.addAction(self.document_menu.menuAction())
        self.menuBar.addAction(self.function_menu.menuAction())
        self.menuBar.addAction(self.traversal_menu.menuAction())
        self.menuBar.addAction(self.view_menu.menuAction())
        self.menuBar.addAction(self.Setting_menu.menuAction())
        self.menuBar.addAction(self.database_menu.menuAction())
        self.menuBar.addAction(self.help_menu.menuAction())
        self.init_database()
        # calculate menu
        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        self.vertical_widget.setVisible(False)
        self.view_widget.setText_init("MatLego", "Lego of Materials")

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.document_menu.setTitle(_translate("MainWindow", "File"))
        self.function_menu.setTitle(_translate("MainWindow", "Function"))
        self.traversal_menu.setTitle(_translate("MainWindow", "Traversal"))
        self.view_menu.setTitle(_translate("MainWindow", "View"))
        self.database_menu.setTitle(_translate("MainWindow", "Database"))
        self.help_menu.setTitle(_translate("MainWindow", "Help"))
        self.Setting_menu.setTitle(_translate("MainWindow", "Settings"))
        # actionOpen
        self.actionNew.setText(_translate("MainWindow", "New"))
        self.actionNew.setShortcut('Ctrl+N')
        self.actionOpen.setText(_translate("MainWindow", "Open"))
        self.actionOpen.setShortcut('Ctrl+O')
        # self.actionOpen.setStatusTip("打开")  # 状态栏显示
        # actionsave
        self.actionSave.setText(_translate("MainWindow", "Save"))
        # actionExit
        self.actionExit.setText(_translate("MainWindow", "Exit"))
        # actionRotate
        self.actionRotate.setText(_translate("MainWindow", "Rotate"))
        # actionCut
        self.actionCut.setText(_translate("MainWindow", "Exfoliate"))
        # actionGuide
        self.actionGuide.setText(_translate("MainWindow", "Guide"))
        self.actionGuide.setEnabled(False)
        # Setting
        self.action_setting_atom_par.setText(_translate("MainWindow", "Atom parameter"))
        self.action_setting_displayed_par.setText(_translate("MainWindow", "Displayed parameters"))
        # actionStack
        self.actionStack2.setText(_translate("MainWindow", "Stack(2)"))

        self.actionStack1.setText(_translate("MainWindow", "Stack(1)"))
        # actionTraversal
        self.actionTraversal_rotate.setText(_translate("MainWindow", "Traversal-Rotate"))
        self.actionTraversal_Stack1.setText(_translate("MainWindow", "Traversal-Stack(1)"))
        self.actionTraversal_Stack2.setText(_translate("MainWindow", "Traversal-Stack(2)"))
        self.actionTraversal_make_devices.setText(_translate("MainWindow", "Traversal-Make divices"))
        # TraversalCut
        self.actionTraversal_Cut.setText(_translate("MainWindow", "Traversal-Exfoliate"))
        # actionTwistlittle_angle
        self.actiontwist_little_angle.setText(_translate("MainWindow", "Twist little angle"))
        self.actionMulti_layer_opitmization.setText(_translate("MainWindow", "Multi-layer stack"))
        # actionClassification
        self.actionClassification.setText(_translate("MainWindow", "Classification"))
        # actionImport
        self.actionImport.setText(_translate("MainWindow", "Import"))
        self.actionImport.setShortcut("Ctrl+I")
        # actionExport
        self.actionExport.setText(_translate("MainWindow", "Export"))
        self.actionViewProject.setText(_translate("MainWindow", "Project"))
        self.action_viewfrom_a.setText(_translate("MainWindow", "a-axis"))
        self.action_viewfrom_b.setText(_translate("MainWindow", "b-axis"))
        self.action_viewfrom_c.setText(_translate("MainWindow", "c-axis"))
        self.actionViewCell.setText(_translate("MainWindow", "Cell"))
        self.actionNo_coordinate_system.setText(_translate("MainWindow", "None"))
        self.actionViewCoordinate_System.setText(_translate("MainWindow", "Cartesian"))
        self.actionDaytimemode.setText(_translate("MainWindow", "Daytime mode"))
        self.actionNightmode.setText(_translate("MainWindow", "Night mode"))

        self.actionViewCoordinate_cell.setText(_translate("MainWindow", "Cell Coordinate"))
        self.action3D_coordinate.setText(_translate("MainWindow", "3D coordinate"))
        self.actionFixed_coordinate.setText(_translate("MainWindow", "Fixed coordinate"))
        self.actionshow_atom_index.setText(_translate("MainWindow", "Atom index"))
        self.actionshow_None.setText(_translate("MainWindow", "None"))
        self.actionshow_atom_element.setText(_translate("MainWindow", "Element symbol"))
        # self.actionshow_coordinate_system.setText(_translate("MainWindow", "Coordinate symbol"))
        self.actionViewGrid.setText(_translate("MainWindow", "Grid"))
        self.actionViewDatabase.setText(_translate("MainWindow", "Database"))
        self.actionViewText.setText(_translate("MainWindow", "Text"))
        self.actionOrthogonal.setText(_translate("MainWindow", "Parallel"))
        self.actionPerspective.setText(_translate("MainWindow", "Perspective"))
        self.action_db_connect.setText(_translate("MainWindow", "Connect"))
        self.action_add_data.setText(_translate("MainWindow", "Add CIF"))
        self.action_cif_to_db.setText(_translate("MainWindow", "Create a database"))

    def init_database(self):
        self.datatable.setRowCount(10)  # 添加信息
        self.datatable.setShowGrid(True)
        self.vertical_widget.setMaximumWidth(1200)

class create_database_window(QtWidgets.QWidget):
    Signal_emit_number = QtCore.pyqtSignal(float)
    def __init__(self, parent=None):
        super(create_database_window, self).__init__(parent)
        layoutH = QtWidgets.QHBoxLayout(self)
        progress_label = QtWidgets.QLabel("Progress:")
        self.pbar = QtWidgets.QProgressBar()
        self.pbar.setValue(15)
        layoutH.addWidget(progress_label)
        layoutH.addWidget(self.pbar)
        self.show()

    def pbar_number(self, value):
        self.pbar.setValue(value)

class Append_data_window(QtWidgets.QWidget):
    signal_emit = QtCore.pyqtSignal(int, str)
    def __init__(self, ID_lst, parent=None):
        super(Append_data_window, self).__init__(parent)
        self.ID_lst = ID_lst
        Hlayout = QtWidgets.QHBoxLayout(self)
        Id_label = QtWidgets.QLabel("ID:")
        self.ID_linedit = QtWidgets.QLineEdit()
        Formula_label = QtWidgets.QLabel("Formula:")
        self.Formula_linedit = QtWidgets.QLineEdit()
        Hlayout.addWidget(Id_label)
        Hlayout.addWidget(self.ID_linedit)
        Hlayout.addWidget(Formula_label)
        Hlayout.addWidget(self.Formula_linedit)
        self.button = QtWidgets.QPushButton("OK")
        self.button.clicked.connect(self.determine)
        Hlayout.addWidget(self.button)
        self.show()

    def determine(self):
        try:
            self.ID = int(self.ID_linedit.text())
            if self.ID not in self.ID_lst:
                self.Formula = self.Formula_linedit.text()
                self.signal_emit.emit(self.ID, self.Formula)
                self.close()
            else:
                QtWidgets.QMessageBox.warning(self, 'error',
                                              'ID already exists.')
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error',
                                          'Input error.')


class edit_data_table_window(QtWidgets.QWidget):
    signal_emit_information = QtCore.pyqtSignal(int, list, list)
    signal_emit_delete = QtCore.pyqtSignal(int)
    def __init__(self, ID_lst, header_lst, information_lst, row, parent=None):
        super(edit_data_table_window, self).__init__(parent)
        self.header_lst = header_lst
        self.information_lst = information_lst
        self.row = row
        self.ID_lst = ID_lst
        layout1 = QtWidgets.QHBoxLayout(self)
        self.edittable = QtWidgets.QTableWidget()
        layout1.addWidget(self.edittable)
        self.edittable.setColumnCount(len(self.header_lst))
        self.edittable.setHorizontalHeaderLabels(self.header_lst)
        self.edittable.setRowCount(1)
        self.edittable.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.edittable.showGrid()
        # self.edittable.setE
        if self.header_lst[-1] != 'I-V curve (Theory)':
            for i, info_piece in enumerate(self.information_lst):
                newitem = QtWidgets.QTableWidgetItem(info_piece)
                newitem.setTextAlignment(5 | 5)
                self.edittable.setItem(0, i, newitem)
        else:
            for i in range(len(self.information_lst) - 1):
                newitem = QtWidgets.QTableWidgetItem(self.information_lst[i])
                newitem.setTextAlignment(5 | 5)
                self.edittable.setItem(0, i, newitem)
            self.browsebutton = QtWidgets.QPushButton(self.information_lst[4])
            self.edittable.setCellWidget(0, 4, self.browsebutton)
            self.browsebutton.clicked.connect(self.browse)

            # newitem = QtWidgets.QTableWidgetItem()
        self.button1 = QtWidgets.QPushButton('Edit')
        self.button2 = QtWidgets.QPushButton('Delete')
        layout1.addWidget(self.button1)
        layout1.addWidget(self.button2)
        self.button1.clicked.connect(self.determine)
        self.button2.clicked.connect(self.delete)
        self.resize(1000, 130)
        self.show()

    def browse(self):
        try:
            image_dir, filetype = QtGui.QFileDialog.getOpenFileName(self, '选择文件', '', 'files(*.png , *.jpg)')
            self.edittable.removeCellWidget(0, 4)
            newitem = QtWidgets.QTableWidgetItem(self.information_lst[4] + '+' + image_dir)
            newitem.setTextAlignment(5 | 5)
            self.edittable.setItem(0, 4, newitem)
            # self.browsebutton.setText(image_dir)
        except Exception as e:
            print(e)

    def determine(self):
        try:
            if self.header_lst[-1] == 'Schottky_barrier(eV)':
                self.header_lst = ['Hetero_junction', 'Optimal_Match', 'Layer', 'Binding_energy', 'Schottky_barrier']
            elif self.header_lst[-1] == 'I-V curve (Theory)':
                self.header_lst = ['Device', 'Optimal_MatchD', 'LayerD', 'Schottky_barrierD', 'Image_dir']
            else:
                self.header_lst = ['ID', 'Formula', 'A', 'B', 'C', 'Angle_a', 'Angle_b', 'Angle_c', 'Volume']

            information_lst = []
            for i in range(self.edittable.columnCount()):
                try:
                    information_lst.append(self.edittable.item(0, i).text())
                except Exception as e:
                    print(e)
                    information_lst.append(self.browsebutton.text())
            if self.header_lst != ['ID', 'Formula', 'A', 'B', 'C', 'Angle_a', 'Angle_b', 'Angle_c', 'Volume']:
                self.signal_emit_information.emit(self.row, self.header_lst, information_lst)
                self.close()
            elif int(information_lst[0]) not in self.ID_lst:
                self.signal_emit_information.emit(self.row, self.header_lst, information_lst)
                self.close()
            else:
                QtWidgets.QMessageBox.warning(self, 'Forbidden',
                                              'ID already exists.')
        except Exception as e:
            print(e)

    def delete(self):
        reply = QtWidgets.QMessageBox.question(self, "Confirm delete", "Are you sure you want to delete this line?",
                                       QtWidgets.QMessageBox.Yes| QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.Yes)
        try:
            if reply == QtWidgets.QMessageBox.Yes:
                self.signal_emit_delete.emit(self.row)
                self.close()
        except Exception as e:
            print(e)

class edit_displayed_par_window(QtWidgets.QWidget):
    signal_emit_information = QtCore.pyqtSignal(int, int, int, int)
    signal_emit_fontsize = QtCore.pyqtSignal(int)
    signal_emit_fov_distance = QtCore.pyqtSignal(int, int)
    close_num = 0
    def __init__(self, max_number, fontsize, fov, distance, parent=None):
        super(edit_displayed_par_window, self).__init__(parent)
        self.fontsize = fontsize
        self.max_number = max_number
        self.fov = fov
        self.distance = distance
        Vlayout = QtWidgets.QVBoxLayout(self)
        Formlayout = QtWidgets.QFormLayout()
        button_Hlayout = QtWidgets.QHBoxLayout()
        font_size_Hlayout = QtWidgets.QHBoxLayout()
        max_atom_displaylabel = QtWidgets.QLabel("Maximum number of atoms displayed:")
        self.max_atom_lineeit = QtWidgets.QLineEdit()
        self.max_atom_lineeit.setPlaceholderText(str(self.max_number))
        self.max_atom_lineeit.setText(str(self.max_number))
        ok_button = QtWidgets.QPushButton("Ok")
        cancel_button = QtWidgets.QPushButton("Cancel")
        button_Hlayout.addWidget(ok_button)
        button_Hlayout.addWidget(cancel_button)
        ok_button.clicked.connect(self.determine)
        cancel_button.clicked.connect(self.close_)
        font_size_displaylabel = QtWidgets.QLabel("Atom index/Element font size:")
        self.fontsize_spinbox = QtWidgets.QSpinBox()
        self.fontsize_spinbox.setMinimum(1)
        self.fontsize_spinbox.setValue(self.fontsize)
        self.fontsize_spinbox.valueChanged.connect(self.setfont_size)
        font_size_Hlayout.addWidget(font_size_displaylabel)
        font_size_Hlayout.addWidget(self.fontsize_spinbox)
        fov_Hlayout = QtWidgets.QHBoxLayout()
        self.fov_slider_widget = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.fov_slider_widget.setRange(1, 179)  # 3
        self.fov_slider_widget.valueChanged.connect(self.on_change_func_fov)
        self.fov_slider_widget.setValue(self.fov)
        self.fov_label = QtWidgets.QLabel("FOV: {}°".format(self.fov))
        fov_Hlayout.addWidget(self.fov_label)
        fov_Hlayout.addWidget(self.fov_slider_widget)
        Formlayout.addRow(max_atom_displaylabel, self.max_atom_lineeit)
        Vlayout.addLayout(Formlayout)
        Vlayout.addLayout(font_size_Hlayout)
        Vlayout.addLayout(fov_Hlayout)
        Vlayout.addLayout(button_Hlayout)
        self.setWindowTitle("Set displayed parameter")
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        self.show()

    def closeEvent(self, event):
        try:
            print("in")
            self.close_()
        except Exception as e:
            print(e)

    def on_change_func_fov(self):
        try:
            fov = self.fov_slider_widget.value()
            self.fov_label.setText("FOV: {}°".format(fov))
            width = self.distance * np.tan(self.fov / 2 / 180 * np.pi)
            distance = width / np.tan(fov / 2 / 180 * np.pi)
            self.signal_emit_fov_distance.emit(fov, distance)
        except Exception as e:
            print(e)

    def close_(self):
        try:
            self.signal_emit_information.emit(self.max_number, self.fontsize, self.fov, self.distance)
            self.close()
        except Exception as e:
            print(e)

    def setfont_size(self):
        try:
            self.signal_emit_fontsize.emit(self.fontsize_spinbox.value())
        except Exception as e:
            print(e)

    def determine(self):
        try:
            self.max_number = int(float(self.max_atom_lineeit.text()))
            self.fov = self.fov_slider_widget.value()
            width = self.distance * np.tan(self.fov / 2 / 180 * np.pi)
            self.distance = width / np.tan(self.fov / 2 / 180 * np.pi)
            self.fontsize = self.fontsize_spinbox.value()
            self.signal_emit_information.emit(self.max_number, self.fontsize, self.fov, self.distance)
            self.close()
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error',
                                          'Input error!')
        else:
            QtWidgets.QMessageBox.information(self, 'Information:', 'Edit successfully')

class periodic_element_table_window(QtWidgets.QWidget):
    """ The periodic table of element. """
    signal_emit_atom_par = QtCore.pyqtSignal(list, list)
    def __init__(self, color_lst, radii_lst, parent=None):
        super(periodic_element_table_window, self).__init__(parent)
        self.color_lst = color_lst
        self.radii_lst = radii_lst
        self.setWindowTitle("Edit atom's parameter")
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        Vlayout = QtWidgets.QVBoxLayout(self)
        Gridlayout = QtWidgets.QGridLayout()
        Vlayout.addLayout(Gridlayout)
        Hlayout = QtWidgets.QHBoxLayout()
        Vlayout.addLayout(Hlayout)
        button_cancel = QtWidgets.QPushButton("Cancel")
        button_cancel.clicked.connect(self.close)
        button_ok = QtWidgets.QPushButton("OK")
        button_ok.clicked.connect(self.determine)
        Hlayout.addWidget(button_cancel)
        Hlayout.addWidget(button_ok)
        H_button = QtWidgets.QPushButton("H")
        H_button.setFont(QtGui.QFont('Microsoft YaHei', 10))
        H_button.setStyleSheet(
            "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
            "QPushButton{background-color:rgb(198,224,205)}"  # 按键背景色
            "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
            "QPushButton{border-radius:6px}"  # 圆角半径
            "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
        )
        H_button.clicked.connect(self.determine_element_formula)
        He_button = QtWidgets.QPushButton('He')
        He_button.setFont(QtGui.QFont('Microsoft YaHei', 10))
        He_button.setStyleSheet(
            "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
            "QPushButton{background-color:rgb(231, 239, 227)}"  # 按键背景色
            "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
            "QPushButton{border-radius:6px}"  # 圆角半径
            "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
        )
        He_button.clicked.connect(self.determine_element_formula)
        Li_button = QtWidgets.QPushButton('Li')
        Li_button.setFont(QtGui.QFont('Microsoft YaHei', 10))
        Li_button.setStyleSheet(
            "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
            "QPushButton{background-color:rgb(242, 213, 231)}"  # 按键背景色
            "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
            "QPushButton{border-radius:6px}"  # 圆角半径
            "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
        )

        Li_button.clicked.connect(self.determine_element_formula)
        Be_button = QtWidgets.QPushButton('Be')
        Be_button.setFont(QtGui.QFont('Microsoft YaHei', 10))
        Be_button.setStyleSheet(
            "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
            "QPushButton{background-color:rgb(182, 177, 239)}"  # 按键背景色
            "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
            "QPushButton{border-radius:6px}"  # 圆角半径
            "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
        )

        Be_button.clicked.connect(self.determine_element_formula)
        Na_button = QtWidgets.QPushButton('Na')
        Na_button.clicked.connect(self.determine_element_formula)
        Na_button.setFont(QtGui.QFont('Microsoft YaHei', 10))
        Na_button.setStyleSheet(
            "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
            "QPushButton{background-color:rgb(242, 213, 231)}"  # 按键背景色
            "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
            "QPushButton{border-radius:6px}"  # 圆角半径
            "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
        )
        Mg_button = QtWidgets.QPushButton('Mg')
        Mg_button.clicked.connect(self.determine_element_formula)
        Mg_button.setFont(QtGui.QFont('Microsoft YaHei', 10))
        Mg_button.setStyleSheet(
            "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
            "QPushButton{background-color:rgb(182, 177, 239)}"  # 按键背景色
            "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
            "QPushButton{border-radius:6px}"  # 圆角半径
            "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
        )
        Cs_button = QtWidgets.QPushButton('Cs')
        Cs_button.clicked.connect(self.determine_element_formula)
        Cs_button.setFont(QtGui.QFont('Microsoft YaHei', 10))

        Cs_button.setStyleSheet(
            "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
            "QPushButton{background-color:rgb(242, 213, 231)}"  # 按键背景色
            "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
            "QPushButton{border-radius:6px}"  # 圆角半径
            "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
        )
        Ba_button = QtWidgets.QPushButton('Ba')
        Ba_button.clicked.connect(self.determine_element_formula)
        Ba_button.setFont(QtGui.QFont('Microsoft YaHei', 10))
        Ba_button.setStyleSheet(
            "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
            "QPushButton{background-color:rgb(182, 177, 239)}"  # 按键背景色
            "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
            "QPushButton{border-radius:6px}"  # 圆角半径
            "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
        )
        Fr_button = QtWidgets.QPushButton('Fr')
        Fr_button.clicked.connect(self.determine_element_formula)
        Fr_button.setFont(QtGui.QFont('Microsoft YaHei', 10))

        Fr_button.setStyleSheet(
            "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
            "QPushButton{background-color:rgb(242, 213, 231)}"  # 按键背景色
            "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
            "QPushButton{border-radius:6px}"  # 圆角半径
            "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
        )
        Ra_button = QtWidgets.QPushButton('Ra')
        Ra_button.clicked.connect(self.determine_element_formula)
        Ra_button.setFont(QtGui.QFont('Microsoft YaHei', 10))
        Ra_button.setStyleSheet(
            "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
            "QPushButton{background-color:rgb(182, 177, 239)}"  # 按键背景色
            "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
            "QPushButton{border-radius:6px}"  # 圆角半径
            "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
        )

        La_Lu_button = QtWidgets.QPushButton('La-Lu')
        La_Lu_button.setFont(QtGui.QFont('Microsoft YaHei', 10))
        La_Lu_button.setStyleSheet(
            "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
            "QPushButton{background-color:rgb(252, 239, 121)}"  
            "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
            "QPushButton{border-radius:6px}"  # 圆角半径
            "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
        )
        La_Lu_button.clicked.connect(self.determine_element_formula)
        Ac_Lr_button = QtWidgets.QPushButton('Ac-Lr')
        Ac_Lr_button.setFont(QtGui.QFont('Microsoft YaHei', 10))
        Ac_Lr_button.setStyleSheet(
            "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
            "QPushButton{background-color:rgb(252, 239, 121)}"  
            "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
            "QPushButton{border-radius:6px}"  # 圆角半径
            "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
        )
        Ac_Lr_button.clicked.connect(self.determine_element_formula)
        Gridlayout.addWidget(H_button, 1, 1)
        Gridlayout.addWidget(He_button, 1, 18)
        Gridlayout.addWidget(Li_button, 2, 1)
        Gridlayout.addWidget(Be_button, 2, 2)
        Gridlayout.addWidget(Na_button, 3, 1)
        Gridlayout.addWidget(Mg_button, 3, 2)
        Gridlayout.addWidget(Cs_button, 6, 1)
        Gridlayout.addWidget(Ba_button, 6, 2)
        Gridlayout.addWidget(La_Lu_button, 6, 3)
        Gridlayout.addWidget(Fr_button, 7, 1)
        Gridlayout.addWidget(Ra_button, 7, 2)
        Gridlayout.addWidget(Ac_Lr_button, 7, 3)
        IA_label = QtWidgets.QLabel('IA')
        IA_label.setFont(QtGui.QFont('Times New Roman', 8))
        IA_label.setAlignment(QtCore.Qt.AlignCenter)
        VIIIA_label = QtWidgets.QLabel('VIIIA')
        VIIIA_label.setFont(QtGui.QFont('Times New Roman', 8))
        VIIIA_label.setAlignment(QtCore.Qt.AlignCenter)
        IIA_label = QtWidgets.QLabel('IIA')
        IIA_label.setFont(QtGui.QFont('Times New Roman', 8))
        IIA_label.setAlignment(QtCore.Qt.AlignCenter)
        IIIA_label = QtWidgets.QLabel('IIIA')
        IIIA_label.setFont(QtGui.QFont('Times New Roman', 8))
        IIIA_label.setAlignment(QtCore.Qt.AlignCenter)
        IVA_lable = QtWidgets.QLabel('IVA')
        IVA_lable.setFont(QtGui.QFont('Times New Roman', 8))
        IVA_lable.setAlignment(QtCore.Qt.AlignCenter)
        VA_label = QtWidgets.QLabel('VA')
        VA_label.setFont(QtGui.QFont('Times New Roman', 8))
        VA_label.setAlignment(QtCore.Qt.AlignCenter)
        VIA_label = QtWidgets.QLabel('VIA')
        VIA_label.setFont(QtGui.QFont('Times New Roman', 8))
        VIA_label.setAlignment(QtCore.Qt.AlignCenter)
        VIIA_label = QtWidgets.QLabel('VIIA')
        VIIA_label.setFont(QtGui.QFont('Times New Roman', 8))
        VIIA_label.setAlignment(QtCore.Qt.AlignCenter)
        IIIB_label = QtWidgets.QLabel('IIIB')
        IIIB_label.setFont(QtGui.QFont('Times New Roman', 8))
        IIIB_label.setAlignment(QtCore.Qt.AlignCenter)
        IVB_label = QtWidgets.QLabel('IVB')
        IVB_label.setFont(QtGui.QFont('Times New Roman', 8))
        IVB_label.setAlignment(QtCore.Qt.AlignCenter)
        VB_label = QtWidgets.QLabel('VB')
        VB_label.setFont(QtGui.QFont('Times New Roman', 8))
        VB_label.setAlignment(QtCore.Qt.AlignCenter)
        VIB_label = QtWidgets.QLabel('VIB')
        VIB_label.setFont(QtGui.QFont('Times New Roman', 8))
        VIB_label.setAlignment(QtCore.Qt.AlignCenter)
        VIIB_label = QtWidgets.QLabel('VIIB')
        VIIB_label.setFont(QtGui.QFont('Times New Roman', 8))
        VIIB_label.setAlignment(QtCore.Qt.AlignCenter)
        VIII_label = QtWidgets.QLabel('VIII')
        VIII_label.setFont(QtGui.QFont('Times New Roman', 8))
        VIII_label.setAlignment(QtCore.Qt.AlignCenter)
        IB_label = QtWidgets.QLabel('IB')
        IB_label.setFont(QtGui.QFont('Times New Roman', 8))
        IB_label.setAlignment(QtCore.Qt.AlignCenter)
        IIB_label = QtWidgets.QLabel('IIB')
        IIB_label.setFont(QtGui.QFont('Times New Roman', 8))
        IIB_label.setAlignment(QtCore.Qt.AlignCenter)
        Gridlayout.addWidget(IA_label, 0, 1)
        Gridlayout.addWidget(VIIIA_label, 0, 18)
        Gridlayout.addWidget(IIA_label, 1, 2)
        Gridlayout.addWidget(IIIA_label, 1, 13)
        Gridlayout.addWidget(IVA_lable, 1, 14)
        Gridlayout.addWidget(VA_label, 1, 15)
        Gridlayout.addWidget(VIA_label, 1, 16)
        Gridlayout.addWidget(VIIA_label, 1, 17)
        Gridlayout.addWidget(IIIB_label, 3, 3)
        Gridlayout.addWidget(IVB_label, 3, 4)
        Gridlayout.addWidget(VB_label, 3, 5)
        Gridlayout.addWidget(VIB_label, 3, 6)
        Gridlayout.addWidget(VIIB_label, 3, 7)
        Gridlayout.addWidget(VIII_label, 3, 8, 1, 3)
        Gridlayout.addWidget(IB_label, 3, 11)
        Gridlayout.addWidget(IIB_label, 3, 12)
        for i in range(7):
            label = QtWidgets.QLabel(str(i + 1))
            label.setAlignment(QtCore.Qt.AlignCenter)
            Gridlayout.addWidget(label, i + 1, 0)
        for i in range(5, 11):
            button = QtWidgets.QPushButton(crys_data.Element_formula[i])
            button.setFont(QtGui.QFont('Microsoft YaHei', 10))
            if i == 5:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(237, 140, 87)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif 6 <= i <= 8:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(179, 239, 151)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif i == 9:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(141, 224, 103)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            else:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(231, 239, 227)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            button.clicked.connect(self.determine_element_formula)
            Gridlayout.addWidget(button, 2, 8 + i)
        for i in range(13, 19):
            button = QtWidgets.QPushButton(crys_data.Element_formula[i])
            button.setFont(QtGui.QFont('Microsoft YaHei', 10))
            if i == 13:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(88, 181, 239)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif i == 14:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(237, 140, 87)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif 15 <= i <= 16:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(179, 239, 151)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif i == 17:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(141, 224, 103)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            else:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(231, 239, 227)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            button.clicked.connect(self.determine_element_formula)
            Gridlayout.addWidget(button, 3, i)
        for i in range(19, 37):
            button = QtWidgets.QPushButton(crys_data.Element_formula[i])
            button.setFont(QtGui.QFont('Microsoft YaHei', 10))
            if 21 <= i <= 30:
                button.setStyleSheet(
                "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                "QPushButton{background-color:rgb(159, 231, 249)}"  # 按键背景色(金属)
                "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                "QPushButton{border-radius:6px}"  # 圆角半径
                "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif i == 19:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(242, 213, 231)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif i == 20:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(182, 177, 239)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif i == 31:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(88, 181, 239)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif 32 <= i <= 33:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(237, 140, 87)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif i == 34:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(179, 239, 151)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif i == 35:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(141, 224, 103)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )

            else:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(231, 239, 227)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            button.clicked.connect(self.determine_element_formula)
            Gridlayout.addWidget(button, 4, i - 18)
        for i in range(37, 55):
            button = QtWidgets.QPushButton(crys_data.Element_formula[i])
            button.setFont(QtGui.QFont('Microsoft YaHei', 10))
            if 39 <= i <= 48:
                button.setStyleSheet(
                "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                "QPushButton{background-color:rgb(159, 231, 249)}"  # (金属)
                "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                "QPushButton{border-radius:6px}"  # 圆角半径
                "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif i == 37:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(242, 213, 231)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif i == 38:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(182, 177, 239)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )

            elif 49 <= i <= 50:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(88, 181, 239)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif 51 <= i <= 52:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(237, 140, 87)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif i == 53:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(141, 224, 103)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            else:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(231, 239, 227)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            button.clicked.connect(self.determine_element_formula)
            Gridlayout.addWidget(button, 5, i - 36)
        for i in range(72, 87):
            button = QtWidgets.QPushButton(crys_data.Element_formula[i])
            button.setFont(QtGui.QFont('Microsoft YaHei', 10))
            if i <= 80:
                button.setStyleSheet(
                "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                "QPushButton{background-color:rgb(159, 231, 249)}"  # 按键背景色(金属)
                "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                "QPushButton{border-radius:6px}"  # 圆角半径
                "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif 81 <= i <= 83:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(88, 181, 239)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif i == 84:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(237, 140, 87)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif i == 85:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(141, 224, 103)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            else:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(231, 239, 227)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            button.clicked.connect(self.determine_element_formula)
            Gridlayout.addWidget(button, 6, i - 68)
        for i in range(104, 119):
            button = QtWidgets.QPushButton(crys_data.Element_formula[i])
            button.setFont(QtGui.QFont('Microsoft YaHei', 10))
            if i <= 112:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(159, 231, 249)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif 113 <= i <= 116:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(88, 181, 239)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            elif i == 117:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(141, 224, 103)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )
            else:
                button.setStyleSheet(
                    "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                    "QPushButton{background-color:rgb(231, 239, 227)}"  # 按键背景色
                    "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                    "QPushButton{border-radius:6px}"  # 圆角半径
                    "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
                )

            button.clicked.connect(self.determine_element_formula)
            Gridlayout.addWidget(button, 7, i - 100)
        for i in range(57, 72):
            button = QtWidgets.QPushButton(crys_data.Element_formula[i])
            button.setFont(QtGui.QFont('Microsoft YaHei', 10))
            button.setStyleSheet(
                "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                "QPushButton{background-color:rgb(252, 239, 121)}"  # 按键背景色
                "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                "QPushButton{border-radius:6px}"  # 圆角半径
                "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
            )
            button.clicked.connect(self.determine_element_formula)
            Gridlayout.addWidget(button, 8, i - 53)
        for i in range(89, 104):
            button = QtWidgets.QPushButton(crys_data.Element_formula[i])
            button.setFont(QtGui.QFont('Microsoft YaHei', 10))
            button.setStyleSheet(
                "QPushButton{color:rgb(101,153,26)}"  # 按键前景色
                "QPushButton{background-color:rgb(252, 239, 121)}"  # 按键背景色
                "QPushButton:hover{color:blue}"  # 光标移动到上面后的前景色
                "QPushButton{border-radius:6px}"  # 圆角半径
                "QPushButton:pressed{background-color:rgb(180,180,180);border: None;}"  # 按下时的样式
            )
            button.clicked.connect(self.determine_element_formula)
            Gridlayout.addWidget(button, 9, i - 85)
        self.show()

    def determine_element_formula(self):
        try:
            sender = self.sender()
            print(sender.text())
            self.number = crys_data.Element_formula.index(sender.text())
            color = self.color_lst[self.number]
            colortrans = (color[0] * 255, color[1] * 255, color[2] * 255)
            radii = self.radii_lst[self.number]
            try:
                self.edit_information_window.close()
            except:
                pass
            self.edit_information_window = edit_atom_para(radii, colortrans)
            self.edit_information_window.signal_emit_change.connect(self.change)
        except Exception as e:
            print(e)

    def change(self, lst):
        try:
            Rgb = list(lst[1].getRgb())
            self.color_lst[self.number] = (Rgb[0]/255, Rgb[1] / 255, Rgb[2]/255, .5)
        except Exception as e:
            print(e)
        try:
            self.radii_lst[self.number] = lst[0]
        except Exception as e:
            print(e)

    def determine(self):
        try:
            self.signal_emit_atom_par.emit(self.color_lst, self.radii_lst)
            QtWidgets.QMessageBox.information(self, 'Information:', 'Edit successfully')
            self.close()
        except Exception as e:
            print(e)

class edit_atom_para(QtWidgets.QWidget):
    signal_emit_change = QtCore.pyqtSignal(list)
    def __init__(self, radii, color, parent=None):
        super(edit_atom_para, self).__init__(parent)
        self.radii = radii
        self.color = color
        self.params = [
            {'name': 'Basic parameter data types', 'type': 'group', 'children': [
                {'name': 'Atom radii', 'type': 'float', 'value': self.radii, 'step': .1},
                {'name': 'Color', 'type': 'color', 'value': self.color, 'tip': "This is a color button"}
            ]
             }]
        self.p = pg.parametertree.Parameter.create(name='params', type='group', children=self.params)
        self.p.sigTreeStateChanged.connect(self.change)
        self.t = pg.parametertree.ParameterTree()
        self.t.setParameters(self.p, showTop=False)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.t)
        for child in self.p.children():
            child.sigValueChanging.connect(self.valueChanging)
            for ch2 in child.children():
                ch2.sigValueChanging.connect(self.valueChanging)
        self.button_cancel = QtWidgets.QPushButton("Cancel")
        Hlayout = QtWidgets.QHBoxLayout()
        self.button_ok = QtWidgets.QPushButton("Ok")
        self.button_cancel.clicked.connect(self.close)
        self.button_ok.clicked.connect(self.determine)
        Hlayout.addWidget(self.button_cancel)
        Hlayout.addWidget(self.button_ok)
        layout.addLayout(Hlayout)
        self.setWindowTitle("Edit atom's parameter")
        self.setWindowIcon(QtGui.QIcon("Main.png"))

        self.show()
        self.resize(600, 400)

    def change(self, param, changes):
        try:
            print("tree changes:")
            for param, change, data in changes:
                path = self.p.childPath(param)
                if path is not None:
                    childName = '.'.join(path)
                else:
                    childName = param.name()

            if childName == "Basic parameter data types.Atom radii":
                self.radii = data
            else:
                self.color = data
        except Exception as e:
            print(e)

    def valueChanging(self, param, value):
        print("Value changing (not finalized): %s %s" % (param, value))

    def save(self):
        global state
        state = self.p.saveState()

    def restore(self):
        global state
        add = self.p['Save/Restore functionality', 'Restore State', 'Add missing items']
        rem = self.p['Save/Restore functionality', 'Restore State', 'Remove extra items']
        self.p.restoreState(state, addChildren=add, removeChildren=rem)

    def determine(self):
        l = [self.radii, self.color]
        self.signal_emit_change.emit(l)
        self.close()



class add_data_window(QtWidgets.QWidget):
    signal_emit_ID = QtCore.pyqtSignal(str, str, int)
    def __init__(self, cif_dir, database_dir, parent=None):
        super(add_data_window, self).__init__(parent)
        self.cif_dir = cif_dir
        self.database_dir = database_dir
        layout1 = QtWidgets.QHBoxLayout(self)
        label = QtWidgets.QLabel("ID：")
        self.ID_lineedit = QtWidgets.QLineEdit()
        self.button1 = QtWidgets.QPushButton('Ok')
        self.button1.clicked.connect(self.determine)
        layout1.addWidget(label)
        layout1.addWidget(self.ID_lineedit)
        layout1.addWidget(self.button1)
        self.show()

    def determine(self):
        try:
            ID = int(self.ID_lineedit.text())
            self.signal_emit_ID.emit(self.cif_dir, self.database_dir, ID)
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, 'error', 'ID input error!')


class database_to_cif(QtWidgets.QWidget):
    signal_emit_formula = QtCore.pyqtSignal(str, str, str)
    def __init__(self, database_dir, parent=None):
        super(database_to_cif, self).__init__(parent)
        self.database_dir = database_dir
        layout = QtWidgets.QHBoxLayout(self)
        label = QtWidgets.QLabel("Formula：")
        self.button = QtWidgets.QPushButton('Ok')
        self.button.clicked.connect(self.determine)
        self.formula_lineedit = QtWidgets.QLineEdit()
        layout.addWidget(label)
        layout.addWidget(self.formula_lineedit)
        layout.addWidget(self.button)
        self.show()

    def determine(self):
        cif_dir = QtGui.QFileDialog.getExistingDirectory(self, 'Select cif files directory')
        self.signal_emit_formula.emit(self.database_dir, cif_dir + '/' + str(self.formula_lineedit.text()) + '.cif',
                                      str(self.formula_lineedit.text()))
        self.close()


class create_2D_devices_window(QtWidgets.QWidget):
    """ The layout of create 2D devices"""
    signal_emit_information = QtCore.pyqtSignal(float, list, list, float, float, float, float, float, float, float)

    def __init__(self, parent=None):
        super(create_2D_devices_window, self).__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)  # 垂直布局
        layout0 = QtWidgets.QHBoxLayout()
        Tab = QtWidgets.QTabWidget()
        # Tab.setFont(QtGui.QFont('Times New Roman', 10))
        self.supercell_par_tab = QtWidgets.QWidget()
        self.error_limit_tab = QtWidgets.QWidget()
        Tab.addTab(self.error_limit_tab, "Error limit")
        Tab.addTab(self.supercell_par_tab, "Supercell parameter")
        normalstrain_label = QtWidgets.QLabel("Normal strain limit(%):")
        self.Normalstrain_lineedit = QtWidgets.QLineEdit()
        shearstrain_label = QtWidgets.QLabel("Shear strain limit(%):")
        self.Shearstrain_lineedit = QtWidgets.QLineEdit()
        areaerror_label = QtWidgets.QLabel("Area error limit(%):")
        self.areaerror_lineedit = QtWidgets.QLineEdit()
        Formlayout = QtWidgets.QFormLayout(self.error_limit_tab)
        Formlayout.addRow(normalstrain_label, self.Normalstrain_lineedit)
        Formlayout.addRow(shearstrain_label, self.Shearstrain_lineedit)
        Formlayout.addRow(areaerror_label, self.areaerror_lineedit)
        Gridlayout = QtWidgets.QGridLayout(self.supercell_par_tab)
        vacuummdis_label = QtWidgets.QLabel("Vacuum distance(Å):")
        self.vacuumdis_lineedit = QtWidgets.QLineEdit()
        maxlength_label = QtWidgets.QLabel("Max(a,b) length of supercell(Å):")
        self.maxlength_lineedit = QtWidgets.QLineEdit()
        aberror_label = QtWidgets.QLabel("Supercell a,b relative error(%)")
        self.aberror_lineedit = QtWidgets.QLineEdit()
        mingama_label = QtWidgets.QLabel("Minimum sin(γ)：")
        self.mingama_lineedit = QtWidgets.QLineEdit()
        self.mingama_lineedit.setPlaceholderText('0')
        maxgama_label = QtWidgets.QLabel("Maximum sin(γ)：")
        self.maxgama_lineedit = QtWidgets.QLineEdit()
        self.maxgama_lineedit.setPlaceholderText('1')
        Gridlayout.addWidget(vacuummdis_label, 0, 0, 1, 3)
        Gridlayout.addWidget(self.vacuumdis_lineedit, 0, 3, 1, 3)
        Gridlayout.addWidget(mingama_label, 0, 6, 1, 3)
        Gridlayout.addWidget(self.mingama_lineedit, 0, 9, 1, 3)
        Gridlayout.addWidget(maxlength_label, 2, 0, 1, 3)
        Gridlayout.addWidget(self.maxlength_lineedit, 2, 3, 1, 3)
        Gridlayout.addWidget(maxgama_label, 2, 6, 1, 3)
        Gridlayout.addWidget(self.maxgama_lineedit, 2, 9, 1, 3)
        Gridlayout.addWidget(aberror_label, 3, 0, 1, 3)
        Gridlayout.addWidget(self.aberror_lineedit, 3, 3, 1, 3)
        self.tablewidget = QtWidgets.QTableWidget()
        self.resize(700, 700)
        self.tablewidget.setColumnCount(2)
        self.tablewidget.setHorizontalHeaderLabels(['Conductivity (from bottom to top)', 'layer distance(Å)'])
        self.tablewidget.resizeColumnsToContents()
        self.text_widgt = QtWidgets.QTextEdit()
        self.text_widgt.setReadOnly(True)
        self.text_widgt.resize(24, 20)
        self.text_widgt.setStyleSheet("color: rgb(255, 255, 255);background-color: rgb(35, 38, 41);")
        Numlayer_button = QtWidgets.QPushButton("Ok")
        layernumber_label = QtWidgets.QLabel("Number of layers：")
        self.layernumber_Spinbox = QtWidgets.QSpinBox()
        self.layernumber_Spinbox.setMinimum(2)
        layout0.addWidget(layernumber_label)
        layout0.addWidget(self.layernumber_Spinbox)
        layout0.addWidget(Numlayer_button)
        Numlayer_button.clicked.connect(self.tablechange)
        layout2 = QtWidgets.QHBoxLayout()
        self.pbar = QtWidgets.QProgressBar()
        self.button = QtWidgets.QPushButton('Start')
        layout2.addWidget(self.pbar)
        layout2.addWidget(self.button)
        layout.addLayout(layout0)
        layout.addWidget(Tab)
        layout.addLayout(layout2)
        layout.addWidget(self.tablewidget)
        layout.addWidget(self.text_widgt)
        self.button.clicked.connect(self.determine)
        self.setWindowTitle("Make 2D devices")
        self.setWindowIcon(QtGui.QIcon("Main.png"))

        self.show()

    def tablechange(self):
        try:
            if self.layernumber_Spinbox.text() == "":
                self.tablewidget.setRowCount(0)
            else:
                num = int(self.layernumber_Spinbox.text())
                self.tablewidget.setRowCount(num)
                self.l = []
                label = QtWidgets.QLabel("0")
                self.tablewidget.setCellWidget(0, 1, label)
                for i in range(num):
                    comBox = QtWidgets.QComboBox()
                    comBox.addItem("Conductor")
                    comBox.addItem("Semiconductor")
                    comBox.addItem("Insulator")
                    self.tablewidget.setCellWidget(i, 0, comBox)
                    self.l.append(comBox)
        except Exception as e:
            print(e)

    def determine(self):
        try:
            vacuum_dis = float(self.vacuumdis_lineedit.text())
            conductivity_lst = []
            layer_dis_lst = []
            for i in range(self.tablewidget.rowCount()):
                conductivity_lst.append(self.l[i].currentText())
            for i in range(1, self.tablewidget.rowCount()):
                layer_dis = float(self.tablewidget.item(i, 1).text())
                layer_dis_lst.append(layer_dis)
            max_singama = float(self.maxgama_lineedit.text())
            min_singama = float(self.mingama_lineedit.text())
            ab_relative_error = float(self.aberror_lineedit.text()) / 100
            Normal_strain_limit = float(self.Normalstrain_lineedit.text()) / 100
            Shear_strain_limit = float(self.Shearstrain_lineedit.text()) / 100

            area_error = float(self.areaerror_lineedit.text()) / 100
            max_len = float(self.maxlength_lineedit.text())
            self.signal_emit_information.emit(vacuum_dis, layer_dis_lst, conductivity_lst, max_singama, min_singama,
                                              ab_relative_error, Normal_strain_limit, Shear_strain_limit,
                                              area_error, max_len)
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')
            print(e)


class tree_Edit_dic_window(QtWidgets.QWidget):
    """ The layout of the self.tree set dic window. """
    signal_edit_text = QtCore.pyqtSignal(str)
    def __init__(self, text, parent=None):
        super(tree_Edit_dic_window, self).__init__(parent)
        self.text = text
        layout0 = QtWidgets.QHBoxLayout(self)  # 水平布局
        label = QtWidgets.QLabel("Set Text ：")
        self.settext_lineedit = QtWidgets.QLineEdit(self.text)
        self.determine_button = QtWidgets.QPushButton('OK')
        layout0.addWidget(label)
        layout0.addWidget(self.settext_lineedit)
        layout0.addWidget(self.determine_button)
        self.determine_button.clicked.connect(self.determine)
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        self.setWindowTitle("Edit text")
        self.show()

    def determine(self):
        try:
            text = self.settext_lineedit.text()
            if text != "":
                self.signal_edit_text.emit(text)
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')
            self.close()
        except Exception as e:
            print(e)


class Create_folder_window(QtWidgets.QWidget):
    """The window of create a folder"""
    signal_emit_information = QtCore.pyqtSignal(str)
    def __init__(self, parent=None):
        super(Create_folder_window, self).__init__(parent)
        self.setWindowTitle("Create a folder")
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        Hlayout = QtWidgets.QHBoxLayout(self)
        label = QtWidgets.QLabel("Name:")
        self.name_lineedit = QtWidgets.QLineEdit()
        self.button = QtWidgets.QPushButton("OK")
        self.button.clicked.connect(self.determine)
        Hlayout.addWidget(label)
        Hlayout.addWidget(self.name_lineedit)
        Hlayout.addWidget(self.button)
        self.show()

    def determine(self):
        try:
            name = self.name_lineedit.text()
            if name != '':
                self.signal_emit_information.emit(name)
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')
            self.close()
        except:
            pass

class tree_set_text_window(QtWidgets.QWidget):
    """ The layout of the self.tree set text window. """
    signal_edit_text = QtCore.pyqtSignal(str)
    signal_edit_text1 = QtCore.pyqtSignal(str)
    def __init__(self, parent=None):
        super(tree_set_text_window, self).__init__(parent)
        layout0 = QtWidgets.QHBoxLayout(self)  # 水平布局
        label = QtWidgets.QLabel("Set Text ：")
        self.treetext = QtWidgets.QComboBox()
        self.treetext.addItem("bulk")
        self.treetext.addItem("layer")
        self.treetext.addItem("stack")
        settext_label = QtWidgets.QLabel("Set Text ：")
        self.treetext1 = QtWidgets.QComboBox()
        self.treetext1.addItem("Conductor")
        self.treetext1.addItem("Semiconductor")
        self.treetext1.addItem("Insulator")
        self.treetext1.addItem("")
        self.determine_button = QtWidgets.QPushButton('OK')
        self.determine_button1 = QtWidgets.QPushButton('OK')
        layout0.addWidget(label)
        layout0.addWidget(self.treetext)
        layout0.addWidget(self.determine_button)
        layout0.addWidget(settext_label)
        layout0.addWidget(self.treetext1)
        layout0.addWidget(self.determine_button1)
        self.determine_button.clicked.connect(self.determine)
        self.determine_button1.clicked.connect(self.determine1)
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        self.setWindowTitle("Edit text")
        self.show()

    def determine1(self):
        try:
            text = self.treetext1.currentText()
            self.signal_edit_text1.emit(text)
        except Exception as e:
            print(e)

    def determine(self):
        try:
            text = self.treetext.currentText()
            self.signal_edit_text.emit(text)
        except Exception as e:
            print(e)

class Stack_main_one_window(QtWidgets.QWidget):
    signal_emit = QtCore.pyqtSignal(float, float, float, float)
    def __init__(self, parent=None):
        super(Stack_main_one_window, self).__init__(parent)
        Vlayout = QtWidgets.QVBoxLayout(self)
        self.setWindowTitle("A stacking direction is known(a1, a2).")
        self.setWindowIcon(QtGui.QIcon("Main.png"))

        Tab = QtWidgets.QTabWidget()
        self.supercell_par_tab = QtWidgets.QWidget()
        self.error_limit_tab = QtWidgets.QWidget()
        Tab.addTab(self.error_limit_tab, "Error limit")
        Tab.addTab(self.supercell_par_tab, "Supercell parameter")
        error_limit_Formlayout = QtWidgets.QFormLayout(self.error_limit_tab)
        supercell_par_Formlayout = QtWidgets.QFormLayout(self.supercell_par_tab)
        Normalstrain_label = QtWidgets.QLabel("Normal strain limit(%)：")
        self.Normalstrain_lineedit = QtWidgets.QLineEdit()
        shear_strain_label = QtWidgets.QLabel("Shear strain limit(%)：")
        self.Shearstrain_lineedit = QtWidgets.QLineEdit()
        error_limit_Formlayout.addRow(Normalstrain_label, self.Normalstrain_lineedit)
        error_limit_Formlayout.addRow(shear_strain_label, self.Shearstrain_lineedit)
        layer_dis_label = QtWidgets.QLabel("layer distance(Å)：")
        self.layerdis_lineedit = QtWidgets.QLineEdit()
        vaccum_dis_label = QtWidgets.QLabel("Vacuum distance(Å)：")
        self.vacuumdis_lineedit = QtWidgets.QLineEdit()
        supercell_par_Formlayout.addRow(layer_dis_label, self.layerdis_lineedit)
        supercell_par_Formlayout.addRow(vaccum_dis_label, self.vacuumdis_lineedit)
        self.Start_button = QtWidgets.QPushButton('Start')
        self.Start_button.clicked.connect(self.determine)
        self.Cancel_button = QtWidgets.QPushButton('Cancel')
        self.Cancel_button.clicked.connect(self.close)
        button_Hlayout = QtWidgets.QHBoxLayout()
        button_Hlayout.addWidget(self.Start_button)
        button_Hlayout.addWidget(self.Cancel_button)
        Vlayout.addWidget(Tab)
        Vlayout.addLayout(button_Hlayout)
        self.setWindowTitle("Stack (a1 and a2 have same orientation)")
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        self.show()

    def determine(self):
        try:
            self.Normal_strain_limit_var = float(self.Normalstrain_lineedit.text()) / 100
            self.Shear_strain_limit_var = float(self.Shearstrain_lineedit.text()) / 100
            self.vacuum_distance_var = float(self.vacuumdis_lineedit.text())
            self.layer_distance_var = float(self.layerdis_lineedit.text())
            self.signal_emit.emit(self.Normal_strain_limit_var, self.Shear_strain_limit_var, self.layer_distance_var, self.vacuum_distance_var)
            self.close()
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')
            print(e)

class traversal_stack_window_1(QtWidgets.QWidget):
    """ The layout of the traversal_stack window. """
    signaltraversalstack = QtCore.pyqtSignal(int)
    signal_start_pbar = QtCore.pyqtSignal(float)
    layerdistance_var = ...
    vacuumdistance_var = ...
    Normal_strain_limit_var = ...
    Shear_strain_limit_var = ...

    def __init__(self, parent=None):
        super(traversal_stack_window_1, self).__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)
        Tab = QtWidgets.QTabWidget()
        self.supercell_par_tab = QtWidgets.QWidget()
        self.error_limit_tab = QtWidgets.QWidget()
        Tab.addTab(self.error_limit_tab, "Error limit")
        Tab.addTab(self.supercell_par_tab, "Supercell parameter")
        error_Formlayout = QtWidgets.QFormLayout(self.error_limit_tab)
        supercell_par_formlayout = QtWidgets.QFormLayout(self.supercell_par_tab)

        Hlayout1 = QtWidgets.QHBoxLayout()
        vacuumdis_label = QtWidgets.QLabel("Vacuum distance(Å)：")
        self.vacuumdis_lineedit = QtWidgets.QLineEdit()
        layerdis_label = QtWidgets.QLabel("Layer distance(Å)：")
        self.layerdis_lineedit = QtWidgets.QLineEdit()
        Normalstrain_label = QtWidgets.QLabel("Normal strain limit(%)：")
        self.Normalstrain_lineedit = QtWidgets.QLineEdit()
        Shearstrain_label = QtWidgets.QLabel("Shear strain limit(%)：")
        self.Shearstrain_lineedit = QtWidgets.QLineEdit()
        self.determine_button = QtWidgets.QPushButton('Start')
        self.determine_button.clicked.connect(self.determine)
        self.pbar = QtWidgets.QProgressBar()
        self.Text = QtWidgets.QTextEdit()
        supercell_par_formlayout.addRow(vacuumdis_label, self.vacuumdis_lineedit)
        supercell_par_formlayout.addRow(layerdis_label, self.layerdis_lineedit)
        error_Formlayout.addRow(Normalstrain_label, self.Normalstrain_lineedit)
        error_Formlayout.addRow(Shearstrain_label, self.Shearstrain_lineedit)
        Hlayout1.addWidget(self.pbar)
        Hlayout1.addWidget(self.determine_button)
        layout.addWidget(Tab)
        layout.addLayout(Hlayout1)
        layout.addWidget(self.Text)
        self.signal_start_pbar.connect(self.start_pbar)
        self.setWindowTitle("Traversal stack")
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        self.show()

    def start_pbar(self, num):
        self.pbar.setValue(num)


    def determine(self):
        try:
            self.layerdistance_var = float(self.layerdis_lineedit.text())
            self.vacuumdistance_var = float(self.vacuumdis_lineedit.text())
            self.Normal_strain_limit_var = float(self.Normalstrain_lineedit.text())
            self.Shear_strain_limit_var = float(self.Shearstrain_lineedit.text())
            self.signaltraversalstack.emit(1)

        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')



class traversal_stack_window_2(QtWidgets.QWidget):
    """ The layout of the traversal_stack window. """
    signaltraversalstack = QtCore.pyqtSignal(int)
    signal_start_pbar = QtCore.pyqtSignal(float)
    vacuumdistance_var = ...
    layerdistance_var = ...
    Normal_strain_limit_var = ...
    def __init__(self, parent=None):
        super(traversal_stack_window_2, self).__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)

        Tab = QtWidgets.QTabWidget()
        self.supercell_par_tab = QtWidgets.QWidget()
        self.error_limit_tab = QtWidgets.QWidget()
        Tab.addTab(self.error_limit_tab, "Error limit")
        Tab.addTab(self.supercell_par_tab, "Supercell parameter")

        supercell_par_Formlayout = QtWidgets.QFormLayout(self.supercell_par_tab)
        error_limit_Formlayout = QtWidgets.QFormLayout(self.error_limit_tab)

        Hlayout2 = QtWidgets.QHBoxLayout()
        Hlayout3 = QtWidgets.QHBoxLayout()

        vacuumdis_label = QtWidgets.QLabel("Vacuum distance(Å)：")
        self.vacuumdis_lineedit = QtWidgets.QLineEdit()
        layerdis_label = QtWidgets.QLabel("Layer distance(Å)：")
        self.layerdis_lineedit = QtWidgets.QLineEdit()
        supercell_par_Formlayout.addRow(vacuumdis_label, self.vacuumdis_lineedit)
        supercell_par_Formlayout.addRow(layerdis_label, self.layerdis_lineedit)
        Normalstrain_label = QtWidgets.QLabel("Normal strain limit(%)：")
        self.Normalstrain_lineedit = QtWidgets.QLineEdit()
        error_limit_Formlayout.addRow(Normalstrain_label, self.Normalstrain_lineedit)
        self.determine_button = QtWidgets.QPushButton('Start')
        self.determine_button.clicked.connect(self.determine)
        self.pbar = QtWidgets.QProgressBar()
        self.Text = QtWidgets.QTextEdit()
        Hlayout2.addWidget(self.pbar)
        Hlayout2.addWidget(self.determine_button)
        Hlayout3.addWidget(self.Text)
        layout.addWidget(Tab)
        layout.addLayout(Hlayout2)
        layout.addLayout(Hlayout3)
        self.signal_start_pbar.connect(self.start_pbar)
        self.setWindowTitle("Traversal stack")
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        self.show()

    def start_pbar(self, num):
        self.pbar.setValue(num)

    def determine(self):
        try:
            self.layerdistance_var = float(self.layerdis_lineedit.text())
            self.vacuumdistance_var = float(self.vacuumdis_lineedit.text())
            self.Normal_strain_limit_var = float(self.Normalstrain_lineedit.text())
            self.signaltraversalstack.emit(2)
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')


class traversal_rotate_window(QtWidgets.QWidget):
    """ The layout of the traversal_stack window. """
    signaltraversalstack = QtCore.pyqtSignal()
    signal_start_pbar = QtCore.pyqtSignal(float)
    def __init__(self, parent=None):
        super(traversal_rotate_window, self).__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)
        Hlayout4 = QtWidgets.QHBoxLayout()  # 水平布局

        Tab = QtWidgets.QTabWidget()
        # Tab.setFont(QtGui.QFont('Times New Roman', 10))
        self.supercell_par_tab = QtWidgets.QWidget()
        self.error_limit_tab = QtWidgets.QWidget()
        Tab.addTab(self.error_limit_tab, "Error limit")
        Tab.addTab(self.supercell_par_tab, "Supercell parameter")
        vacuumdis_label = QtWidgets.QLabel("Vacuum distance(Å)：")
        self.vacuumdis_lineedit = QtWidgets.QLineEdit()
        Maxlength_label = QtWidgets.QLabel("Max(a, b) length of supercell(Å): ")
        self.maxlength_lineedit = QtWidgets.QLineEdit()
        mingama_label = QtWidgets.QLabel("Minimum sin(γ)：")
        self.mingama_lineedit = QtWidgets.QLineEdit()
        self.mingama_lineedit.setPlaceholderText('0')
        Gridlayout = QtWidgets.QGridLayout(self.supercell_par_tab)
        Gridlayout.addWidget(vacuumdis_label, 0, 0, 1, 3)
        Gridlayout.addWidget(self.vacuumdis_lineedit, 0, 3, 1, 3)
        Gridlayout.addWidget(Maxlength_label, 1, 0, 1, 3)
        Gridlayout.addWidget(self.maxlength_lineedit, 1, 3, 1, 3)
        Gridlayout.addWidget(mingama_label, 1, 6, 1, 3)
        Gridlayout.addWidget(self.mingama_lineedit, 1, 9, 1, 3)
        layerdis_label = QtWidgets.QLabel("Layer distance(Å)：")
        self.layerdis_lineedit = QtWidgets.QLineEdit()
        aberror_label = QtWidgets.QLabel("Supercell a,b relative error(%)：")
        self.aberror_lineedit = QtWidgets.QLineEdit()
        maxgama_label = QtWidgets.QLabel("Maximum sin(γ)：")
        self.maxgama_lineedit = QtWidgets.QLineEdit()
        self.maxgama_lineedit.setPlaceholderText('1')
        Gridlayout.addWidget(layerdis_label, 0, 6, 1, 3)
        Gridlayout.addWidget(self.layerdis_lineedit, 0, 9, 1, 3)
        Gridlayout.addWidget(aberror_label, 2, 0, 1, 3)
        Gridlayout.addWidget(self.aberror_lineedit, 2, 3, 1, 3)
        Gridlayout.addWidget(maxgama_label, 2, 6, 1, 3)
        Gridlayout.addWidget(self.maxgama_lineedit, 2, 9, 1, 3)
        Normalstrain_label = QtWidgets.QLabel("Normal strain limit(%)：")
        self.normalstrain_lineedit = QtWidgets.QLineEdit()
        shearstrain_label = QtWidgets.QLabel("Shear strain limit(%)：")
        self.shearstrain_lineedit = QtWidgets.QLineEdit()
        arealimit_label = QtWidgets.QLabel("Area error limit(%)：")
        self.areaerror_lineedit = QtWidgets.QLineEdit()
        Formlayout = QtWidgets.QFormLayout(self.error_limit_tab)
        Formlayout.addRow(Normalstrain_label, self.normalstrain_lineedit)
        Formlayout.addRow(shearstrain_label, self.shearstrain_lineedit)
        Formlayout.addRow(arealimit_label, self.areaerror_lineedit)
        self.determine_button = QtWidgets.QPushButton('Start')
        self.determine_button.clicked.connect(self.determine)
        self.pbar = QtWidgets.QProgressBar()
        self.Text = QtWidgets.QTextEdit()
        Hlayout4.addWidget(self.pbar)
        Hlayout4.addWidget(self.determine_button)
        layout.addWidget(Tab)
        layout.addLayout(Hlayout4)
        layout.addWidget(self.Text)
        self.signal_start_pbar.connect(self.start_pbar)
        self.setWindowTitle("Traversal-Rotate")
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        self.show()

    def start_pbar(self, num):
        self.pbar.setValue(num)

    def determine(self):
        try:
            self.layerdistance_var = float(self.layerdis_lineedit.text())
            self.vacuumdistance_var = float(self.vacuumdis_lineedit.text())
            self.normal_strain_var = float(self.normalstrain_lineedit.text())
            self.shear_strain_var = float(self.shearstrain_lineedit.text())
            self.max_length_var = float(self.maxlength_lineedit.text())
            self.area_tolerance_var = float(self.areaerror_lineedit.text())
            self.minimum_gamma_var = float(self.mingama_lineedit.text())
            self.maximum_gamma_var = float(self.maxgama_lineedit.text())
            self.ab_relative_error_var = float(self.aberror_lineedit.text())
            if self.layerdistance_var > 0 and self.vacuumdistance_var >=0 and 100 > self.normal_strain_var >= 0 \
                    and self.max_length_var > 0 and 100 > self.area_tolerance_var >= 0 and \
                    1 >= self.maximum_gamma_var >= self.minimum_gamma_var >= 0 and self.ab_relative_error_var > 0:

                self.signaltraversalstack.emit()
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')



class scatter_plot_move_layer_window(QtWidgets.QWidget):
    x_move = 0
    y_move = 0
    z_move = 0
    emit_information_signal = QtCore.pyqtSignal(float, float, float)
    def __init__(self, crys1, crys2, color_lst, radii_lst, parent=None):
        super(scatter_plot_move_layer_window, self).__init__(parent)
        pg.setConfigOption('background', (35, 38, 41))
        self.radii_lst = radii_lst
        self.color_lst = color_lst
        self.resize(1000, 2000)
        self.crys1 = crys1  # down
        self.crys2 = crys2  # up
        layout = QtWidgets.QVBoxLayout(self)  # 垂直布局
        self.determine_button = QtWidgets.QPushButton('OK')
        self.view = pg.GraphicsLayoutWidget()
        self.p1 = self.view.addPlot()
        self.vb = self.p1.vb
        self.s1 = pg.ScatterPlotItem(pen=pg.mkPen(None))
        self.s1.setPxMode(False)
        layout.addWidget(self.view)
        self.sliderx = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        cell_par = self.crys1.get_cell_lengths_and_angles()
        x_range = np.floor(cell_par[0] * 100) + 1
        y_range = np.floor(cell_par[1] * np.sin(cell_par[5] / 180 * np.pi) * 100) + 1
        crys2_pos_lst = list(self.crys2.get_positions())
        self.z_value = 1000000000000
        z_max = -10
        for pos in crys2_pos_lst:
            if pos[2] < self.z_value:
                self.z_value = pos[2]
            if pos[2] > z_max:
                z_max = pos[2]
        houdu = z_max - self.z_value
        z_range = np.floor((cell_par[2] - houdu) * 100)
        self.sliderx.setRange(0, x_range)  # 3
        self.sliderx.valueChanged.connect(self.on_change_funcx)
        self.label_sliderx = QtWidgets.QLabel('x:0(Å)', self)
        self.label_sliderx.setFont(QtGui.QFont('Arial Black', 10))
        self.slidery = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.slidery.setRange(0, y_range)  # 3
        self.slidery.valueChanged.connect(self.on_change_funcy)
        self.label_slidery = QtWidgets.QLabel('y:0(Å)', self)
        self.label_slidery.setFont(QtGui.QFont('Arial Black', 10))
        self.sliderz = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.sliderz.setRange(0, z_range)  # 3
        self.sliderz.valueChanged.connect(self.on_change_funcz)
        self.sliderz.setValue(np.floor(self.z_value * 100))
        self.label_sliderz = QtWidgets.QLabel('z:{}(Å)'.format(np.floor(self.z_value * 100)/100), self)
        self.label_sliderz.setFont(QtGui.QFont('Arial Black', 10))
        self.x_move = float(self.sliderx.value() / 100)
        self.y_move = float(self.slidery.value() / 100)
        self.z_move = float(self.sliderz.value() / 100)
        h_layout1 = QtWidgets.QHBoxLayout()
        h_layout1.addWidget(self.label_sliderx)
        h_layout1.addWidget(self.sliderx)
        h_layout2 = QtWidgets.QHBoxLayout()
        h_layout2.addWidget(self.label_slidery)
        h_layout2.addWidget(self.slidery)
        h_layout3 = QtWidgets.QHBoxLayout()
        h_layout3.addWidget(self.label_sliderz)
        h_layout3.addWidget(self.sliderz)
        h_layout3.addWidget(self.determine_button)
        layout.addLayout(h_layout1)
        layout.addLayout(h_layout2)
        layout.addLayout(h_layout3)
        self.determine_button.clicked.connect(self.determine)
        self.setWindowTitle("Atoms Plot move layer")
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        self.plot()

    def plot(self):               #  orijin
        try:
            self.p1.clear()
            self.s1 = pg.ScatterPlotItem(pen=pg.mkPen(None))
            self.s1.setPxMode(False)
            self.lattice_enlarge()
            self.xd = [ori_pos_piece[0] for ori_pos_piece in self.pos_lst_orijin]
            self.yd = [ori_pos_piece[1] for ori_pos_piece in self.pos_lst_orijin]
            self.xu = [upori_pos_piece[0] for upori_pos_piece in self.pos_lst_orijin_up]
            self.yu = [upori_pos_piece[1] for upori_pos_piece in self.pos_lst_orijin_up]
            self.lattice_plot()
            self.p1.addItem(self.s1)
            self.show()
        except Exception as e:
            print(e)


    def lattice_plot(self):
        try:
            array_bx = np.array(self.x_axisb)
            array_by = np.array(self.y_axisb)
            x = np.append(self.xd, self.xu)
            y = np.append(self.yd, self.yu)
            pos = np.array([x, y])
            spots = [{'pos': pos[:, i], 'data': 1, 'brush': pg.mkBrush(255, 250, 250, 120),
                      'size':self.radii_lst_crys1[i] * 2} for i in range(len(self.xd))] +\
                     [{'pos': pos[:, i], 'data': 1, 'brush': pg.mkBrush(255, 20, 147, 120),
                       'size': self.radii_lst_crys2[i - len(self.xd)] * 2}
                      for i in range(len(self.xd), len(x))]
            self.p1.plot(array_bx, array_by)
            self.s1.addPoints(spots)
        except Exception as e:
            print(e)

    def on_change_funcx(self):
        try:
            self.label_sliderx.setText("x:" + str(self.sliderx.value()/100) + "(Å)")
            self.x_move = float(self.sliderx.value()/100)
            self.new_x_up = list(map(lambda x: x + self.x_move, self.xu))
            self.new_y_up = list(map(lambda y: y + self.y_move, self.yu))
            self.new_x = np.append(self.xd, self.new_x_up)       # 所有x_down点 + x_up点
            self.new_y = np.append(self.yd, self.new_y_up)
            pos = np.array([self.new_x, self.new_y])
            self.new_spot = [{'pos': pos[:, i], 'data': 1, 'brush': pg.mkBrush(255, 250, 250, 120),
              'size': self.radii_lst_crys1[i] * 2} for i in range(len(self.xd))] + \
            [{'pos': pos[:, i], 'data': 1, 'brush': pg.mkBrush(255, 20, 147, 120),
              'size': self.radii_lst_crys2[i - len(self.xd)] * 2}
             for i in range(len(self.xd), len(self.new_x))]
            self.repaint()
        except Exception as e:
            print(e)

    def on_change_funcy(self):
        try:
            self.label_slidery.setText("y:" + str(self.slidery.value()/100) + "(Å)")
            self.y_move = float(self.slidery.value()/100)
            self.new_x_up = list(map(lambda x: x + self.x_move, self.xu))
            self.new_y_up = list(map(lambda y: y + self.y_move, self.yu))
            self.new_x = np.append(self.xd, self.new_x_up)       # 所有x_down点 + x_up点
            self.new_y = np.append(self.yd, self.new_y_up)
            pos = np.array([self.new_x, self.new_y])
            self.new_spot = [{'pos': pos[:, i], 'data': 1, 'brush': pg.mkBrush(255, 250, 250, 120),
                              'size': self.radii_lst_crys1[i] * 2} for i in range(len(self.xd))] + \
                            [{'pos': pos[:, i], 'data': 1, 'brush': pg.mkBrush(255, 20, 147, 120),
                              'size': self.radii_lst_crys2[i - len(self.xd)] * 2}
                             for i in range(len(self.xd), len(self.new_x))]
            self.repaint()
        except Exception as e:
            print(e)

    def on_change_funcz(self):
        try:
            self.label_sliderz.setText("z:" + str(self.sliderz.value() / 100) + "(Å)")
            self.z_move = float(self.sliderz.value() / 100)
        except Exception as e:
            print(e)

    def repaint(self):
        try:
            array_bx = np.array(self.x_axisb)
            array_by = np.array(self.y_axisb)
            self.p1.clear()
            self.p1.addItem(self.s1)
            self.p1.plot(array_bx, array_by)
            self.s1.setPoints(self.new_spot)
        except Exception as e:
            print(e)

    def lattice_enlarge(self):
        try:
            print("in lattice enlarge")
            cell1 = list(self.crys1.get_cell())
            crys1_pos_3 = list(self.crys1.get_positions())
            z2_max = max(_[2] for _ in crys1_pos_3)
            crys1_pos_2 = [np.array([pos[0], pos[1]])
                           for pos in crys1_pos_3 if abs(pos[2] - z2_max) < .5]
            crys1_atomic_number = self.crys1.get_atomic_numbers()
            self.radii_lst_crys1 = [self.radii_lst[crys1_atomic_number[i]]
                               for i, pos in enumerate(crys1_pos_3) if abs(pos[2] - z2_max) < .5]
            crys2_pos_3 = list(self.crys2.get_positions())
            crys2_pos_2 = [np.array([pos[0], pos[1]])
                           for pos in crys2_pos_3 if abs(pos[2] - self.z_value) < .5]
            crys2_atomic_number = self.crys2.get_atomic_numbers()
            self.radii_lst_crys2 = [self.radii_lst[crys2_atomic_number[i]]
                               for i, pos in enumerate(crys2_pos_3) if abs(pos[2] - self.z_value) < .5]
            vectora1 = np.array(cell1[0][:2])
            self.vectora1 = vectora1
            vectorb1 = np.array(cell1[1][:2])
            self.vectorb1 = vectorb1
            self.x_axisb = [0, self.vectora1[0], self.vectora1[0] + self.vectorb1[0], self.vectorb1[0], 0]
            self.y_axisb = [0, self.vectora1[1], self.vectora1[1] + self.vectorb1[1], self.vectorb1[1], 0]
            self.pos_lst_orijin = []
            self.pos_lst_orijin_up = []
            gedian_lst_down = []
            for i in range(0, 2):
                for j in range(0, 2):
                    pos = vectora1 * i + vectorb1 * j
                    gedian_lst_down.append(pos)
            self.radii_lst_crys1 *= len(gedian_lst_down)
            for pingyi_pos in gedian_lst_down:
                for pos1 in crys1_pos_2:
                    self.pos_lst_orijin.append(pos1 + pingyi_pos)
            gedian_lst_up = []
            for i in range(-1, 1):
                for j in range(-1, 1):
                    pos = vectora1 * i + vectorb1 * j
                    gedian_lst_up.append(pos)
            for pingyi_pos in gedian_lst_up:
                for pos2 in crys2_pos_2:
                    self.pos_lst_orijin_up.append(pos2 + pingyi_pos)
            self.radii_lst_crys2 *= len(gedian_lst_up)
        except Exception as e:
            print(e)

    def determine(self):
        try:
            self.emit_information_signal.emit(self.x_move, self.y_move, self.z_move - self.z_value)
            self.close()
        except Exception as e:
            print(e)


class scatter_plot_window(QtWidgets.QWidget):
    len_max = 15
    normal_strain_limit = .02
    shear_strain_limit = .01
    area_error_limit = .045
    error_limit = .0224
    emit_information_signal = QtCore.pyqtSignal(float, list, list, float, float, float, float, float, list, list)
    def __init__(self, crys1, crys2, parent=None):
        super(scatter_plot_window, self).__init__(parent)
        pg.setConfigOption('background', (35, 38, 41))
        self.resize(1000, 2000)
        self.crys1 = crys1      # down
        self.crys2 = crys2      # up
        cell1 = list(self.crys1.get_cell())
        self.vectora1 = np.array(cell1[0][:2])
        self.vectorb1 = np.array(cell1[1][:2])
        cell2 = list(self.crys2.get_cell())
        self.vectora2 = np.array(cell2[0][:2])
        self.vectorb2 = np.array(cell2[1][:2])
        layout = QtWidgets.QVBoxLayout(self)  # 垂直布局
        Vsplitter = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        self.set_par_button = QtWidgets.QPushButton('Set parameter')
        self.Show_all_Checkbox = QtWidgets.QCheckBox("Show all results")
        self.Show_all_Checkbox.clicked.connect(self.show_all)
        self.determine_button = QtWidgets.QPushButton('OK')
        Tab = QtWidgets.QTabWidget()
        self.supercell_par_tab = QtWidgets.QWidget()
        self.parmameter_tab = QtWidgets.QWidget()
        Tab.addTab(self.parmameter_tab, "Parameter")
        Tab.addTab(self.supercell_par_tab, "Supercell parameter")
        supercell_par_formlayout = QtWidgets.QFormLayout(self.supercell_par_tab)
        parameter_gridlayout = QtWidgets.QGridLayout(self.parmameter_tab)
        self.normalstrain_lineedit = QtWidgets.QLineEdit()
        self.normalstrain_lineedit.setPlaceholderText(str(self.normal_strain_limit * 100))
        self.normalstrain_lineedit.setText(str(self.normal_strain_limit * 100))
        self.shearstrain_lineedit = QtWidgets.QLineEdit()
        self.shearstrain_lineedit.setPlaceholderText(str(self.shear_strain_limit * 100))
        self.shearstrain_lineedit.setText(str(self.shear_strain_limit * 100))
        self.areaerror_lineedit = QtWidgets.QLineEdit()
        self.areaerror_lineedit.setText(str(self.area_error_limit * 100))
        self.areaerror_lineedit.setPlaceholderText(str(self.area_error_limit * 100))
        area_label = QtWidgets.QLabel("Area error limit(%):")
        area_label.setBuddy(self.areaerror_lineedit)
        normal_strain_label = QtWidgets.QLabel("Normal strain limit(%)：")
        normal_strain_label.setBuddy(self.normalstrain_lineedit)
        shear_strain_label = QtWidgets.QLabel("Shear strain limit(%)：")
        shear_strain_label.setBuddy(self.shearstrain_lineedit)
        self.maxlength_lineedit = QtWidgets.QLineEdit()
        self.maxlength_lineedit.setText(str(self.len_max))
        self.maxlength_lineedit.setPlaceholderText(str(self.len_max))
        max_length_label = QtWidgets.QLabel("Max length(Å): ")
        max_length_label.setBuddy(self.maxlength_lineedit)
        parameter_gridlayout.addWidget(normal_strain_label, 0, 0, 1, 3)
        parameter_gridlayout.addWidget(self.normalstrain_lineedit, 0, 3, 1, 3)
        parameter_gridlayout.addWidget(shear_strain_label, 0, 6, 1, 3)
        parameter_gridlayout.addWidget(self.shearstrain_lineedit, 0, 9, 1, 3)
        parameter_gridlayout.addWidget(area_label, 1, 0, 1, 3)
        parameter_gridlayout.addWidget(self.areaerror_lineedit, 1, 3, 1, 3)
        parameter_gridlayout.addWidget(max_length_label, 1, 6, 1, 3)
        parameter_gridlayout.addWidget(self.maxlength_lineedit, 1, 9, 1, 2)
        parameter_gridlayout.addWidget(self.set_par_button, 1, 11, 1, 1)
        self.layerdis_lineedit = QtWidgets.QLineEdit()
        layerdis_label = QtWidgets.QLabel("layer distance(Å)：")
        self.vacuumdis_lineedit = QtWidgets.QLineEdit()
        vacuumdis_label = QtWidgets.QLabel("Vacuum layer distance(Å)：")
        supercell_par_formlayout.addRow(layerdis_label, self.layerdis_lineedit)
        supercell_par_formlayout.addRow(vacuumdis_label, self.vacuumdis_lineedit)
        self.view = pg.GraphicsLayoutWidget()
        self.p1 = self.view.addPlot()
        self.vb = self.p1.vb
        self.s1 = pg.ScatterPlotItem(size=10, pen=pg.mkPen(None))
        pg.setConfigOption('background', (35, 38, 41))
        text_showall_Hsplitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.show_all_table_widget = QtWidgets.QTableWidget()
        self.show_all_table_widget.setColumnCount(7)
        self.show_all_table_widget.setHorizontalHeaderLabels(['a', 'b', 'γ(°)', 'Area', 'Normal strain', 'Shear strain',
                                                      'Area relative error'])
        self.show_all_table_widget.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.show_all_table_widget.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.show_all_table_widget.itemClicked.connect(self.handleItemClick_on_table)

        self.text_widgt = QtWidgets.QTextEdit()
        self.text_widgt.setReadOnly(True)
        self.text_widgt.resize(24, 20)
        self.text_widgt.setStyleSheet("color: rgb(255, 255, 255);background-color: rgb(35, 38, 41);")
        text_showall_Hsplitter.addWidget(self.text_widgt)
        text_showall_Hsplitter.addWidget(self.show_all_table_widget)
        Vsplitter.addWidget(self.view)
        Vsplitter.addWidget(text_showall_Hsplitter)
        self.show_all_table_widget.setHidden(True)
        self.slider = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.slider.setFixedSize(700, 50)  # 2
        self.slider.setRange(0, 1800)  # 3
        self.slider.valueChanged.connect(self.on_change_func)  # 5
        self.label_slider = QtWidgets.QLabel('0', self)
        self.label_slider.setFont(QtGui.QFont('Arial Black', 20))
        rotate_Hlayout = QtWidgets.QHBoxLayout()
        rotate_Hlayout.addWidget(self.slider)
        rotate_Hlayout.addWidget(self.label_slider)
        rotate_Hlayout.addWidget(self.Show_all_Checkbox)
        rotate_Hlayout.addWidget(self.determine_button)
        Vsplitter.addWidget(Tab)
        layout.addWidget(Vsplitter)
        layout.addLayout(rotate_Hlayout)
        self.set_par_button.clicked.connect(self.set_parameter)
        self.determine_button.clicked.connect(self.determine)
        self.setWindowTitle("Scatter Plot")
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        self.plot()

    def on_change_func(self):
        self.theta = float(self.slider.value()/10) / 180 * np.pi
        self.label_slider.setText(str(self.slider.value()/10) + "°")
        self.rotate_matrix = np.array([[np.cos(self.theta), np.sin(self.theta)], [-np.sin(self.theta), np.cos(self.theta)]])
        self.new_x_up = []
        self.new_y_up = []
        for pos in zip(self.xu, self.yu):
            new_pos = np.array(pos) @ self.rotate_matrix
            self.new_x_up.append(new_pos[0])
            self.new_x_up.append(new_pos[1])
        self.new_x = np.append(self.xd, self.new_x_up)       # 所有x_down点 + x_up点
        self.new_y = np.append(self.yd, self.new_y_up)
        pos = np.array([self.new_x, self.new_y])
        self.new_spot = [{'pos': pos[:, i], 'data': 1, 'brush':pg.mkBrush(255, 250, 250, 120)} for i in range(len(self.xd))] + \
                [{'pos': pos[:, i], 'data': 1, 'brush': pg.mkBrush(255, 20, 147, 120)} for i in range(len(self.xd), len(self.new_x))]
        self.find_min()

    def plot(self):
        self.p1.clear()
        self.s1 = pg.ScatterPlotItem(size=10, pen=pg.mkPen(None))
        self.lattice_enlarge()
        self.xd = [ori_pos_piece[0] for ori_pos_piece in self.pos_lst_orijin]
        self.yd = [ori_pos_piece[1] for ori_pos_piece in self.pos_lst_orijin]
        self.xu = [upori_pos_piece[0] for upori_pos_piece in self.pos_lst_orijin_up]
        self.yu = [upori_pos_piece[1] for upori_pos_piece in self.pos_lst_orijin_up]
        self.lattice_plot()
        self.p1.addItem(self.s1)
        self.show()

    def lattice_enlarge(self):
        cell1 = list(self.crys1.get_cell())
        vectora1 = np.array(cell1[0][:2])
        vectorb1 = np.array(cell1[1][:2])
        a_len1 = np.linalg.norm(vectora1)
        b_len1 = np.linalg.norm(vectorb1)
        if a_len1 < b_len1:
            n1 = int(self.len_max / a_len1)
        else:
            n1 = int(self.len_max / b_len1)
        cell2 = list(self.crys2.get_cell())
        vectora2 = np.array(cell2[0][:2])
        vectorb2 = np.array(cell2[1][:2])
        a_len2 = np.linalg.norm(vectora2)
        b_len2 = np.linalg.norm(vectorb2)
        if a_len2 < b_len2:
            n2 = int(self.len_max / a_len2)
        else:
            n2 = int(self.len_max / b_len2)
        self.pos_lst_orijin = []
        self.pos_lst_orijin_up = []
        for i in range(-2 * n2, 2 * n2):
            for j in range(-2 * n2, 2 * n2):
                pos = vectora2 * i + vectorb2 * j
                if np.linalg.norm(pos) <= self.len_max:
                    self.pos_lst_orijin_up.append(pos)
        for i in range(-2 * n1, 2 * n1):
            for j in range(2 * n1):
                pos = vectora1 * i + vectorb1 * j
                if np.linalg.norm(pos) <= self.len_max:
                    self.pos_lst_orijin.append(pos)

    def lattice_plot(self):
        try:
            x = np.append(self.xd, self.xu)
            y = np.append(self.yd, self.yu)
            pos = np.array([x, y])
            spots = [{'pos': pos[:, i], 'data': 1, 'brush': pg.mkBrush(255, 250, 250, 120)} for i in range(len(self.xd))] +\
                     [{'pos': pos[:, i], 'data': 1, 'brush': pg.mkBrush(255, 20, 147, 120)} for i in range(len(self.xd), len(x))]
            self.s1.addPoints(spots)
        except Exception as e:
            print(e)

    def set_parameter(self):
        try:
            self.len_max = float(self.maxlength_lineedit.text())
            self.normal_strain_limit = float(self.normalstrain_lineedit.text()) / 100
            self.normalstrain_lineedit.setPlaceholderText(str(self.normal_strain_limit * 100))
            self.shear_strain_limit = float(self.shearstrain_lineedit.text()) / 100
            self.shearstrain_lineedit.setPlaceholderText(str(self.shear_strain_limit * 100))
            self.error_limit = np.sqrt(self.normal_strain_limit ** 2 + self.shear_strain_limit ** 2)
            self.area_error_limit = float(self.areaerror_lineedit.text()) / 100
            self.areaerror_lineedit.setPlaceholderText(str(self.area_error_limit * 100))
            self.maxlength_lineedit.setPlaceholderText(str(self.len_max))
            self.label_slider.setText("0")
            self.slider.setValue(0)
            self.theta = 0
            self.plot()
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, 'error', "Input error!")
            print(e)

    def show_all(self):
        try:
            if self.Show_all_Checkbox.isChecked():
                self.all_lst = self.show_all_peidui_lst()
                if self.all_lst is not None:
                    print("all_lst = ", self.all_lst)
                    self.show_all_table_widget.setHidden(False)
                    self.text_widgt.setHidden(True)
                    self.show_all_table_widget.setRowCount(0)
                    self.show_all_table_widget.setRowCount(len(self.all_lst))
                    for i in range(len(self.all_lst)):
                        for j in range(2, len(self.all_lst[0])):
                            text = str(self.all_lst[i][j])
                            newitem = QtWidgets.QTableWidgetItem(text)
                            newitem.setTextAlignment(5 | 5)
                            self.show_all_table_widget.setItem(i, j - 2, newitem)
            else:
                self.show_all_table_widget.setHidden(True)
                self.text_widgt.setHidden(False)
        except Exception as e:
            print(e)


    def determine(self):
        try:
            layer_dis = float(self.layerdis_lineedit.text())
            vacuum_dis = float(self.vacuumdis_lineedit.text())
            self.emit_information_signal.emit(self.theta, list(self.a1), list(self.b1), layer_dis,
                                              vacuum_dis,
                                              self.normal_strain_limit, self.shear_strain_limit,
                                              self.area_error_limit,
                                              list(self.a2), list(self.b2))
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, 'error', "Input error!")
            print(e)

    def handleItemClick_on_table(self):
        """ To handle item click on result tablewidget."""
        try:
            index_lst = []
            for item in self.show_all_table_widget.selectedItems():
                index_lst.append(item.row())
            print("index_lst = ", index_lst)
            index = index_lst[-1]
            print("index = ", index_lst[-1])
            a_lst, b_lst = self.all_lst[index][0], self.all_lst[index][1]
            self.a1 = deepcopy(a_lst[0])
            self.b1 = deepcopy(b_lst[0])
            self.a2 = deepcopy(a_lst[1])
            self.b2 = deepcopy(b_lst[1])
            self.plot_cell()
        except Exception as e:
             print(e)

    def show_all_peidui_lst(self):
        try:
            # par_lst = [a_lst, b_lst, a, b, gama, area, normal_strain, shear_strain, area_strain]
            all_lst = []
            for i in range(len(self.peidui_lst)):
                for j in range(i + 1, len(self.peidui_lst)):
                    par_lst = [0 for _ in range(9)]
                    n_vector1 = np.cross(self.peidui_lst[i][0], self.peidui_lst[j][0])
                    square2 = abs(np.cross(self.peidui_lst[i][1], self.peidui_lst[j][1]))
                    square1 = abs(n_vector1)
                    if 1e-6 < square1 and abs(square1 - square2) < self.area_error_limit * square1:
                        shear_strain = abs(np.arccos((self.peidui_lst[i][0] @ self.peidui_lst[j][0]) /
                                                     np.linalg.norm(self.peidui_lst[i][0]) / np.linalg.norm(
                            self.peidui_lst[j][0])) -
                                           np.arccos((self.peidui_lst[i][1] @ self.peidui_lst[j][1]) /
                                                     np.linalg.norm(self.peidui_lst[i][1]) / np.linalg.norm(
                                               self.peidui_lst[j][1])))
                        if shear_strain < self.shear_strain_limit:
                            if n_vector1 > 0:  # 保证右手系
                                a_lst = self.peidui_lst[i]
                                b_lst = self.peidui_lst[j]
                            else:
                                b_lst = self.peidui_lst[i]
                                a_lst = self.peidui_lst[j]
                            par_lst[0], par_lst[1], par_lst[7] = a_lst, b_lst, round(shear_strain, 5)
                            par_lst[5] = round(square1, 2)
                            par_lst[8] = round(abs(square1 - square2) / square1, 5)
                            par_lst[2], par_lst[3] = round(np.linalg.norm(a_lst[0]), 2), \
                                                     round(np.linalg.norm(b_lst[0]), 2)
                            gama = round(np.arccos((a_lst[0] @ b_lst[0]) / np.linalg.norm(a_lst[0]) / \
                                   np.linalg.norm(b_lst[0])) / np.pi * 180, 2)
                            par_lst[4] = gama
                            norm_straina = round(abs(np.linalg.norm(a_lst[0]) - np.linalg.norm(a_lst[1])) / \
                                           np.linalg.norm(a_lst[0]), 5)
                            norm_strainb = round(abs(np.linalg.norm(b_lst[0]) - np.linalg.norm(b_lst[1])) / \
                                           np.linalg.norm(b_lst[0]), 5)
                            par_lst[6] = [norm_straina, norm_strainb]
                            all_lst.append(par_lst)
            if len(all_lst) != 0:
                all_lst.sort(key=itemgetter(5))
                return all_lst
        except Exception as e:
            print(e)

    def find_min(self):
        try:
            self.peidui_lst = []
            self.new_vectora2 = self.vectora2 @ self.rotate_matrix
            self.new_vectorb2 = self.vectorb2 @ self.rotate_matrix
            solve_matrix = np.linalg.inv(np.array([self.new_vectora2, self.new_vectorb2]))
            for _ in zip(self.xd, self.yd):
                pos1 = np.array(_)
                r = pos1 @ solve_matrix
                r1 = np.floor(r[0])
                r2 = np.floor(r[1])
                pos2 = self.new_vectora2 * r1 + self.new_vectorb2 * r2
                if np.linalg.norm(pos1 - pos2) < self.error_limit * np.linalg.norm(pos1):
                    if abs(pos1 @ pos1 - pos1 @ pos2) < self.normal_strain_limit * (pos1 @ pos1):
                        self.peidui_lst.append([pos1, pos2])
                pos2 = self.new_vectora2 * (r1 + 1) + self.new_vectorb2 * r2
                if np.linalg.norm(pos1 - pos2) < self.error_limit * np.linalg.norm(pos1):
                    if abs(pos1 @ pos1 - pos1 @ pos2) < self.normal_strain_limit * (pos1 @ pos1):
                        self.peidui_lst.append([pos1, pos2])
                pos2 = self.new_vectora2 * (r1 + 1) + self.new_vectorb2 * (r2 + 1)
                if np.linalg.norm(pos1 - pos2) < self.error_limit * np.linalg.norm(pos1):
                    if abs(pos1 @ pos1 - pos1 @ pos2) < self.normal_strain_limit * (pos1 @ pos1):
                        self.peidui_lst.append([pos1, pos2])
                pos2 = self.new_vectora2 * r1 + self.new_vectorb2 * (r2 + 1)
                if np.linalg.norm(pos1 - pos2) < self.error_limit * np.linalg.norm(pos1):
                    if abs(pos1 @ pos1 - pos1 @ pos2) < self.normal_strain_limit * (pos1 @ pos1):
                        self.peidui_lst.append([pos1, pos2])
            if len(self.peidui_lst) < 2:
                self.text = ""
                self.p1.clear()
                self.p1.addItem(self.s1)
                self.s1.setPoints(self.new_spot)
                self.text_widgt.setText(self.text)
                if self.Show_all_Checkbox.isChecked():
                    self.show_all_table_widget.setRowCount(0)
            else:
                self.min = 1000000
                self.a_lst = ...
                self.b_lsg = ...
                for index in range(len(self.peidui_lst)):
                    for j in range(index + 1, len(self.peidui_lst)):
                        n_vector1 = np.cross(self.peidui_lst[index][0], self.peidui_lst[j][0])
                        square2 = abs(np.cross(self.peidui_lst[index][1], self.peidui_lst[j][1]))
                        square1 = abs(n_vector1)
                        if 1e-6 < square1 < self.min and abs(square1 - square2) < self.area_error_limit * square1:
                            shear_strain = abs(np.arccos((self.peidui_lst[index][0] @ self.peidui_lst[j][0]) /
                                      np.linalg.norm(self.peidui_lst[index][0]) / np.linalg.norm(self.peidui_lst[j][0])) -
                                              np.arccos((self.peidui_lst[index][1] @ self.peidui_lst[j][1]) /
                                      np.linalg.norm(self.peidui_lst[index][1]) / np.linalg.norm(self.peidui_lst[j][1])))
                            if shear_strain < self.shear_strain_limit:
                                if n_vector1 > 0:     # 保证右手系
                                    self.a_lst = self.peidui_lst[index]
                                    self.b_lst = self.peidui_lst[j]
                                else:
                                    self.b_lst = self.peidui_lst[index]
                                    self.a_lst = self.peidui_lst[j]
                                self.min = square1
                                print("self.min = ", self.min)
                if self.min == 1000000:
                    self.text = ""
                    self.p1.clear()
                    self.p1.addItem(self.s1)
                    self.s1.setPoints(self.new_spot)
                    self.text_widgt.setText(self.text)
                    if self.Show_all_Checkbox.isChecked():
                        self.show_all_table_widget.setRowCount(0)
                else:
                    self.a1 = self.a_lst[0]
                    self.b1 = self.b_lst[0]
                    self.a2 = self.a_lst[1]
                    self.b2 = self.b_lst[1]
                    Area_limit_error = abs(abs(np.cross(self.a1, self.b1)) - abs(np.cross(self.a_lst[1], self.b_lst[1]))) / \
                                       abs(np.cross(self.a1, self.b1))
                    delta_gama = abs(np.arccos((self.a1 @ self.b1) / np.linalg.norm(self.a1) / np.linalg.norm(self.b1))\
                                 / np.pi * 180 - np.arccos((self.a2 @ self.b2) / np.linalg.norm(self.a2) /
                                                           np.linalg.norm(self.b2))\
                                 / np.pi * 180)

                    self.text = "Optimal match(xy coordinate):  a = {}, b = {}".format(self.a1, self.b1) + "\n" +\
                                '-'*16 + '\n' + \
                                "Minimum square = {}".format(self.min) + '\n' + \
                                "Area relative error limit = {}".format(Area_limit_error) + '\n' + \
                                "Δγ(°) = {}".format(delta_gama)

                    self.text_widgt.setText(self.text)
                    self.plot_cell()
                    if self.Show_all_Checkbox.isChecked():
                        self.show_all()

        except Exception as e:
            print(e)

    def plot_cell(self):
        try:
            self.x_axisb = [0, self.a1[0], self.a1[0] + self.b1[0], self.b1[0], 0]
            self.y_axisb = [0, self.a1[1], self.a1[1] + self.b1[1], self.b1[1], 0]
            array_bx = np.array(self.x_axisb)
            array_by = np.array(self.y_axisb)
            self.p1.clear()
            self.p1.addItem(self.s1)
            self.p1.plot(array_bx, array_by)
            self.s1.setPoints(self.new_spot)
        except Exception as e:
            print(e)

class zhexian_plot_window(QtWidgets.QWidget):
    """ The layout and function of zhexian_plot window. """
    signal_determine = QtCore.pyqtSignal(int, int, float, float)
    def __init__(self, lista, listb, obj1, obj2, axis_lst, parent=None):
        pg.setConfigOption('background', (35, 38, 41))
        super(zhexian_plot_window, self).__init__(parent)
        self.lista = lista
        self.listb = listb
        self.obj1 = obj1
        self.obj2 = obj2
        self.axis_lst = axis_lst
        layout = QtWidgets.QVBoxLayout(self)    # 垂直布局
        Hlayout = QtWidgets.QHBoxLayout()
        splitter2 = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        self.determine_button = QtWidgets.QPushButton('OK')
        self.layerdis_lineedit = QtWidgets.QLineEdit()
        self.Vacuumdis_lineedit = QtWidgets.QLineEdit()

        layerdis_label = QtWidgets.QLabel("Layer distance(Å):")
        layerdis_label.setBuddy(self.layerdis_lineedit)
        vacuumdis_label = QtWidgets.QLabel("Vacuum distance(Å):")
        vacuumdis_label.setBuddy(self.Vacuumdis_lineedit)
        self.lastClicked = []
        self.lastClicked2 = []
        self.view = pg.GraphicsLayoutWidget()  ## GraphicsView with GraphicsLayout inserted by default
        self.plot_window1 = self.view.addPlot(title="Supercell(a) - Normal Strain(a) plot")
        self.plot_window2 = self.view.addPlot(title="Supercell(b) - Normal Strain(b) plot")
        self.plot_window1.setLabel('left', "Supercell(a)", units='Å')
        self.plot_window1.setLabel('bottom', "Normal strain(a)", units='')
        self.plot_window2.setLabel('left', "Supercell(b)", units='Å')
        self.plot_window2.setLabel('bottom', "Normal strain(b)", units='')
        self.text_widgt = QtWidgets.QTextEdit()
        self.text_widgt.setReadOnly(True)
        # self.text_widgt.resize(24, 20)
        self.text_widgt.setStyleSheet("color: rgb(255, 255, 255);background-color: rgb(35, 38, 41);")
        splitter2.addWidget(self.view)
        splitter2.addWidget(self.text_widgt)
        layout.addWidget(splitter2)
        Hlayout.addWidget(layerdis_label)
        Hlayout.addWidget(self.layerdis_lineedit)
        Hlayout.addWidget(vacuumdis_label)
        Hlayout.addWidget(self.Vacuumdis_lineedit)
        Hlayout.addWidget(self.determine_button)
        layout.addLayout(Hlayout)
        self.determine_button.clicked.connect(self.determine)
        self.setWindowTitle("Normal strain a(b) - Supercell a(b) Plot")
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        self.plot(self.axis_lst)
        self.set_text_init(self.obj1, self.obj2)
        self.resize(1000, 1000)

    def plot(self, plotlist):        # 画折线图
        self.list = plotlist
        self.x_axisa = [_[0] for _ in self.list[0]]
        self.y_axisa = [_[1] for _ in self.list[0]]
        array_ax = np.array(self.x_axisa)
        array_ay = np.array(self.y_axisa)
        self.plot_window1.plot(array_ax, array_ay)
        self.scene1 = pg.ScatterPlotItem(size=10, pen=pg.mkPen(None), brush=pg.mkBrush(0, 255, 255, 255))
        self.scene1.setData(x=self.x_axisa, y=self.y_axisa)
        self.plot_window1.addItem(self.scene1)
        self.scene1.sigClicked.connect(self.clicked)
        self.plot1 = pg.SignalProxy(self.plot_window1.scene().sigMouseMoved, rateLimit=60,
                                    slot=self.MouseMoveEvent1)
        self.plot2 = pg.SignalProxy(self.plot_window2.scene().sigMouseMoved, rateLimit=60,
                                    slot=self.MouseMoveEvent2)
        self.scene2 = pg.ScatterPlotItem(size=10, pen=pg.mkPen(None), brush=pg.mkBrush(255, 255, 0, 255))
        self.x_axisb = [_[0] for _ in self.list[1]]
        self.y_axisb = [_[1] for _ in self.list[1]]
        array_bx = np.array(self.x_axisb)
        array_by = np.array(self.y_axisb)
        self.plot_window2.plot(array_bx, array_by)
        self.scene2.setData(x=np.array(self.x_axisb), y=np.array(self.y_axisb))
        self.plot_window2.addItem(self.scene2)
        self.scene2.sigClicked.connect(self.clicked2)
        self.show()         # 绘图

    def set_text_init(self, crys1, crys2):
        crys1_formula = crys1.get_chemical_formula(mode='hill')
        crys2_formula = crys2.get_chemical_formula(mode='hill')
        crys1_info = crys1.get_cell_lengths_and_angles()
        crys2_info = crys2.get_cell_lengths_and_angles()
        info = crys1_formula + ' stack ' + crys2_formula + '\n' + '-'*32 + '\n' + \
               crys1_formula + ': a=' + str(crys1_info[0]) + '    b = ' + str(crys1_info[1]) + \
               '    gamma=' + str(crys1_info[5]) + '\n' + crys2_formula + ': a=' + str(crys2_info[0]) + \
               '    b = ' + str(crys2_info[1]) + '    gamma=' + str(crys1_info[5])
        self.text_widgt.setText(info)

    def judge_point1(self, float1, float2):
        length_min = 1000         # 判断用于点击的点的index
        point1_index = ...
        for index, x_axisa in enumerate(self.x_axisa):
            length = np.sqrt((x_axisa * 10000 - float1 * 10000) ** 2 + (self.y_axisa[index] - float2) ** 2)
            if length < length_min:
                length_min = length
                point1_index = index
        self.findnum1(point1_index)

    def judge_point2(self, float1, float2):
        length_min = 1000         # 判断用于点击的点的index
        point2_index = ...
        for index, x_axisb in enumerate(self.x_axisb):
            length = np.sqrt((x_axisb * 10000 - float1 * 10000) ** 2 + (self.y_axisb[index] - float2) ** 2)
            if length < length_min:
                length_min = length
                point2_index = index
        self.findnum2(point2_index)

    def findnum1(self, number1):
        """ client chooses the Normal_straina - supercella dot.(GLOBAL THE NUMBER)"""
        self.zhexianplot1_info = ""
        self.point1_index = number1
        Normal_straina = str(self.lista[self.point1_index][2])
        supercella = str(self.lista[self.point1_index][3])
        self.zhexianplot1_info = 'Normal strain(a):' + Normal_straina + '    Supercell(a):' + supercella + \
                                 '\n' + '-' * 32
        try:
            if self.point2_index is not None:
                self.text_widgt.setText(self.zhexianplot1_info + '\n' + '\n' + self.zhexianplot2_info)
                obj1_num = len(self.obj1.get_positions())  # client选择的Atoms对象
                obj2_num = len(self.obj2.get_positions())
                crys1_cell_lst = self.obj1.get_cell_lengths_and_angles()
                crys2_cell_lst = self.obj2.get_cell_lengths_and_angles()
                obj_strain_lcma = self.axis_lst[0][self.point1_index]  # client选择的strain-lcma列表
                obj_strain_lcmb = self.axis_lst[1][self.point2_index]  # client选择的strain-lcmb列表
                multi_a_crys1 = round(obj_strain_lcma[1] / crys1_cell_lst[0])
                multi_a_crys2 = round(obj_strain_lcma[1] / crys2_cell_lst[0])
                multi_b_crys1 = round(obj_strain_lcmb[1] / crys1_cell_lst[1])
                multi_b_crys2 = round(obj_strain_lcmb[1] / crys2_cell_lst[1])
                atom_num_sum = int(obj1_num * multi_a_crys1 * multi_b_crys1 + obj2_num * multi_a_crys2 * multi_b_crys2)
                at_num_text = "Number of atoms in supercell:" + str(atom_num_sum)
                self.text_widgt.setText(
                    at_num_text + '\n' + self.zhexianplot1_info + '\n' + '\n' + self.zhexianplot2_info)
            else:
                self.text_widgt.setText(self.zhexianplot1_info)
        except Exception as e:
            print(e)

    def findnum2(self, number2):
        """ client chooses the strainsum - supercellb dot.(GLOBAL THE INDEX)"""
        self.zhexianplot2_info = ""
        self.point2_index = number2
        Normal_strainb = str(self.listb[self.point2_index][2])
        supercellb = str(self.listb[self.point2_index][3])
        self.zhexianplot2_info = 'Normal strain(b):' + Normal_strainb + '    Supercell(b):' + supercellb + \
                                 '\n' + '-' * 32
        try:
            if self.point1_index is not None:
                self.text_widgt.setText(
                    self.zhexianplot1_info + '\n' + '\n' + self.zhexianplot2_info)
                obj1_num = len(self.obj1.get_positions())  # client选择的Atoms对象
                obj2_num = len(self.obj2.get_positions())
                crys1_cell_lst = self.obj1.get_cell_lengths_and_angles()
                crys2_cell_lst = self.obj2.get_cell_lengths_and_angles()
                obj_strain_lcma = self.axis_lst[0][self.point1_index]  # client选择的strain-lcma列表
                obj_strain_lcmb = self.axis_lst[1][self.point2_index]  # client选择的strain-lcmb列表
                multi_a_crys1 = round(obj_strain_lcma[1] / crys1_cell_lst[0])
                multi_a_crys2 = round(obj_strain_lcma[1] / crys2_cell_lst[0])
                multi_b_crys1 = round(obj_strain_lcmb[1] / crys1_cell_lst[1])
                multi_b_crys2 = round(obj_strain_lcmb[1] / crys2_cell_lst[1])
                atom_num_sum = int(obj1_num * multi_a_crys1 * multi_b_crys1 + obj2_num * multi_a_crys2 * multi_b_crys2)
                at_num_text = "Number of atoms in supercell:" + str(atom_num_sum)
                self.text_widgt.setText(
                    at_num_text + '\n' + self.zhexianplot1_info + '\n' + '\n' + self.zhexianplot2_info)
            else:
                self.text_widgt.setText(self.zhexianplot2_info)
        except Exception as e:
            print(e)

    def MouseMoveEvent1(self, event):
        global mousePoint1
        pos = event[0]
        vb = self.plot_window1.vb
        mousePoint1 = vb.mapSceneToView(pos)

    def MouseMoveEvent2(self, event):
        global mousePoint2
        pos = event[0]
        vb = self.plot_window2.vb
        mousePoint2 = vb.mapSceneToView(pos)

    def clicked(self, event, points):
        global mousePoint1
        try:
            for p in self.lastClicked:
                p.resetPen()
            for p in points:
                p.setPen('r', width=2)
                self.judge_point1(mousePoint1.x(), mousePoint1.y())
            self.lastClicked = points
        except Exception as e:
            print(e)

    def clicked2(self, event, points):
        global mousePoint2
        try:
            for p in self.lastClicked2:
                p.resetPen()
            for p in points:
                p.setPen('r', width=2)
                self.judge_point2(mousePoint2.x(), mousePoint2.y())
            self.lastClicked2 = points
        except Exception as e:
            print(e)

    def determine(self):
        try:
            num = 1
            layer_dis = float(self.layerdis_lineedit.text())
            num = 3
            vacuum_dis = float(self.Vacuumdis_lineedit.text())
            num = 2
            self.point1_index = int(self.point1_index)
            self.point2_index = int(self.point2_index)
            self.signal_determine.emit(self.point1_index, self.point2_index, layer_dis, vacuum_dis)  # 选定,与self.layertransform连接
            self.close()
        except Exception as e:
            if num == 1:
                QtWidgets.QMessageBox.warning(self, 'error', 'Layer distance input error!')
            elif num == 2:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please click point to choose!')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Vacuum distance input error!')
            print(e)

class multi_layer_optimize_window(QtWidgets.QWidget):
    signal_emit_infomation = QtCore.pyqtSignal(float, list, list, float, float, float, float, float, float, float, float)
    def __init__(self, obj_lst, parent=None):
        super(multi_layer_optimize_window, self).__init__(parent)
        Vlayout = QtWidgets.QVBoxLayout(self)
        Hlayout4 = QtWidgets.QHBoxLayout()
        self.btn1 = QtWidgets.QRadioButton("a1,a2;b1,b2 have same orientation")
        self.btn1.setChecked(False)
        self.btn1.toggled.connect(lambda: self.btnstate(self.btn1))
        Vlayout.addWidget(self.btn1)
        self.btn2 = QtWidgets.QRadioButton("a1,a2 have same orientation")
        self.btn2.setChecked(False)
        self.btn2.toggled.connect(lambda: self.btnstate(self.btn2))
        Vlayout.addWidget(self.btn2)
        self.btn3 = QtWidgets.QRadioButton("Cell orientation is uncertain")
        self.btn3.setChecked(False)
        self.btn3.toggled.connect(lambda: self.btnstate(self.btn3))
        Vlayout.addWidget(self.btn3)
        Tab = QtWidgets.QTabWidget()
        self.supercell_par_tab = QtWidgets.QWidget()
        self.error_limit_tab = QtWidgets.QWidget()
        Tab.addTab(self.error_limit_tab, "Error limit")
        Tab.addTab(self.supercell_par_tab, "Supercell parameter")
        error_limit_Formlayout = QtWidgets.QFormLayout(self.error_limit_tab)
        supercell_par_Formlayout = QtWidgets.QFormLayout(self.supercell_par_tab)
        self.area_label = QtWidgets.QLabel("Area relative error(%)")
        self.area_label.setHidden(True)
        self.areaerror_lineedit = QtWidgets.QLineEdit()
        self.areaerror_lineedit.setHidden(True)
        self.shear_strain_label = QtWidgets.QLabel("Shear strain limit(%)")
        self.shear_strain_label.setHidden(True)
        self.shear_strain_lineedit = QtWidgets.QLineEdit()
        self.shear_strain_lineedit.setHidden(True)
        self.normal_strain_label = QtWidgets.QLabel("Normal strain limit(%)")
        self.normal_strain_label.setHidden(True)
        self.normal_strain_lineedit = QtWidgets.QLineEdit()
        self.normal_strain_lineedit.setHidden(True)
        error_limit_Formlayout.addRow(self.normal_strain_label, self.normal_strain_lineedit)
        error_limit_Formlayout.addRow(self.shear_strain_label, self.shear_strain_lineedit)
        error_limit_Formlayout.addRow(self.area_label, self.areaerror_lineedit)
        self.vacuum_label = QtWidgets.QLabel("Vacuum distance(Å):")
        self.vacuum_label.setHidden(True)
        self.vacuum_dis_lineedit = QtWidgets.QLineEdit()
        self.vacuum_dis_lineedit.setHidden(True)
        self.max_singama_label = QtWidgets.QLabel("Maximum sin(γ):")
        self.max_singama_lineedit = QtWidgets.QLineEdit()
        self.max_singama_lineedit.setPlaceholderText('1')
        self.min_singama_label = QtWidgets.QLabel("Minimum sin(γ):")
        self.min_singama_lineedit = QtWidgets.QLineEdit()
        self.min_singama_lineedit.setPlaceholderText('0')
        self.max_singama_label.setHidden(True)
        self.max_singama_lineedit.setHidden(True)
        self.min_singama_label.setHidden(True)
        self.min_singama_lineedit.setHidden(True)
        self.ab_label = QtWidgets.QLabel("Supercell a,b relative error(%):")
        self.ab_error_lineedit = QtWidgets.QLineEdit()
        self.max_length_label = QtWidgets.QLabel("Maximum side length of supercell(Å)")
        self.max_length_lineedit = QtWidgets.QLineEdit()
        self.ab_label.setHidden(True)
        self.ab_error_lineedit.setHidden(True)
        self.max_length_label.setHidden(True)
        self.max_length_lineedit.setHidden(True)
        supercell_par_Formlayout.addRow(self.vacuum_label, self.vacuum_dis_lineedit)
        supercell_par_Formlayout.addRow(self.max_singama_label, self.max_singama_lineedit)
        supercell_par_Formlayout.addRow(self.min_singama_label, self.min_singama_lineedit)
        supercell_par_Formlayout.addRow(self.ab_label, self.ab_error_lineedit)
        supercell_par_Formlayout.addRow(self.max_length_label, self.max_length_lineedit)
        self.Start_pushbutton = QtWidgets.QPushButton("Start")
        self.Cancel_pushbutton = QtWidgets.QPushButton("Cancel")
        self.Cancel_pushbutton.clicked.connect(self.close)
        self.Start_pushbutton.setHidden(True)
        self.Start_pushbutton.clicked.connect(self.determine)
        Hlayout4.addWidget(self.Start_pushbutton)
        Hlayout4.addWidget(self.Cancel_pushbutton)
        self.obj_lst = obj_lst
        self.tablewidget = QtWidgets.QTableWidget()
        self.tablewidget.setHidden(True)
        # self.resize(530, 400)
        self.tablewidget.setColumnCount(1)
        self.tablewidget.setHorizontalHeaderLabels(['Hetero-junction'])
        Vlayout.addWidget(Tab)
        Vlayout.addLayout(Hlayout4)
        Vlayout.addWidget(self.tablewidget)
        self.setTable()
        self.setWindowTitle("Multi-layer stack")
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        self.show()

    def setTable(self):
        try:
            number = len(self.obj_lst)
            self.tablewidget.setRowCount(2*number - 1)
            for index, obj in enumerate(self.obj_lst):
                text = obj.get_chemical_formula(mode='hill')
                newitem = QtWidgets.QTableWidgetItem(text)
                newitem.setTextAlignment(5 | 5)
                newitem.setFlags(QtCore.Qt.ItemIsEnabled)
                self.tablewidget.setItem(2*index, 0, newitem)
            c = []
            for index in range(2*number - 1):
                if index % 2 == 0:
                    c.append("{} layer".format(int((index+2)/2)))
                else:
                    c.append("layer distance:")
            self.tablewidget.setVerticalHeaderLabels(c)
        except Exception as e:
            print(e)



    def btnstate(self, btn):
        try:
            self.resize(700, 700)
            self.choosebutton = btn
            self.tablewidget.setHidden(False)
            self.vacuum_label.setHidden(False)
            self.vacuum_dis_lineedit.setHidden(False)
            self.normal_strain_label.setHidden(False)
            self.normal_strain_lineedit.setHidden(False)
            self.Start_pushbutton.setHidden(False)
            self.shear_strain_label.setHidden(False)
            self.shear_strain_lineedit.setHidden(False)

            if self.choosebutton.text() == "Cell orientation is uncertain":
                self.area_label.setHidden(False)
                self.areaerror_lineedit.setHidden(False)
                self.shear_strain_lineedit.setHidden(False)
                self.shear_strain_label.setHidden(False)
                self.max_singama_label.setHidden(False)
                self.max_singama_lineedit.setHidden(False)
                self.min_singama_label.setHidden(False)
                self.min_singama_lineedit.setHidden(False)
                self.ab_label.setHidden(False)
                self.ab_error_lineedit.setHidden(False)
                self.max_length_label.setHidden(False)
                self.max_length_lineedit.setHidden(False)
            elif self.choosebutton.text() == "a1,a2;b1,b2 have same orientation":
                self.area_label.setHidden(True)
                self.areaerror_lineedit.setHidden(True)
                self.shear_strain_label.setHidden(True)
                self.shear_strain_lineedit.setHidden(True)
                self.max_singama_label.setHidden(True)
                self.max_singama_lineedit.setHidden(True)
                self.min_singama_label.setHidden(True)
                self.min_singama_lineedit.setHidden(True)
                self.ab_label.setHidden(True)
                self.ab_error_lineedit.setHidden(True)
                self.max_length_label.setHidden(True)
                self.max_length_lineedit.setHidden(True)
            else:
                self.area_label.setHidden(True)
                self.areaerror_lineedit.setHidden(True)
                self.max_singama_label.setHidden(True)
                self.max_singama_lineedit.setHidden(True)
                self.min_singama_label.setHidden(True)
                self.min_singama_lineedit.setHidden(True)
                self.ab_label.setHidden(True)
                self.ab_error_lineedit.setHidden(True)
                self.max_length_label.setHidden(True)
                self.max_length_lineedit.setHidden(True)
        except Exception as e:
            print(e)


    def determine(self):
        try:
            num = 3
            if self.choosebutton.text() == "a1,a2;b1,b2 have same orientation":
                num = 2
            elif self.choosebutton.text() == "a1,a2 have same orientation":
                num = 1
            elif self.choosebutton.text() == "Cell orientation is uncertain":
                num = 0
            if num == 3:
                QtWidgets.QMessageBox.warning(self, 'error', 'Choose stack style!')
            else:
                layer_dis_lst = []
                for i in range(len(self.obj_lst) - 1):
                    layer_dis_lst.append(float(self.tablewidget.item(2*i+1, 0).text()))
                vacuum_dis = float(self.vacuum_dis_lineedit.text())
                Normal_strain = float(self.normal_strain_lineedit.text()) / 100
                try:
                    Shear_strain = float(self.shear_strain_lineedit.text()) / 100
                except:
                    Shear_strain = 1000
                try:
                    area_error = float(self.areaerror_lineedit.text()) / 100
                except:
                    area_error = 100
                try:
                    max_gama = float(self.max_singama_lineedit.text())
                except:
                    max_gama = 1
                try:
                    min_gama = float(self.min_singama_lineedit.text())
                except:
                    min_gama = 0
                try:
                    ab_relative_error = float(self.ab_error_lineedit.text()) / 100
                except:
                    ab_relative_error = 100000
                try:
                    max_len = float(self.max_length_lineedit.text())
                except:
                    max_len = 0
                self.signal_emit_infomation.emit(num, self.obj_lst, layer_dis_lst, vacuum_dis,
                                                 area_error,  Normal_strain, Shear_strain,
                                                 max_gama, min_gama, ab_relative_error, max_len)
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error', 'Input Error!')


# class exit_message_window(QtWidgets.QWidget):
#     def __init__(self, parent=None):
#         super(exit_message_window, self).__init__(parent)
#         Hlayout = QtWidgets.QHBoxLayout(self)
#         yes_button = QtWidgets.QMessageBox.Yes
#         No_button = QtWidgets.QMessageBox.Yes
#         Hlayout.addWidget(yes_button)
#         Hlayout.addWidget(No_button)
#         self.show()

class Twisted_little_window(QtWidgets.QWidget):
    signal_emit_information = QtCore.pyqtSignal(float, float, float, float, float, float, float, float, float)
    def __init__(self, parent=None):
        super(Twisted_little_window, self).__init__(parent)
        Vlayout = QtWidgets.QVBoxLayout(self)
        Hlayout = QtWidgets.QHBoxLayout()
        Tab = QtWidgets.QTabWidget()
        # Tab.setFont(QtGui.QFont('Times New Roman', 10))
        self.supercell_par_tab = QtWidgets.QWidget()
        self.error_limit_tab = QtWidgets.QWidget()
        Tab.addTab(self.error_limit_tab, "Error limit")
        Tab.addTab(self.supercell_par_tab, "Supercell parameter")
        Normalstrain_label = QtWidgets.QLabel("Normal strain limit(%)：")
        self.Normal_strain_lineedit = QtWidgets.QLineEdit()
        shearstrain_label = QtWidgets.QLabel("Shear strain limit(%)：")
        self.Shear_strain_lineedit = QtWidgets.QLineEdit()
        areaerror_label = QtWidgets.QLabel("Area error limit(%)：")
        self.area_error_lineedit = QtWidgets.QLineEdit()
        Formlayout = QtWidgets.QFormLayout(self.error_limit_tab)
        Formlayout.addRow(Normalstrain_label, self.Normal_strain_lineedit)
        Formlayout.addRow(shearstrain_label, self.Shear_strain_lineedit)
        Formlayout.addRow(areaerror_label, self.area_error_lineedit)
        vacuumdis_label = QtWidgets.QLabel("Vacuum distance(Å) ：")
        self.vacuum_dis_lineedit = QtWidgets.QLineEdit()
        layerdis_label = QtWidgets.QLabel("Layer distance(Å)：")
        self.layer_dis_lineedit = QtWidgets.QLineEdit()
        mingama_label = QtWidgets.QLabel("Minimum sin(γ)：")
        self.min_gamma_lineedit = QtWidgets.QLineEdit()
        self.min_gamma_lineedit.setPlaceholderText('0')
        maxgama_label = QtWidgets.QLabel("Maximum sin(γ)：")
        self.max_singama_lineedit = QtWidgets.QLineEdit()
        self.max_singama_lineedit.setPlaceholderText('1')
        aberror_label = QtWidgets.QLabel("Supercell a,b relative error(%):")
        self.ab_relative_error_lineedit = QtWidgets.QLineEdit()
        rotateangle_label = QtWidgets.QLabel("Rotate angle(°)：")
        self.angle_widget = QtWidgets.QLineEdit()
        Gridlayout = QtWidgets.QGridLayout(self.supercell_par_tab)
        Gridlayout.addWidget(vacuumdis_label, 0, 0, 1, 3)
        Gridlayout.addWidget(self.vacuum_dis_lineedit, 0, 3, 1, 3)
        Gridlayout.addWidget(layerdis_label, 0, 6, 1, 3)
        Gridlayout.addWidget(self.layer_dis_lineedit, 0, 9, 1, 3)
        Gridlayout.addWidget(aberror_label, 1, 0, 1, 3)
        Gridlayout.addWidget(self.ab_relative_error_lineedit, 1, 3, 1, 3)
        Gridlayout.addWidget(rotateangle_label, 1, 6, 1, 3)
        Gridlayout.addWidget(self.angle_widget, 1, 9, 1, 3)
        Gridlayout.addWidget(mingama_label, 2, 0, 1, 3)
        Gridlayout.addWidget(self.min_gamma_lineedit, 2, 3, 1, 3)
        Gridlayout.addWidget(maxgama_label, 2, 6, 1, 3)
        Gridlayout.addWidget(self.max_singama_lineedit, 2, 9, 1, 3)
        self.determine_button = QtWidgets.QPushButton('Start')
        self.cancel_button = QtWidgets.QPushButton('Cancel')
        self.cancel_button.clicked.connect(self.close)
        self.determine_button.clicked.connect(self.determine)
        Hlayout.addWidget(self.determine_button)
        Hlayout.addWidget(self.cancel_button)
        Vlayout.addWidget(Tab)
        Vlayout.addLayout(Hlayout)
        self.setWindowTitle("Twist little angle")
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        self.show()

    def determine(self):
        try:
            self.layerdistance_var = float(self.layer_dis_lineedit.text())
            self.vacuumdistance_var = float(self.vacuum_dis_lineedit.text())
            self.Normal_strain_limit_var = float(self.Normal_strain_lineedit.text()) / 100
            self.Shear_strain_limit_var = float(self.Shear_strain_lineedit.text()) / 100
            self.area_tolerance_var = float(self.area_error_lineedit.text()) / 100
            self.minimum_gamma_var = float(self.min_gamma_lineedit.text())
            self.maximum_gamma_var = float(self.max_singama_lineedit.text())
            self.ab_relative_error_var = float(self.ab_relative_error_lineedit.text()) / 100
            self.angle_var = float(self.angle_widget.text())
            if self.layerdistance_var > 0 and self.vacuumdistance_var >= 0 and 100 > self.Normal_strain_limit_var >= 0 \
                    and 100 > self.area_tolerance_var >= 0 and 1 >= self.maximum_gamma_var \
                    >= self.minimum_gamma_var >= 0 and self.ab_relative_error_var > 0:

                self.signal_emit_information.emit(self.layerdistance_var, self.vacuumdistance_var,
                                                  self.Normal_strain_limit_var, self.Shear_strain_limit_var,
                                                  self.area_tolerance_var,
                                                  self.minimum_gamma_var, self.angle_var,
                                                  self.maximum_gamma_var, self.ab_relative_error_var)
                self.close()
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')


class Classification_window(QtWidgets.QWidget):
    signal_emit_dirs = QtCore.pyqtSignal(str, str)
    input_directory = ""
    output_directory = ""
    def __init__(self, parent=None):
        super(Classification_window, self).__init__(parent)
        Vlayout = QtWidgets.QVBoxLayout(self)
        Hlayout1 = QtWidgets.QHBoxLayout()
        Hlayout2 = QtWidgets.QHBoxLayout()
        Hlayout3 = QtWidgets.QHBoxLayout()
        cifdir_label = QtWidgets.QLabel("CIF Directory:")
        self.cif_dir_lineedit = QtWidgets.QLineEdit()
        self.button1 = QtWidgets.QPushButton("Browse")
        self.button1.clicked.connect(self.CIF_dir)
        self.cif_dir_lineedit.setEnabled(False)
        self.cif_dir_lineedit.setStyleSheet("color: rgb(255, 255, 255);background-color: rgb(35, 38, 41);")
        self.setWindowTitle("Screen out the block 2D material CIF file")
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        Hlayout1.addWidget(cifdir_label)
        Hlayout1.addWidget(self.cif_dir_lineedit)
        Hlayout1.addWidget(self.button1)
        exportaddr_label = QtWidgets.QLabel("Export address:")
        self.export_addr_lineedit = QtWidgets.QLineEdit()
        self.export_addr_lineedit.setStyleSheet("color: rgb(255, 255, 255);background-color: rgb(35, 38, 41);")
        self.export_addr_lineedit.setEnabled(False)
        self.button2 = QtWidgets.QPushButton("Browse")
        self.button2.clicked.connect(self.export_dir)
        Hlayout2.addWidget(exportaddr_label)
        Hlayout2.addWidget(self.export_addr_lineedit)
        Hlayout2.addWidget(self.button2)
        self.button = QtWidgets.QPushButton("Start")
        self.button.clicked.connect(self.determine)
        self.pbar = QtWidgets.QProgressBar()
        Hlayout3.addWidget(self.button)
        Hlayout3.addWidget(self.pbar)
        Vlayout.addLayout(Hlayout1)
        Vlayout.addLayout(Hlayout2)
        Vlayout.addLayout(Hlayout3)
        self.show()

    def CIF_dir(self):
        self.input_directory = QtGui.QFileDialog.getExistingDirectory(self, 'Select Classify directory')
        self.cif_dir_lineedit.setText(self.input_directory)

    def export_dir(self):
        self.output_directory = QtGui.QFileDialog.getExistingDirectory(self, 'Select Export directory')
        self.export_addr_lineedit.setText(self.output_directory)

    def determine(self):
        try:
            if self.input_directory and self.output_directory:
                self.signal_emit_dirs.emit(self.input_directory, self.output_directory)
                # self.close()
            elif not self.input_directory:
                QtWidgets.QMessageBox.warning(self, 'error', 'Choose CIF directory!')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Choose export directory!')
        except Exception as e:
            print(e)


class Change_gama_window(QtWidgets.QWidget):
    """ To change crys' γ"""

    signal_gama_error = QtCore.pyqtSignal(list)
    def __init__(self, parent=None):
        super(Change_gama_window, self).__init__(parent)
        Hlayout = QtWidgets.QHBoxLayout(self)
        objectgama_label = QtWidgets.QLabel("Change γ to (degree)：")
        self.objective_gama_lineedit = QtWidgets.QLineEdit()
        error_label = QtWidgets.QLabel("Error limit(°): ")
        self.error_limit_lineedit = QtWidgets.QLineEdit()
        self.button = QtWidgets.QPushButton("Ok")
        Hlayout.addWidget(objectgama_label)
        Hlayout.addWidget(self.objective_gama_lineedit)
        Hlayout.addWidget(error_label)
        Hlayout.addWidget(self.error_limit_lineedit)
        Hlayout.addWidget(self.button)
        self.button.clicked.connect(self.determine)
        self.setWindowTitle("To change γ")
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        self.show()

    def determine(self):
        try:
            gamadegree = float(self.objective_gama_lineedit.text())
            errorlimit = float(self.error_limit_lineedit.text())
            if 0<gamadegree<180 and errorlimit> 0:
                self.signal_gama_error.emit([gamadegree, errorlimit])
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')



class add_vacuum_layer_window(QtWidgets.QWidget):
    """ The layout of the add vacuum layer window. """
    sinal_vacuum = QtCore.pyqtSignal(float)
    def __init__(self, parent=None):
        super(add_vacuum_layer_window, self).__init__(parent)
        self.setWindowTitle("Add vacuum layer")
        self.setWindowIcon(QtGui.QIcon("Main.png"))

        layout = QtWidgets.QHBoxLayout(self)
        vacuumdis_label = QtWidgets.QLabel("vacuum layer distance(Å)：")
        self.vacuum_dis_lineedit = QtWidgets.QLineEdit()
        self.determinebutton = QtWidgets.QPushButton("OK")
        layout.addWidget(vacuumdis_label)
        layout.addWidget(self.vacuum_dis_lineedit)
        layout.addWidget(self.determinebutton)
        self.determinebutton.clicked.connect(self.determine)
        self.show()

    def determine(self):
        try:
            vacuum_distance = float(self.vacuum_dis_lineedit.text())
            self.sinal_vacuum.emit(vacuum_distance)
            self.close()
        except Exception as e:
            print(e)


class cut_window(QtWidgets.QWidget):
    """ The layout of the cut window. """
    signal_custom_cut_exf = QtCore.pyqtSignal(list, int)
    num = 0
    def __init__(self, cell, parent=None):
        super(cut_window, self).__init__(parent)
        self.cell = cell
        self.resize(550, 500)
        zonglayout = QtWidgets.QVBoxLayout(self)
        layout = QtWidgets.QHBoxLayout()
        self.btn1 = QtWidgets.QRadioButton("xyz")
        self.btn1.setChecked(False)
        self.btn1.toggled.connect(lambda: self.btnstate(self.btn1))
        layout.addWidget(self.btn1)
        # self.btn2 = QtWidgets.QRadioButton("Customize")
        self.btn2 = QtWidgets.QRadioButton("uvw")
        self.btn2.setChecked(False)
        self.btn2.toggled.connect(lambda: self.btnstate(self.btn2))
        layout.addWidget(self.btn2)
        self.btn3 = QtWidgets.QRadioButton("hkl")
        self.btn3.setChecked(False)
        self.btn3.toggled.connect(lambda: self.btnstate(self.btn3))
        layout.addWidget(self.btn3)
        self.setWindowTitle("Cut Style")
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        self.determinebutton = QtWidgets.QPushButton("OK")
        layout.addWidget(self.determinebutton)
        self.determinebutton.clicked.connect(self.determine)
        zonglayout.addLayout(layout)
        self.tablewidget = QtWidgets.QTableWidget()
        self.tablewidget.setRowCount(3)
        self.tablewidget.setColumnCount(3)
        zonglayout.addWidget(self.tablewidget)
        self.show()

    def btnstate(self, btn):
        try:
            self.tablewidget.setRowCount(3)
            self.tablewidget.setColumnCount(3)
            self.choosebutton = btn
            if self.choosebutton.text() == "xyz":
                self.tablewidget.setRowCount(3)
                self.tablewidget.setHorizontalHeaderLabels(['x(Å)', 'y(Å)', 'z(Å)'])
                self.tablewidget.setVerticalHeaderLabels(['a-axis', 'b-axis', 'c-axis'])
            elif self.choosebutton.text() == 'uvw':
                self.tablewidget.setRowCount(3)
                self.tablewidget.setHorizontalHeaderLabels(['u', 'v', 'w'])
                self.tablewidget.setVerticalHeaderLabels(['a', 'b', 'c'])
            elif self.choosebutton.text() == 'hkl':
                self.tablewidget.setHorizontalHeaderLabels(['h', 'k', 'l'])
                self.tablewidget.setRowCount(1)
                self.tablewidget.setVerticalHeaderLabels(['Integer'])
        except Exception as e:
            print(e)

    def determine(self):
        try:
            if self.choosebutton.text() == "xyz":
                cell_par = []
                for i in range(3):
                    l = []
                    for j in range(3):
                        l.append(float(self.tablewidget.item(i, j).text()))
                    cell_par.append(l)
                self.signal_custom_cut_exf.emit(cell_par, 1)
            elif self.choosebutton.text() == 'uvw':
                cell_par = []
                for i in range(3):
                    l = []
                    for j in range(3):
                        l.append(float(self.tablewidget.item(i, j).text()))
                    cell_par.append(l)
                self.signal_custom_cut_exf.emit(cell_par, 2)
            elif self.choosebutton.text() == 'hkl':
                cell_par = []
                for j in range(3):
                    cell_par.append(float(self.tablewidget.item(0, j).text()))
                self.signal_custom_cut_exf.emit(cell_par, 3)
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error', 'Input Error')


class New_file_window(QtWidgets.QWidget):
    """ The layout of the New file window. """
    signal_determine = QtCore.pyqtSignal(str)
    def __init__(self, parent=None):
        print("in")
        super(New_file_window, self).__init__(parent)
        self.setWindowTitle('New cifFile')
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        layout = QtWidgets.QVBoxLayout(self)    # 垂直布局
        self.tablewidget = QtWidgets.QTableWidget()
        Hlayout1 = QtWidgets.QHBoxLayout()
        Hlayout2 =QtWidgets.QHBoxLayout()
        Hlayout2_5 = QtWidgets.QHBoxLayout()
        Hlayout0 = QtWidgets.QHBoxLayout()
        labela = QtWidgets.QLabel("a length(Å): ")
        self.a_lineedit = QtWidgets.QLineEdit()
        labelb = QtWidgets.QLabel("b length(Å): ")
        self.b_lineedit = QtWidgets.QLineEdit()
        labelc = QtWidgets.QLabel("c length(Å): ")
        self.c_lineedit = QtWidgets.QLineEdit()
        labelalpha = QtWidgets.QLabel("α(degree): ")
        self.alpha_lineedit = QtWidgets.QLineEdit()
        labelbetta = QtWidgets.QLabel("β(degree): ")
        self.beta_lineedit = QtWidgets.QLineEdit()
        labelgamma = QtWidgets.QLabel("γ(degree): ")
        self.gama_lineedit = QtWidgets.QLineEdit()
        labelatom_number = QtWidgets.QLabel("Number of crystal cell atoms: ")
        self.atom_number_lineedit = QtWidgets.QLineEdit()
        self.determine_button0 = QtWidgets.QPushButton('OK')
        self.determine_button0.clicked.connect(self.settable_number)
        Hlayout0.addWidget(labelatom_number)
        Hlayout0.addWidget(self.atom_number_lineedit)
        Hlayout0.addWidget(self.determine_button0)
        Hlayout1.addWidget(labela)
        Hlayout1.addWidget(self.a_lineedit)
        Hlayout1.addWidget(labelalpha)
        Hlayout1.addWidget(self.alpha_lineedit)
        Hlayout2.addWidget(labelb)
        Hlayout2.addWidget(self.b_lineedit)
        Hlayout2.addWidget(labelbetta)
        Hlayout2.addWidget(self.beta_lineedit)
        Hlayout2_5.addWidget(labelc)
        Hlayout2_5.addWidget(self.c_lineedit)
        Hlayout2_5.addWidget(labelgamma)
        Hlayout2_5.addWidget(self.gama_lineedit)
        layout.addLayout(Hlayout0)
        layout.addLayout(Hlayout1)
        layout.addLayout(Hlayout2)
        layout.addLayout(Hlayout2_5)
        self.determine_button = QtWidgets.QPushButton('OK')
        self.determine_button.clicked.connect(self.determine)
        Hlayout3 = QtWidgets.QHBoxLayout(self)
        self.name_label = QtWidgets.QLabel("File Name: ")
        self.name_lineedit = QtWidgets.QLineEdit()
        Hlayout3.addWidget(self.name_label)
        Hlayout3.addWidget(self.name_lineedit)
        Hlayout3.addWidget(self.determine_button)
        layout4 = QtWidgets.QHBoxLayout()
        self.btn1 = QtWidgets.QRadioButton("xyz")
        self.btn1.setChecked(False)
        self.btn1.toggled.connect(lambda: self.btnstate(self.btn1))
        layout4.addWidget(self.btn1)
        self.btn2 = QtWidgets.QRadioButton("uvw")
        self.btn2.setChecked(False)
        self.btn2.toggled.connect(lambda: self.btnstate(self.btn2))
        layout4.addWidget(self.btn2)
        layout.addLayout(layout4)
        layout.addWidget(self.tablewidget)
        layout.addLayout(Hlayout3)
        self.show()

    def btnstate(self, btn):
        try:
            self.tablewidget.setColumnCount(4)
            self.choosebutton = btn
            if self.choosebutton.text() == "xyz":
                self.tablewidget.setHorizontalHeaderLabels(["Atom(Formula or atomic number)", 'x(Å)', 'y(Å)', 'z(Å)'])
                self.tablewidget.resizeColumnsToContents()
            elif self.choosebutton.text() == 'uvw':
                self.tablewidget.setHorizontalHeaderLabels(["Atom(Formula or atomic number)", 'u', 'v', 'w'])
                self.tablewidget.resizeColumnsToContents()
        except Exception as e:
            print(e)

    def settable_number(self):
        try:
            number = int(self.atom_number_lineedit.text())
            self.tablewidget.setRowCount(number)
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')
            print(e)

    def determine(self):
        try:
            test = 0
            alpha = float(self.alpha_lineedit.text())
            betta = float(self.beta_lineedit.text())
            gama = float(self.gama_lineedit.text())
            a = float(self.a_lineedit.text())
            b = float(self.b_lineedit.text())
            c = float(self.c_lineedit.text())
            self.par_lst_matrix = self.make_par_lst_matrix(alpha, betta, gama, a, b, c)  # 生成晶胞向量
            test = 1
            self.atom_lst = []
            for row in range(self.tablewidget.rowCount()):
                lst = []
                for column in range(self.tablewidget.columnCount()):
                    if self.tablewidget.item(row, column):
                        if column == 0:
                            lst.append(self.tablewidget.item(row, column).text())
                        else:
                            lst.append(float(self.tablewidget.item(row, column).text()))
                    else:
                        break
                if len(lst) == 4:
                    self.atom_lst.append(lst)
            if len(self.atom_lst) == self.tablewidget.rowCount():
                self.par_lst_species = [_[0] for _ in self.atom_lst]
                if self.choosebutton.text() == 'xyz':
                    self.par_lst_coords = [_[1:4] for _ in self.atom_lst]
                    self.par_lst_coords = [_[1:4] for _ in self.atom_lst]
                else:
                    a = np.array(self.par_lst_matrix[0])
                    b = np.array(self.par_lst_matrix[1])
                    c = np.array(self.par_lst_matrix[2])
                    self.par_lst_coords = [list(atom[1] * a + atom[2] * b + atom[3] * c) for atom in self.atom_lst]
                self.cifname = self.name_lineedit.text()
                if self.cifname == '':
                    QtWidgets.QMessageBox.warning(self, 'error', 'Name Input error!')
                else:
                    self.generate()
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Number input error!')
        except Exception as e:
            print(e)
            if test == 0:
                QtWidgets.QMessageBox.warning(self, 'error', 'Cell parameter Input error!')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Table Pos Input error!')


    def coords_transform(self, A1, B1, C1):
        eps = .0001
        A = A1 / 180 * np.pi
        B = B1 / 180 * np.pi
        C = C1 / 180 * np.pi
        z = np.sqrt(1 - np.cos(B) ** 2 - (np.cos(A) * np.sin(C)) ** 2)
        a_vec = [1, 0, 0]
        b_vec = [np.cos(C), np.sin(C), 0]
        c_vec = [np.cos(B), np.cos(A) * np.sin(C), z]
        par_lst = [a_vec, b_vec, c_vec]

        for vec in par_lst:
            for i in range(3):
                if abs(vec[i]) < eps:
                    vec[i] = 0
        return par_lst

    def num_multiple(self, lst, mult):
        lst_new = []
        for i in lst:
            lst_new.append(i * mult)
        return lst_new

    def make_par_lst_matrix(self, A, B, C, a, b, c):
        self.par_lst = self.coords_transform(A, B, C)
        a_vec = self.num_multiple(self.par_lst[0], a)
        b_vec = self.num_multiple(self.par_lst[1], b)
        c_vec = self.num_multiple(self.par_lst[2], c)
        return [a_vec, b_vec, c_vec]
        # self.par_lst_matrix = [a_vec, b_vec, c_vec]

    def generate(self):
        try:
            self.directory = QtGui.QFileDialog.getExistingDirectory(self, 'Select directory')
            structure = Structure(self.par_lst_matrix, self.par_lst_species, self.par_lst_coords)
            slab = CifWriter(structure, write_magmoms=True)
            slab.write_file(self.directory + '/' + self.cifname + '.cif')
            self.signal_determine.emit(self.directory + '/' + self.cifname + '.cif')
        except Exception as e:
            print(e)

class Create_supercell_window(QtWidgets.QWidget):      # 扩胞
    """ The layout of the creat supercell window. """
    signala = QtCore.pyqtSignal(int)
    signalb = QtCore.pyqtSignal(int)
    signalc = QtCore.pyqtSignal(int)
    def __init__(self, parent=None):
        super(Create_supercell_window, self).__init__(parent)
        self.setWindowTitle('add cell')
        self.setWindowIcon(QtGui.QIcon("Main.png"))

        self.resize(300, 100)
        layout = QtWidgets.QVBoxLayout(self)    # 垂直布局
        Hlayout = QtWidgets.QHBoxLayout()
        self.numberSpinBoxa = QtWidgets.QSpinBox()
        self.numberSpinBoxb = QtWidgets.QSpinBox()
        self.numberSpinBoxc = QtWidgets.QSpinBox()
        self.numberSpinBoxa.setMinimum(1)
        self.numberSpinBoxb.setMinimum(1)
        self.numberSpinBoxc.setMinimum(1)
        self.numberSpinBoxa.valueChanged.connect(self.valuea)
        self.numberSpinBoxb.valueChanged.connect(self.valueb)
        self.numberSpinBoxc.valueChanged.connect(self.valuec)
        labela = QtWidgets.QLabel("a-axis: ")
        labelb = QtWidgets.QLabel("b-axis: ")
        labelc = QtWidgets.QLabel("c-axis: ")
        Hlayout.addWidget(labela)
        Hlayout.addWidget(self.numberSpinBoxa)
        Hlayout.addWidget(labelb)
        Hlayout.addWidget(self.numberSpinBoxb)
        Hlayout.addWidget(labelc)
        Hlayout.addWidget(self.numberSpinBoxc)
        self.setUnitCell = QtWidgets.QPushButton("Set UnitCell")
        self.determine = QtWidgets.QPushButton("OK")
        layout.addLayout(Hlayout)
        Hlayout1 = QtWidgets.QHBoxLayout()
        Hlayout1.addWidget(self.setUnitCell)
        Hlayout1.addWidget(self.determine)
        layout.addLayout(Hlayout1)
        self.setWindowTitle("Create supercell")
        self.setWindowIcon(QtGui.QIcon("Main.png"))
        self.show()

    def valuea(self):
        self.signala.emit(self.numberSpinBoxa.value())

    def valueb(self):
        self.signalb.emit(self.numberSpinBoxb.value())

    def valuec(self):
        self.signalc.emit(self.numberSpinBoxc.value())

class set_cell_window(QtWidgets.QWidget):
    signal_emit_information = QtCore.pyqtSignal(int, float, float, float, float, float, float)
    def __init__(self, cell_par):
        super(set_cell_window, self).__init__()
        self.cell_par = cell_par
        self.setWindowTitle("Set Cell")
        self.setWindowIcon(QtGui.QIcon("Main1.png"))
        Vlayout = QtWidgets.QVBoxLayout(self)
        Hlayout0 = QtWidgets.QHBoxLayout()
        self.btn1 = QtWidgets.QRadioButton("Atoms move.")
        self.btn1.setChecked(False)
        self.btn1.toggled.connect(lambda: self.btnstate(self.btn1))
        Hlayout0.addWidget(self.btn1)
        self.btn2 = QtWidgets.QRadioButton("Atoms don't move.")
        self.btn2.setChecked(False)
        self.btn2.toggled.connect(lambda: self.btnstate(self.btn2))
        Hlayout0.addWidget(self.btn2)
        button = QtWidgets.QPushButton("Ok")
        button.clicked.connect(self.determine)
        Hlayout0.addWidget(button)
        Hlayout1 = QtWidgets.QHBoxLayout()
        Hlayout2 = QtWidgets.QHBoxLayout()
        Hlayout3 = QtWidgets.QHBoxLayout()
        labela = QtWidgets.QLabel("a(Å): ")
        self.a_lineedit = QtWidgets.QLineEdit(str(self.cell_par[0]))
        labelalpha = QtWidgets.QLabel("α(degree): ")
        self.alpha_lineedit = QtWidgets.QLineEdit(str(self.cell_par[3]))
        labelb = QtWidgets.QLabel("b(Å): ")
        self.b_lineedit = QtWidgets.QLineEdit(str(self.cell_par[1]))
        labelbetta = QtWidgets.QLabel("β(degree): ")
        self.betta_lineedit = QtWidgets.QLineEdit(str(self.cell_par[4]))
        labelc = QtWidgets.QLabel("c(Å): ")
        self.c_lineedit = QtWidgets.QLineEdit(str(self.cell_par[2]))
        labelgama = QtWidgets.QLabel("γ(degree): ")
        self.gama_lineedit = QtWidgets.QLineEdit(str(self.cell_par[5]))
        Hlayout1.addWidget(labela)
        Hlayout1.addWidget(self.a_lineedit)
        Hlayout1.addWidget(labelalpha)
        Hlayout1.addWidget(self.alpha_lineedit)
        Hlayout2.addWidget(labelb)
        Hlayout2.addWidget(self.b_lineedit)
        Hlayout2.addWidget(labelbetta)
        Hlayout2.addWidget(self.betta_lineedit)
        Hlayout3.addWidget(labelc)
        Hlayout3.addWidget(self.c_lineedit)
        Hlayout3.addWidget(labelgama)
        Hlayout3.addWidget(self.gama_lineedit)
        Vlayout.addLayout(Hlayout0)
        Vlayout.addLayout(Hlayout1)
        Vlayout.addLayout(Hlayout2)
        Vlayout.addLayout(Hlayout3)
        self.show()

    def btnstate(self, btn):
        try:
            self.choosebutton = btn
            if self.choosebutton.text() == "Atoms move.":
                self.move_num = 1
            elif self.choosebutton.text() == "Atoms don't move.":
                self.move_num = 0
        except Exception as e:
            print(e)

    def determine(self):
        try:
            a = float(self.a_lineedit.text())
            b = float(self.b_lineedit.text())
            c = float(self.c_lineedit.text())
            alpha = float(self.alpha_lineedit.text())
            betta = float(self.betta_lineedit.text())
            gama = float(self.gama_lineedit.text())
            self.signal_emit_information.emit(self.move_num, a, b, c, alpha, betta, gama)
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')

class add_atomwindow(QtWidgets.QWidget):
    signal_emit_information = QtCore.pyqtSignal(int, list)
    def __init__(self):
        super(add_atomwindow, self).__init__()
        self.setWindowTitle("Add an atom")
        self.setWindowIcon(QtGui.QIcon("Main1.png"))
        self.resize(600, 300)
        Vlayout = QtWidgets.QVBoxLayout(self)
        Hlayout0 = QtWidgets.QHBoxLayout()
        label_atomicnumber = QtWidgets.QLabel("Atomic number: ")
        self.atomicnumber_lineedit = QtWidgets.QLineEdit()
        self.element_button = QtWidgets.QPushButton("...")
        self.element_button.clicked.connect(self.element_table_show)
        self.button = QtWidgets.QPushButton("Ok")
        Hlayout0.addWidget(label_atomicnumber)
        Hlayout0.addWidget(self.atomicnumber_lineedit)
        Hlayout0.addWidget(self.element_button)
        Hlayout0.addWidget(self.button)
        Vlayout.addLayout(Hlayout0)
        self.tablewidget = QtWidgets.QTableWidget()
        self.tablewidget.setRowCount(1)
        self.tablewidget.setColumnCount(3)
        self.tablewidget.setHorizontalHeaderLabels(['x(Å)', 'y(Å)', 'z(Å)'])
        self.tablewidget.setVerticalHeaderLabels(['Position'])
        Vlayout.addWidget(self.tablewidget)
        self.show()
        self.button.clicked.connect(self.determine)

    def element_table_show(self):
        try:
            self.periodic_window = add_atom_periodic_window()
            self.periodic_window.signal_atomic_number.connect(self.setAtomic_number)
        except Exception as e:
            print(e)

    def setAtomic_number(self, atomic_number):
        try:
            self.atomicnumber_lineedit.setText(str(atomic_number))
        except:
            pass

    def determine(self):
        try:
            atomic_number = int(self.atomicnumber_lineedit.text())
            pos = []
            for j in range(3):
                pos.append(float(self.tablewidget.item(0, j).text()))
            self.signal_emit_information.emit(atomic_number, pos)
            self.close()
        except:
            QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')

class add_atom_periodic_window(periodic_element_table_window):
    signal_atomic_number = QtCore.pyqtSignal(int)
    def __init__(self):
        super(add_atom_periodic_window, self).__init__([], [])
        self.setWindowTitle("Add an atom")

    def determine_element_formula(self):
        try:
            sender = self.sender()
            print(sender.text())
            self.number = crys_data.Element_formula.index(sender.text())
        except Exception as e:
            print(e)

    def determine(self):
        try:
            self.signal_atomic_number.emit(self.number)
            self.close()
        except:
            pass

class Replace_Atom_window(periodic_element_table_window):
    signal_atomic_number = QtCore.pyqtSignal(int)
    def __init__(self):
        super(Replace_Atom_window, self).__init__([], [])
        self.setWindowTitle("Replace atoms")

    def determine_element_formula(self):
        try:
            sender = self.sender()
            print(sender.text())
            self.number = crys_data.Element_formula.index(sender.text())
        except Exception as e:
            print(e)

    def determine(self):
        try:
            self.signal_atomic_number.emit(self.number)
            self.close()
        except:
            pass

"""
Main Function of Matlego 1.5, depend mainly on ase and pyqtgraph.
MatLego:1.5
Finish time: 2019/3/5
Main function: to build up model for hetero-junction, especially 2D electron Devices.

"""


class mywindow(QtWidgets.QMainWindow, Ui_MainWindow):  # 实现代码与界面的分离
    """  The fuction of the Mainwindow. """
    datanumber = 1     # database
    row_lst = []  # database
    true_id_lst = []      # database
    dic_Formula_Atoms = {}  # key = formula value = Atoms_object
    dic_Hetero_junction = {}    # key = formula value = Object Hetero_junction
    cell_view_num = 0        # whether draw cell or not
    coordinatesystem_view_num = 2  # Draw cartestian coordinate first
    coordinatesystem_fixed_num = 0
    grid_itemnum = 0         # 是否画网格
    plot_num = 0             # Ball-stick(1), Ball-filling(0)
    max_number_show = 1500
    show_atom_index_judge = False
    show_atom_element_judge = False
    exit_box_show = True
    up_down = False
    remove_atom_index_num = None
    def __init__(self):
        super(mywindow, self).__init__()
        self.key_double_click = True
        self.setupUi(self)
        # self.exit_message_window = exit_message_window()
        self.view_widget.installEventFilter(self)
        self.setWindowIcon(QtGui.QIcon("Main1.png"))
        self.setWindowTitle("MatLego")
        self.actionNew.triggered.connect(self.newfile)
        self.actionNewToolBar.triggered.connect(self.newfile)
        self.actionOpen.triggered.connect(self.openfile)
        self.actionOpenToolBar.triggered.connect(self.openfile)
        self.actionExit_toolbar.triggered.connect(self.close_)
        self.actionExit.triggered.connect(self.close_)
        self.actionCut.triggered.connect(self.cut_customize)
        self.actionRotate.triggered.connect(self.Rotate)
        # 批量操作
        self.actionTraversal_rotate.triggered.connect(self.traversalrotate)
        self.actionTraversal_Stack1.triggered.connect(self.traversalstack1)
        self.actionTraversal_Stack2.triggered.connect(self.traversalstack2)
        self.actionTraversal_make_devices.triggered.connect(self.traversal_devices)
        self.actionTraversal_Cut.triggered.connect(self.traversalcut)
        self.actionCutToolBar.triggered.connect(self.judgetool)
        self.actionExport.triggered.connect(self.Export)
        self.actionImport.triggered.connect(self.Import)
        self.actionStack2.triggered.connect(self.Stack_main_two)
        self.actionStack1.triggered.connect(self.Stack_main_one)     # one stack direction is known
        self.actionStack2Toolbar.triggered.connect(self.Stack_main_two)  # two stack direction is known
        self.actiontwist_little_angle.triggered.connect(self.twist_little_angle)
        self.actionMulti_layer_opitmization.triggered.connect(self.multi_layer_optimize)
        self.actionClassification.triggered.connect(self.classification_main)
        # two kinds
        self.actiondrag_rectangle.triggered.connect(self.dragrectangle)
        self.actionNormalchoose.triggered.connect(self.normalchoose)
        self.actionSetlayer.triggered.connect(self.setlayer)         # tool setlayer
        self.actionMovelayer.triggered.connect(self.movelayer)       # tool movelayer
        self.actionleft_screen.triggered.connect(self.left_screen)
        self.actionright_screen.triggered.connect(self.right_screen)
        self.actionup_screen.triggered.connect(self.up_screen)
        self.actiondown_screen.triggered.connect(self.down_screen)
        self.actionmiddle_screen.triggered.connect(self.middle_screen)
        self.actionleft_screen.triggered.connect(self.left_screen)
        self.actionright_screen.triggered.connect(self.right_screen)
        self.actionup_screen.triggered.connect(self.up_screen)
        self.actiondown_screen.triggered.connect(self.down_screen)
        self.actionmiddle_screen.triggered.connect(self.middle_screen)
        self.actiontranslate_choose.triggered.connect(self.translate_choose)
        self.action_viewfrom_a.triggered.connect(self.viewfroma)
        self.action_viewfrom_b.triggered.connect(self.viewfromb)
        self.action_viewfrom_c.triggered.connect(self.viewfromc)
        self.actionSave.triggered.connect(self.Save_cif)
        self.actionViewCell.triggered.connect(self.view_cell)

        self.actionNo_coordinate_system.triggered.connect(self.No_coordinate_system)
        self.actionViewCoordinate_System.triggered.connect(self.viewCoordinate_system)
        self.actionViewCoordinate_cell.triggered.connect(self.viewCoordinatcell_system)
        self.action3D_coordinate.triggered.connect(self.threeD_coordinate)
        self.actionFixed_coordinate.triggered.connect(self.fixed_coordinate)
        self.actionViewGrid.triggered.connect(self.viewgrid)
        self.actionshow_None.triggered.connect(self.show_None)
        self.actionshow_atom_index.triggered.connect(self.show_atom_index)
        self.actionshow_atom_element.triggered.connect(self.show_atom_element)
        # self.actionshow_coordinate_system.triggered.connect(self.show_coordinate_label)
        # right click，right_menubar
        self.view_widget.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.view_widget.customContextMenuRequested['QPoint'].connect(self.rightMenuShow_view_widget)
        self.actionViewDatabase.triggered.connect(self.viewdatabase)
        self.actionViewProject.triggered.connect(self.viewproject)
        self.actionViewText.triggered.connect(self.viewtext)
        self.project_tree.clicked.connect(self.onTreeClicked)
        self.project_tree.doubleClicked.connect(self.doubleclickedontree)
        self.project_tree.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.project_tree.customContextMenuRequested['QPoint'].connect(self.rightMenuShow_tree)
        self.datatable.itemDoubleClicked.connect(self.double_clicked_on_datatable)
        self.datatable.itemClicked.connect(self.clicked_on_datatable)
        self.Device_tablewidget.itemClicked.connect(self.clicked_on_tab3wg2)
        self.object_tablewidget.itemClicked.connect(self.handleItemClick_on_object)  # object box click
        self.object_tablewidget.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.object_tablewidget.customContextMenuRequested['QPoint'].connect(self.rightMenuShow_Object)
        self.calculate_distance_button.clicked.connect(self.calculatedistance)
        self.calculate_degree_button.clicked.connect(self.calculatedegree)
        self.calculate_vector_button.clicked.connect(self.calculatevector)
        self.filterbutton.clicked.connect(self.filterfrom_to)
        self.ballfilling_button.clicked.connect(self.ballfillingtype)
        self.ball_stick_button.clicked.connect(self.ball_and_stick)
        self.stick_button.clicked.connect(self.Stick)
        self.searchLineEdit.textChanged.connect(self.search)
        self.searchdatabase_LineEdit.textChanged.connect(self.searchdatabase)
        self.btn1.toggled.connect(lambda: self.btnstate(self.btn1))
        self.btn2.toggled.connect(lambda: self.btnstate(self.btn2))
        self.btn3.toggled.connect(lambda: self.btnstate(self.btn3))
        # Two viewports
        self.actionPerspective.triggered.connect(self.viewport_perspective)
        self.actionOrthogonal.triggered.connect(self.viewport_orthogonal)
        # database
        self.action_db_connect.triggered.connect(self.connect_db)
        self.action_cif_to_db.triggered.connect(self.database_cif_to_db_main)    # create a database
        self.action_add_data.triggered.connect(self.use_cif_add_data)                 # Add data to database
        self.action_setting_atom_par.triggered.connect(self.edit_atom_para)
        self.action_setting_displayed_par.triggered.connect(self.edit_displayed_atom)
        self.Atoms_color = crys_data.atom_color
        self.Atom_radii = crys_data.atom_radii_lst
        self.datatable.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.datatable.customContextMenuRequested['QPoint'].connect(self.rightMenuShow_database)
        self.view_widget.opts['distance'] = 50
        self.view_widget.opts['fov'] = 150


    def closeEvent(self, event):
        # Ask the user if they actually want to quit using custom message box
        try:
            print("in")
            qm = QtWidgets.QMessageBox.question(self, 'Question', 'Are you sure you want to exit MatLego?')
            if qm == QtWidgets.QMessageBox.Yes:
                super().close()
            else:
                event.ignore()
        except Exception as e:
            print(e)

    def close_(self):
        try:
            qm = QtWidgets.QMessageBox.question(self, 'Question',
                                                'Are you sure you want to exit MatLego?')
            if qm == QtWidgets.QMessageBox.Yes:
                self.close()
        except Exception as e:
            print(e)

    def edit_displayed_atom(self):
        try:
            self.edit_display_par_window.close()
        except:
             pass
        try:
            if self.view_widget.num_mouse_track == 0:
                fov = self.view_widget.opts['fov']
                distance = self.view_widget.opts['distance']
                self.edit_display_par_window = edit_displayed_par_window(self.max_number_show,
                                                                         self.view_widget.set_font_size, fov, distance)
                self.edit_display_par_window.signal_emit_information.connect(self.determine_edit_display_par)
                self.edit_display_par_window.signal_emit_fontsize.connect(self.edit_fontsize)
                self.edit_display_par_window.signal_emit_fov_distance.connect(self.edit_fov_distance)
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose drag mode.')
        except Exception as e:
            print(e)

    def edit_fov_distance(self, fov, distance):
        try:
            self.view_widget.opts['fov'] = fov
            self.view_widget.opts['distance'] = distance
            self.view_widget.update()
        except Exception as e:
            print(e)

    def edit_fontsize(self, fontsize):
        try:
            self.view_widget.set_font_size = fontsize
            self.view_widget.update()
        except Exception as e:
            print(e)

    def determine_edit_display_par(self, max_number, fontsize, fov, distance):
        try:
            self.max_number_show = max_number
            self.view_widget.set_font_size = fontsize
            self.view_widget.opts['fov'] = fov
            self.view_widget.opts['distance'] = distance
            self.view_widget.update()
        except Exception as e:
            print(e)

    def edit_atom_para(self):
        try:
            self.edit_atom_para_window.close()
        except:
            pass
        try:
            self.edit_atom_para_window = periodic_element_table_window(self.Atoms_color, self.Atom_radii)
            self.edit_atom_para_window.signal_emit_atom_par.connect(self.edit_color_radii_lst)
        except Exception as e:
            print(e)

    def edit_color_radii_lst(self, color_lst, radii_lst):
        try:
            self.Atoms_color = color_lst
            self.Atom_radii = radii_lst
            self.plot(self.Atomsobject)
        except Exception as e:
            print(e)

    def rightMenuShow_database(self):
        self.row_num = None
        self.column_num = None
        for i in self.datatable.selectionModel().selection().indexes():
            self.row_num = i.row()
            self.column_num = i.column()
        rightMenu = QtWidgets.QMenu(self.menuBar)
        self.Edit_Action = QtWidgets.QAction(self)
        self.Edit_Action.setText("Edit")
        self.Edit_Action.triggered.connect(self.beforeEdit)
        self.copy_Action = QtWidgets.QAction(self)
        self.copy_Action.setText("Copy")
        self.copy_Action.triggered.connect(self.copy_database)
        self.copyline_Action = QtWidgets.QAction(self)
        self.copyline_Action.setText("Copy a row")
        self.copyline_Action.triggered.connect(self.copy_linedatabase)
        self.paste_Action = QtWidgets.QAction(self)
        self.paste_Action.setText("Paste")
        self.paste_Action.triggered.connect(self.paste_database)
        self.append_Action = QtWidgets.QAction(self)
        self.append_Action.setText("Append data")
        self.append_Action.triggered.connect(self.append_data)
        self.add_line_Action = QtWidgets.QAction(self)
        self.add_line_Action.setText("Add a row")
        self.add_line_Action.triggered.connect(self.add_line)
        rightMenu.addAction(self.Edit_Action)
        rightMenu.addSeparator()
        rightMenu.addAction(self.copy_Action)
        rightMenu.addAction(self.copyline_Action)
        rightMenu.addAction(self.paste_Action)
        rightMenu.addSeparator()
        rightMenu.addAction(self.append_Action)
        rightMenu.addAction(self.add_line_Action)
        rightMenu.exec_(QtGui.QCursor.pos())

    def copy_linedatabase(self):
        try:
            self.copy_text = ""
            self.informatino_lst = []
            for column in range(self.datatable.columnCount()):
                self.informatino_lst.append(self.datatable.item(self.row_num, column).text())
        except Exception as e:
            print(e)

    def paste_database(self):
        try:
            item = self.datatable.item(self.row_num, self.column_num)
            if self.copy_text:
                self.double_clicked_on_datatable(item, paste=True)
            elif self.informatino_lst:
                self.double_clicked_on_datatable(item, pasteline=True)

        except Exception as e:
            print(e)
    def copy_database(self):
        try:
            self.informatino_lst = []
            self.copy_text = self.datatable.item(self.row_num, self.column_num).text()
        except Exception as e:
            print(e)

    def add_line(self):
        try:
            if self.datanumber == 1:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose Hetero-junction or Device to add line.')
            elif self.datanumber == 2:    # hetero-junction
                data = self.database_class.info_data_base
                if self.row_lst:
                    row = self.row_lst[self.row_num]
                else:
                    row = self.row_num
                Currentdata = list(data[row])
                informatino_lst = Currentdata[1:10]
                header_lst = ['ID', 'Formula', 'A', 'B', 'C', 'Angle_a', 'Angle_b', 'Angle_c', 'Volume']
                self.after_append_data_window(informatino_lst[0], informatino_lst[1])
                self.edit_database(self.row_num + 1, header_lst, informatino_lst, addline=True)
                self.searchdatabase()
        except Exception as e:
            print(e)

    def beforeEdit(self):
        try:
            item = self.datatable.item(self.row_num, self.column_num)
            self.double_clicked_on_datatable(item)
        except Exception as e:
            print(e)

    def append_data(self):
        try:
            self.append_data_window.close()
        except:
            pass

        try:
            self.database_class.exhibit_all_info(fenkuai=False)
            ID_lst = list(set([info[1] for info in self.database_class.info_data_base]))
            self.append_data_window = Append_data_window(ID_lst)
            self.append_data_window.signal_emit.connect(self.after_append_data_window)
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error',
                                          'Connect database please.')
    def after_append_data_window(self, ID, formula):
        try:
            info = self.database_class.info_data_base
            TrueID = max(_[0] for _ in info) + 1
            self.database_class.insert_row_into_test(TrueID, ID, formula)
            self.data_base_init()
        except Exception as e:
            print(e)

    def clicked_on_tab3wg2(self, item):
        try:
            print("text = ", item.text())
            if item.column() == 4:
                if item.text() != "":
                    l = item.text().split('+')
                    address_lst = []
                    for address in l:
                        if os.path.exists(address):
                            address_lst.append(address)
                    self.show_image_table.setColumnCount(1)
                    self.show_image_table.setRowCount(len(address_lst))
                    self.show_image_table.setIconSize(QtCore.QSize(300, 200))
                    self.show_image_table.setColumnWidth(0, 300)
                    for index, address_ in enumerate(address_lst):
                        newitem = QtWidgets.QTableWidgetItem()
                        try:
                            Icon = QtGui.QIcon(address_)
                            newitem.setIcon(QtGui.QIcon(Icon))
                        except Exception as e:
                            print(e)
                        newitem.setFlags(QtCore.Qt.ItemIsEnabled)
                        self.show_image_table.setItem(index, 0, newitem)
                        self.show_image_table.setRowHeight(index, 200)
        except Exception as e:
            print(e)

    def clicked_on_datatable(self, item):
        try:
            print("text = ", item.text())
            if self.datanumber == 3 and item.column() == 4:
                if item.text() != "":
                    l = item.text().split('+')
                    address_lst = []
                    for address in l:
                        if os.path.exists(address):
                            address_lst.append(address)
                    self.show_image_table.setColumnCount(1)
                    self.show_image_table.setRowCount(len(address_lst))
                    self.show_image_table.setIconSize(QtCore.QSize(300, 200))
                    self.show_image_table.setColumnWidth(0, 300)
                    for i, address in enumerate(address_lst):
                        newitem = QtWidgets.QTableWidgetItem()
                        try:
                            Icon = QtGui.QIcon(address)
                            newitem.setIcon(QtGui.QIcon(Icon))
                        except Exception as e:
                            print(e)
                        newitem.setFlags(QtCore.Qt.ItemIsEnabled)
                        self.show_image_table.setItem(i, 0, newitem)
                        self.show_image_table.setRowHeight(i, 200)
        except Exception as e:
            print(e)

    def double_clicked_on_datatable(self, item, paste=False, pasteline=False):
        try:
            self.database_class.exhibit_all_info(fenkuai=False)
            ID_lst = list(set([_[1] for _ in self.database_class.info_data_base]))
            if self.true_id_lst:
                if self.Reflect_id_lst:
                    True_ID = self.Reflect_id_lst[item.row()]
                else:
                    True_ID = self.true_id_lst[item.row()]
            else:
                True_ID = self.database_class.info_data_base[item.row()]

            for data in self.database_class.info_data_base:
                if data[0] == True_ID:
                    break

            ID_lst.remove(data[1])
            information_lst = []
            header_lst = []
            for i in range(self.datatable.columnCount()):
                information_lst.append(self.datatable.item(item.row(), i).text())
                header = self.datatable.horizontalHeaderItem(i).text()
                header_lst.append(header)
            if paste:
                information_lst[item.column()] = self.copy_text
            elif pasteline:
                if len(information_lst) == len(self.informatino_lst):
                    information_lst = deepcopy(self.informatino_lst)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error',
                                                  'Different shape.')
            try:
                self.edit_information_window.close()
            except:
                pass
            self.edit_information_window = edit_data_table_window(ID_lst, header_lst, information_lst, item.row())
            self.edit_information_window.signal_emit_information.connect(self.edit_database)
            self.edit_information_window.signal_emit_delete.connect(self.delete_database)
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error',
                                          'Connect database please.')

    def delete_database(self, row):
        try:
            if self.datanumber != 1:
                if self.true_id_lst:
                    True_ID = self.true_id_lst[row]
                else:
                    self.database_class.exhibit_all_info(fenkuai=False)
                    True_ID = self.database_class.info_data_base[row][0]
                self.database_class.delete_row(True_ID, ID=True)
                try:
                    if self.true_id_lst:
                        self.true_id_lst.remove(True_ID)
                except Exception as e:
                    print(e)
                self.data_base_init()
            else:
                if self.true_id_lst:
                    true_ID = self.Reflect_id_lst[row]
                else:
                    self.database_class.exhibit_all_info(fenkuai=True)
                    true_ID = self.database_class.info_data_base[row][0][0]
                self.database_class.exhibit_all_info(fenkuai=False)
                datafenkuai_false = self.database_class.info_data_base
                for data in datafenkuai_false:
                    if data[0] == true_ID:
                        break
                for True_ID in (_ for _ in datafenkuai_false if _[1] == data[1]):
                    self.database_class.delete_row(True_ID, ID=True)
                    try:
                        if self.true_id_lst:
                            self.true_id_lst.remove(True_ID)
                    except Exception as e:
                        print(e)
                self.data_base_init()
        except Exception as e:
            print(e)

    def edit_database(self, row, header_lst, information_lst, addline=False):
        try:
            if self.datanumber != 1:
                if self.true_id_lst:
                    true_ID = self.true_id_lst[row]
                else:
                    self.database_class.exhibit_all_info(fenkuai=False)
                    true_ID = self.database_class.info_data_base[row][0]
                for i in range(len(information_lst)):
                    try:
                        information_lst[i] = "'" + information_lst[i] + "'"
                        self.database_class.insert_data(true_ID, header_lst[i], information_lst[i], ID=True)
                    except Exception as e:
                        print(e)
            else:            # delete crystal information
                if self.true_id_lst:
                    true_ID = self.Reflect_id_lst[row]
                else:
                    self.database_class.exhibit_all_info(fenkuai=True)
                    true_ID = self.database_class.info_data_base[row][0][0]
                self.database_class.exhibit_all_info(fenkuai=False)
                datafenkuai_false = self.database_class.info_data_base
                for data in datafenkuai_false:
                    if data[0] == true_ID:
                        break
                for true_ID in (_ for _ in datafenkuai_false if _[1] == data[1]):
                    for i in range(len(information_lst)):
                        try:
                            information = "'" + information_lst[i] + "'"
                            self.database_class.insert_data(true_ID, header_lst[i], information, ID=True)
                        except Exception as e:
                            print(e)
            if addline == False:
                QtWidgets.QMessageBox.information(self, 'Information:', 'Edit successfully')
            self.data_base_init()
        except Exception as e:
            print(e)


    def connect_db(self):
        try:
            self.database_dir, filetype = QtWidgets.QFileDialog.getOpenFileName(self, '选择文件', '', 'files(*.db)')
            self.data_base_init()
            if self.database_dir:
                QtWidgets.QMessageBox.information(self, 'Information:', 'Connect successfully')
        except Exception as e:
            print(e)

    def data_base_init(self):
        try:
            self.database_class = Database(self.database_dir)
            self.database_class.exhibit_all_info()
            print(self.database_class.info_data_base)
            self.database_class.exhibit_all_info(fenkuai=True)
            print(self.database_class.info_data_base)
            if self.datanumber == 1:
                self.btn1.setChecked(True)
                self.databaseshow(1)
            elif self.datanumber == 2:
                self.btn2.setChecked(True)
                self.databaseshow(2)
            elif self.datanumber == 3:
                self.btn3.setChecked(True)
                self.databaseshow(3)
            self.plot(self.Atomsobject)
        except Exception as e:
            print(e)

    def use_cif_add_data(self):
        cif_dir, filetype = QtWidgets.QFileDialog.getOpenFileName(self, '选择文件', '', 'files(*.cif)')
        try:
            self.database_class.exhibit_all_info()
            info = self.database_class.info_data_base
            self.insert_into_db(cif_dir, self.database_dir, info[-1][0] + 1)
            # self.database_cif_to_db(self.database_dir, cif_dir)
            self.data_base_init()
        except Exception as e:
            print(e)

    def database_cif_to_db_main(self):
        try:
            cif_dir = QtGui.QFileDialog.getExistingDirectory(self, 'Select cif files directory')
            r_paths = os.listdir(cif_dir)
            paths = [cif_dir + '/' + x for x in r_paths if x != '.DS_Store']  # '/Users/mac/Desktop/ground'
            # create a db
            db_dir = cif_dir + '/default.db'
            self.set_up_db(db_dir)
            # emit_zong_num
            zong_num = len(paths)
            try:
                self.create_data_base_window.close()
            except:
                pass
            self.create_data_base_window = create_database_window()
            self.create_data_base_window.Signal_emit_number.connect(self.create_data_base_window.pbar_number)
            i = 1
            num = 0
            while paths:
                try:
                    dir = paths.pop(0)
                    self.insert_into_db(dir, db_dir, i)
                    i += 1
                    num += 1
                    QtWidgets.QApplication.processEvents()
                    self.create_data_base_window.Signal_emit_number.emit(num / zong_num * 100)
                except Exception as e:
                    print(e)
            self.create_data_base_window.close()
            QtWidgets.QMessageBox.information(self, 'Information:',
                                               'The address of the database is ' + '\n' + db_dir)
        except Exception as e:
            print(e)

    def insert_into_db(self, cif_dir, db_dir, i):  # 录入数据库材料的cif文件路径， 自己建立的材料ID
        '''
        由于把 本地某个cif文件 的 物理性质 储存进相应的表格。
        :param cif_dir: cif文件所在路径
        :param db_dir: 数据库的路径
        :param i: 材料ID(默认从一递增，也可以自己输入)
        :return: None
        '''
        try:
            f = open(cif_dir, 'r')  # dir = '/Users/mac/Desktop/1.cif'
            content = f.read()
            Atomsobject = read(cif_dir)
            r_formula = Atomsobject.get_chemical_formula(mode='hill')
            formula = "'" + r_formula + "'"
            findall('_symmetry_space_group_name_H-M(.*?)\n', content)[0].strip()
            a = float(findall('_cell_length_a(.*?)\n', content)[0].strip())
            b = float(findall('_cell_length_b(.*?)\n', content)[0].strip())
            c = float(findall('_cell_length_c(.*?)\n', content)[0].strip())
            angle_alpha = float(findall('_cell_angle_alpha(.*?)\n', content)[0].strip())
            angle_beta = float(findall('_cell_angle_beta(.*?)\n', content)[0].strip())
            angle_gamma = float(findall('_cell_angle_gamma(.*?)\n', content)[0].strip())
            volume = float(findall('_cell_volume(.*?)\n', content)[0].strip())
            conn = sqlite3.connect(db_dir)
            cn = conn.cursor()

            cn.execute(
                "INSERT INTO test (True_ID, ID, Formula, A, B, C, Angle_a, Angle_b, Angle_C, Volume) VALUES (%d, %d, %s, %g, %g, %g, %g, %g, %g, %g)" % (
                    i, i, formula, a, b, c, angle_alpha, angle_beta, angle_gamma, volume))

            conn.commit()
            conn.close()
            f.close()
        except Exception as e:
            print(e)

    def set_up_db(self, dir):  # dir: 数据库的路径 -- /Users/mac/Desktop/crystal.db
        '''
        用于建立 数据库 和 储存物理性质的表格 。
        :param dir: 函数所建立数据库的路径
        :return: None
        '''
        conn = sqlite3.connect(dir)
        cn = conn.cursor()
        cn.execute('''CREATE TABLE test
                                                         (True_ID  INT PRIMARY KEY    NOT NULL,
                                                          ID                  INT     NOT NULL,
                                                          Formula             TEXT    NOT NULL,
                                                          A                   TEXT    ,
                                                          B                   TEXT    ,
                                                          C                   TEXT    ,
                                                          Angle_a             TEXT    ,
                                                          Angle_b             TEXT    ,
                                                          Angle_c             TEXT    ,
                                                          Volume              TEXT    ,

                                                          Hetero_junction     TEXT,
                                                          Optimal_Match       TEXT,
                                                          Layer               TEXT,
                                                          Binding_energy      TEXT,
                                                          Schottky_barrier    TEXT,

                                                          Device              Text, 
                                                          Optimal_MatchD      Text, 
                                                          LayerD              Text,   
                                                          Schottky_barrierD   TEXT,  
                                                          Image_dir           Text)''')
        conn.commit()
        conn.close()

    def viewport_perspective(self):
        try:
            distance = self.view_widget.opts['distance']
            fov = self.view_widget.opts['fov']
            width = distance * np.tan(fov / 2 / 180 * np.pi)
            self.view_widget.opts['fov'] = 150
            self.view_widget.opts['distance'] = int(width / np.tan(150 / 2 / 180 * np.pi))
            self.plot(self.Atomsobject)
        except Exception as e:
            print(e)

    def viewport_orthogonal(self):
        try:
            distance = self.view_widget.opts['distance']
            fov = self.view_widget.opts['fov']
            width = distance * np.tan(fov / 2 / 180 * np.pi)
            self.view_widget.opts['fov'] = 60
            self.view_widget.opts['distance'] = int(width / np.tan(60 / 2 / 180 * np.pi))
            self.plot(self.Atomsobject)
        except Exception as e:
            print(e)

    def Rotate(self):
        """To stack along different orientations."""
        try:
            formula_lst = []
            itemlst = []
            for item in self.project_tree.selectedItems():
                formula_lst.append(item.text(1))
                itemlst.append(item)
            if len(formula_lst) != 2:
                QtWidgets.QMessageBox.warning(self, 'error',
                                              'Please choose 2 crystal in Projectbox' + '\n' + '(Press Ctrl or Shift).')
            elif itemlst[0].text(0) == 'bulk' or itemlst[1].text(0) == 'bulk':
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose stack or layer(doubleclick to change)')
            else:
                self.rotate_obj1 = self.dic_Formula_Atoms[formula_lst[0]]  # Atoms1，对象
                self.rotate_obj2 = self.dic_Formula_Atoms[formula_lst[1]]  # Atoms2，对象
                self.scatter_plot = scatter_plot_window(self.rotate_obj1, self.rotate_obj2)
                self.scatter_plot.emit_information_signal.connect(self.plot_after_scatter)
        except Exception as e:
            print(e)

    def plot_after_scatter(self, degree, a1_vector, b1_vector, layer_distance, vaccum_distance, Normal_strain_limit,
                           shear_strain_limit, area_error, a2_vector, b2_vector):
        try:
            Hetero_junction_object = Hetero_junction()
            Hetero_junction_object.Rotate_angle = degree/np.pi*180
            Hetero_junction_object.Normal_strain_limit = Normal_strain_limit
            Hetero_junction_object.Shear_strain_limit = shear_strain_limit
            Hetero_junction_object.Area_error_limit = area_error

            cell1 = self.rotate_obj1.get_cell()
            cell1_a, cell1_b = cell1[0][:2], cell1[1][:2]
            A = np.array([cell1_a, cell1_b]).T
            r1a = np.linalg.solve(A, np.array(a1_vector))
            r1b = np.linalg.solve(A, np.array(b1_vector))
            k11 = list(map(lambda x: round(x), r1a))
            k12 = list(map(lambda x: round(x), r1b))

            optimal_match_layer1 = [k11, k12]
            r1a = np.append(r1a, 0)
            r1b = np.append(r1b, 0)
            layer_down = self.deal_with_rotate(cut(self.rotate_obj1, r1a, r1b, [0, 0, 1], origo=[0, 0, 0]))
            # layer_up
            cell2_before_rotate = self.rotate_obj2.get_cell()
            rotate_matrix = np.array(
                [[np.cos(degree), np.sin(degree)], [-np.sin(degree), np.cos(degree)]])
            new_cell2_a = cell2_before_rotate[0][:2] @ rotate_matrix
            new_cell2_b = cell2_before_rotate[1][:2] @ rotate_matrix
            A = np.array([new_cell2_a, new_cell2_b]).T
            r2a = np.linalg.solve(A, np.array(a2_vector))
            r2b = np.linalg.solve(A, np.array(b2_vector))  # 求解线性方程组，简单,直角坐标系下----用晶胞坐标系表示
            k21, k22 = list(r2a), list(r2b)
            for i in range(len(k21)):
                k21[i] = round(k21[i])
                k22[i] = round(k22[i])
            optimal_match_layer2 = [k21, k22]
            Hetero_junction_object.Optimal_orientation = [optimal_match_layer1, optimal_match_layer2]

            r2a = np.append(r2a, 0)
            r2b = np.append(r2b, 0)
            layer_up = cut(self.rotate_obj2, r2a, r2b, [0, 0, 1], origo=[0, 0, 0])
            self.layer_up = self.deal_with_rotate(layer_up)
            layer_up_cell = self.layer_up.get_cell_lengths_and_angles()
            # start stacking
            layer_down_cell = layer_down.get_cell_lengths_and_angles()
            layer_up_cell[0], layer_up_cell[1], layer_up_cell[5] = layer_down_cell[0], \
                                                                   layer_down_cell[1], layer_down_cell[5]
            self.layer_up.set_cell(layer_up_cell, scale_atoms=True)
            self.layer_up.translate(np.array([0, 0, layer_distance + layer_down.get_cell()[2][2]]))
            layer_down_cell[2] += self.layer_up.get_cell_lengths_and_angles()[2] + layer_distance + .01
            layer_down.extend(self.layer_up)
            layer_down.set_cell(layer_down_cell)
            self.Atomsobject = deepcopy(layer_down)
            layer_down = self.set_vacuum_layer(vaccum_distance, addproject=False)
            # plot
            self.plot(layer_down, clear=True, globalAtomsobject=False, dictionary=True)
            Hetero_junction_object.Name = self.dirkey
            Hetero_junction_object.Atomsobject = layer_down
            self.add_hetero_junction(Hetero_junction_object)
            try:
                parent = self.project_tree.selectedItems()[0].parent()
                stackchild = QtWidgets.QTreeWidgetItem(parent)
                stackchild.setText(0, "stack")
                stackchild.setText(1, self.dirkey)
                if self.project_tree.selectedItems()[0].text(2) != "":
                    text1 = self.project_tree.selectedItems()[0].text(2)
                else:
                    text1 = self.project_tree.selectedItems()[0].text(1)
                if self.project_tree.selectedItems()[1].text(2) != "":
                    text2 = self.project_tree.selectedItems()[1].text(2)
                else:
                    text2 = self.project_tree.selectedItems()[1].text(1)
                stackchild.setText(2, text1 + '-' + text2)
            except Exception as e:
                print(e)
                stackchild.setText(2, self.project_tree.selectedItems()[0].text(2))
            finally:
                self.scatter_plot.close()
        except Exception as e:
            print(e)


    def middle_screen(self):
        """ Put the model in the middle of the screen."""
        try:
            height = self.view_widget.height()
            width = self.view_widget.width()
            self.view_widget.opts['viewport']  = (-width, -height,
                                                 3 * width, 3 * height)
            self.view_widget.update()
        except Exception as e:
            print(e)

    def translate_choose(self):
        """ To translate Crystal in the screen"""
        try:
            self.view_widget.num_mouse_track = 2
        except Exception as e:
            print(e)

    def down_screen(self):
        """ Move the model to the bottom of the screen."""
        try:
            height = self.view_widget.height()
            viewport = list(self.view_widget.getViewport())
            viewport[1] -= 30
            viewport[3] = 2 * height - viewport[1]
            viewport = tuple(viewport)
            self.view_widget.opts['viewport'] = viewport
            self.view_widget.update()
        except Exception as e:
            print(e)


    def up_screen(self):
        """ Move the model to the top of the screen."""
        try:
            height = self.view_widget.height()
            viewport = list(self.view_widget.getViewport())
            viewport[1] += 30
            viewport[3] = 2 * height - viewport[1]
            viewport = tuple(viewport)
            self.view_widget.opts['viewport'] = viewport
            self.view_widget.update()
        except Exception as e:
            print(e)

    def right_screen(self):
        """ Move the model to the right of the screen."""
        try:
            width = self.view_widget.width()
            viewport = list(self.view_widget.getViewport())
            viewport[0] += 30
            viewport[2] = 2 * width - viewport[0]
            viewport = tuple(viewport)
            self.view_widget.opts['viewport'] = viewport
            self.view_widget.update()
        except Exception as e:
            print(e)

    def left_screen(self):
        """ Move the model to the left of the screen."""
        try:
            width = self.view_widget.width()
            viewport = list(self.view_widget.getViewport())
            viewport[0] -= 30
            viewport[2] = 2 * width - viewport[0]
            viewport = tuple(viewport)
            self.view_widget.opts['viewport'] = viewport
            self.view_widget.update()
        except Exception as e:
            print(e)

    def Stick(self):
        """ Stick mode"""
        try:
            current = self.project_tree.currentItem()
            if current is not None and current.parent() is not None and self.Atomsobject is not None:
                if len(self.project_tree.selectedItems()) == 1:  # 用户只选择一个的情况
                    self.plot_num = 2
                    self.plot(self.Atomsobject)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)



    def ball_and_stick(self):
        """ Set the atomic model to the ball-stick model."""
        try:
            current = self.project_tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.project_tree.selectedItems()) == 1 and self.Atomsobject is not None:  # 用户只选择一个的情况
                    self.plot_num = 1
                    self.plot(self.Atomsobject)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)
            print(self.view_widget.getViewport())

    def ballfillingtype(self):
        """ Set the atomic model to the ball-filling model."""
        try:
            current = self.project_tree.currentItem()
            if current is not None and current.parent() is not None and self.Atomsobject is not None:
                if len(self.project_tree.selectedItems()) == 1:  # 用户只选择一个的情况
                    self.plot_num = 0
                    self.plot(self.Atomsobject)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def viewtext(self):
        """ To (not) view the text widget."""
        try:
            if self.tab.isHidden():
                self.tab.setVisible(True)
            else:
                self.tab.setVisible(False)
        except Exception as e:
            print(e)

    def viewproject(self):
        """ To (not) view the project widget."""
        try:
            if self.tab_tree.isHidden():
                self.tab_tree.setVisible(True)
            else:
                self.tab_tree.setVisible(False)
        except Exception as e:
            print(e)

    def viewdatabase(self):
        """ To (not) view the database widget."""
        try:
            if self.vertical_widget.isHidden():
                self.vertical_widget.setVisible(True)
            else:
                self.vertical_widget.setVisible(False)
        except Exception as e:
            print(e)

    def btnstate(self, btn):
        try:
            self.choosebutton = btn
            if self.choosebutton.text() == "Crystal":
                self.databaseshow(1)
            elif self.choosebutton.text() == 'Hetero-junction':
                self.databaseshow(2)
            elif self.choosebutton.text() == 'Device':
                self.databaseshow(3)
        except Exception as e:
            print(e)

    def databaseshow(self, number):
        """ To show the database (Crystal, Hetero-junction, Device)."""
        try:
            if number == 1:        # Crystal
                self.datanumber = 1
                self.datatable.setColumnCount(9)
                self.datatable.setHorizontalHeaderLabels(['ID', 'Formula', 'a(Å)', 'b(Å)', 'c(Å)',
                                                          'α(°)', 'β(°)', 'γ(°)', 'Volume'])
                self.database_class.exhibit_all_info(fenkuai=True)
                data = self.database_class.info_data_base
                if not self.true_id_lst:    # search里没有信息,true_id_lst为空
                    data_info = [k[0] for k in data]
                    self.datatable.setRowCount(len(data_info))
                    for index, data_piece in enumerate(data_info):
                        for k in range(1, 10):
                            newitem = QtWidgets.QTableWidgetItem(str(list(data_piece)[k]))
                            newitem.setTextAlignment(5 | 5)
                            self.datatable.setItem(index, k - 1, newitem)
                else:
                    self.database_class.exhibit_all_info(fenkuai=False)
                    data_fenkuaiFalse = self.database_class.info_data_base
                    true_id_lst = [data_piece[0][0] for data_piece in data]

                    # 交集
                    self.Reflect_id_lst = list(set(true_id_lst).intersection(set(self.true_id_lst)))
                    self.datatable.setRowCount(len(self.Reflect_id_lst))
                    id_lst = [info[0] for info in data_fenkuaiFalse]
                    print("id_lst = ", id_lst)
                    for index, reflect_id in enumerate(self.Reflect_id_lst):
                        for k in range(1, 10):
                            newitem = QtWidgets.QTableWidgetItem(str(data_fenkuaiFalse[id_lst.index(reflect_id)][k]))
                            newitem.setTextAlignment(5 | 5)
                            self.datatable.setItem(index, k - 1, newitem)
                self.datatable.setShowGrid(True)
            elif number == 2:      # Hetero-junction
                self.datanumber = 2
                try:
                    self.datatable.setColumnCount(5)
                    self.datatable.setHorizontalHeaderLabels(['Hetero-junction', 'Optimal Match', '#Layer',
                                                              'Binding energy(eV)', 'Schottky_barrier(eV)'])
                    self.database_class.exhibit_all_info(fenkuai=False)
                    data_info = self.database_class.info_data_base
                    self.Reflect_id_lst = []
                    if not self.true_id_lst:
                        self.datatable.setRowCount(len(data_info))  # 添加信息
                        for index, data_piece in enumerate(data_info):
                            for k in range(10, 15):
                                newitem = QtWidgets.QTableWidgetItem(str(list(data_piece)[k]))
                                newitem.setTextAlignment(5 | 5)
                                self.datatable.setItem(index, k - 10, newitem)
                    else:
                        self.datatable.setRowCount(len(self.true_id_lst))
                        id_lst = [_[0] for _ in data_info]

                        for index, true_id in enumerate(self.true_id_lst):
                            for k in range(10, 15):
                                newitem = QtWidgets.QTableWidgetItem(str(data_info[id_lst.index(true_id)][k]))
                                newitem.setTextAlignment(5 | 5)
                                self.datatable.setItem(index, k - 10, newitem)
                    self.datatable.setShowGrid(True)
                except Exception as e:
                    print(e)
            elif number == 3:
                self.Reflect_id_lst = []
                self.datanumber = 3
                try:
                    self.datatable.setColumnCount(5)
                    self.datatable.setHorizontalHeaderLabels(['Device', 'Optimal Match', '#Layer',
                                                              'Schottky barrier (eV)', 'I-V curve (Theory)'])
                    self.database_class.exhibit_all_info(fenkuai=False)
                    data_info = self.database_class.info_data_base
                    if not self.true_id_lst:
                        self.datatable.setRowCount(len(data_info))  # 添加信息
                        for index, data_piece in enumerate(data_info):
                            for k in range(15, 20):
                                newitem = QtWidgets.QTableWidgetItem(str(list(data_piece)[k]))
                                newitem.setTextAlignment(5 | 5)
                                self.datatable.setItem(index, k - 15, newitem)
                    else:
                        id_lst = [data_piece[0] for data_piece in data_info]
                        self.datatable.setRowCount(len(self.true_id_lst))
                        for index, true_id in enumerate(self.true_id_lst):
                            for k in range(15, 20):
                                newitem = QtWidgets.QTableWidgetItem(str(data_info[id_lst.index(true_id)][k]))
                                newitem.setTextAlignment(5 | 5)
                                self.datatable.setItem(index, k - 15, newitem)
                    self.datatable.setShowGrid(True)
                except Exception as e:
                    print(e)
            self.show_image_table.setHidden(False)
        except Exception as e:
            print(e)

    def searchdatabase(self):
        """ database searcher"""
        try:
            text = self.searchdatabase_LineEdit.text()
            if len(text):
                self.true_id_lst = []
                self.row_lst = []
                textlength = len(text)
                self.database_class.exhibit_all_info(fenkuai=False)
                data_info = self.database_class.info_data_base
                for row, data in enumerate(data_info):
                    num = 0
                    itemtext = data[2]
                    for letter in text:
                        if letter not in itemtext:
                            break
                        else:
                            num += 1
                    if num == textlength:
                        self.true_id_lst.append(data[0])
                        self.row_lst.append(row)
                self.databaseshow(self.datanumber)
            else:
                self.true_id_lst = []
                self.row_lst = []
                self.databaseshow(self.datanumber)
        except Exception as e:
            print(e)

    def search(self):
        """ search function in file Tab."""
        global searchlist
        try:
            for item in searchlist:
                item.setHidden(False)
        except Exception as e:
            print(e)
        try:
            text = self.searchLineEdit.text()
            textlength = len(text)
            if textlength != 0:
                it = QtWidgets.QTreeWidgetItemIterator(self.project_tree)
                searchlist = []
                while it.value():
                    if it.value().text(0) != 'bulk' and it.value().text(0) != 'layer' and it.value().text(0) != 'stack':
                        num = 0
                        for letter in text:
                            if letter not in it.value().text(0):
                                break
                            else:
                                num += 1
                        if num != textlength:
                            searchlist.append(it.value())
                    it += 1
                for item in searchlist:
                    item.setHidden(True)
        except Exception as e:
            print(e)

    def filterfrom_to(self):
        try:
            if self.filterLineEdit1.text() == "" and self.filterLineEdit2.text() == "":
                QtWidgets.QMessageBox.warning(self, 'error', "Please input the number of Atoms")
            elif self.filterLineEdit1.text() == "":
                fromnum = 0
                try:
                    tonum = int(self.filterLineEdit2.text())
                except:
                    QtWidgets.QMessageBox.warning(self, 'error', "Please input an integer")
            elif self.filterLineEdit2.text() == "":
                tonum = 100000000000
                try:
                    fromnum = int(self.filterLineEdit1.text())
                except:
                    QtWidgets.QMessageBox.warning(self, 'error', "Please input an integer")
            else:
                try:
                    fromnum = int(self.filterLineEdit1.text())
                    tonum = int(self.filterLineEdit2.text())
                except Exception as e:
                    print(e)
                    QtWidgets.QMessageBox.warning(self, 'error', "Please input an integer")
            it = QtWidgets.QTreeWidgetItemIterator(self.project_tree)
            num = 1
            l = []
            while it.value():
                if it.value().text(1) != "":
                    crys1 = self.dic_Formula_Atoms[it.value().text(1)]
                    atomnum = len(crys1.get_positions())
                    if atomnum > tonum or atomnum < fromnum:
                        l.append(it.value())
                num += 1
                it += 1
            try:
                self.view_widget.clear()
                self.clear_text()

                for child in l:
                    # self.clear_painter()

                    self.dic_Formula_Atoms.pop(child.text(1))  # 对象字典去除该对象
                    if child.parent().childCount() > 1:
                        child.parent().removeChild(child)
                    else:
                        self.project_tree.takeTopLevelItem(self.project_tree.indexOfTopLevelItem(child.parent()))
            except Exception as e:
                print(e)
        except Exception as e:
            print(e)

    def normalchoose(self):
        """ To change self.view_widget.num1 to 0, and thus change the PaintEvent of self.view_widget."""
        try:
            global x_rectangle_dic
            try:  # 消除之前画的图
                for obj in x_rectangle_dic.values():
                    self.view_widget.removeItem(obj)
            except Exception as e:
                print(e)
            self.view_widget.num_mouse_track = 0  # 改变原来的view_widget 函数,从而重写glviewwidget
        except Exception as e:
            print(e)

    def dragrectangle(self):
        try:
            global x
            try:
                self.view_widget.removeItem(x)
            except:
                pass
            current = self.project_tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.project_tree.selectedItems()) == 1 and self.Atomsobject is not None:  # 用户只选择一个的情况
                    self.view_widget.update()
                    self.view_widget.signal_release.connect(self.release)
                    self.view_widget.num_mouse_track = 1  # 改变原来的view_widget 函数,从而重写glviewwidget
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def release(self, initpos, currentpos):
        try:
            global x_rectangle_dic
            try:  # 消除之前画的图
                for obj in x_rectangle_dic.values():
                    self.view_widget.removeItem(obj)
            except Exception as e:
                print(e)
            x_rectangle_dic = {}  # 字典
            initpos = list(initpos)
            currentpos = list(currentpos)
            # initpos[1] = self.view_widget.height() - initpos[1]
            # currentpos[1] = self.view_widget.height() - currentpos[1]
            x_min = min(initpos[0], currentpos[0])
            x_max = max(initpos[0], currentpos[0])
            y_min = min(initpos[1], currentpos[1])
            y_max = max(initpos[1], currentpos[1])
            for index, pos in enumerate(self.global_positions_lst):  # 获得位置在屏幕上的投影的屏幕坐标
                new_pos = pos + self.translate_pos
                screentuple = gluProject(new_pos[0], new_pos[1], new_pos[2])
                if x_min < screentuple[0] < x_max and y_min < self.view_widget.height() - screentuple[1] < y_max:
                    try:
                        Atomic_number = self.Atomsobject.get_atomic_numbers()[index]

                        md = pg.opengl.MeshData.sphere(rows=10, cols=20)
                        x = pg.opengl.GLMeshItem(meshdata=md, smooth=True,
                                          color=self.Atoms_color[self.Atoms_num_lst[index]], shader='shaded',
                                          drawEdges=True)
                        translate_pos = self.global_positions_lst[index] + self.translate_pos
                        x.translate(translate_pos[0], translate_pos[1], translate_pos[2])
                        if self.plot_num == 0:
                            x.scale(self.Atom_radii[Atomic_number] + .01,
                                    self.Atom_radii[Atomic_number] + .01,
                                    self.Atom_radii[Atomic_number] + .01)
                        elif self.plot_num == 1:
                            x.scale(self.Atom_radii[Atomic_number]/2 + .005,
                                    self.Atom_radii[Atomic_number]/2 + .005,
                                    self.Atom_radii[Atomic_number]/2 + .005)
                        elif self.plot_num == 2:
                            x.scale(.1, .1, .1)
                        self.view_widget.addItem(x)
                        x_rectangle_dic[index] = x
                    except Exception as e:
                        print(e)
        except Exception as e:
            print(e)

    def setlayer(self):
        """ To set layer chosen by drag rectangle mode."""
        global x_rectangle_dic
        try:
            if len(x_rectangle_dic) == 0:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please use drag rectangle to select atoms')
            else:
                qm = QtWidgets.QMessageBox.question(self, 'Question',
                                                    'Would yo like to set {} atoms you choose as a layer?'.format(
                                                        str(len(x_rectangle_dic))))
                if qm == QtWidgets.QMessageBox.Yes:
                    max_height = 0
                    min_height = 10000
                    self.view_widget.num_mouse_track = 0
                    choose_index_lst = []
                    for index in x_rectangle_dic.keys():
                        print("index = ", index)
                        choose_index_lst.append(index)
                        if self.global_positions_lst[index][2] > max_height:
                            max_height = self.global_positions_lst[index][2]
                        if self.global_positions_lst[index][2] < min_height:       # 取最低的原子作为origo
                            min_height = self.global_positions_lst[index][2]
                    height = max_height - min_height
                    choose_Atomsobject = deepcopy(self.Atomsobject)  # 用户选的原子
                    choose_index_lst.sort()
                    for index in range(len(self.global_positions_lst) - 1, -1, -1):  # 所有的索引
                        if index not in choose_index_lst:
                            choose_Atomsobject.pop(index)
                    cell_par = choose_Atomsobject.get_cell_lengths_and_angles()
                    if cell_par[3] == 90 and cell_par[4] == 90:
                        cell_par[2] = height + .01
                        choose_Atomsobject.set_cell(cell_par)
                        choose_Atomsobject.translate(np.array([0, 0, -min_height]))

                    else:
                        cell_par[3], cell_par[4], cell_par[2] = 90, 90, height + .01
                        choose_Atomsobject.set_cell(cell_par)
                        choose_Atomsobject.translate(np.array([0, 0, -min_height]))
                        pos_lst = choose_Atomsobject.get_positions()
                        zhijiao_system = choose_Atomsobject.get_cell()
                        A = zhijiao_system.T
                        new_pos_lst = []
                        for pos in pos_lst:
                            b = pos.T
                            r = np.linalg.solve(A, b)  # 求解线性方程组，直角坐标系下----用晶胞坐标系表示
                            while r[0] < 0:
                                pos += zhijiao_system[0]
                                k = pos.T
                                r = np.linalg.solve(A, k)
                            while r[1] < 0:
                                pos += zhijiao_system[1]
                                k = pos.T
                                r = np.linalg.solve(A, k)
                            new_pos_lst.append(pos)
                        choose_Atomsobject.set_positions(new_pos_lst)
                    self.plot(choose_Atomsobject, clear=True, dictionary=True, globalAtomsobject=True, object=True)
                    text3column = self.judgeconductivity(choose_Atomsobject)
                    childx = QtWidgets.QTreeWidgetItem(self.project_tree.currentItem().parent())
                    childx.setText(1, self.dirkey)
                    childx.setText(0, 'layer')
                    childx.setText(3, text3column)
                    x_rectangle_dic = {}
                    self.Atomsobject = None
        except Exception as e:
            print(e)

    def movelayer(self):
        """ To move layer."""
        try:
            self.up_down = True
            self.supercell1_cell = self.Atomsobject.get_cell_lengths_and_angles()
            self.view_widget.num_mouse_track = 0
            global x_rectangle_dic
            if len(x_rectangle_dic) == 0:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please use drag rectangle to select atoms')
            else:
                atom_num = len(self.global_positions_lst)
                client_lst = list(x_rectangle_dic.keys())
                client_lst.sort()
                orijin1 = deepcopy(self.Atomsobject)  # 用户选的原子
                orijin2 = deepcopy(self.Atomsobject)
                for index in (_ for _ in range(atom_num - 1, -1, -1)):  # 所有的索引
                    if index not in client_lst:
                        orijin1.pop(index)
                    elif index in client_lst:
                        orijin2.pop(index)
                self.crys1_supercell = deepcopy(orijin2)
                self.crys2_supercell = deepcopy(orijin1)  # 用户选的在上层
                self.scatter_plot_move_layer_window = scatter_plot_move_layer_window(
                    self.crys1_supercell, self.crys2_supercell, self.Atoms_color, self.Atom_radii)
                self.scatter_plot_move_layer_window.emit_information_signal.connect(self.after_move_layer_window)
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error', 'Please use drag rectangle to select atoms')

    def after_move_layer_window(self, x, y, z):
        crys2_pos = list(self.crys2_supercell.get_positions())
        new_pos = list(map(lambda x:x + np.array([x, y, z]), crys2_pos))
        try:  # 阔胞后切割
            self.crys2_supercell.set_positions(new_pos)
            self.crys1_supercell.extend(self.crys2_supercell)  # 把crys2加在crys1上
            self.crys1_supercell.set_cell(self.supercell1_cell)
            self.crys1_supercell = cut(self.crys1_supercell, a=(1, 0, 0), b=(0, 1, 0), c=(0, 0, 1), origo=(0, 0, 0))
        except Exception as e:
            print(e)
        # 画出切割后的超胞
        try:
            self.plot(self.crys1_supercell, clear=True, globalAtomsobject=True, dictionary=True)
            if len(self.project_tree.selectedItems()) == 2:
                self.Hetero_junction_object.Atomsobject = deepcopy(self.crys1_supercell)
                self.Hetero_junction_object.Name = self.dirkey
                self.add_hetero_junction(self.Hetero_junction_object)
                parent = self.project_tree.selectedItems()[0].parent()
                stackchild = QtWidgets.QTreeWidgetItem(parent)
                stackchild.setText(0, "stack")
                stackchild.setText(1, self.dirkey)
                if self.project_tree.selectedItems()[0].text(2) != "":
                    text1 = self.project_tree.selectedItems()[0].text(2)
                else:
                    text1 = self.project_tree.selectedItems()[0].text(1)
                if self.project_tree.selectedItems()[1].text(2) != "":
                    text2 = self.project_tree.selectedItems()[1].text(2)
                else:
                    text2 = self.project_tree.selectedItems()[1].text(1)
                stackchild.setText(2, text1 + '-' + text2)
            elif len(self.project_tree.selectedItems()) == 1:
                parent = self.project_tree.selectedItems()[0].parent()
                stackchild = QtWidgets.QTreeWidgetItem(parent)
                stackchild.setText(0, self.project_tree.selectedItems()[0].text(0))
                stackchild.setText(1, self.dirkey)
                stackchild.setText(2, self.project_tree.selectedItems()[0].text(2))
                stackchild.setText(3, self.project_tree.selectedItems()[0].text(3))
        except Exception as e:
            print(e)
            stackchild.setText(2, self.project_tree.selectedItems()[0].text(2))

    def traversalcut(self):
        try:
            it = QtWidgets.QTreeWidgetItemIterator(self.project_tree)
            l = []
            s = []
            while it.value():
                v = it.value().text(0)
                if v == "bulk":
                    l.append(self.dic_Formula_Atoms[it.value().text(1)])
                    s.append(it.value())
                it += 1
            for index, object in enumerate(l):
                try:
                    self.Atomsobject = object
                    self.cut_exf_layer(s[index].parent(), traversal=True)
                except Exception as e:
                    print(e)
        except Exception as e:
            print(e)

    def traversalrotate(self):
        try:
            self.traversalrotate_window.close()
        except:
            pass
        try:
            self.traversalrotate_window = traversal_rotate_window()
            self.traversalrotate_window.signaltraversalstack.connect(self.aftertraversalrotate_window)
        except Exception as e:
            print(e)

    def rotate_all_angles(self, obj1, obj2):
        try:
            rotate_obj1 = deepcopy(obj1)
            rotate_obj2 = deepcopy(obj2)
            cell1 = list(rotate_obj1.get_cell())
            self.vectora1 = np.array(cell1[0][:2])
            self.vectorb1 = np.array(cell1[1][:2])
            self.matrix = np.linalg.inv(np.array([self.vectora1, self.vectorb1]))
            cell2 = list(rotate_obj2.get_cell())
            vectora2 = np.array(cell2[0][:2])
            vectorb2 = np.array(cell2[1][:2])
            self.pos_lst_orijin_up = []
            number = int(self.len_max / min(np.linalg.norm(vectora2), np.linalg.norm(vectorb2)))
            for i in range(- number, number):
                for j in range(-number,  number):
                    pos = vectora2 * i + vectorb2 * j
                    self.pos_lst_orijin_up.append(pos)
            theta = 0
            unit = self.relative_length_limit
            compare_square = 100000000000
            a_down = ...
            b_down = ...
            a_up = ...
            b_up = ...
            theta_true = ...
            cell_2_par = rotate_obj2.get_cell_lengths_and_angles()
            if abs(cell_2_par[0] - cell_2_par[1]) < 10e-6:
                if abs(cell_2_par[5] - 60) < 1e-6 or abs(cell_2_par[5] - 120) < 1e-6:
                    division = 6
                elif abs(cell_2_par[5] - 90) < 1e-6:
                    division = 4
            else:
                division = 1
            while theta < np.pi / division:
                rotate_matrix = np.array(
                    [[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
                self.new_spot_up_lst = list(map(lambda x:x @ rotate_matrix, self.pos_lst_orijin_up))
                a1, b1, a2, b2, min_square = self.find_min_square_orientation()
                if a1 is not None:
                    if min_square < compare_square:
                        compare_square = min_square
                        a_down = a1
                        b_down = b1
                        a_up = a2
                        b_up = b2
                        theta_true = theta
                theta += unit
            if compare_square == 100000000000:
                return "fail"
            else:
                rotate_hetero_junction_object = Hetero_junction()
                rotate_hetero_junction_object.Normal_strain_limit = self.normal_strain_limit
                rotate_hetero_junction_object.Rotate_angle = theta_true / np.pi * 180
                cell1 = rotate_obj1.get_cell()
                cell1_a = cell1[0][:2]
                cell1_b = cell1[1][:2]
                A = np.array([cell1_a, cell1_b]).T
                r1a = np.linalg.solve(A, np.array(a_down))
                r1b = np.linalg.solve(A, np.array(b_down))
                k11 = list(map(lambda x: round(x), r1a))
                k12 = list(map(lambda x: round(x), r1b))
                optimal_match_layer1 = [k11, k12]
                r1a = np.append(r1a, 0)
                r1b = np.append(r1b, 0)
                self.layer_down = self.deal_with_rotate(cut(rotate_obj1, r1a, r1b, [0, 0, 1], origo=[0, 0, 0]))
                # layer_up
                cell2_before_rotate = rotate_obj2.get_cell()
                cell2_a = cell2_before_rotate[0][:2]
                cell2_b = cell2_before_rotate[1][:2]
                rotate_matrix = np.array(
                    [[np.cos(theta_true), np.sin(theta_true)], [-np.sin(theta_true), np.cos(theta_true)]])
                new_cell2_a = cell2_a @ rotate_matrix
                new_cell2_b = cell2_b @ rotate_matrix
                A = np.array([new_cell2_a, new_cell2_b]).T
                r2a = np.linalg.solve(A, np.array(a_up))
                r2b = np.linalg.solve(A, np.array(b_up))  # 求解线性方程组，简单,直角坐标系下----用晶胞坐标系表示
                k21 = list(map(lambda x: round(x), (_ for _ in r2a)))
                k22 = list(map(lambda x: round(x), (_ for _ in r2b)))
                optimal_match_layer2 = [k21, k22]
                rotate_hetero_junction_object.Optimal_orientation = [optimal_match_layer1, optimal_match_layer2]
                r2a = np.append(r2a, 0)
                r2b = np.append(r2b, 0)
                self.layer_up = self.deal_with_rotate(cut(rotate_obj2, r2a, r2b, [0, 0, 1], origo=[0, 0, 0]))
                # start stacking
                layer_up_cell = self.layer_up.get_cell_lengths_and_angles()
                layer_down_cell = self.layer_down.get_cell_lengths_and_angles()
                layer_up_cell[0], layer_up_cell[1], layer_up_cell[5] = \
                    layer_down_cell[0], layer_down_cell[1], layer_down_cell[5]
                self.layer_up.set_cell(layer_up_cell, scale_atoms=True)
                self.layer_up.translate(np.array([0, 0, self.layerlength + self.layer_down.get_cell()[2][2]]))
                layer_down_cell[2] += self.layer_up.get_cell_lengths_and_angles()[2] + self.layerlength
                self.layer_down.extend(self.layer_up)
                self.layer_down.set_cell(layer_down_cell)
                self.Atomsobject = deepcopy(self.layer_down)
                self.layer_down = self.set_vacuum_layer(self.vacuum_distance, addproject=False)
                self.plot(self.layer_down, plot=False, object=False, Hetero_tab=False, clear=True, globalAtomsobject=False, dictionary=True)
                rotate_hetero_junction_object.Atomsobject = deepcopy(self.layer_down)
                rotate_hetero_junction_object.Name = self.dirkey
                self.add_hetero_junction(rotate_hetero_junction_object)
                key = list(self.dic_Formula_Atoms.keys())[
                    list(self.dic_Formula_Atoms.values()).index(self.orijincrys1)]
                it = QtWidgets.QTreeWidgetItemIterator(self.project_tree)
                while it.value():
                    if it.value().text(1) == key:
                        try:
                            childx = QtWidgets.QTreeWidgetItem(it.value().parent())
                            childx.setText(1, self.dirkey)
                            childx.setText(0, 'stack')
                            if self.it_value1.text(2) != "":
                                text1 = self.it_value1.text(2)
                            else:
                                text1 = self.it_value1.text(1)
                            if self.it_value2.text(2) != "":
                                text2 = self.it_value2.text(2)
                            else:
                                text2 = self.it_value2.text(1)
                            childx.setText(2, text1 + '-' + text2)
                        except Exception as e:
                            print(e)
                        break
                    it += 1
                return theta_true
        except Exception as e:
            print(e)

    def aftertraversalrotate_window(self):
        """Two sides don't know."""
        try:
            self.len_max = self.traversalrotate_window.max_length_var
            self.layerlength = self.traversalrotate_window.layerdistance_var
            self.vacuum_distance = self.traversalrotate_window.vacuumdistance_var
            self.normal_strain_limit = self.traversalrotate_window.normal_strain_var / 100

            self.shear_strain_limit = self.traversalrotate_window.shear_strain_var / 100
            self.relative_length_limit = np.sqrt(self.shear_strain_limit ** 2 + self.normal_strain_limit ** 2)

            self.areatolerance = self.traversalrotate_window.area_tolerance_var / 100
            self.minimum_gama = self.traversalrotate_window.minimum_gamma_var
            self.maximum_gama = self.traversalrotate_window.maximum_gamma_var
            self.ab_relative_error = self.traversalrotate_window.ab_relative_error_var / 100
            l = []  # 用来储存用户要求堆垛的晶胞,layer
            s = []  # stack
            objl = []  # Project栏的it.value,layer
            objs = []  # Project栏的it.value,stack
            it = QtWidgets.QTreeWidgetItemIterator(self.project_tree)
            while it.value():
                v = it.value().text(0)
                if v == "layer":
                    l.append(self.dic_Formula_Atoms[it.value().text(1)])
                    objl.append(it.value())

                it += 1
            zongnum = len(l) * (len(l) - 1) / 2
            failnum = 0
            num = 0
            failmessage = ""
            successmessage = ""
            for i in range(len(l)):
                for j in range(i + 1, len(l)):
                    try:
                        self.it_value1 = objl[i]  # l[i],l[j]的原子对象对应的it.value()
                        self.it_value2 = objl[j]
                        self.obj1 = deepcopy(l[i])
                        self.orijincrys1 = deepcopy(self.obj1)  # Atoms1，对象
                        self.obj2 = deepcopy(l[j])
                        fail = self.rotate_all_angles(self.obj1, self.obj2)
                        if fail != "fail":
                            num += 1
                            successmessage = successmessage + l[i].get_chemical_formula(mode='hill') + "-" + l[
                                j].get_chemical_formula(mode='hill') + "\t" + "Rotate {}°".format(fail / np.pi * 180) + "\n"
                        else:
                            failnum += 1
                            failmessage = failmessage + l[i].get_chemical_formula(mode='hill') \
                                          + "-" + l[j].get_chemical_formula(mode='hill') + "\n"
                        self.traversalrotate_window.signal_start_pbar.emit((num + failnum) / zongnum * 100)
                    except Exception as e:
                        failnum += 1
                        failmessage = failmessage + l[i].get_chemical_formula(mode='hill') + "-" + l[
                            j].get_chemical_formula(mode='hill') + "\n"
                        print(e)
            message1 = "There are {} stack situations in total".format(int(zongnum)) \
                       + "\n" + "{} succeeded".format(num) + "\n" + "{} failed".format(failnum)
            message2 = "\n" + "-" * 20 + "\n" + "failed situations:" + "\n" + failmessage
            message3 = "\n" + "-" * 20 + "\n" + "succeeded situations:" + "\n" + successmessage
            message = message1 + message2 + message3
            self.traversalrotate_window.Text.setText(message)
        except Exception as e:
            print(e)


    def find_min_square_orientation(self):
        try:
            peidui_lst = []
            for pos2 in self.new_spot_up_lst:
                if pos2[1] > -self.relative_length_limit * self.len_max:
                    r = pos2 @ self.matrix
                    r1 = np.floor(r[0])
                    r2 = np.floor(r[1])
                    pos1 = self.vectora1 * r1 + self.vectorb1 * r2
                    if np.linalg.norm(pos1 - pos2) < self.relative_length_limit * np.linalg.norm(pos1):
                        if abs(pos1 @ pos1 - pos1 @ pos2) < self.normal_strain_limit * (pos1 @ pos1):
                            peidui_lst.append([pos1, pos2])
                    pos1 = self.vectora1 * (r1 + 1) + self.vectorb1 * r2
                    if np.linalg.norm(pos1 - pos2) < self.relative_length_limit * np.linalg.norm(pos1):
                        if abs(pos1 @ pos1 - pos1 @ pos2) < self.normal_strain_limit * (pos1 @ pos1):
                            peidui_lst.append([pos1, pos2])
                    pos1 = self.vectora1 * (r1 + 1) + self.vectorb1 * (r2 + 1)
                    if np.linalg.norm(pos1 - pos2) < self.relative_length_limit * np.linalg.norm(pos1):
                        if abs(pos1 @ pos1 - pos1 @ pos2) < self.normal_strain_limit * (pos1 @ pos1):
                            peidui_lst.append([pos1, pos2])
                    pos1 = self.vectora1 * r1 + self.vectorb1 * (r2 + 1)
                    if np.linalg.norm(pos1 - pos2) < self.relative_length_limit * np.linalg.norm(pos1):
                        if abs(pos1 @ pos1 - pos1 @ pos2) < self.normal_strain_limit * (pos1 @ pos1):
                            peidui_lst.append([pos1, pos2])
            if len(peidui_lst) < 2:
                return None, None, None, None, None
            else:
                min_square = 100000000000
                a1 = ...
                b1 = ...
                for i in range(len(peidui_lst)):
                    for j in range(i + 1, len(peidui_lst)):
                         if abs(np.linalg.norm(peidui_lst[i][0]) - np.linalg.norm(peidui_lst[j][0])) / \
                                np.linalg.norm(peidui_lst[i][0]) < self.ab_relative_error:
                             n_vector1 = np.cross(peidui_lst[i][0], peidui_lst[j][0])
                             n_vector2 = np.cross(peidui_lst[i][1], peidui_lst[j][1])
                             square1 = abs(n_vector1)
                             square2 = abs(n_vector2)
                             sin_gama = square1 / np.linalg.norm(peidui_lst[i][0]) / np.linalg.norm(peidui_lst[j][0])
                             if abs(square1 - square2) < self.areatolerance * square1 and \
                                        self.maximum_gama > sin_gama > self.minimum_gama:
                                if square1 < min_square:
                                    shear_strain = abs(np.arccos((peidui_lst[i][0] @ peidui_lst[j][0]) /
                                                                 np.linalg.norm(peidui_lst[i][0]) / np.linalg.norm(
                                        peidui_lst[j][0])) -
                                                       np.arccos((peidui_lst[i][1] @ peidui_lst[j][1]) /
                                                                 np.linalg.norm(peidui_lst[i][1]) / np.linalg.norm(
                                                           peidui_lst[j][1])))
                                    if shear_strain < self.shear_strain_limit:
                                        if n_vector1 > 0:     # 保证右手系
                                            a1 = peidui_lst[i][0]
                                            a2 = peidui_lst[i][1]
                                            b1 = peidui_lst[j][0]
                                            b2 = peidui_lst[j][1]

                                        else:
                                            b1 = peidui_lst[i][0]
                                            b2 = peidui_lst[i][1]
                                            a1 = peidui_lst[j][0]
                                            a2 = peidui_lst[j][1]
                                        min_square = square1
                if min_square == 100000000000:
                    return None, None, None, None, None
                else:
                    return a1, b1, a2, b2, min_square
        except Exception as e:
            print(e)

    def traversal_devices(self):
        try:
            self.create_devices_window.close()
        except:
            pass
        try:
            self.create_devices_window = create_2D_devices_window()
            self.create_devices_window.signal_emit_information.connect(self.after_createdevices_window)
        except Exception as e:
            print(e)

    def after_createdevices_window(self, vacuum_dis, layer_dis_lst, conductivity_lst, max_singama, min_singama,
                                   ab_relative_error, Normal_strain_limit, Shear_strain_limit, area_error, max_len):
        try:
            con = []
            semi = []
            insul = []
            it = QtWidgets.QTreeWidgetItemIterator(self.project_tree)
            while it.value():
                if it.value().text(3) == 'Conductor':  # 用project栏heterostructure一栏判断
                    con.append(it.value())
                elif it.value().text(3) == 'Semiconductor':
                    semi.append(it.value())
                elif it.value().text(3) == 'Insulator':
                    insul.append(it.value())
                it += 1
            sum_lst = []
            sum_lst.append(con)
            sum_lst.append(semi)
            sum_lst.append(insul)
            shunxu_lst = []
            for conductivity in conductivity_lst:
                for conduct in sum_lst:
                    if conduct[0].text(3) == conductivity:
                        shunxu_lst.append(conduct)  # 按用户需求排序
                        break
            if len(shunxu_lst) == len(conductivity_lst):
                Atoms = []
                zongnum = 1
                for shunxu in shunxu_lst:
                    s = []
                    for j in shunxu:
                        s.append(self.dic_Formula_Atoms[j.text(1)])
                    zongnum *= len(s)
                    Atoms.append(s)

                Atom_all_condition = self.make_Atom_lst(Atoms, [])
                self.device_root = QtWidgets.QTreeWidgetItem(self.project_tree)
                fail_num = 0
                success_num = 0
                fail_condition = ""
                success_condition = ""
                for index, obj_lst in enumerate(Atom_all_condition):
                    success_or_not, formula = self.after_multi_layer_window(3, obj_lst, layer_dis_lst, vacuum_dis, area_error,
                                                                            Normal_strain_limit, Shear_strain_limit,
                                                                            max_singama, min_singama,
                                                                            ab_relative_error, max_len,
                                                                            traversal_stack=True)
                    if success_or_not == 'fail':
                        fail_num += 1
                        fail_condition = fail_condition + formula + '\n'
                    else:
                        success_num += 1
                        success_condition = success_condition + formula + '\n'
                    self.create_devices_window.pbar.setValue((index + 1) / zongnum * 100)
                self.device_root.setText(0, "Device")
                self.project_tree.expandAll()

                text_that_set = "There are {} conditions in total".format(zongnum) + '\n' + str(fail_num) + \
                                "fail" + "\t" + str(success_num)  + 'succeed' + "\n" + "-" * 32 + "\n" +\
                                "successful condition:" + "\n" + success_condition + "\n"+ "-"*32 + "\n" + \
                                "fail condition:" + "\n" + fail_condition
                self.create_devices_window.text_widgt.setText(text_that_set)
        except Exception as e:
            print(e)

    def make_Atom_lst(self, Atom_lst, make_lst):
        Atom_layer = Atom_lst.pop(0)
        c = []
        for Atom in Atom_layer:
            if make_lst:
                for already_lst in make_lst:
                    t = deepcopy(already_lst)
                    t.append(Atom)
                    c.append(t)
            else:
                c.append([Atom])
        if len(Atom_lst) == 0:
            return c
        else:
            return self.make_Atom_lst(Atom_lst, c)

    def traversalstack1(self):
        try:
            self.traversalstack_window_1.close()
        except:
            pass
        try:
            self.traversalstack_window_1 = traversal_stack_window_1()
            self.traversalstack_window_1.signaltraversalstack.connect(self.aftertraversalstack_window)
        except Exception as e:
            print(e)

    def traversalstack2(self):
        """ To traversal stack multiple crys. """
        try:
            self.traversalstack_window_2.close()
        except:
            pass
        try:
            self.traversalstack_window_2 = traversal_stack_window_2()
            self.traversalstack_window_2.signaltraversalstack.connect(self.aftertraversalstack_window)
        except Exception as e:
            print(e)


    def aftertraversalstack_window(self, num):
        """ To traversal stack multiple crys. """
        try:
            l = []  # 用来储存用户要求堆垛的晶胞,layer
            objl = []  # Project栏的it.value,layer
            it = QtWidgets.QTreeWidgetItemIterator(self.project_tree)
            while it.value():
                v = it.value().text(0)
                if v == "layer":
                    l.append(self.dic_Formula_Atoms[it.value().text(1)])
                    objl.append(it.value())
                it += 1
            if num == 2:      # two-side know
                self.layerlength = self.traversalstack_window_2.layerdistance_var
                self.vacuum_distance = self.traversalstack_window_2.vacuumdistance_var
                self.normal_strain_limit = self.traversalstack_window_2.Normal_strain_limit_var / 100
                zongnum = len(l) * (len(l) - 1) / 2
                failnum = 0
                num = 0
                failmessage = ""
                successmessage = ""
                for i in range(len(l)):
                    for j in range(i + 1, len(l)):
                        try:
                            self.compare_cell(l[i], l[j], traversalstack=True)
                            self.orijincrys1 = deepcopy(l[i])         # 没有转换之前的
                            self.orijincrys2 = deepcopy(l[j])         # 没有转换之前的
                            self.it_value1 = objl[i]         # l[i],l[j]的原子对象对应的it.value()
                            self.it_value2 = objl[j]
                            if self.crys1 is not None and self.crys2 is not None:
                                self.layer_transform(traversalstack=True)
                                num += 1
                                successmessage = successmessage + l[i].get_chemical_formula(mode='hill') + "-" + l[j].get_chemical_formula(mode='hill') + "\n"
                            else:
                                failnum += 1
                                failmessage = failmessage + l[i].get_chemical_formula(mode='hill')  \
                                              + "-" + l[j].get_chemical_formula(mode='hill')  + "\n"
                            self.traversalstack_window_2.signal_start_pbar.emit((num + failnum) / zongnum * 100)
                        except Exception as e:
                            failnum += 1
                            failmessage = failmessage + l[i].get_chemical_formula(mode='hill') + "-" + l[j].get_chemical_formula(mode='hill') + "\n"
                            print(e)
                message1 = "There are {} stack situations in total".format(int(zongnum)) \
                           + "\n" + "{} succeeded".format(num) + "\n" + "{} failed".format(failnum)
                message2 = "\n" + "-" * 20 + "\n" + "failed situations:" + "\n" + failmessage
                message3 = "\n" + "-" * 20 + "\n" + "succeeded situations:" + "\n" + successmessage
                message = message1 + message2 + message3
                self.traversalstack_window_2.Text.setText(message)
                self.Atomsobject = None
            else:    # one-sid-known
                layerlength = self.traversalstack_window_1.layerdistance_var
                vacuum_distance = self.traversalstack_window_1.vacuumdistance_var
                normal_strain_limit = self.traversalstack_window_1.Normal_strain_limit_var / 100
                shear_strain_limit = self.traversalstack_window_1.Shear_strain_limit_var / 100
                zongnum = len(l) * (len(l) - 1) / 2
                failnum = 0
                num = 0
                failmessage = ""
                successmessage = ""
                for i in range(len(l)):
                    for j in range(i + 1, len(l)):
                        try:
                            self.it_value1 = objl[i]  # l[i],l[j]的原子对象对应的it.value()
                            self.it_value2 = objl[j]
                            self.obj1 = deepcopy(l[i])
                            self.orijincrys1 = deepcopy(self.obj1)  # Atoms1，对象
                            self.obj2 = deepcopy(l[j])
                            fail = self.stack_one(normal_strain_limit, shear_strain_limit, layerlength,
                                                  vacuum_distance, traversal=True)
                            if fail != "fail":
                                num += 1
                                successmessage = successmessage + l[i].get_chemical_formula(mode='hill') + "-" + l[
                                    j].get_chemical_formula(mode='hill') + "\n"
                            else:
                                failnum += 1
                                failmessage = failmessage + l[i].get_chemical_formula(mode='hill') \
                                              + "-" + l[j].get_chemical_formula(mode='hill') + "\n"
                            self.traversalstack_window_1.signal_start_pbar.emit((num + failnum) / zongnum * 100)
                        except Exception as e:
                            failnum += 1
                            failmessage = failmessage + l[i].get_chemical_formula(mode='hill') + "-" + l[
                                j].get_chemical_formula(mode='hill') + "\n"
                            print(e)
                message1 = "There are {} stack situations in total".format(int(zongnum)) \
                           + "\n" + "{} succeeded".format(num) + "\n" + "{} failed".format(failnum)
                message2 = "\n" + "-" * 20 + "\n" + "failed situations:" + "\n" + failmessage
                message3 = "\n" + "-" * 20 + "\n" + "succeeded situations:" + "\n" + successmessage
                message = message1 + message2 + message3
                self.traversalstack_window_1.Text.setText(message)
                self.Atomsobject = None
        except Exception as e:
            print(e)

    def calculatevector(self):
        """   To calculate the vector of two atoms."""
        try:
            global index_lst
            if len(index_lst) == 2:
                pos1 = self.Atomsobject.get_positions()[index_lst[0]]
                pos2 = self.Atomsobject.get_positions()[index_lst[1]]
                vector = pos2 - pos1
                self.objectTextwidget.setText("Vector = " + '\n' + str(vector))
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Need 2 Atoms')
        except:
            QtWidgets.QMessageBox.warning(self, 'error', 'Choose 2 Atoms in Object Box')

    def calculatedistance(self):
        """   To calculate the distance of two atoms."""
        try:
            global index_lst
            if len(index_lst) == 2:
                pos1 = self.Atomsobject.get_positions()[index_lst[0]]
                pos2 = self.Atomsobject.get_positions()[index_lst[1]]
                distance = np.linalg.norm(pos1 - pos2)
                self.objectTextwidget.setText("distance(Å) = " + str(distance))
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Need 2 Atoms')
        except:
            QtWidgets.QMessageBox.warning(self, 'error', 'Choose 2 Atoms in Object Box')

    def calculatedegree(self):
        """   To calculate the degree of three atoms' angle."""
        try:
            global index_lst
            if len(index_lst) == 3:
                pos1 = self.Atomsobject.get_positions()[index_lst[0]]
                pos2 = self.Atomsobject.get_positions()[index_lst[1]]
                pos3 = self.Atomsobject.get_positions()[index_lst[2]]
                vector1 = pos1 - pos2
                vector2 = pos3 - pos2
                dc = vector1 @ vector2
                mod_vector1 = np.linalg.norm(vector1)
                mod_vector2 = np.linalg.norm(vector2)
                degree = np.arccos(dc / (mod_vector1 * mod_vector2)) / np.pi * 180
                self.objectTextwidget.setText("degree = " + str(degree) + '°')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Need 3 Atoms')
        except:
            QtWidgets.QMessageBox.warning(self, 'error', 'Choose 3 Atoms in Object Box')

    def handleItemClick_on_object(self):
        """  To handle Atom clicked in the objectbox. """
        current = self.project_tree.currentItem()
        try:
            if current is not None and current.parent() is not None:
                if len(self.project_tree.selectedItems()) == 1 and self.Atomsobject is not None:  # 用户只选择一个的情况
                    try:
                        global x
                        self.view_widget.removeItem(x)  # 移走双击选中的原子
                        x = None
                        self.remove_atom_index_num = None
                    except Exception as e:
                        print(e)
                    global index_lst
                    index_lst = [item.row() for i, item in enumerate(self.object_tablewidget.selectedItems()) if i % 3 == 1]
                    self.tab.setCurrentWidget(self.Atom_tab)
                    try:
                        global obj_lst_table
                        for Atom in obj_lst_table:
                            self.view_widget.removeItem(Atom)
                    except Exception as e:
                        print(e)
                    try:
                        obj_lst_table = []
                        text = ''
                        for index in index_lst:
                            Atomic_number = self.Atomsobject.get_atomic_numbers()[index]
                            md = pg.opengl.MeshData.sphere(rows=10, cols=20)
                            Atom = pg.opengl.GLMeshItem(meshdata=md, smooth=True,
                                              color=self.Atoms_color[self.Atoms_num_lst[index]],
                                              shader='shaded', drawEdges=True)
                            new_transpos = self.global_positions_lst[index] + self.translate_pos
                            Atom.translate(new_transpos[0], new_transpos[1], new_transpos[2])
                            if self.plot_num == 0:
                                Atom.scale(self.Atom_radii[Atomic_number] + .01,
                                        self.Atom_radii[Atomic_number] + .01,
                                        self.Atom_radii[Atomic_number] + .01)
                            elif self.plot_num == 1:
                                Atom.scale(self.Atom_radii[Atomic_number]/2 + .005,
                                        self.Atom_radii[Atomic_number]/2 + .005,
                                        self.Atom_radii[Atomic_number]/2 + .005)
                            elif self.plot_num == 2:
                                Atom.scale(0.1, 0.1, 0.1)

                            self.view_widget.addItem(Atom)
                            j = crys_data.Element_formula[Atomic_number] + "(index:" + str(index) + ")" + "\n" + "Atomic number:" + str(
                                Atomic_number) + "\n" + "position:" + str(self.global_positions_lst[index]) + '\n' + "-" * 32 + '\n'
                            text += j
                            obj_lst_table.append(Atom)
                        self.Atom_textwidget.setText(text)
                    except Exception as e:
                        print(e)

                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def show_None(self):
        """ To delete atom text."""

        try:
            self.show_atom_element_judge = False
            self.show_atom_index_judge = False
            self.plot(self.Atomsobject)

        except Exception as e:
            print(e)

    def show_atom_index(self):
        """ To plot atom index."""
        try:
            self.show_atom_index_judge = True
            self.show_atom_element_judge = False

            self.plot(self.Atomsobject)
        except Exception as e:
            print(e)

    def show_atom_element(self):
        """ To show atom element."""
        try:
            self.show_atom_index_judge = False
            self.show_atom_element_judge = True
            self.plot(self.Atomsobject)
        except Exception as e:
            print(e)

    def viewgrid(self):
        """ Whether plot grid or not."""
        try:
            self.grid_itemnum += 1
            self.plot(self.Atomsobject, dictionary=False)
        except Exception as e:
            print(e)

    def viewCoordinatcell_system(self):  # 晶胞坐标
        """ Whether plot cell coordinate system or not."""
        try:
            self.coordinatesystem_view_num = 1
            self.plot(self.Atomsobject, dictionary=False)
        except Exception as e:
            print(e)


    def No_coordinate_system(self):
        """ No coordiante system"""
        try:
            self.coordinatesystem_view_num = 0
            self.plot(self.Atomsobject, dictionary=False)
        except Exception as e:
            print(e)

    def threeD_coordinate(self):
        """ ThreeD coordinate """
        try:
            self.coordinatesystem_fixed_num = 0
            self.plot(self.Atomsobject)
        except Exception as e:
            print(e)

    def fixed_coordinate(self):
        """ fixed coordinate """
        try:
            self.coordinatesystem_fixed_num = 1
            self.plot(self.Atomsobject)
        except Exception as e:
            print(e)

    def viewCoordinate_system(self):  # 笛卡尔坐标（直角坐标）
        """ Whether plot cartestian coordinate system or not."""
        try:
            self.coordinatesystem_view_num = 2
            self.plot(self.Atomsobject, dictionary=False)
        except Exception as e:
            print(e)

    def view_cell(self):
        """ Whether see cell or not."""
        try:
            self.cell_view_num += 1
            if self.cell_view_num % 2 == 1:
                self.actionViewCell.setIconVisibleInMenu(False)
            else:
                self.actionViewCell.setIconVisibleInMenu(True)
            self.plot(self.Atomsobject, dictionary=False)
        except:
            ...

    def viewfroma(self):
        """ To view from a axis."""
        try:
            cell_par = self.Atomsobject.get_cell()
            if cell_par[0][0] != 0:  # x != 0
                tantheta = cell_par[0][1] / cell_par[0][0]
                if cell_par[0][1] >= 0:  # 一二象限情况
                    theta = np.arctan(tantheta) / np.pi * 180
                else:  # 三四象限情况
                    theta = np.arctan(tantheta) / np.pi * 180 + 180
            elif cell_par[0][1] == 0:  # x = 0,y = 0的情况
                theta = -90
            elif cell_par[0][1] > 0:  # x = 0, y > 0的情况
                theta = 90
            elif cell_par[0][1] < 0:  # x = 0, y < 0的情况
                theta = -90  # 方位角

            if cell_par[0][0] == 0 and cell_par[0][1] == 0:
                elevation = 90
            else:
                tanelevation = cell_par[0][2] / np.sqrt(cell_par[0][0] ** 2 + cell_par[0][1] ** 2)
                if cell_par[0][2] >= 0:
                    elevation = np.arctan(tanelevation) / np.pi * 180
                else:
                    elevation = np.arctan(tanelevation) / np.pi * 180 + 180
            self.view_widget.setCameraPosition(elevation=elevation, azimuth=theta)
        except:
            ...

    def viewfromb(self):
        """ To view from b axis."""
        try:
            cell_par = self.Atomsobject.get_cell()
            if cell_par[1][0] != 0:  # x != 0
                tantheta = cell_par[1][1] / cell_par[1][0]
                if cell_par[1][1] >= 0:  # 一二象限情况
                    theta = np.arctan(tantheta) / np.pi * 180
                else:  # 三四象限情况
                    theta = np.arctan(tantheta) / np.pi * 180 + 180
                if tantheta < 0 and abs(theta - 0) < .0001:
                    theta = 180
            elif cell_par[1][1] == 0:  # y = 0的情况
                theta = -90
            elif cell_par[1][1] > 0:  # y > 0的情况
                theta = 90
            elif cell_par[1][1] < 0:  # y < 0的情况
                theta = -90  # 方位角
            if cell_par[1][0] == 0 and cell_par[1][1] == 0:
                elevation = 90
            else:
                tanelevation = cell_par[1][2] / np.sqrt(cell_par[1][0] ** 2 + cell_par[1][1] ** 2)
                if cell_par[1][2] >= 0:
                    elevation = np.arctan(tanelevation) / np.pi * 180
                else:
                    elevation = np.arctan(tanelevation) / np.pi * 180 + 180
            self.view_widget.setCameraPosition(elevation=elevation, azimuth=theta)
        except:
            ...

    def viewfromc(self):
        """ To view from c axis."""
        try:
            cell_par = self.Atomsobject.get_cell()
            if cell_par[2][0] != 0:  # x = 0
                tantheta = cell_par[2][1] / cell_par[2][0]
                if cell_par[2][1] >= 0:  # 一二象限情况
                    theta = np.arctan(tantheta) / np.pi * 180
                else:  # 三四象限情况
                    theta = np.arctan(tantheta) / np.pi * 180 + 180
                if tantheta < 0 and abs(theta - 0) < .0001:
                    theta = 180
            elif cell_par[2][1] == 0:  # y = 0的情况
                theta = -90
            elif cell_par[2][1] > 0:  # y > 0的情况
                theta = 90
            elif cell_par[2][1] < 0:  # y < 0的情况
                theta = -90  # 方位角

            if cell_par[2][0] == 0 and cell_par[2][1] == 0:
                elevation = 90
            else:
                tanelevation = cell_par[2][2] / np.sqrt(cell_par[2][0] ** 2 + cell_par[2][1] ** 2)
                if cell_par[2][2] >= 0:
                    elevation = np.arctan(tanelevation) / np.pi * 180
                else:
                    elevation = np.arctan(tanelevation) / np.pi * 180 + 180
            self.view_widget.setCameraPosition(elevation=elevation, azimuth=theta)
        except:
            ...

    def newfile(self):
        """ To creat new cif file"""
        try:
            self.new_window.close()
        except:
            pass
        try:
            self.new_window = New_file_window()
            self.new_window.signal_determine.connect(self.opennewfile)
        except:
            ...

    def opennewfile(self, str):
        """ To open file."""
        try:
            self.fileName_choose = str
            self.openfile(new=True)
        except:
            pass

    def Edit_text(self):
        """ To react to the action_Edit_text on tree."""
        current = self.project_tree.currentItem()
        try:
            try:
                self.setTexttree_window.close()
            except:
                 pass
            try:
                self.setExportDic_window.close()
            except:
                pass
            if len(self.project_tree.selectedItems()) == 1:
                if current.text(1) != "":
                    self.setTexttree_window = tree_set_text_window()
                    self.setTexttree_window.signal_edit_text.connect(self.edit_tree_text)
                    self.setTexttree_window.signal_edit_text1.connect(self.edit_tree_text_conduct)
                else:
                    self.setExportDic_window = tree_Edit_dic_window(current.text(0))
                    self.setExportDic_window.signal_edit_text.connect(self.edit_dic_name)
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except:
            ...

    def doubleclickedontree(self):
        """ To react to the action_Edit_text on tree."""
        current = self.project_tree.currentItem()
        try:
            try:
                self.setTexttree_window.close()
            except:
                 pass
            if current.text(1) != "":
                self.setTexttree_window = tree_set_text_window()
                self.setTexttree_window.signal_edit_text.connect(self.edit_tree_text)
                self.setTexttree_window.signal_edit_text1.connect(self.edit_tree_text_conduct)
        except:
            ...

    def edit_dic_name(self, text):
        current = self.project_tree.currentItem()
        try:
            current.setText(0, text)
        except:
            ...

    def edit_tree_text_conduct(self, text):
        """ To edit the text of tree. (Conductor, Semiconductor, Insulator)"""
        current = self.project_tree.currentItem()
        try:
            current.setText(3, text)
        except:
            ...

    def edit_tree_text(self, text):
        """ To edit the text of tree. (bulk, layer, stack)"""
        current = self.project_tree.currentItem()
        try:
            current.setText(0, text)
            if text == 'bulk' or text == 'stack':      # 如果是bulk或stack没有conductivity
                current.setText(3, "")
            elif text == 'layer':
                if current.text(3) == '':
                    Text3column = self.judgeconductivity(self.dic_Formula_Atoms[current.text(1)])
                    current.setText(3, Text3column)
        except:
            ...

    # 点击树
    def onTreeClicked(self):
        """ Reaction of click on self.treewidget """
        try:
            try:  # 移去drag rectangle选中的Item
                global x_rectangle_dic
                for obj in x_rectangle_dic.values():
                    self.view_widget.removeItem(obj)
                x_rectangle_dic = {}
            except:
                ...
            current = self.project_tree.currentItem()
            if current.text(1) == "":
                self.view_widget.clear()
                self.clear_text()
                self.Atomsobject = None
            elif current.text(1) not in self.dic_Hetero_junction.keys():
                self.plot(self.dic_Formula_Atoms[current.text(1)], dictionary=False, globalAtomsobject=True)  # 不加对象
            else:
                self.plot(self.dic_Formula_Atoms[current.text(1)], dictionary=False, globalAtomsobject=True,
                          Hetero_junction=self.dic_Hetero_junction[current.text(1)])
        except:
            ...

    def rightMenuShow_Object(self):
        try:
            rightMenu = QtWidgets.QMenu(self.menuBar)
            self.replace_atom_Action = QtWidgets.QAction(self)
            self.replace_atom_Action.setText("Replace Atom")
            self.replace_atom_Action.triggered.connect(self.replace_atom)
            rightMenu.addAction(self.replace_atom_Action)
            self.remove_atom_Action = QtWidgets.QAction(self)
            self.remove_atom_Action.setText("Remove Atom")
            self.remove_atom_Action.triggered.connect(self.remove_atom)
            rightMenu.addAction(self.remove_atom_Action)
            rightMenu.addSeparator()
            self.actionset_unit_cell = QtWidgets.QAction(self)
            self.actionset_unit_cell.setText("Set Unit Cell")
            self.actionset_unit_cell.triggered.connect(self.set_unit_cell)
            rightMenu.addAction(self.actionset_unit_cell)
            rightMenu.exec_(QtGui.QCursor.pos())
        except Exception as e:
            print(e)


    def rightMenuShow_tree(self):
        try:
            rightMenu = QtWidgets.QMenu(self.menuBar)
            self.delete_action = QtWidgets.QAction(self)
            self.delete_action.setText("Delete")
            self.delete_action.triggered.connect(self.deleteobject)
            rightMenu.addAction(self.delete_action)
            rightMenu.addSeparator()
            self.Edit_text_action = QtWidgets.QAction(self)
            self.Edit_text_action.setText("Edit text")
            self.Edit_text_action.triggered.connect(self.Edit_text)
            rightMenu.addAction(self.Edit_text_action)
            self.Create_folder_action = QtWidgets.QAction(self)
            self.Create_folder_action.setText("Create a folder")
            self.Create_folder_action.triggered.connect(self.Create_folder)
            rightMenu.addAction(self.Create_folder_action)
            rightMenu.exec_(QtGui.QCursor.pos())
        except Exception as e:
            print(e)

    def Create_folder(self):
        """ To create a folder to export."""
        try:
            self.create_folderwindow.close()
        except:
             pass
        try:
            self.create_folderwindow = Create_folder_window()
            self.create_folderwindow.signal_emit_information.connect(self.after_create_folder_window)
        except:
            ...

    def after_create_folder_window(self, name):
        try:
            root = QtWidgets.QTreeWidgetItem(self.project_tree)
            root.setText(0, name)
        except:
            ...

    def rightMenuShow_view_widget(self):
        """ Mouse right button click on self.view_widget to show a menubar."""
        try:
            if self.view_widget.num_mouse_track == 0:
                rightMenu = QtWidgets.QMenu()
                self.add_vacuum_layer = QtWidgets.QAction(self)
                self.add_vacuum_layer.setText("Add vacuum layer")
                self.add_vacuum_layer.triggered.connect(self.addvacuumlayer)
                self.add_vacuum_layer.setEnabled(True)
                rightMenu.addAction(self.add_vacuum_layer)
                self.add_cell_Action = QtWidgets.QAction(self)
                self.add_cell_Action.setText("Create supercell")
                self.add_cell_Action.triggered.connect(self.plot_add_cell)
                rightMenu.addAction(self.add_cell_Action)
                rightMenu.addSeparator()
                self.remove_atom_Action = QtWidgets.QAction(self)
                self.remove_atom_Action.setText("Remove Atom")
                self.remove_atom_Action.triggered.connect(self.remove_atom)
                self.remove_atom_Action.setEnabled(True)
                rightMenu.addAction(self.remove_atom_Action)
                self.replace_atom_Action = QtWidgets.QAction(self)
                self.replace_atom_Action.setText("Replace Atom")
                self.replace_atom_Action.triggered.connect(self.replace_atom)
                self.replace_atom_Action.setEnabled(True)
                rightMenu.addAction(self.replace_atom_Action)
                self.paste_atom_Action = QtWidgets.QAction(self)
                self.paste_atom_Action.setText("Add Atom")
                self.paste_atom_Action.triggered.connect(self.paste_atom)
                self.paste_atom_Action.setEnabled(True)
                rightMenu.addAction(self.paste_atom_Action)
                rightMenu.addSeparator()
                self.setcell_Action = QtWidgets.QAction(self)
                self.setcell_Action.setText("Set cell")
                self.setcell_Action.triggered.connect(self.setcell_atom_no_move)
                rightMenu.addAction(self.setcell_Action)

                rightMenu.addSeparator()
                self.actionset_acute_angle = QtWidgets.QAction(self)
                self.actionset_acute_angle.setText("γ to sharp angle")
                self.actionset_acute_angle.triggered.connect(self.to_acute_angle)
                self.actionset_acute_angle.setEnabled(True)
                rightMenu.addAction(self.actionset_acute_angle)
                self.actionset_abuse_angle = QtWidgets.QAction(self)
                self.actionset_abuse_angle.setText("γ to blunt angle")
                self.actionset_abuse_angle.triggered.connect(self.to_blunt_angle)
                self.actionset_abuse_angle.setEnabled(True)
                rightMenu.addAction(self.actionset_abuse_angle)
                self.changegama = QtWidgets.QAction(self)
                self.changegama.setText("Change γ degree")
                self.changegama.triggered.connect(self.change_gama)
                self.changegama.setEnabled(False)
                rightMenu.addAction(self.changegama)
                rightMenu.exec_(QtGui.QCursor.pos())
        except Exception as e:
            print(e)

    def change_gama(self):
        """ To change gama according to the clients."""
        current = self.project_tree.currentItem()
        try:
            if current is not None and current.parent() is not None:
                if len(self.project_tree.selectedItems()) == 1 and self.Atomsobject is not None:  # 用户只选择一个的情况
                    try:
                        self.change_gama_window.close()
                    except:
                         pass
                    try:
                        self.change_gama_window = Change_gama_window()
                        # self.change_gama_window.signal_gama_error.connect(self.trans_gama)
                    except Exception as e:
                        print(e)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:  # parent()是self.tree的情况
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except:
            ...

    def set_unit_cell(self):
        try:
            global obj_lst_table
            if len(obj_lst_table) != 3:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose 3 atoms in Objectbox')
            else:
                global index_lst
                pos1 = self.Atomsobject.get_positions()[index_lst[0]]
                pos2 = self.Atomsobject.get_positions()[index_lst[1]]
                pos3 = self.Atomsobject.get_positions()[index_lst[2]]
                a_vector = pos2 - pos1
                b_vector = pos3 - pos1
                a_vector[2] = 0
                b_vector[2] = 0
                A = self.Atomsobject.get_cell().T
                r1 = np.linalg.solve(A, a_vector)  # Crystal coordinate
                r2 = np.linalg.solve(A, b_vector)
                layer = cut(self.Atomsobject, a=r1.T,
                            b=r2.T, c=(0, 0, 1), origo=(pos1[0], pos1[1], 0))
                Atoms = self.deal_with_rotate(layer)
                self.plot(Atoms)
                current = self.project_tree.currentItem()
                childx = QtWidgets.QTreeWidgetItem(current.parent())
                childx.setText(1, self.dirkey)
                childx.setText(0, 'layer')
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error', 'Please choose 3 atoms in Objectbox')

    def to_acute_angle(self):
        try:
            current = self.project_tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.project_tree.selectedItems()) == 1 and self.Atomsobject is not None:  # 用户只选择一个的情况
                    Atomsobject = deepcopy(self.Atomsobject)
                    if Atomsobject.get_cell_lengths_and_angles()[5] > 90:
                        Crystype = CrystalOperationsNum1(Atomsobject)
                        Atomsobject = Crystype.gamaToAcuteangle()
                        self.plot(Atomsobject, clear=True, globalAtomsobject=False, dictionary=True)
                        childx = QtWidgets.QTreeWidgetItem(self.project_tree.currentItem().parent())
                        childx.setText(1, self.dirkey)
                        childx.setText(0, current.text(0))
                        childx.setText(3, current.text(3))
                        childx.setText(2, current.text(2))
                        self.Atomsobject = None
                    else:
                        QtWidgets.QMessageBox.information(self, 'Information:', 'γ is acute angle already')
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def to_blunt_angle(self):
        try:
            current = self.project_tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.project_tree.selectedItems()) == 1 and self.Atomsobject is not None:  # 用户只选择一个的情况
                    Atomsobject = deepcopy(self.Atomsobject)
                    if Atomsobject.get_cell_lengths_and_angles()[5] < 90:
                        Crystype = CrystalOperationsNum1(Atomsobject)
                        Atomsobject = Crystype.gamaToObtuseangle()
                        self.plot(Atomsobject, clear=True, globalAtomsobject=False, dictionary=True)
                        childx = QtWidgets.QTreeWidgetItem(self.project_tree.currentItem().parent())
                        childx.setText(1, self.dirkey)
                        childx.setText(0, current.text(0))
                        childx.setText(3, current.text(3))
                        childx.setText(2, current.text(2))
                        self.Atomsobject = None
                    else:
                        QtWidgets.QMessageBox.information(self, 'Information:', 'γ is blunt angle already')
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except:
            ...

    def deal_with_rotate(self, Atomsobject):
        """ after rotate cut no atoms left"""
        try:
            cell_par = Atomsobject.get_cell_lengths_and_angles()
            Atomsobject.set_cell(cell_par, scale_atoms=True)
            return Atomsobject
        except Exception as e:
            print(e)

    def addvacuumlayer(self):
        """ To add vacuum layer."""
        current = self.project_tree.currentItem()
        try:
            if current is not None and current.parent() is not None:
                if len(self.project_tree.selectedItems()) == 1 and self.Atomsobject is not None:
                    try:
                        self.vacuum_window.close()
                    except:
                         pass
                    self.vacuum_window = add_vacuum_layer_window()
                    self.vacuum_window.sinal_vacuum.connect(self.set_vacuum_layer)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def set_vacuum_layer(self, vacuum_dis, addproject=True):  # c轴方向加上真空层距离，非垂直距离
        """
        To create and draw a crystal with vacuum layer
        :param vacuum_layer: the vacuum layer distance
        create a Atoms object with vacuum layer
        """
        try:
            current = self.project_tree.currentItem()
            Atomsobject = deepcopy(self.Atomsobject)
            c_length = Atomsobject.get_cell_lengths_and_angles()[2]
            cell_par = Atomsobject.get_cell()
            cell_par[2] *= (c_length + vacuum_dis) / c_length
            Atomsobject.set_cell(cell_par)
            Atomsobject.translate(cell_par[2] / (c_length + vacuum_dis) * vacuum_dis / 2)
            if addproject == True:
                self.plot(Atomsobject, clear=True, globalAtomsobject=False, dictionary=True)
                childx = QtWidgets.QTreeWidgetItem(self.project_tree.currentItem().parent())
                childx.setText(1, self.dirkey)
                childx.setText(0, current.text(0))
                childx.setText(3, current.text(3))
                childx.setText(2, current.text(2))
            else:
                return Atomsobject
        except Exception as e:
            print(e)

    def deleteobject(self):
        """
        To delete the Atomsobject and Items which client chooses in the self.tree widget
        """
        try:
            if self.exit_box_show == True:
                cb = QtWidgets.QCheckBox("Don't ask me again.")
                cb.setFont(QtGui.QFont("宋体", 10))
                msg_box = QtWidgets.QMessageBox(QtWidgets.QMessageBox.Question, self.tr("Confirm delete"),
                                      self.tr("Are you sure you want to delete files you choose?"))
                msg_box.addButton(QtWidgets.QMessageBox.Yes)
                msg_box.addButton(QtWidgets.QMessageBox.No)
                msg_box.setCheckBox(cb)
                msg_box.setWindowIcon(QtGui.QIcon("Main.png"))
                reply = msg_box.exec_()
                try:
                    if cb.isChecked():
                        self.exit_box_show = False
                        if reply == QtWidgets.QMessageBox.Yes:
                            self.after_deleteobject()
                    else:
                        if reply == QtWidgets.QMessageBox.Yes:
                            self.after_deleteobject()
                except:
                    pass
            else:
                self.after_deleteobject()
        except:
            pass

    def after_deleteobject(self):
        try:
            current = self.project_tree.currentItem()
            try:
                if self.show_atom_element_judge:
                    self.show_atom_element_judge = False
                    self.plot(self.Atomsobject)
                    self.show_atom_element_judge = True
                if self.show_atom_index_judge:
                    self.show_atom_index_judge = False
                    self.plot(self.Atomsobject)
                    self.show_atom_index_judge = True
            except:
                pass
            if len(self.project_tree.selectedItems()) == 1:  # 用户只选择一个的情况
                if current is not None and current.parent() is not None:  # parent()不是self.tree的情况
                    self.view_widget.clear()
                    self.clear_text()
                    self.dic_Formula_Atoms.pop(current.text(1))  # 对象字典去除该对象
                    try:
                        self.dic_Hetero_junction.pop(current.text(1))
                    except Exception as e:
                        print(e)

                    if current.parent().childCount() > 1:
                        current.parent().removeChild(current)
                    else:
                        self.project_tree.takeTopLevelItem(self.project_tree.indexOfTopLevelItem(current.parent()))
                elif current == None:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
                elif current.parent() == None:  # parent()是self.tree的情况
                    for i in range(current.childCount()):
                        self.dic_Formula_Atoms.pop(current.child(i).text(1))
                        try:
                            self.dic_Hetero_junction.pop(current.child(i).text(1))
                        except Exception as e:
                            print(e)
                    self.project_tree.takeTopLevelItem(self.project_tree.indexOfTopLevelItem(current))
                    self.clear_text()
                    # self.clear_painter()
                    self.view_widget.clear()
            elif len(self.project_tree.selectedItems()) > 1:
                for item in self.project_tree.selectedItems():
                    if item is not None and item.parent() is not None:  # parent()不是self.tree的情况
                        # self.clear_painter()
                        self.view_widget.clear()

                        self.clear_text()
                        self.dic_Formula_Atoms.pop(item.text(1))  # 对象字典去除该对象
                        try:
                            self.dic_Hetero_junction.pop(item.text(1))
                        except Exception as e:
                            print(e)
                        if item.parent().childCount() > 1:
                            item.parent().removeChild(item)
                        else:
                            self.project_tree.takeTopLevelItem(self.project_tree.indexOfTopLevelItem(item.parent()))
                    elif item.parent() == None:  # parent()是self.tree的情况
                        for i in range(item.childCount()):
                            self.dic_Formula_Atoms.pop(item.child(i).text(1))
                            try:
                                self.dic_Hetero_junction.pop(item.child(i).text(1))
                            except Exception as e:
                                print(e)
                        self.project_tree.takeTopLevelItem(self.project_tree.indexOfTopLevelItem(item))
                        self.clear_text()
                        self.view_widget.clear()
                        # self.clear_painter()
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            self.Atomsobject = None
        except:
            pass

    def setcell_atom_no_move(self):
        try:
            current = self.project_tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.project_tree.selectedItems()) == 1 and self.Atomsobject is not None:  # 用户只选择一个的情况
                    try:
                        self.setcell_window.close()
                    except:
                        pass
                    cell_par = self.Atomsobject.get_cell_lengths_and_angles()
                    self.setcell_window = set_cell_window(cell_par)
                    self.setcell_window.signal_emit_information.connect(self.after_setcellwindow)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def after_setcellwindow(self, number, a, b, c, alpha, betta, gama):
        try:
            current = self.project_tree.currentItem()
            Atomsobject = deepcopy(self.Atomsobject)
            cell_par = Atomsobject.get_cell_lengths_and_angles()
            cell_par[0] = a
            cell_par[1] = b
            cell_par[2] = c
            cell_par[3] = alpha
            cell_par[4] = betta
            cell_par[5] = gama
            Atomsobject.set_cell(cell_par, scale_atoms=number)
            self.plot(Atomsobject, clear=True, dictionary=True)
            child1 = QtWidgets.QTreeWidgetItem(current.parent())
            child1.setText(0, current.text(0))
            child1.setText(1, self.dirkey)
            child1.setText(2, current.text(2))
            child1.setText(3, current.text(3))
            self.setcell_window.close()
        except Exception as e:
            print(e)

    def remove_atom(self):
        """To remove the Atoms which is chosen by the client."""
        try:
            global index_lst
            current = self.project_tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.project_tree.selectedItems()) == 1 and self.Atomsobject is not None:  # 用户只选择一个的情况
                    Atomsobject = deepcopy(self.Atomsobject)
                    try:
                        if self.remove_atom_index_num:
                            Atomsobject.pop(self.remove_atom_index_num)
                            self.remove_atom_index_num = None
                            self.plot(Atomsobject, clear=True, dictionary=True)
                            child1 = QtWidgets.QTreeWidgetItem(current.parent())
                            child1.setText(0, current.text(0))
                            child1.setText(1, self.dirkey)
                            child1.setText(2, current.text(2))
                            child1.setText(3, current.text(3))
                            self.Atomsobject = None
                        elif index_lst:
                            for _ in sorted(index_lst, reverse=True):
                                Atomsobject.pop(_)
                            index_lst = []
                            self.plot(Atomsobject, clear=True, dictionary=True)
                            child1 = QtWidgets.QTreeWidgetItem(current.parent())
                            child1.setText(0, current.text(0))
                            child1.setText(1, self.dirkey)
                            child1.setText(2, current.text(2))
                            child1.setText(3, current.text(3))
                            self.Atomsobject = None
                    except:
                        QtWidgets.QMessageBox.warning(self, 'error', 'Please doubleclick an atom in screen or '
                                                                   'choose atoms in object box.')

                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def replace_atom(self):
        """To replace the Atoms which is chosen by the client."""
        try:
            current = self.project_tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.project_tree.selectedItems()) == 1 and self.Atomsobject is not None:  # 用户只选择一个的情况
                    try:
                        self.replace_window.close()
                    except:
                         pass
                    global index_lst
                    try:
                        if self.remove_atom_index_num or index_lst:
                            self.replace_window = Replace_Atom_window()
                            self.replace_window.signal_atomic_number.connect(self.afterreplace_window)
                    except Exception as e:
                        QtWidgets.QMessageBox.warning(self, 'error', 'Please doubleclick an atom in screen or '
                                                                     'choose atoms in object box.')
                        print(e)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)


    def afterreplace_window(self, atomicnumber):
        """ To replace an Atom."""
        try:
            global index_lst
            current = self.project_tree.currentItem()
            Atomsobject = deepcopy(self.Atomsobject)
            if self.remove_atom_index_num:
                Atomsobject.numbers[self.remove_atom_index_num] = atomicnumber
                self.remove_atom_index_num = None
                self.plot(Atomsobject, clear=True, dictionary=True)
                child1 = QtWidgets.QTreeWidgetItem(current.parent())
                child1.setText(0, current.text(0))
                child1.setText(1, self.dirkey)
                child1.setText(2, current.text(2))
                child1.setText(3, current.text(3))
                self.Atomsobject = None   # Can't operate
            elif index_lst:
                for index in index_lst:
                    Atomsobject.numbers[index] = atomicnumber
                    self.plot(Atomsobject, clear=True, dictionary=True)
                    child1 = QtWidgets.QTreeWidgetItem(current.parent())
                    child1.setText(0, current.text(0))
                    child1.setText(1, self.dirkey)
                    child1.setText(2, current.text(2))
                    child1.setText(3, current.text(3))
                index_lst = []
                self.Atomsobject = None
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, 'error', 'Please click an atom or choose an atom in Object box.')
            print(e)

    def paste_atom(self):
        """To paste the Atoms and vector which is chosen by the client."""
        try:
            current = self.project_tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.project_tree.selectedItems()) == 1 and self.Atomsobject is not None:  # 用户只选择一个的情况
                    try:
                        self.paste_window.close()
                    except:
                         pass
                    self.paste_window = add_atomwindow()
                    self.paste_window.signal_emit_information.connect(self.afterpaste_window)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def afterpaste_window(self, atomic_number, pos):
        try:
            current = self.project_tree.currentItem()
            Atomsobject = deepcopy(self.Atomsobject)
            new_atom = Atom(atomic_number, position=pos)
            Atomsobject.append(new_atom)
            self.plot(Atomsobject, clear=True, dictionary=True)
            child1 = QtWidgets.QTreeWidgetItem(current.parent())
            child1.setText(0, current.text(0))
            child1.setText(1, self.dirkey)
            child1.setText(2, current.text(2))
            child1.setText(3, current.text(3))
            self.Atomsobject = None
        except Exception as e:
            print(e)

    # 导出cif文件

    def Save_cif(self):
        """To save one cif file which client chooses in the self.tree widget."""
        try:
            current = self.project_tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.project_tree.selectedItems()) == 1 and self.Atomsobject is not None:  # 用户只选择一个的情况
                    output_file = QtWidgets.QFileDialog.getSaveFileName(self, "Save as cif", "{}.cif".format(
                        self.Atomsobject.get_chemical_formula(mode='hill')))
                    f = AseAtomsAdaptor.get_structure(self.Atomsobject)
                    print(AseAtomsAdaptor.get_structure(self.Atomsobject))
                    f1 = str(f) + '\n'
                    j_pv_lst = findall('abc\s\s\s:(.*?)\n', f1)[0]  # abc   :  19.257300  19.569178  21.133988
                    j1_pv_lst = j_pv_lst.split(' ')  # abc   :   6.419100   6.523059   7.044663
                    while '' in j1_pv_lst:
                        j1_pv_lst.remove('')
                    par_lst_matrix = self.Atomsobject.get_cell()
                    material = findall('Full\sFormula\s(.*?)\n', f1)[0].lstrip('(').rstrip(')')
                    elements = material.split(' ')
                    zmb_lst = [chr(i) for i in range(97, 123)] + [chr(i) for i in range(65, 91)]
                    szb_lst = [str(i) for i in range(0, 10)]
                    element_lst = []
                    number_lst = []
                    for element in elements:  # 这个循环之后，得到原子与对应的原子数两个列表
                        letter_lst = list(element)
                        symbol_lst = []
                        num_lst = []
                        element_lst1 = []
                        number_lst1 = []
                        for letter in letter_lst:
                            if letter in szb_lst:
                                num_lst.append(letter)
                                # print(num_lst)
                            if letter in zmb_lst:
                                symbol_lst.append(letter)
                            element1 = ''.join(symbol_lst)
                            element_lst1.append(element1)
                            number1 = ''.join(num_lst)
                            number_lst1.append(number1)
                        ys = 'a'  # 元素
                        gs = '0'  # 个数
                        for i in element_lst1:
                            if len(i) >= len(ys):
                                ys = i
                        element_lst.append(i)
                        for i in number_lst1:
                            if len(i) >= len(gs):
                                gs = i
                        number_lst.append(i)

                    par_lst_species = []  # 用于Cifwrite参数(species)的
                    for i, element in enumerate(element_lst):
                        num = int(number_lst[i])
                        for j in range(num):
                            par_lst_species.append(element)
                    # 每个原子的坐标
                    ord_lst = []  # 最终Cifwriter所需要的coords参数
                    ord_lst2 = []  # 储存的形式为
                    for element in element_lst:
                        ord_lst1 = findall(element + '\s\s\s\s(.*?)\n', f1)
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
                    slab.write_file(output_file[0])
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    # 扩胞
    def plot_add_cell(self):  # add_cellAction 连接
        """To plot crys which has added cell on time."""
        try:
            # self.key_double_click = False
            current = self.project_tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.project_tree.selectedItems()) == 1 and self.Atomsobject is not None:  # 用户只选择一个的情况
                    self.plot_cell = False
                    self.c = Create_supercell_window()
                    # self.c.close.triggered.connect(self.key_double_click_setTrue)
                    self.repeat_a = 1
                    self.repeat_b = 1
                    self.repeat_c = 1
                    self.c.setUnitCell.clicked.connect(self.setUnitCell)
                    self.c.determine.clicked.connect(self.add_celloutput)
                    self.c.signala.connect(self.repeata)
                    self.c.signalb.connect(self.repeatb)
                    self.c.signalc.connect(self.repeatc)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def add_celloutput(self):
        """To create crys after create supercell."""
        try:
            root = QtWidgets.QTreeWidgetItem(self.project_tree)
            self.plot(self.after_repeat, dictionary=True)
            root.setText(0, self.dirkey)
            child = QtWidgets.QTreeWidgetItem(root)
            child.setText(0, "bulk")
            child.setText(1, self.dirkey)
            self.c.close()
            self.project_tree.expandAll()
        except Exception as e:
            print(e)

    def setUnitCell(self):
        self.plot_cell = True
        self.key_double_click = True
        self.repeat_add_cell()

    def repeata(self, int):
        self.repeat_a = int
        self.repeat_add_cell()

    def repeatb(self, int):
        self.repeat_b = int
        self.repeat_add_cell()

    def repeatc(self, int):
        self.repeat_c = int
        self.repeat_add_cell()

    def repeat_add_cell(self):
        try:
            if self.plot_cell == False:
                a = self.Atomsobject.get_cell()
                after_repeat = self.Atomsobject.repeat((self.repeat_a, self.repeat_b, self.repeat_c))
                after_repeat.set_cell(a)
                self.plot(after_repeat, dictionary=False, globalAtomsobject=False, object=False)  # 中间变量，不用加入object
            elif self.plot_cell == True:
                print("self.repeat_a = ", self.repeat_a)
                self.after_repeat = self.Atomsobject.repeat((self.repeat_a, self.repeat_b, self.repeat_c))
                self.plot(self.after_repeat, dictionary=False, globalAtomsobject=False, object=False)  # 中间变量，不用加入object
                self.plot_cell = False
        except Exception as e:
            print(e)

    def eventFilter(self, watched, event):
        """ One Fuctions:
         MouseDblClick to choose an Atom.
        :param watched: the widget where the event happens.
        :param event: ...
        """
        if watched == self.view_widget:
            if event.type() == QtCore.QEvent.MouseButtonDblClick and self.key_double_click == True and \
                    self.view_widget.num_mouse_track != 1:
                try:
                    mousex = event.pos().x()
                    mousey = event.pos().y()
                    try:
                        length_min = 25
                        choose_position = ...
                        for positions in self.global_positions_lst:  # 获得位置在屏幕上的投影的屏幕坐标
                            screentuple = gluProject(positions[0] + self.translate_pos[0], positions[1] +
                                                     self.translate_pos[1], positions[2] + self.translate_pos[2])
                            screenlen = np.sqrt((screentuple[0] - mousex) ** 2 +
                                                ((self.view_widget.height() - screentuple[1]) - mousey) ** 2)
                            if screenlen < length_min:
                                length_min = screenlen
                                choose_position = positions
                        if length_min != 25:
                            if self.Atomsobject:
                                choose_position_index = self.global_positions_lst.index(choose_position)
                                self.tab.setCurrentWidget(self.Atom_tab)
                                self.clickatoms_showmessage(choose_position_index)
                            else:
                                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
                        else:
                            global x
                            try:
                                self.view_widget.removeItem(x)
                                self.tab.setCurrentWidget(self.Crystal_tab)
                                self.Atom_textwidget.setText("")
                                self.objectTextwidget.setText("")
                            except Exception as e:
                                print(e)

                            try:
                                global obj_lst_table
                                for obj in obj_lst_table:
                                    self.view_widget.removeItem(obj)
                                    self.tab.setCurrentWidget(self.Crystal_tab)
                                    self.Atom_textwidget.setText("")
                                    self.objectTextwidget.setText("")
                            except Exception as e:
                                print(e)

                            try:
                                global x_rectangle_dic
                                for obj in x_rectangle_dic.values():
                                    self.view_widget.removeItem(obj)
                                    self.tab.setCurrentWidget(self.Crystal_tab)
                                    self.Atom_textwidget.setText("")
                                    self.objectTextwidget.setText("")
                                x_rectangle_dic = {}
                            except Exception as e:
                                print(e)
                    except Exception as e:
                        print(e)
                except Exception as e:
                    print(e)
        return pg.opengl.GLViewWidget.eventFilter(self, watched, event)

    # 双击原子后显示的原子信息
    def clickatoms_showmessage(self, index, clear=True):
        """
        After Mouse dblClicked, to show the atom's message
        :param index: the index of atom which was chosen by the client
        :param clear: To clear the x chosen before.(choose only one Atoms by mouse dblClick)
        """
        try:
            Atomic_number = self.Atomsobject.get_atomic_numbers()[index]
        except Exception as e:
            print(e)
        if clear == True:
            global x, obj_lst_table, index_lst
            try:
                index_lst = []
                self.view_widget.removeItem(x)
            except Exception as e:
                print(e)
            try:
                for Atom in obj_lst_table:
                    self.view_widget.removeItem(Atom)
                obj_lst_table = []
            except Exception as e:
                print(e)
        try:
            md = pg.opengl.MeshData.sphere(rows=10, cols=20)
            x = pg.opengl.GLMeshItem(meshdata=md, smooth=True, color=self.Atoms_color[self.Atoms_num_lst[index]],
                              shader='shaded', drawEdges=True)
            x.translate(self.global_positions_lst[index][0], self.global_positions_lst[index][1], self.global_positions_lst[index][2])
            x.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])
            if self.plot_num == 0:
                x.scale(self.Atom_radii[Atomic_number] + .01,
                        self.Atom_radii[Atomic_number] + .01,
                        self.Atom_radii[Atomic_number] + .01)
            elif self.plot_num == 1:
                x.scale(self.Atom_radii[Atomic_number]/2 + .005,
                        self.Atom_radii[Atomic_number]/2 + .005,
                        self.Atom_radii[Atomic_number]/2 + .005)
            elif self.plot_num == 2:
                x.scale(0.1, 0.1, 0.1)

            self.view_widget.addItem(x)
            zhijiao_system = self.Atomsobject.get_cell()
            A = zhijiao_system.T
            b = np.array(self.global_positions_lst[index]).T
            r = np.linalg.solve(A, b)  # 求解线性方程组，直角坐标系下----用晶胞坐标系表示
            text = crys_data.Element_formula[Atomic_number] + "(index:" + str(index) + ")" + "\n" + "Atomic number:" + str(
                Atomic_number) + "\n" + "[xyz]:" + str(self.global_positions_lst[index]) + \
                   "\n" + "[uvw]:" + str(r)
            self.Atom_textwidget.setText(text)
            self.remove_atom_index_num = index
        except:
            pass

    def openfile(self, new=False):
        """
        To open a file, draw it and put it in the self.tree widget
        :param new:create two forms, open file and new file
        """
        if new == False:
            self.fileName_choose, filetype = QtWidgets.QFileDialog.getOpenFileName(self, '选择文件', '', 'files(*.cif , *.vasp)')
        if self.fileName_choose:
            try:
                print("in")
                self.Atomsobject = deepcopy(read(self.fileName_choose, index=None, format=None, parallel=True))
                root = QtWidgets.QTreeWidgetItem(self.project_tree)
                self.plot(self.Atomsobject, dictionary=True)
                root.setText(0, self.dirkey)
                child = QtWidgets.QTreeWidgetItem(root)
                child.setText(0, "bulk")
                child.setText(1, self.dirkey)
                self.project_tree.expandAll()
                self.view_widget.setText_init("", "")
            except:
                pass

    def clear_text(self):
        """ To clear the textwidget and objectbox."""
        try:
            self.object_tablewidget.setRowCount(0)  # 清空object——box的信息
        except:
            pass

        try:
            self.objectTextwidget.setText("")
        except Exception as e:
            print(e)
        try:
            self.Crystal_textwidget.setText("")
        except Exception as e:
            print(e)
        try:
            self.Atom_textwidget.setText("")
        except Exception as e:
            print(e)

    def add_hetero_junction(self, Hetero_junction_class):
        try:
            Change_key = Hetero_junction_class.Name
            self.dic_Hetero_junction[Change_key] = Hetero_junction_class  # 在self.dic_formula_Atoms里添加Atomsobject
            print(self.dic_Hetero_junction)
        except Exception as e:
            print(e)

    def plot(self, Atomsobject, plot=True, clear=True, dictionary=True, globalAtomsobject=True,
             object=True, Hetero_tab=True, Hetero_junction=False):  # plot 会创造对象self.Atomsobject,目前有bug
        """
        Normally, to plot an ase.Atoms object.
        :param Atomsobject:the ase.Atoms object to be plotted
        :param clear: to clear the painter if clear == True
        :param dictionary: add the Atomsobject to self.dic_Formula_Atoms(dictionary)
        :param globalAtomsobject: whether create self.Atomsobject or not.
        :param object: whether add atoms to self.tablewidget or not
        """
        # glBegin
        Number = len(Atomsobject.get_positions())
        Atomsobject = deepcopy(sort(Atomsobject))
        if clear == True:
            self.view_widget.clear()
            self.clear_text()
        if globalAtomsobject == True:  # 若Atomsobject =True创建全局变量
            self.Atomsobject = deepcopy(Atomsobject)
        if dictionary == True:
            self.dirkey = Atomsobject.get_chemical_formula(mode='hill')
            orijin_dirkey = self.dirkey
            num = 1
            while True:
                if not self.dic_Formula_Atoms.get(self.dirkey, 0):
                    self.dirkey = orijin_dirkey + '({})'.format(num)
                    num += 1
                else:
                    break
            self.dic_Formula_Atoms[self.dirkey] = Atomsobject
        if object == True:
            if Number < self.max_number_show:
                self.add_object(Atomsobject)
        if Hetero_tab == True:
            if Number < self.max_number_show:
                self.search_database_tab(Atomsobject)
        if plot == True:
            try:
                self.view_widget.setText_init("", "")
            except Exception as e:
                print(e)
            self.plot_Atoms(Atomsobject)
            # 画坐标轴
            self.plot_coordinate_system()
            # 画晶胞,画线
            self.plot_cell_line()
            if self.grid_itemnum % 2 == 0:
                self.grid = pg.opengl.GLGridItem()  # 画网格
                self.grid.scale(1, 1, 1)
                self.grid.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])
                self.view_widget.addItem(self.grid)
            else:
                ...
            text = self.all_info(Atomsobject)
            self.tab.setCurrentWidget(self.Crystal_tab)
            if not Hetero_junction:
                self.Crystal_textwidget.setText(text)
            else:
                try:
                    Normal_strain_limit = Hetero_junction.Normal_strain_limit
                    Normal_strain_limit *= 100
                except:
                    pass
                try:
                    Shear_strain_limit = Hetero_junction.Shear_strain_limit
                    Shear_strain_limit *= 100
                except:
                    pass

                try:
                    Area_error_limit = Hetero_junction.Area_error_limit
                    Area_error_limit *= 100
                except:
                    pass
                text = 'Rotate angle(°):' + str(Hetero_junction.Rotate_angle)  + '\t' +  \
                       'Optimal match:' + str(Hetero_junction.Optimal_orientation) + '\n' + \
                       'Normal strain limit(%):' + str(Normal_strain_limit) + '\t' + \
                       'Shear strain limit(%):' + str(Shear_strain_limit) + '\t' + \
                       'Area error limit(%):' + str(Area_error_limit)  + \
                       '\n' + '-' * 60 + \
                        '\n' + text
                self.Crystal_textwidget.setText(text)

    def search_database_tab(self, Atomsobject):
        try:
            formula = Atomsobject.get_chemical_formula(mode='hill')
            self.database_class.exhibit_all_info(fenkuai=False)
            data_info = self.database_class.info_data_base
            tab3_info = [_ for _ in data_info if _[2] == formula]
            self.Hetero_tablewidget.setRowCount(len(tab3_info))
            self.Device_tablewidget.setRowCount(len(tab3_info))
            for index, tab3_infopiece in enumerate(tab3_info):
                for k in range(10, 15):
                    newitem = QtWidgets.QTableWidgetItem(str(list(tab3_infopiece)[k]))
                    newitem.setTextAlignment(5 | 5)
                    self.Hetero_tablewidget.setItem(index, k - 10, newitem)
                for k in range(15, 20):
                    newitem = QtWidgets.QTableWidgetItem(str(list(tab3_infopiece)[k]))
                    newitem.setTextAlignment(5 | 5)
                    self.Device_tablewidget.setItem(index, k - 15, newitem)
            self.Device_tablewidget.setShowGrid(True)
            self.Hetero_tablewidget.setShowGrid(True)
        except Exception as e:
            print(e)

    def plot_Atoms(self, Atomsobject):
        """ To plot Atoms in two ways (Ball and stick/Space filling) """
        try:
            cell = Atomsobject.get_cell()
            self.translate_pos = -(cell[0] + cell[1] + cell[2]) / 2
            positions = Atomsobject.get_positions()
            self.Atoms_num_lst = Atomsobject.get_atomic_numbers()
            self.global_positions_lst = positions.tolist()
            self.cell_array = Atomsobject.get_cell()
            if len(positions) < self.max_number_show:
                k = deepcopy(positions)
                if self.show_atom_index_judge == True:
                    text_lst = [str(i) for i in range(len(k))]
                    k = map(lambda x: x + self.translate_pos, k)
                    self.view_widget.setAtomstext(text_lst, k)
                if self.show_atom_element_judge == True:
                    text_lst = [crys_data.Element_formula[self.Atoms_num_lst[i]] for i in range(len(k))]
                    k = map(lambda x: x + self.translate_pos, k)
                    self.view_widget.setAtomstext(text_lst, k)
                self.key_double_click = True
                if self.plot_num == 0:
                    md = pg.opengl.MeshData.sphere(rows=10, cols=20)
                    for i, pos in enumerate(positions):
                        new_trans_pos = pos + self.translate_pos
                        Atom = pg.opengl.GLMeshItem(meshdata=md, smooth=True,
                                                         color=self.Atoms_color[self.Atoms_num_lst[i]],
                                                         shader='shaded',
                                                         drawFaces=True)
                        Atom.translate(new_trans_pos[0], new_trans_pos[1], new_trans_pos[2])
                        Atom.scale(self.Atom_radii[self.Atoms_num_lst[i]],
                                        self.Atom_radii[self.Atoms_num_lst[i]],
                                        self.Atom_radii[self.Atoms_num_lst[i]])
                        self.view_widget.addItem(Atom)
                elif self.plot_num == 1:
                    md = pg.opengl.MeshData.sphere(rows=10, cols=20)
                    for i, pos in enumerate(positions):
                        Atom = pg.opengl.GLMeshItem(meshdata=md, smooth=True,
                                                    color=self.Atoms_color[self.Atoms_num_lst[i]],
                                                    shader='shaded',
                                                    drawFaces=True)
                        new_trans_pos = pos + self.translate_pos
                        Atom.translate(new_trans_pos[0], new_trans_pos[1], new_trans_pos[2])
                        radii = self.Atom_radii[self.Atoms_num_lst[i]] / 2
                        Atom.scale(radii, radii, radii)
                        self.view_widget.addItem(Atom)
                    all_vector = Atomsobject.get_all_distances(vector=True)
                    all_dis = Atomsobject.get_all_distances(vector=False)
                    order = len(positions)
                    vander_wals_matrix = np.diag([crys_data.vander_wals_radii[self.Atoms_num_lst[i]] for i in range(order)])
                    vander_wals_matrix = all_dis + np.ones((order, order)) * 1.3 - \
                                         np.transpose(np.ones((order, order)) @ vander_wals_matrix) - np.ones((order, order)) @ vander_wals_matrix + \
                                         np.eye(order) * 1000
                    dis_or_not_matrix = (vander_wals_matrix > 0)
                    pos_i, pos_j = np.where(dis_or_not_matrix == 0)
                    for k in range(len(pos_i)):
                        i = pos_i[k]
                        j = pos_j[k]
                        if crys_data.metal_or_not[self.Atoms_num_lst[i]] == 0 or \
                                crys_data.metal_or_not[self.Atoms_num_lst[j]] == 0:
                            # 不画金属键
                            if abs(crys_data.electronegativity[self.Atoms_num_lst[i]] -
                                crys_data.electronegativity[self.Atoms_num_lst[j]]) < 1.5:
                                    # 不画离子键
                                    dis = all_dis[i][j]
                                    cylinder = pg.opengl.MeshData.cylinder(rows=10, cols=20, radius=[.1, .1],
                                                                           length=dis / 2, offset=False)
                                    valence_band = pg.opengl.GLMeshItem(meshdata=cylinder, smooth=True,
                                                                  color=self.Atoms_color[self.Atoms_num_lst[i]],
                                                                  shader='shaded',
                                                                  drawFaces=True)

                                    pos_vec = all_vector[i][j]   # pos[j] - pos[i]
                                    unit_pos_vec = pos_vec / np.linalg.norm(pos_vec)
                                    rotate_axis = unit_pos_vec + np.array([0, 0, 1])
                                    valence_band.rotate(180, rotate_axis[0], rotate_axis[1], rotate_axis[2])
                                    new_trans_pos = positions[i] + self.translate_pos
                                    valence_band.translate(new_trans_pos[0], new_trans_pos[1], new_trans_pos[2])
                                    self.view_widget.addItem(valence_band)
                elif self.plot_num == 2:        # Stick mode
                    md = pg.opengl.MeshData.sphere(rows=10, cols=20)
                    for i, pos in enumerate(positions):
                        Atom = pg.opengl.GLMeshItem(meshdata=md, smooth=True,
                                                    color=self.Atoms_color[self.Atoms_num_lst[i]],
                                                    shader='shaded',
                                                    drawFaces=True)
                        new_trans_pos = pos + self.translate_pos
                        Atom.translate(new_trans_pos[0], new_trans_pos[1], new_trans_pos[2])
                        radii = 0.1
                        Atom.scale(radii, radii, radii)
                        self.view_widget.addItem(Atom)
                    all_vector = Atomsobject.get_all_distances(vector=True)
                    all_dis = Atomsobject.get_all_distances(vector=False)
                    order = len(positions)
                    vander_wals_matrix = np.diag([crys_data.vander_wals_radii[self.Atoms_num_lst[i]] for i in range(order)])
                    vander_wals_matrix = all_dis + np.ones((order, order)) * 1.3 - \
                                         np.transpose(np.ones((order, order)) @ vander_wals_matrix) - np.ones(
                        (order, order)) @ vander_wals_matrix + \
                                         np.eye(order) * 1000
                    dis_or_not_matrix = (vander_wals_matrix > 0)
                    pos_i, pos_j = np.where(dis_or_not_matrix == 0)
                    for k in range(len(pos_i)):
                        i = pos_i[k]
                        j = pos_j[k]
                        if crys_data.metal_or_not[self.Atoms_num_lst[i]] == 0 or \
                                crys_data.metal_or_not[self.Atoms_num_lst[j]] == 0:
                            # 不画金属键
                            if abs(crys_data.electronegativity[self.Atoms_num_lst[i]] -
                                   crys_data.electronegativity[self.Atoms_num_lst[j]]) < 1.5:
                                # 不画离子键
                                dis = all_dis[i][j]
                                cylinder = pg.opengl.MeshData.cylinder(rows=10, cols=20, radius=[.1, .1],
                                                                       length=dis / 2, offset=False)
                                valence_band = pg.opengl.GLMeshItem(meshdata=cylinder, smooth=True,
                                                                    color=self.Atoms_color[self.Atoms_num_lst[i]],
                                                                    shader='shaded',
                                                                    drawFaces=True)

                                pos_vec = all_vector[i][j]  # pos[j] - pos[i]
                                unit_pos_vec = pos_vec / np.linalg.norm(pos_vec)
                                rotate_axis = unit_pos_vec + np.array([0, 0, 1])
                                valence_band.rotate(180, rotate_axis[0], rotate_axis[1], rotate_axis[2])
                                new_trans_pos = positions[i] + self.translate_pos
                                valence_band.translate(new_trans_pos[0], new_trans_pos[1], new_trans_pos[2])
                                self.view_widget.addItem(valence_band)


            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Too many Atoms to show!')
                self.key_double_click = False
        except Exception as e:
            print(e)

    def add_object(self, Atomsobject):
        try:
            Atomsobject_atomic_numbers_lst = Atomsobject.get_atomic_numbers()
            f = len(Atomsobject_atomic_numbers_lst)
            self.object_tablewidget.setRowCount(f)  # 添加信息
            self.object_tablewidget.setVerticalHeaderLabels([str(i) for i in range(f)])
            for i in range(f):
                newitem = QtWidgets.QTableWidgetItem(crys_data.Element_formula[Atomsobject_atomic_numbers_lst[i]])
                newitem.setTextAlignment(5 | 5)
                self.object_tablewidget.setItem(i, 0, newitem)
            for i in range(f):
                newitem = QtWidgets.QTableWidgetItem("")
                color = self.Atoms_color[Atomsobject_atomic_numbers_lst[i]]
                colortrans = (color[0] * 255, color[1] * 255, color[2] * 255)
                self.object_tablewidget.setItem(i, 1, newitem)
                self.object_tablewidget.item(i, 1).setBackground(QtGui.QBrush(
                    QtGui.QColor(colortrans[0], colortrans[1], colortrans[2])))

            for i in range(f):
                newitem = QtWidgets.QTableWidgetItem(
                    str(self.Atom_radii[Atomsobject_atomic_numbers_lst[i]]))
                newitem.setTextAlignment(5 | 5)
                self.object_tablewidget.setItem(i, 2, newitem)
            self.object_tablewidget.setShowGrid(True)
        except Exception as e:
            print(e)

    def plot_coordinate_system(self):
        """ To plot coordinate system."""
        try:
            if self.coordinatesystem_view_num == 2:  # 直角坐标系
                if self.coordinatesystem_fixed_num == 0:
                    plt = MyaxisItem(QtGui.QVector3D(5, 5, 5))
                else:
                    self.view_widget.setcoordinate_text(self.translate_pos, "Cartestian")
                    # plt = MyaxisItem(QtGui.QVector3D(5, 5, 5), fixed=True)
            elif self.coordinatesystem_view_num == 1:  # 1的时候画晶体坐标系
                if self.coordinatesystem_fixed_num == 0:
                    plt = MyaxisItem(cell_coordinate=list(self.cell_array))
                else:
                    self.view_widget.setcoordinate_text(self.translate_pos, "Cell", cell_array=list(self.cell_array))


                    # plt = MyaxisItem(cell_coordinate=list(self.cell_array), fixed=True)

            plt.translate(self.translate_pos[0] - 5, self.translate_pos[1] - 5, self.translate_pos[2])
            self.view_widget.addItem(plt)
        except Exception as e:
            print(e)

    def plot_cell_line(self):
        try:
            if self.cell_view_num % 2 == 0:
                cell_array_line = []
                for i in range(3):
                    if self.cell_array[i][0] != 0:
                        vector1_x = np.array([0, self.cell_array[i][0]])
                        vector1_y = vector1_x / self.cell_array[i][0] * self.cell_array[i][1]
                        vector1_z = vector1_x / self.cell_array[i][0] * self.cell_array[i][2]
                        pts = np.vstack([vector1_x, vector1_y, vector1_z]).transpose()
                        cell_array_line.append(pts)
                        plt = pg.opengl.GLLinePlotItem(pos=pts, color=pg.glColor((1, 2)), width=5 / 10, antialias=True)
                        plt.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])
                        self.view_widget.addItem(plt)
                    elif self.cell_array[i][2] != 0:
                        vector1_z = np.array([0, self.cell_array[i][2]])
                        vector1_y = vector1_z / self.cell_array[i][2] * self.cell_array[i][1]
                        vector1_x = vector1_z / self.cell_array[i][2] * self.cell_array[i][0]
                        pts = np.vstack([vector1_x, vector1_y, vector1_z]).transpose()
                        plt = pg.opengl.GLLinePlotItem(pos=pts, color=pg.glColor((1, 2)), width=5 / 10, antialias=True)
                        plt.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])
                        self.view_widget.addItem(plt)
                        cell_array_line.append(pts)
                    elif self.cell_array[i][1] != 0:
                        vector1_y = np.array([0, self.cell_array[i][1]])
                        vector1_z = vector1_y / self.cell_array[i][1] * self.cell_array[i][2]
                        vector1_x = vector1_y / self.cell_array[i][1] * self.cell_array[i][0]
                        pts = np.vstack([vector1_x, vector1_y, vector1_z]).transpose()
                        plt = pg.opengl.GLLinePlotItem(pos=pts, color=pg.glColor((1, 2)), width=5 / 10, antialias=True)
                        plt.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])
                        self.view_widget.addItem(plt)
                        cell_array_line.append(pts)
                # 画线
                for i in range(3):
                    lst = [0, 1, 2]
                    lst.remove(i)
                    x = cell_array_line[i] + self.cell_array[lst[0]]
                    plt = pg.opengl.GLLinePlotItem(pos=x, color=pg.glColor((1, 2)), width=5 / 10, antialias=True)
                    plt.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])
                    self.view_widget.addItem(plt)
                    y = cell_array_line[i] + self.cell_array[lst[1]]
                    plt = pg.opengl.GLLinePlotItem(pos=y, color=pg.glColor((1, 2)), width=5 / 10, antialias=True)
                    plt.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])
                    self.view_widget.addItem(plt)
                    z = cell_array_line[i] + self.cell_array[lst[0]] + self.cell_array[lst[1]]
                    plt = pg.opengl.GLLinePlotItem(pos=z, color=pg.glColor((1, 2)), width=5 / 10, antialias=True)
                    plt.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])
                    self.view_widget.addItem(plt)
            else:
                pass
        except Exception as e:
            print(e)

    def all_info(self, Atomsobject):
        """
        To get the information of ase.Atoms object.
        :param Atomsobject: ase.Atoms object
        :return: crysinfo
        """
        try:
            crys = Atomsobject
            f = str(AseAtomsAdaptor.get_structure(crys))
            lst = f.split('Sites {}'.format(num))
            f1 = lst[0]
            f2 = lst[1]
            spacegroup = str(get_spacegroup(crys))
            space_group = 'Space group: ' + findall('\s\s(.*)\n', spacegroup)[0]
            crys_info = f1 + '\n' + space_group
            volume = str(crys.get_volume())
            unit_volume = 'Unit cell volume: ' + str(volume)
            crys_info = crys_info + '\n' + unit_volume + ' Å^3'
            f2 = 'Cell params:' + f2
            crys_info = crys_info + '\n\n' + f2
            return crys_info
        except Exception as e:
            print(e)

    def judgetool(self):
        """ To cut self.Atomsobject.(defalut form)"""
        try:
            current = self.project_tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.project_tree.selectedItems()) == 1 and self.Atomsobject is not None:  # 用户只选择一个的情况
                    self.cut_exf_layer(current.parent())
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def cut_customize(self):
        """ Has two form: default and cutomize"""
        try:
            current = self.project_tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.project_tree.selectedItems()) == 1 and self.Atomsobject is not None:  # 用户只选择一个的情况
                    try:
                        self.cutwindow.close()
                    except:
                        pass
                    cell = self.Atomsobject.get_cell()
                    self.cutwindow = cut_window(cell)
                    self.cutwindow.signal_custom_cut_exf.connect(self.custom_cut)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def custom_cut(self, cell_par, number):
        """ To cut according to the client."""
        try:
            customize_cut_object = deepcopy(self.Atomsobject)
            if number == 1:         # "xyz"
                A = customize_cut_object.get_cell().T
                cell_coordinat_system_par = [np.linalg.solve(A, np.array(_).T).T for _ in cell_par]
                layer = cut(customize_cut_object, a=cell_coordinat_system_par[0],
                                b=cell_coordinat_system_par[1],
                                c=cell_coordinat_system_par[2], origo=(0, 0, 0))
            elif number == 2:         # "uvw"
                layer = cut(customize_cut_object, a=cell_par[0],
                                            b=cell_par[1],
                                            c=cell_par[2], origo=(0, 0, 0))
            else:
                cell = customize_cut_object.get_cell()
                a = cell[0]
                b = cell[1]
                c = cell[2]
                h = cell_par[0]
                k = cell_par[1]
                l = cell_par[2]
                cell_par = [a, b, c]
                lst_face_par = [h, k, l]
                if 0 in lst_face_par:
                    num = lst_face_par.count(0)
                    if num == 3:
                        QtWidgets.QMessageBox.warning(self, 'error', "Input error!")
                        return None
                    elif num == 2:
                        r_lst = []
                        for index, par in enumerate(lst_face_par):
                            if par == 0:
                                c = [0, 0, 0]
                                c[index] = 1
                                r_lst.append(np.array(c))
                        r1 = r_lst[0]
                        r2 = r_lst[1]
                    else:
                        index = lst_face_par.index(0)
                        r1 = [0, 0, 0]
                        r1[index] = 1
                        r1 = np.array(r1)
                        c = [cell_par[k] / par for k, par in enumerate(lst_face_par) if k != index]
                        orientation2 = (c[0] - c[1]).T
                        A = customize_cut_object.get_cell().T
                        r2 = np.linalg.solve(A, orientation2).T
                        r3 = np.cross(r1, r2) / \
                         np.linalg.norm(np.cross(r1, r2)) / 20
                    layer = cut(customize_cut_object, a=r1 * 10,
                                b=r2 * 10,
                                c=r3, origo=(0, 0, 0))
                else:
                    orientation1 = b/k - a/h
                    orientation2 = b/k - c/l
                    A = customize_cut_object.get_cell().T
                    b1 = np.array(orientation1).T
                    r1 = np.linalg.solve(A, b1).T  # 求解线性方程组，直角坐标系下----用晶胞坐标系表示
                    b2 = np.array(orientation2).T
                    r2 = np.linalg.solve(A, b2).T
                    r3 = np.cross(r1, r2) / \
                                   np.linalg.norm(np.cross(r1, r2)) / 20
                    layer = cut(customize_cut_object, a=r1 * 10 * h,
                                b=r2 * 10 * k,
                                c=r3, origo=(0, 0, 0))
            Atoms = self.deal_with_rotate(layer)
            self.plot(Atoms)
            current = self.project_tree.currentItem()
            childx = QtWidgets.QTreeWidgetItem(current.parent())
            childx.setText(1, self.dirkey)
            childx.setText(0, 'layer')
        except Exception as e:
            print(e)

    def cut_exf_layer(self, parent, traversal=False):
        """ Default cut, remove the exfoliable layer distance. (vanderwals1 + vanderwals2 -1.3 < dis())"""
        try:
            after_add_cell_self_ATO = self.Atomsobject.repeat((1, 1, 2))
            pos_lst = after_add_cell_self_ATO.get_positions()
            all_dis = after_add_cell_self_ATO.get_all_distances(vector=False)
            order = len(pos_lst)
            Atomic_number_lst = after_add_cell_self_ATO.get_atomic_numbers()
            vander_wals_matrix = np.diag([crys_data.vander_wals_radii[Atomic_number_lst[i]] for i in range(order)])
            vander_wals_matrix = all_dis + np.ones((order, order)) * 1.3 - \
                                 np.transpose(np.ones((order, order)) @ vander_wals_matrix) - np.ones(
                (order, order)) @ vander_wals_matrix
            dis_or_not_matrix = (vander_wals_matrix > 0)
            gouzaolist = [_ + [i] for i, _ in enumerate(pos_lst)]
            gouzaolist.sort(key=lambda x:x[2])# 根据z轴由小到大排序
            min_z = gouzaolist[0][2]
            height = 0
            exfoliat_height = 0
            index_lst = [gouzaolist[0][3]]
            for i in range(len(gouzaolist) - 1):
                if not dis_or_not_matrix[gouzaolist[i][3]][gouzaolist[i + 1][3]]:  # valence bond
                    height += (gouzaolist[i+1][2] - gouzaolist[i][2])
                    index_lst.append(gouzaolist[i+1][3])
                elif (gouzaolist[i+1][2] - gouzaolist[i][2]) / \
                        all_dis[gouzaolist[i][3]][gouzaolist[i + 1][3]] < .5:
                    height += (gouzaolist[i + 1][2] - gouzaolist[i][2])
                    index_lst.append(gouzaolist[i + 1][3])
                else:
                    exfoliat_height = gouzaolist[i+1][2] - gouzaolist[i][2]
                    break
            if not exfoliat_height:
                if traversal == False:
                    QtWidgets.QMessageBox.warning(self, 'error', "Can't exfoliate.")
            else:
                for index in range(len(gouzaolist) - 1, -1, -1):
                    if index not in index_lst:
                        after_add_cell_self_ATO.pop(index)
                cell_par = after_add_cell_self_ATO.get_cell_lengths_and_angles()
                if cell_par[3] == 90 and cell_par[4] == 90:
                    cell_par[2] = height + .01
                    after_add_cell_self_ATO.set_cell(cell_par)
                    after_add_cell_self_ATO.translate(np.array([0, 0, -min_z]))
                else:
                    cell_par[2:5] = [height + .01, 90, 90]
                    after_add_cell_self_ATO.set_cell(cell_par)
                    after_add_cell_self_ATO.translate(np.array([0, 0, -min_z]))
                    pos_lst = after_add_cell_self_ATO.get_positions()
                    zhijiao_system = after_add_cell_self_ATO.get_cell()
                    A = zhijiao_system.T
                    new_pos_lst = []
                    for pos in pos_lst:
                        b = pos.T
                        r = np.linalg.solve(A, b)  # 求解线性方程组，直角坐标系下----用晶胞坐标系表示
                        while r[0] < 0:
                            pos += zhijiao_system[0]
                            k = pos.T
                            r = np.linalg.solve(A, k)
                        while r[1] < 0:
                            pos += zhijiao_system[1]
                            k = pos.T
                            r = np.linalg.solve(A, k)
                        new_pos_lst.append(pos)
                    after_add_cell_self_ATO.set_positions(new_pos_lst)
                self.plot(after_add_cell_self_ATO, dictionary=True, clear=True, globalAtomsobject=False)
                self.Atomsobject = None
                Text3column = self.judgeconductivity(after_add_cell_self_ATO)
                childx = QtWidgets.QTreeWidgetItem(parent)
                childx.setText(1, self.dirkey)
                childx.setText(0, 'layer')
                childx.setText(3, Text3column)
        except Exception as e:
            print(e)

    def Stack_main_two(self):
        """ To stack two layer according to the client."""
        try:
            itemlst = [_ for _ in self.project_tree.selectedItems()]
            formula_lst = [_.text(1) for _ in itemlst]
            if len(formula_lst) != 2:
                QtWidgets.QMessageBox.warning(self, 'error',
                                              'Please choose 2 crystal in Projectbox' + '\n' + '(Press Ctrl or Shift).')
            elif itemlst[0].text(0) == 'bulk' or itemlst[1].text(0) == 'bulk':
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose stack or layer(doubleclick to change)')
            else:
                obj1 = self.dic_Formula_Atoms[formula_lst[0]]  # Atoms1，对象
                obj2 = self.dic_Formula_Atoms[formula_lst[1]]  # Atoms2，对象
                self.compare_cell(obj1, obj2)
                try:
                    if self.crys1 == None or self.crys2 == None:
                        QtWidgets.QMessageBox.warning(self, 'error', "Can't stack, for the γ difference is too big. ")
                    else:
                        self.zhexianplot = zhexian_plot_window(self.lista, self.listb, obj1, obj2, self.axis_lst)
                        self.zhexianplot.signal_determine.connect(self.after_zhexian_plot)
                except Exception as e:
                    print(e)
        except Exception as e:
            print(e)

    def Stack_main_one(self):
        try:
            formula_lst = []
            itemlst = []
            for item in self.project_tree.selectedItems():
                formula_lst.append(item.text(1))
                itemlst.append(item)
            if len(formula_lst) != 2:
                QtWidgets.QMessageBox.warning(self, 'error',
                                              'Please choose 2 crystal in Projectbox' + '\n' + '(Press Ctrl or Shift).')
            elif itemlst[0].text(0) == 'bulk' or itemlst[1].text(0) == 'bulk':
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose stack or layer(doubleclick to change)')
            else:
                try:
                    self.stack_one_window.close()
                except:
                     pass
                self.stack_one_window = Stack_main_one_window()
                self.stack_one_window.signal_emit.connect(self.stack_one)
                self.orijincrys1 = deepcopy(self.dic_Formula_Atoms[formula_lst[0]])  # Atoms1，对象
                obj2 = self.dic_Formula_Atoms[formula_lst[1]]  # Atoms2，对象
                crystype1 = CrystalOperationsNum1(self.orijincrys1)
                self.obj1 = deepcopy(crystype1.gamaToAcuteangle())  # 在切的时候统一转成锐角
                crystype2 = CrystalOperationsNum1(obj2)
                self.obj2 = deepcopy(crystype2.gamaToAcuteangle())  # 在切的时候统一转成锐角
        except Exception as e:
            print(e)

    def stack_one(self, Normal_strain_limit, Shear_strain_limit, layer_dis, vacuum_dis, traversal=False):
        try:
            cell1 = self.obj1.get_cell()
            cell2 = self.obj2.get_cell()
            cell1_par = self.obj1.get_cell_lengths_and_angles()
            cell2_par = self.obj2.get_cell_lengths_and_angles()
            self.cell1 = np.array([cell1[0][:2], cell1[1][:2]])
            self.cell2 = np.array([cell2[0][:2], cell2[1][:2]])
            a1 = cell1_par[0]
            a2 = cell2_par[0]
            atimes_lst = self.multigaosifunction_with_start_value([a1, a2], [0, 0], Normal_strain_limit)
            a1_times = atimes_lst[0]
            a2_times = atimes_lst[1]
            b1 = cell1_par[1]
            b2 = cell2_par[1]
            b_vertical1 = b1 * np.sin(cell1_par[5] / 180 * np.pi)
            b_vertical2 = b2 * np.sin(cell2_par[5] / 180 * np.pi)
            btimes_lst = self.multigaosifunction_with_start_value([b_vertical1, b_vertical2], [0, 0], Normal_strain_limit)
            b1_btimes = btimes_lst[0]
            b2_btimes = btimes_lst[1]
            init_b1 = deepcopy(self.cell1[1] * b1_btimes)
            init_b2 = deepcopy(self.cell2[1] * b2_btimes)
            init_a1 = deepcopy(self.cell1[0])
            init_a2 = deepcopy(self.cell2[0])
            number1 = - int(b1 * np.cos(cell1_par[5] / 180 * np.pi) * b1_btimes / a1) - 1
            number2 = - int(b2 * np.cos(cell2_par[5] / 180 * np.pi) * b2_btimes / a2) - 1
            number1, number2 = self.find_a_var(init_b1, init_b2,
                                               init_a1, init_a2, Normal_strain_limit, Shear_strain_limit, number1, number2)

            if (number1 is None or number2 is None) and traversal == False:
                QtWidgets.QMessageBox.warning(self, 'error',
                                              "Don't find.")
            elif (number1 is None or number2 is None) and traversal == True:
                return "fail"
            else:
                a1 = (a1_times, 0, 0)
                b1 = (number1, b1_btimes, 0)
                a2 = (a2_times, 0, 0)
                b2 = (number2, b2_btimes, 0)
                self.Atomsobject = cut(self.obj1, a=a1, b=b1, c=[0, 0, 1], origo=(0, 0, 0))
                crys1_cell = self.Atomsobject.get_cell_lengths_and_angles()
                crys2_new = cut(self.obj2, a=a2, b=b2, c=[0, 0, 1], origo=(0, 0, 0))
                crys2_cell = crys2_new.get_cell_lengths_and_angles()
                crys2_cell[0], crys2_cell[1], crys2_cell[5] = crys1_cell[0], crys1_cell[1], crys1_cell[5]
                crys2_new.set_cell(crys2_cell, scale_atoms=True)
                cell1_new = self.Atomsobject.get_cell()
                crys2_new.translate(np.array([0, 0, layer_dis + cell1_new[2][2]]))
                crys1_cell[2] = crys1_cell[2] + crys2_cell[2] + layer_dis
                self.Atomsobject.extend(crys2_new)
                self.Atomsobject.set_cell(crys1_cell)
                Atomsobject = self.set_vacuum_layer(vacuum_dis, addproject=False)
                self.plot(Atomsobject, plot=False, object=False, clear=False, globalAtomsobject=False, dictionary=True,
                          Hetero_tab=False)
                key = list(self.dic_Formula_Atoms.keys())[
                    list(self.dic_Formula_Atoms.values()).index(self.orijincrys1)]
                Hetero_junction_object = Hetero_junction()
                Hetero_junction_object.Optimal_orientation = [[[a1_times, 0], [number1, b1_btimes]],
                                                              [[a2_times, 0], [number2, b2_btimes]]]
                Hetero_junction_object.Rotate_angle = 0
                Hetero_junction_object.Atomsobject = Atomsobject
                Hetero_junction_object.Name = self.dirkey
                Hetero_junction_object.Normal_strain_limit = Normal_strain_limit
                Hetero_junction_object.Shear_strain_limit = Shear_strain_limit
                Hetero_junction_object.Area_error_limit = np.sqrt(Normal_strain_limit ** 2 +
                                                                  Shear_strain_limit ** 2) * 2
                self.add_hetero_junction(Hetero_junction_object)
                it = QtWidgets.QTreeWidgetItemIterator(self.project_tree)
            if traversal==True:
                while it.value():
                    if it.value().text(1) == key:
                        try:
                            childx = QtWidgets.QTreeWidgetItem(it.value().parent())
                            childx.setText(1, self.dirkey)
                            childx.setText(0, 'stack')
                            if self.it_value1.text(2) != "":
                                text1 = self.it_value1.text(2)
                            else:
                                text1 = self.it_value1.text(1)
                            if self.it_value2.text(2) != "":
                                text2 = self.it_value2.text(2)
                            else:
                                text2 = self.it_value2.text(1)
                            childx.setText(2, text1 + '-' + text2)
                        except Exception as e:
                            print(e)
                        break
                    it += 1
            else:
                try:
                    parent = self.project_tree.selectedItems()[0].parent()
                    stackchild = QtWidgets.QTreeWidgetItem(parent)
                    stackchild.setText(0, "stack")
                    stackchild.setText(1, self.dirkey)
                    if self.project_tree.selectedItems()[0].text(2) != "":
                        text1 = self.project_tree.selectedItems()[0].text(2)
                    else:
                        text1 = self.project_tree.selectedItems()[0].text(1)
                    if self.project_tree.selectedItems()[1].text(2) != "":
                        text2 = self.project_tree.selectedItems()[1].text(2)
                    else:
                        text2 = self.project_tree.selectedItems()[1].text(1)
                    stackchild.setText(2, text1 + '-' + text2)
                except Exception as e:
                    print(e)
                    stackchild.setText(2, self.project_tree.selectedItems()[0].text(2))
        except Exception as e:
            print(e)

    def find_a_var(self, b1, b2, a1, a2, Normal_strain_limit, Shear_strain_limit, number1, number2):
        try:
            len_error = np.sqrt(Normal_strain_limit ** 2 + Shear_strain_limit ** 2)
            while True:
                b1_new = deepcopy(b1 + number1 * a1)
                b2_new = deepcopy(b2 + number2 * a2)
                if np.linalg.norm(b1_new - b2_new) < len_error * np.linalg.norm(b1_new) and \
                        abs(b1_new @ b1_new - b1_new @ b2_new) < \
                        Normal_strain_limit * (b1_new @ b1_new) and \
                        abs(np.arctan(b1_new[1] / b1_new[0]) - np.arctan(b2_new[1] / b2_new[0])) < Shear_strain_limit:
                    return number1, number2
                elif b1_new[0] > b2_new[0]:
                    number2 += 1
                else:
                    number1 += 1
        except Exception as e:
            print(e)


    def after_zhexian_plot(self, num1, num2, length, vacuum_dis):
        """ globalize layerlength"""
        try:
            self.point1_index = num1
            self.point2_index = num2
            self.layerlength = length
            self.vacuum_distance = vacuum_dis
            self.layer_transform()
        except Exception as e:
            print(e)

    def multigaosifunction_with_start_value(self, a_lst, startvalue_lst, error_limit, times_ornot=True):
        """ With Start Value:gaosifunction"""
        a0_orijin = a_lst[0]
        a_init = deepcopy(a_lst)
        timesnum = 0
        while True:
            count_num = 0
            timesnum += 1
            a0 = startvalue_lst[0] + a0_orijin * timesnum
            times_lst = [timesnum]
            length_lst = [a0]
            for i in range(1, len(a_init)):
                timesi = round((a0 - startvalue_lst[i]) / a_init[i])
                if abs(a0 - timesi * a_init[i] - startvalue_lst[i]) > a0 * error_limit:
                    break
                else:
                    count_num += 1
                    times_lst.append(int(timesi))
                    length_lst.append(timesi * a_init[i] + startvalue_lst[i])
            if count_num == len(a_init) - 1:
                if times_ornot == True:
                    return times_lst
                else:
                    return length_lst

    def chase_lst(self, float1, float2):
        lst = [Fraction(1, 1)]
        orijin1 = float1
        orijin2 = float2
        strain = abs(float1 - float2) / float1
        num1 = 1
        num2 = 1
        num = 0
        zongnum = int((strain - .005) / .001)
        while num < zongnum and strain > .005 and num1 < 50 and num2 < 50:
            if float1 > float2:
                float2 += orijin2
                num2 += 1
                strain_new = abs(float1 - float2) / float1
                if strain - strain_new > .001:
                    lst.append(Fraction(num2, num1))
                    strain = strain_new
                    num += 1
                    print(Fraction(num2, num1))
            else:
                float1 += orijin1
                num1 += 1
                strain_new = abs(float1 - float2) / float1
                if strain - strain_new > .001:
                    lst.append(Fraction(num2, num1))
                    strain = strain_new
                    num += 1
                    print(Fraction(num2, num1))
        return lst


    def compare_cell(self, crys1, crys2, traversalstack=False):  # 绘图    # 计算两原胞a，b轴,返回应变-超胞a,b长的列表
        """
        To compare two crys.
        :param crys1: the bottom crys
        :param crys2: the top crys
        :param traversalstack: whether is traversalstack form or not
        :return: lista(b) ,which have strain1, strain2, strain_sum, lcm,
        int1/100(crys1'a-axis), int2/100(cry2'a-axis) information
        """
        self.crys1 = deepcopy(crys1)         # globalize crys1,crys2
        self.crys2 = deepcopy(crys2)
        crys1_cell_lst = crys1.get_cell_lengths_and_angles()
        crys2_cell_lst = crys2.get_cell_lengths_and_angles()
        shear_strain = abs(crys1_cell_lst[5] - crys2_cell_lst[5]) / 180 * np.pi
        if shear_strain < .001745:  # 若两晶体的gama角之差小于0.1°，则认为属于一种晶系
            crys1_a_axis_len,  crys2_a_axis_len, crys1_b_axis_len, crys2_b_axis_len = \
                crys1_cell_lst[0], crys2_cell_lst[0], crys1_cell_lst[1], crys2_cell_lst[1]
            # lista:strain1, strain2, strain_sum, lcm, int1(晶体1的a-axis), int2(晶体2的a-axis)
            if traversalstack == False:
                a_ratio_of_elong = self.chase_lst(crys1_a_axis_len, crys2_a_axis_len)
                b_ratio_of_elong = self.chase_lst(crys1_b_axis_len, crys2_b_axis_len)
                # 用追逐法
                self.lista = self.make_listab(crys1_a_axis_len, crys2_a_axis_len, a_ratio_of_elong)
                self.listb = self.make_listab(crys1_b_axis_len, crys2_b_axis_len, b_ratio_of_elong)
                self.lista.sort(key=itemgetter(2))  # lista以strainsum排序
                self.listb.sort(key=itemgetter(2))
                functiona = [a[2:4] for a in self.lista]
                functionb = [b[2:4] for b in self.listb]
                # functiona = [strainsum, lcma]
                # functionb = [strainsum, lcmb]
                self.axis_lst = [functiona, functionb]
            else:     # traversal_stack
                a_times_lst = self.multigaosifunction_with_start_value([crys1_a_axis_len, crys2_a_axis_len], [0, 0],
                                                                       self.normal_strain_limit)
                self.a1_times = a_times_lst[0]
                self.a2_times = a_times_lst[1]
                b_times_lst = self.multigaosifunction_with_start_value([crys1_b_axis_len, crys2_b_axis_len], [0, 0],
                                                                       self.normal_strain_limit)
                self.b1_times = b_times_lst[0]
                self.b2_times = b_times_lst[1]

                self.traversal_Hetero_junction_object = Hetero_junction()
                self.traversal_Hetero_junction_object.Optimal_orientation = [[[self.a1_times, 0],
                                                                             [0, self.b1_times]],
                                                                             [[self.a2_times, 0],
                                                                             [0, self.b2_times]]]
                self.traversal_Hetero_junction_object.Normal_strain_limit = self.normal_strain_limit
                self.traversal_Hetero_junction_object.Shear_strain_limit = 0
                self.traversal_Hetero_junction_object.Area_error_limit = 2 * self.normal_strain_limit
        else:
            self.crys1 = None        # 如果大于0.5°则没法进行layer_transform
            self.crys2 = None

    def make_listab(self, float1, float2, conti_frac_lst):
        """ To replace the function (least_common_multiple_float) in src.function"""
        lcm_strain_sum_lst = []
        for conti_frac in conti_frac_lst:
            if conti_frac == 0:
                continue
            length1 = float1 * conti_frac.denominator
            length2 = float2 * conti_frac.numerator
            if length1 >= length2:
                lcm = length1
                strain1 = 0
                strain2 = (length1 - length2)/length2
                strain_sum = strain1 + strain2
                crys1_length = float1
                crys2_length = lcm / conti_frac.numerator
            else:
                lcm = length2
                strain2 = 0
                strain1 = (length2 - length1) / length1
                strain_sum = strain1 + strain2
                crys2_length = float2
                crys1_length = lcm / conti_frac.denominator
            lst = [strain1, strain2, strain_sum, lcm, crys1_length, crys2_length]
            # 晶体1的应变, 晶体2的应变 晶体1， 2的应变和 最小公倍数 变形后晶体1长度 变形后晶体2长度
            lcm_strain_sum_lst.append(lst)
        return lcm_strain_sum_lst

    def layer_transform(self, traversalstack=False):
        """ To change two crys' a,b axis."""
        try:
            print("in layer_transform")
            if self.crys1 is not None and self.crys2 is not None:
                if traversalstack == False:
                    print(self.crys1, self.crys2)
                    orijin1 = deepcopy(self.crys1)  # 不改变crys1的晶胞状态
                    orijin2 = deepcopy(self.crys2)
                    crys1_cell_lst = orijin1.get_cell_lengths_and_angles().tolist()  # 返回存有[应变， 最小公倍数]的列表
                    crys2_cell_lst = orijin2.get_cell_lengths_and_angles().tolist()
                    if self.point1_index is None or self.point2_index is None:
                        ...
                    else:
                        obj_strain_lcma = self.axis_lst[0][self.point1_index]  # client选择的strain-lcma列表
                        obj_strain_lcmb = self.axis_lst[1][self.point2_index]  # client选择的strain-lcmb列表
                        crys1_cell_lst[0] = self.lista[self.point1_index][4]  # 晶体1的a-axis
                        crys1_cell_lst[1] = self.listb[self.point2_index][4]  # 晶体1的b-axis
                        crys2_cell_lst[0] = self.lista[self.point1_index][5]  # 晶体2的a-axis
                        crys2_cell_lst[1] = self.listb[self.point2_index][5]  # 晶体2的b-axis
                        orijin1.set_cell(crys1_cell_lst)
                        orijin2.set_cell(crys2_cell_lst)  # 按照client需求改变晶胞大小,里面的原子做仿射变换
                        multi_a_crys1 = round(obj_strain_lcma[1] / crys1_cell_lst[0])
                        multi_a_crys2 = round(obj_strain_lcma[1] / crys2_cell_lst[0])
                        multi_b_crys1 = round(obj_strain_lcmb[1] / crys1_cell_lst[1])
                        multi_b_crys2 = round(obj_strain_lcmb[1] / crys2_cell_lst[1])
                        self.Hetero_junction_object = Hetero_junction()
                        self.Hetero_junction_object.Rotate_angle = 0
                        self.Hetero_junction_object.Normal_strain_limit = max(obj_strain_lcma[0], obj_strain_lcmb[0])
                        self.Hetero_junction_object.Optimal_orientation = [[[multi_a_crys1, 0], [0, multi_b_crys1]],
                                                                           [[multi_a_crys2, 0], [0, multi_b_crys2]]]
                        self.Hetero_junction_object.Shear_strain_limit = 0
                        self.Hetero_junction_object.Area_error_limit = 2 * max(obj_strain_lcma[0] , obj_strain_lcmb[0])
                        self.crys1_supercell = orijin1.repeat((multi_a_crys1, multi_b_crys1, 1))
                        self.crys2_supercell = orijin2.repeat((multi_a_crys2, multi_b_crys2, 1))

                        self.point1_index = None
                        self.point2_index = None
                        self.stack(self.layerlength, 0)
                else:
                    self.crys1_supercell = self.crys1.repeat((self.a1_times, self.b1_times, 1))
                    self.crys2_supercell = self.crys2.repeat((self.a2_times, self.b2_times, 1))
                    self.stack(self.layerlength, 1)
        except Exception as e:
            print(e)

    def stack(self, layerlength, traversalstack):  # 把crys2堆垛在crys1上,返回成crys1
        """ To stack crys2 on crys1 with layerlength"""
        try:
            print("in stack")
            self.supercell1_cell = self.crys1_supercell.get_cell_lengths_and_angles()
            self.supercell2_cell = self.crys2_supercell.get_cell_lengths_and_angles()
            self.supercell2_cell[0], self.supercell2_cell[1] = self.supercell1_cell[0], self.supercell1_cell[1]
            self.crys2_supercell.set_cell(self.supercell2_cell, scale_atoms=True)
            supercell1_cell_ndarray = self.crys1.get_cell()
            self.crys2_supercell.translate(np.array([0, 0, layerlength + supercell1_cell_ndarray[2][2]]))
            self.supercell1_cell[2] = self.supercell1_cell[2] + self.supercell2_cell[2] + layerlength + .01
            if traversalstack == 0:  # not - traversalstack
                try:
                    self.crys1_supercell.extend(self.crys2_supercell)
                    self.crys1_supercell.set_cell(self.supercell1_cell)
                    self.Atomsobject = deepcopy(self.crys1_supercell)
                    self.crys1_supercell = self.set_vacuum_layer(self.vacuum_distance, addproject=False)
                except Exception as e:
                    print(e)
                try:
                    self.plot(self.crys1_supercell, clear=True, globalAtomsobject=True, dictionary=True)
                    self.Hetero_junction_object.Atomsobject = deepcopy(self.crys1_supercell)
                    self.Hetero_junction_object.Name = self.dirkey
                    self.add_hetero_junction(self.Hetero_junction_object)
                    parent = self.project_tree.selectedItems()[0].parent()
                    stackchild = QtWidgets.QTreeWidgetItem(parent)
                    stackchild.setText(0, "stack")
                    stackchild.setText(1, self.dirkey)
                    if self.project_tree.selectedItems()[0].text(2) != "":
                        text1 = self.project_tree.selectedItems()[0].text(2)
                    else:
                        text1 = self.project_tree.selectedItems()[0].text(1)
                    if self.project_tree.selectedItems()[1].text(2) != "":
                        text2 = self.project_tree.selectedItems()[1].text(2)
                    else:
                        text2 = self.project_tree.selectedItems()[1].text(1)
                    stackchild.setText(2, text1 + '-' + text2)
                except Exception as e:
                    print(e)
            else:  # 批量处理，用户没法自定义
                self.crys1_supercell.extend(self.crys2_supercell)
                self.crys1_supercell.set_cell(self.supercell1_cell)
                self.Atomsobject = deepcopy(self.crys1_supercell)
                self.crys1_supercell = self.set_vacuum_layer(self.vacuum_distance, addproject=False)
                self.plot(self.crys1_supercell, plot=False, object=False, clear=False, globalAtomsobject=False,
                          dictionary=True, Hetero_tab=False)
                self.traversal_Hetero_junction_object.Atomsobject = deepcopy(self.Atomsobject)
                self.traversal_Hetero_junction_object.Name = self.dirkey
                self.traversal_Hetero_junction_object.Rotate_angle = 0
                self.add_hetero_junction(self.traversal_Hetero_junction_object)
                key = list(self.dic_Formula_Atoms.keys())[list(self.dic_Formula_Atoms.values()).index(self.orijincrys1)]
                it = QtWidgets.QTreeWidgetItemIterator(self.project_tree)
                while it.value():
                    if it.value().text(1) == key:
                        try:
                            childx = QtWidgets.QTreeWidgetItem(it.value().parent())
                            childx.setText(1, self.dirkey)
                            childx.setText(0, 'stack')
                            if self.it_value1.text(2) != "":
                                text1 = self.it_value1.text(2)
                            else:
                                text1 = self.it_value1.text(1)
                            if self.it_value2.text(2) != "":
                                text2 = self.it_value2.text(2)
                            else:
                                text2 = self.it_value2.text(1)
                            childx.setText(2, text1 + '-' + text2)
                        except Exception as e:
                            print(e)
                        break
                    it += 1
        except Exception as e:
            print(e)

    def Export(self):
        """ To export cif files as a dictionary."""
        try:
            output_file = QtWidgets.QFileDialog.getSaveFileName(self, "Save as cif",
                                                      "untitled.dic")
            if output_file[0] == "":
                ...
            else:
                folder = os.path.exists(output_file[0])
                if not folder:  # 判断是否存在文件夹如果不存在则创建为文件夹	,建文件时如果路径不存在会创建这个路径
                    os.makedirs(output_file[0])
                it = QtWidgets.QTreeWidgetItemIterator(self.project_tree)
                while it.value():
                    if it.value().text(1) == '':  # 根目录
                        try:
                            root = output_file[0] + '/{}'.format(it.value().text(0))  # 根目录地址
                            os.makedirs(root)
                        except Exception as e:
                            print(e)
                    else:
                        if it.value().text(1) in self.dic_Hetero_junction.keys():
                            self.export_to_cif(self.dic_Formula_Atoms[it.value().text(1)], root, it.value().text(1),
                                               Hetero_junction=self.dic_Hetero_junction[it.value().text(1)])
                        else:
                            self.export_to_cif(self.dic_Formula_Atoms[it.value().text(1)], root, it.value().text(1))

                    it += 1
        except Exception as e:
            print(e)

    def export_to_cif(self, crys, dir, name, Hetero_junction=False):  # 输出crys.cif文件
        """ To export cif files as a dictionary."""
        f = AseAtomsAdaptor.get_structure(crys)
        f1 = str(f) + '\n'
        j_pv_lst = findall('abc\s\s\s:(.*?)\n', f1)[0]  # abc   :  19.257300  19.569178  21.133988
        j1_pv_lst = j_pv_lst.split(' ')  # abc   :   6.419100   6.523059   7.044663
        while '' in j1_pv_lst:
            j1_pv_lst.remove('')
        par_lst_matrix = crys.get_cell()
        # 物质种类（比如：Lu2 Al4）
        y1 = findall('Full\sFormula\s(.*?)\n', f1)[0]
        y1.lstrip('(').rstrip(')')
        material = findall('Full\sFormula\s(.*?)\n', f1)[0].lstrip('(').rstrip(')')
        elements = material.split(' ')
        zmb_lst = [chr(i) for i in range(97, 123)] + [chr(i) for i in range(65, 91)]
        szb_lst = [str(i) for i in range(0, 10)]
        element_lst = []
        number_lst = []
        for element in elements:  # 这个循环之后，得到原子与对应的原子数两个列表
            letter_lst = list(element)
            symbol_lst = []
            num_lst = []
            element_lst1 = []
            number_lst1 = []
            for letter in letter_lst:
                if letter in szb_lst:
                    num_lst.append(letter)
                if letter in zmb_lst:
                    symbol_lst.append(letter)
                element1 = ''.join(symbol_lst)
                element_lst1.append(element1)
                number1 = ''.join(num_lst)
                number_lst1.append(number1)
            ys = 'a'  # 元素
            gs = '0'  # 个数
            for i in element_lst1:
                if len(i) >= len(ys):
                    ys = i
            element_lst.append(i)
            for i in number_lst1:
                if len(i) >= len(gs):
                    gs = i
            number_lst.append(i)
        par_lst_species = []  # 用于Cifwrite参数(species)的
        for i, element in enumerate(element_lst):
            num = int(number_lst[i])
            for j in range(num):
                par_lst_species.append(element)
        # 每个原子的坐标
        ord_lst = []  # 最终Cifwriter所需要的coords参数
        ord_lst2 = []  # 储存的形式为
        for element in element_lst:
            ord_lst1 = findall(element + '\s\s\s\s(.*?)\n', f1)
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
        slab.write_file(dir + "\{}.cif".format(name))


    def Import(self):  # 两层或一层文件夹
        """ To import dictionary with cif files."""
        try:
            directory = QtGui.QFileDialog.getExistingDirectory(self, 'Select directory')
            names = []
            dir_lst = []
            for dirpath, dirs, files in os.walk(directory):  # 递归遍历当前目录和所有子目录的文件和目录
                for name in files:  # files保存的是所有的文件名
                    if os.path.splitext(name)[1] in ['.cif', '.vasp']:
                        file_path = os.path.join(dirpath, name)  # 加上路径，dirpath是遍历时文件对应的路径
                        names.append(name)
                        dir_lst.append(file_path)
            for i, name in enumerate(names):
                try:
                    self.Atomsobject = deepcopy(read(dir_lst[i]))
                    self.plot(self.Atomsobject, plot=False, object=False, clear=False, dictionary=True,
                              globalAtomsobject=False, Hetero_tab=False)
                    root = QtWidgets.QTreeWidgetItem(self.project_tree)
                    root.setText(0, name)
                    child = QtWidgets.QTreeWidgetItem(root)
                    child.setText(0, "bulk")
                    child.setText(1, self.dirkey)
                except Exception as e:
                    print(e)
            self.project_tree.expandAll()
        except Exception as e:
            print(e)

    def judgeconductivity(self, Atomsobject):
        """ To judge the conductivity of an Atoms object."""
        try:
            Numberlist = Atomsobject.get_atomic_numbers()
            Numberset_lst = list(set(Numberlist))
            Elementlist = [crys_data.Element_formula[number] for number in Numberset_lst]
            num = 0
            conductlist = [crys_data.conductor_Elemental_composition, crys_data.semiconductor_Element_composition,
                           crys_data.Insulator_Element_composition]
            for i in range(3):
                for conduct_element in conductlist[i]:
                    if len(conduct_element) == len(Elementlist):
                        for element in Elementlist:  # Atomsobject 的元素
                            if element not in conduct_element:
                                break
                            if element == Elementlist[-1]:
                                if i == 0:
                                    num = 1
                                if i == 1:
                                    num = 2
                                if i == 2:
                                    num = 3
            if num == 1:
                return "Conductor"
            elif num == 2:
                return "Semiconductor"
            elif num == 3:
                return "Insulator"
            elif num == 0:
                return ""
        except Exception as e:
            print(e)

    def multi_layer_optimize(self):
        try:
            try:
                num = 1
                formula_lst = []
                itemlst = []
                for item in self.project_tree.selectedItems():
                    formula_lst.append(item.text(1))
                    itemlst.append(item)
                for item in itemlst:
                    if item.text(0) == 'bulk':
                        num = 0
                        break
                if len(formula_lst) <= 2:
                    QtWidgets.QMessageBox.warning(self, 'error',
                                                  'Please choose at least 3 crystal in Projectbox' + '\n' +
                                                  '(Press Ctrl or Shift).')
                elif num == 0:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose stack or layer(doubleclick to change)')
                else:
                    try:
                        self.multi_layer_window.close()
                    except:
                        pass
                    obj_lst = []
                    for formula in formula_lst:
                        crystype = CrystalOperationsNum1(self.dic_Formula_Atoms[formula])
                        Atomsobject = crystype.gamaToAcuteangle()   # toacute angle
                        obj_lst.append(Atomsobject)
                    self.multi_layer_window = multi_layer_optimize_window(obj_lst)
                    self.multi_layer_window.signal_emit_infomation.connect(self.after_multi_layer_window)
            except Exception as e:
                print(e)
        except Exception as e:
            print(e)

    def after_multi_layer_window(self, num, obj_lst, layer_dis_lst, vacuum_dis, Area_error_limit,
                                 Normal_strain_limit, Shear_strain_limit, max_gama, min_gama, ab_relative_error, max_len,
                                 traversal_stack=False):
        try:
            if num == 2:      # 2 orientation already know
                init_angle = obj_lst[0].get_cell_lengths_and_angles()[5]
                jingxi = True
                for obj_index in range(1, len(obj_lst)):
                    angle = obj_lst[obj_index].get_cell_lengths_and_angles()[5]
                    if abs(angle - init_angle) > .001:
                        jingxi = False
                        break
                if jingxi == False and traversal_stack == False:
                    QtWidgets.QMessageBox.warning(self, 'error',
                                                  'Different gama!')
                else:
                    a_lst = [_.get_cell_lengths_and_angles()[0] for _ in obj_lst]
                    b_lst = [_.get_cell_lengths_and_angles()[1] for _ in obj_lst]
                    startvalue_lst = [0 for _ in range(len(a_lst))]
                    new_times_a = self.multigaosifunction_with_start_value(a_lst, startvalue_lst, Normal_strain_limit)
                    new_times_b = self.multigaosifunction_with_start_value(b_lst, startvalue_lst, Normal_strain_limit)
                    for obj_index in range(len(obj_lst)):
                        obj_lst[obj_index] *= (new_times_a[obj_index], new_times_b[obj_index], 1)
                    first_layer_cell = obj_lst[0].get_cell_lengths_and_angles()
                    for obj_index in range(1, len(obj_lst)):
                        cell = obj_lst[obj_index].get_cell_lengths_and_angles()
                        cell[0], cell[1], cell[5] = first_layer_cell[0], first_layer_cell[1], first_layer_cell[5]
                        obj_lst[obj_index].set_cell(cell, scale_atoms=True)
                    supercell_init = deepcopy(obj_lst[0])
                    for obj_index in range(1, len(obj_lst)):
                        cell1_new = supercell_init.get_cell()
                        layer_up = deepcopy(obj_lst[obj_index])
                        layer_up.translate(np.array([0, 0, layer_dis_lst[obj_index - 1] + cell1_new[2][2]]))
                        supercell_init.extend(layer_up)
                        super_cell_init = supercell_init.get_cell_lengths_and_angles()
                        super_cell_init[2] += layer_up.get_cell_lengths_and_angles()[2] + layer_dis_lst[obj_index - 1]
                        supercell_init.set_cell(super_cell_init)
                    self.Atomsobject = deepcopy(supercell_init)
                    supercell_init = self.set_vacuum_layer(vacuum_dis, addproject=False)
                    self.plot(supercell_init, clear=True, globalAtomsobject=False, dictionary=True)
                    Hetero_junction_object = Hetero_junction()
                    Hetero_junction_object.Name = self.dirkey
                    Hetero_junction_object.Atomsobject = supercell_init
                    Hetero_junction_object.Normal_strain_limit = Normal_strain_limit
                    Hetero_junction_object.Shear_strain_limit = 0
                    Hetero_junction_object.Area_error_limit = 2 * Normal_strain_limit
                    Hetero_junction_object.Rotate_angle = 0
                    self.add_hetero_junction(Hetero_junction_object)
                    try:
                        parent = self.project_tree.selectedItems()[0].parent()
                        stackchild = QtWidgets.QTreeWidgetItem(parent)
                        stackchild.setText(0, "stack")
                        stackchild.setText(1, self.dirkey)
                        text = self.project_tree.selectedItems()[0].text(1)
                        for item_index in range(1, len(self.project_tree.selectedItems())):
                            text += '-' + self.project_tree.selectedItems()[item_index].text(1)
                        stackchild.setText(2, text)
                    except Exception as e:
                        print(e)
            elif num == 1:   # a1, a2 already know
                a_lst = [obj_piece.get_cell_lengths_and_angles()[0] for obj_piece in obj_lst]
                b_timesin_gama = [_.get_cell_lengths_and_angles()[1] * \
                               np.sin(_.get_cell_lengths_and_angles()[5] / 180 * np.pi) for _ in obj_lst]
                startvalue_lst = [0 for _ in range(len(a_lst))]
                a_length = self.multigaosifunction_with_start_value(a_lst, startvalue_lst, Normal_strain_limit,
                                                                    times_ornot=False)
                b_blength = self.multigaosifunction_with_start_value(b_timesin_gama, startvalue_lst,
                                                                     Normal_strain_limit, times_ornot=False)
                inita = []
                for obj in obj_lst:
                    cell = obj.get_cell()
                    inita.append(cell[1][0] - cell[0][0])
                b_a_length = self.multi_vector_gaosifunction(a_lst, inita, b_blength, Normal_strain_limit, Shear_strain_limit)
                # b_a_length = self.multigaosifunction_with_start_value(a_lst, inita, Normal_strain_limit, times_ornot=False)
                layer_lst = []
                optimal_match = []
                for obj_index, obj in enumerate(obj_lst):
                    # obj = obj_lst[obj_index]
                    a_vector = np.array([a_length[obj_index], 0, 0])
                    b_vector = np.array([b_a_length[obj_index], b_blength[obj_index], 0])
                    A = obj.get_cell().T
                    b1 = a_vector.T
                    b2 = b_vector.T
                    r1 = np.linalg.solve(A, b1)
                    r2 = np.linalg.solve(A, b2)
                    k1, k2 = list(r1), list(r2)
                    optimal_match.append([[round(k1[0]), round(k1[1])], [round(k2[0]), round(k2[1])]])
                    layer = cut(obj, a=r1,
                                b=r2,
                                c=(0, 0, 1), origo=(0, 0, 0))
                    layer_lst.append(layer)
                supercell_init = deepcopy(layer_lst[0])
                layer_down_cell = supercell_init.get_cell_lengths_and_angles()
                for obj_index in range(1, len(layer_lst)):
                    cell1_new = supercell_init.get_cell()
                    layer_up = deepcopy(layer_lst[obj_index])
                    layer_up_cell = layer_up.get_cell_lengths_and_angles()
                    layer_up_cell[0], layer_up_cell[1], layer_up_cell[5] = layer_down_cell[0], \
                                                                           layer_down_cell[1],layer_down_cell[5]
                    layer_up.set_cell(layer_up_cell, scale_atoms=True)
                    layer_up.translate(np.array([0, 0, layer_dis_lst[obj_index - 1] + cell1_new[2][2]]))
                    supercell_init.extend(layer_up)
                    super_cell_init = supercell_init.get_cell_lengths_and_angles()
                    super_cell_init[2] += layer_up.get_cell_lengths_and_angles()[2] + layer_dis_lst[obj_index - 1]
                    supercell_init.set_cell(super_cell_init)
                self.Atomsobject = deepcopy(supercell_init)
                supercell_init = self.set_vacuum_layer(vacuum_dis, addproject=False)

                self.plot(supercell_init, clear=True, globalAtomsobject=False, dictionary=True)
                Hetero_junction_object = Hetero_junction()
                Hetero_junction_object.Name = self.dirkey
                Hetero_junction_object.Atomsobject = supercell_init
                Hetero_junction_object.Optimal_orientation = optimal_match
                Hetero_junction_object.Normal_strain_limit = Normal_strain_limit
                Hetero_junction_object.Shear_strain_limit = Shear_strain_limit
                Hetero_junction_object.Area_error_limit = np.sqrt(Normal_strain_limit ** 2 + Shear_strain_limit ** 2) * 2
                Hetero_junction_object.Rotate_angle = 0
                self.add_hetero_junction(Hetero_junction_object)
                try:
                    parent = self.project_tree.selectedItems()[0].parent()
                    stackchild = QtWidgets.QTreeWidgetItem(parent)
                    stackchild.setText(0, "stack")
                    stackchild.setText(1, self.dirkey)
                    text = self.project_tree.selectedItems()[0].text(1)
                    for item_index in range(1, len(self.project_tree.selectedItems())):
                        text += '-' + self.project_tree.selectedItems()[item_index].text(1)
                    stackchild.setText(2, text)
                except Exception as e:
                    print(e)
            else:        # orientation is unknown
                vector_error_limit = np.sqrt(Normal_strain_limit ** 2 + Shear_strain_limit ** 2)
                angle_lst = []
                for obj_index in range(1, len(obj_lst)):
                    cell_par = obj_lst[obj_index].get_cell_lengths_and_angles()
                    if abs(cell_par[0] - cell_par[1]) < 10e-6:
                        if abs(cell_par[5] - 60) < 1e-6 or abs(cell_par[5] - 120) < 1e-6:
                            division = 6
                        elif abs(cell_par[5] - 90) < 1e-6:
                            division = 4
                        else:
                            division = 1
                    else:
                        division = 1
                    angle_lst.append(180 / division)
                # create 1 layer lattice
                cell_1_layer = obj_lst[0].get_cell()
                cell_1_vectora = np.array([cell_1_layer[0][0], cell_1_layer[0][1]])
                cell_1_vectorb = np.array([cell_1_layer[1][0], cell_1_layer[1][1]])
                layer1_pos = self.layer1_pos(cell_1_vectora, cell_1_vectorb, len_max=max_len)

                twoD_cell_iter = ([obj.get_cell()[0][0:2], obj.get_cell()[1][0:2]] for obj in obj_lst[1:])
                each_layer_peidui_lst = []
                for cell_index, twoD_cell in enumerate(twoD_cell_iter):
                    angle = int(angle_lst[cell_index]) / 180 * np.pi
                    fuhe_yaoqiu_lst = []
                    theta = 0
                    while theta < angle:
                        rotate_angle = theta
                        rotate_matrix = np.array(
                            [[np.cos(rotate_angle), -np.sin(rotate_angle)],
                             [np.sin(rotate_angle), np.cos(rotate_angle)]])
                        new_cell = twoD_cell @ rotate_matrix
                        inv_matrix = np.linalg.inv(new_cell)
                        peidui_lst = []
                        for pos1 in layer1_pos:
                            r = pos1 @ inv_matrix
                            r1 = np.floor(r[0])
                            r2 = np.floor(r[1])
                            pos2 = r1 * new_cell[0] + \
                                   r2 * new_cell[1]
                            if np.linalg.norm(pos1 - pos2) < vector_error_limit * np.linalg.norm(pos1):
                                if abs(pos1 @ pos1 - pos1 @ pos2) < Normal_strain_limit * (pos1 @ pos1):
                                    peidui_lst.append([pos1, pos2])
                            pos2 = (r1 + 1) * new_cell[0] + \
                                   r2 * new_cell[1]
                            if np.linalg.norm(pos1 - pos2) < vector_error_limit * np.linalg.norm(pos1):
                                if abs(pos1 @ pos1 - pos1 @ pos2) < Normal_strain_limit * (pos1 @ pos1):
                                    peidui_lst.append([pos1, pos2])
                            pos2 = r1 * new_cell[0] + \
                                   (r2 + 1) * new_cell[1]
                            if np.linalg.norm(pos1 - pos2) < vector_error_limit * np.linalg.norm(pos1):
                                if abs(pos1 @ pos1 - pos1 @ pos2) < Normal_strain_limit * (pos1 @ pos1):
                                    peidui_lst.append([pos1, pos2])
                            pos2 = (r1 + 1) * new_cell[0] + \
                                   (r2 + 1) * new_cell[1]
                            if np.linalg.norm(pos1 - pos2) < vector_error_limit * np.linalg.norm(pos1):
                                if abs(pos1 @ pos1 - pos1 @ pos2) < Normal_strain_limit * (pos1 @ pos1):
                                    peidui_lst.append([pos1, pos2])
                        if len(peidui_lst) < 2:
                            break
                        else:
                            for i in range(len(peidui_lst)):
                                for j in range(i + 1, len(peidui_lst)):
                                    if abs(np.linalg.norm(peidui_lst[i][0]) - np.linalg.norm(peidui_lst[j][0])) < \
                                            ab_relative_error * np.linalg.norm(peidui_lst[i][0]):

                                        n_vector1 = np.cross(peidui_lst[i][0], peidui_lst[j][0])
                                        n_vector2 = np.cross(peidui_lst[i][1], peidui_lst[j][1])
                                        square1 = abs(n_vector1)
                                        square2 = abs(n_vector2)
                                        sin_gama = square1 / np.linalg.norm(peidui_lst[i][0]) / np.linalg.norm(
                                            peidui_lst[j][0])
                                        if abs(square1 - square2) < Area_error_limit * square1 and \
                                                max_gama > sin_gama > min_gama:
                                            shear_strain = abs(np.arccos((peidui_lst[i][0] @ peidui_lst[j][0]) /
                                                               np.linalg.norm(peidui_lst[i][0]) / np.linalg.norm(
                                                peidui_lst[j][0])) - np.arccos((peidui_lst[i][1] @ peidui_lst[j][1])
                                            / np.linalg.norm(peidui_lst[i][1]) / np.linalg.norm(peidui_lst[j][1])))
                                            if shear_strain < Shear_strain_limit:
                                                if n_vector1 > 0:  # 保证右手系
                                                    a1 = peidui_lst[i][0]
                                                    b1 = peidui_lst[j][0]
                                                    a2 = peidui_lst[i][1]
                                                    b2 = peidui_lst[j][1]
                                                else:
                                                    b1 = peidui_lst[i][0]
                                                    a1 = peidui_lst[j][0]
                                                    b2 = peidui_lst[i][1]
                                                    a2 = peidui_lst[j][1]
                                                fuhe_yaoqiu_lst.append([a1, b1, rotate_angle, a2, b2])
                        theta += vector_error_limit
                    if len(fuhe_yaoqiu_lst) != 0:
                        each_layer_peidui_lst.append(fuhe_yaoqiu_lst)
                    elif traversal_stack == False:
                        QtWidgets.QMessageBox.warning(self, 'error',
                                                      "Can't find!")
                    else:
                        text = obj_lst[0].get_chemical_formula(mode='hill')
                        for obj_index in range(1, len(obj_lst)):
                            text = text + '-' + obj_lst[obj_index].get_chemical_formula(mode='hill')
                        return 'fail', text

                if len(each_layer_peidui_lst) == len(obj_lst) - 1:
                    optimal_choice = self.solve_each_layer_peidui_lst(each_layer_peidui_lst)
                    if optimal_choice is not None:
                        init_cell = np.array([cell_1_vectora, cell_1_vectorb]).T
                        r1 = np.linalg.solve(init_cell, optimal_choice[0][0].T)
                        r2 = np.linalg.solve(init_cell, optimal_choice[0][1].T)
                        k1 = list(map(lambda x: round(x), r1))
                        k2 = list(map(lambda x: round(x), r2))
                        optimal_lst = [[k1, k2]]
                        r1 = np.append(r1, 0)
                        r2 = np.append(r2, 0)
                        layer_init = cut(obj_lst[0], a=r1, b=r2, c=(0, 0, 1), origo=(0, 0, 0))
                        layer_init = deepcopy(self.deal_with_rotate(layer_init))
                        layer_lst = [layer_init]
                        cell_down = layer_init.get_cell_lengths_and_angles()
                        rotate_lst = [0]
                        for i in range(len(optimal_choice)):
                            obj = obj_lst[i + 1]
                            cell2_before_rotate = obj.get_cell()
                            cell2_vectora, cell2_vectorb = cell2_before_rotate[0][:2], \
                                               cell2_before_rotate[1][:2]
                            degree = optimal_choice[i][2]
                            # 顺时针旋转
                            rotate_matrix = np.array(
                                [[np.cos(degree), np.sin(degree)], [-np.sin(degree), np.cos(degree)]])
                            new_a = optimal_choice[i][3] @ rotate_matrix
                            new_b = optimal_choice[i][4] @ rotate_matrix
                            rotate_lst.append(degree / np.pi * 180)
                            A = np.array([cell2_vectora, cell2_vectorb]).T
                            r2a = np.linalg.solve(A, new_a.T)
                            r2b = np.linalg.solve(A, new_b.T)  # 求解线性方程组，简单,直角坐标系下----用晶胞坐标系表示
                            k1 = list(map(lambda x: round(x), r2a))
                            k2 = list(map(lambda x: round(x), r2b))
                            optimal_lst.append([k1, k2])
                            r2a = np.append(r2a, 0)
                            r2b = np.append(r2b, 0)
                            layer_up = cut(obj, a=r2a, b=r2b, c=[0, 0, 1], origo=[0, 0, 0])
                            k = deepcopy(self.deal_with_rotate(layer_up))
                            cell_up = k.get_cell_lengths_and_angles()
                            cell_up[0], cell_up[1], cell_up[5] = cell_down[0], cell_down[1], cell_down[5]
                            k.set_cell(cell_up, scale_atoms=True)
                            layer_lst.append(k)
                        supercell_init = deepcopy(layer_lst[0])

                        for obj_index in range(1, len(layer_lst)):
                            cell1_new = supercell_init.get_cell()
                            layer_up = deepcopy(layer_lst[obj_index])
                            layer_up.translate(np.array([0, 0, layer_dis_lst[obj_index - 1] + cell1_new[2][2]]))
                            supercell_init.extend(layer_up)
                            super_cell_init = supercell_init.get_cell_lengths_and_angles()
                            super_cell_init[2] += layer_up.get_cell_lengths_and_angles()[2] + layer_dis_lst[obj_index - 1]
                            supercell_init.set_cell(super_cell_init)
                        self.Atomsobject = deepcopy(supercell_init)
                        supercell_init = self.set_vacuum_layer(vacuum_dis, addproject=False)

                        if traversal_stack == False:
                            self.plot(supercell_init, clear=True, globalAtomsobject=False, dictionary=True)
                        else:
                            self.plot(supercell_init, plot=False, globalAtomsobject=False, dictionary=True,
                                      object=False, Hetero_tab=False)
                        Hetero_junction_object = Hetero_junction()
                        Hetero_junction_object.Name = self.dirkey
                        Hetero_junction_object.Atomsobject = supercell_init
                        Hetero_junction_object.Normal_strain_limit = Normal_strain_limit
                        Hetero_junction_object.Shear_strain_limit = Shear_strain_limit
                        Hetero_junction_object.Area_error_limit = Area_error_limit
                        Hetero_junction_object.Rotate_angle = rotate_lst
                        Hetero_junction_object.Optimal_orientation = optimal_lst
                        self.add_hetero_junction(Hetero_junction_object)
                        if traversal_stack == False:
                            try:
                                parent = self.project_tree.selectedItems()[0].parent()
                                stackchild = QtWidgets.QTreeWidgetItem(parent)
                                stackchild.setText(0, "stack")
                                stackchild.setText(1, self.dirkey)
                                text = self.project_tree.selectedItems()[0].text(1)
                                for item_index in range(1, len(self.project_tree.selectedItems())):
                                    text += '-' + self.project_tree.selectedItems()[item_index].text(1)
                                stackchild.setText(2, text)
                            except Exception as e:
                                print(e)
                        else:
                            text = obj_lst[0].get_chemical_formula(mode='hill')
                            for obj_index in range(1, len(obj_lst)):
                                text = text + '-' + obj_lst[obj_index].get_chemical_formula(mode='hill')
                            try:
                                child = QtWidgets.QTreeWidgetItem(self.device_root)
                                child.setText(0, "Stack")
                                child.setText(1, self.dirkey)
                                child.setText(2, text)
                            except Exception as e:
                                print(e)
                            return 'successful', text
                    elif traversal_stack == True:
                        text = obj_lst[0].get_chemical_formula(mode='hill')
                        for obj_index in range(1, len(obj_lst)):
                            text = text + '-' + obj_lst[obj_index].get_chemical_formula(mode='hill')
                        return 'fail', text
                elif traversal_stack == True:
                    text = obj_lst[0].get_chemical_formula(mode='hill')
                    for obj_index in range(1, len(obj_lst)):
                        text = text + '-' + obj_lst[obj_index].get_chemical_formula(mode='hill')
                    return 'fail', text
        except Exception as e:
            print(e)

    def multi_vector_gaosifunction(self, a_lst, startvalue_lst, b_blength, Normal_strain_limit, Shear_strain_limit):
        """ To solve multi layer stack(num = 1) b_alength"""
        try:
            error_limit = np.sqrt(Normal_strain_limit ** 2 + Shear_strain_limit ** 2)
            a0_orijin = a_lst[0]
            a_init = deepcopy(a_lst)
            timesnum = 0
            while True:
                count_num = 0
                timesnum += 1
                a0 = startvalue_lst[0] + a0_orijin * timesnum
                pos1 = np.array([a0, b_blength[0]])
                times_lst = [timesnum]
                length_lst = [a0]
                for i in range(1, len(a_init)):
                    timesi = round((a0 - startvalue_lst[i]) / a_init[i])
                    pos2 = np.array([a_init[i] * timesi + startvalue_lst[i], b_blength[i]])

                    if abs(a0 - timesi * a_init[i] - startvalue_lst[i]) > a0 * error_limit:
                        break
                    elif abs(pos1 @ pos1 - pos1 @ pos2) > Normal_strain_limit * (pos1 @ pos1):
                        break
                    elif abs(np.arctan(pos2[1] / pos2[0]) - np.arctan(pos1[1] / pos1[0])) > Shear_strain_limit:
                        break
                    else:
                        count_num += 1
                        times_lst.append(int(timesi))
                        length_lst.append(timesi * a_init[i] + startvalue_lst[i])
                if count_num == len(a_init) - 1:
                    return length_lst

        except Exception as e:
            print(e)
    def solve_each_layer_peidui_lst(self, each_layer_peidui_lst):
        try:
            each_layer_ab = []
            for each_layer in each_layer_peidui_lst:
                layer_lst = [[list(each_layer_piece[0]), list(each_layer_piece[1])] for each_layer_piece in each_layer]
                each_layer_ab.append(layer_lst)
            chongdie_ab_lst = []
            for ab in each_layer_ab[0]:
                num = 0
                for j in range(1, len(each_layer_ab)):
                    if ab in each_layer_ab[j]:
                        num += 1
                    else:
                        break
                if num == len(each_layer_ab) - 1:
                    chongdie_ab_lst.append(ab)
            if chongdie_ab_lst:
                minsquare = 10000000
                optimal_ab = ...
                for chongdie_ab in chongdie_ab_lst:
                    square = abs(np.cross(chongdie_ab[0], chongdie_ab[1]))
                    if square < minsquare:
                        minsquare = square
                        optimal_ab = chongdie_ab
                index_lst = [_.index(optimal_ab) for _ in each_layer_ab]
                return [_[index_lst[i]] for i, _ in enumerate(each_layer_peidui_lst)]
            else:
                return None
        except Exception as e:
            print(e)

    def layer1_pos(self, vector1a, vector1b, len_max):
        try:
            up_lst = []
            number = int(max(len_max / np.linalg.norm(vector1a), len_max / np.linalg.norm(vector1b)))
            for i in range(-1, 2 * number):
                for j in range(-1, 2 * number):
                    up_pos = vector1a * i + vector1b * j
                    if np.linalg.norm(up_pos) < len_max:
                        up_lst.append(up_pos)
            return up_lst
        except Exception as e:
            print(e)

    def make_multi_angle(self, angle_lst, make_lst):
        try:
            angle = angle_lst.pop(0)
            c = []
            for rotate_angle in range(angle):
                if make_lst:
                    for already_lst in make_lst:
                        t = deepcopy(already_lst)
                        t.append(rotate_angle / 180 * np.pi)
                        c.append(t)
                else:
                    c.append([rotate_angle / 180 * np.pi])
            if len(angle_lst) == 0:
                return c
            else:
                return self.make_multi_angle(angle_lst, c)
        except Exception as e:
            print(e)

    def twist_little_angle(self):
        try:
            formula_lst = []
            itemlst = []
            for item in self.project_tree.selectedItems():
                formula_lst.append(item.text(1))
                itemlst.append(item)
            if len(formula_lst) != 2:
                QtWidgets.QMessageBox.warning(self, 'error',
                                              'Please choose 2 crystal in Projectbox' + '\n' + '(Press Ctrl or Shift).')
            elif itemlst[0].text(0) == 'bulk' or itemlst[1].text(0) == 'bulk':
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose stack or layer(doubleclick to change)')
            else:
                try:
                    self.twist_window.close()
                except:
                     pass
                self.twist_obj1 = self.dic_Formula_Atoms[formula_lst[0]]  # Atoms1，对象
                self.twist_obj2 = self.dic_Formula_Atoms[formula_lst[1]]  # Atoms2，对象
                self.twist_window = Twisted_little_window()
                self.twist_window.signal_emit_information.connect(self.after_small_twisted_window)
        except Exception as e:
            print(e)

    def after_small_twisted_window(self, layerdis, vacuumdis, normal_strain_limit, shear_strain_limit,
                                   area_tolerance, minimum_gamma, angle,
                                   maximum_gamma, ab_relative_error):
        try:
            cell1 = self.twist_obj1.get_cell()
            cell2 = self.twist_obj2.get_cell()
            vector1a = np.array(cell1[0][:2])
            vector1b = np.array(cell1[1][:2])
            vector2a = np.array(cell2[0][:2])
            vector2b = np.array(cell2[1][:2])
            theta = angle / 180 * np.pi
            rotate_matrix = np.array([[np.cos(theta), -np.sin(np.sin(theta))],
                                      [np.sin(theta), np.cos(theta)]])
            new_vector2a = vector2a @ rotate_matrix
            new_vector2b = vector2b @ rotate_matrix
            i = 1
            while True:
                a1_vector, b1_vector, a2_vector, b2_vector = self.solve_twist_little_angle(vector1a, vector1b, new_vector2a, new_vector2b, minimum_gamma,
                                                                                           normal_strain_limit, shear_strain_limit, area_tolerance, 50 * i, maximum_gamma, ab_relative_error)
                if a1_vector is None:
                    i += i
                else:
                    break
            Hetero_junction_object = Hetero_junction()
            Hetero_junction_object.Rotate_angle = angle
            Hetero_junction_object.Normal_strain_limit = normal_strain_limit
            cell1 = self.twist_obj1.get_cell()
            cell1_a, cell1_b = cell1[0][:2], cell1[1][:2]
            A = np.array([cell1_a, cell1_b]).T
            r1a = np.linalg.solve(A, np.array(a1_vector))  # 求解线性方程组，简单,直角坐标系下----用晶胞坐标系表示
            r1b = np.linalg.solve(A, np.array(b1_vector))
            k11 = list(map(lambda x: round(x), r1a))
            k12 = list(map(lambda x: round(x), r1b))
            optimal_match_layer1 = [k11, k12]
            r1a = np.append(r1a, 0)
            r1b = np.append(r1b, 0)
            layer_down = self.deal_with_rotate(cut(self.twist_obj1, r1a, r1b, [0, 0, 1], origo=[0, 0, 0]))

            A = np.array([new_vector2a, new_vector2b]).T
            r2a = np.linalg.solve(A, np.array(a2_vector))
            r2b = np.linalg.solve(A, np.array(b2_vector))  # 求解线性方程组，简单,直角坐标系下----用晶胞坐标系表示
            k21 = list(map(lambda x: round(x), r2a))
            k22 = list(map(lambda x: round(x), r2b))
            optimal_match_layer2 = [k21, k22]
            r2b = np.append(r2b, 0)
            r2a = np.append(r2a, 0)
            Hetero_junction_object.Optimal_orientation = [optimal_match_layer1, optimal_match_layer2]
            layer_up = cut(self.twist_obj2, r2a, r2b, [0, 0, 1], origo=[0, 0, 0])
            self.layer_up = self.deal_with_rotate(layer_up)
            # start stacking
            layer_down_cell = layer_down.get_cell_lengths_and_angles()
            self.layer_up.translate(np.array([0, 0, layerdis + layer_down.get_cell()[2][2]]))
            layer_down_cell[2] += self.layer_up.get_cell_lengths_and_angles()[2] + layerdis + .01
            layer_down.extend(self.layer_up)
            layer_down.set_cell(layer_down_cell)
            # add_vacuum _layer
            self.Atomsobject = deepcopy(layer_down)
            layer_down = self.set_vacuum_layer(vacuumdis, addproject=False)

            # plot
            self.plot(layer_down, clear=True, plot=False, globalAtomsobject=False, dictionary=True, object=False,
                      Hetero_tab=False)
            Hetero_junction_object.Name = self.dirkey
            Hetero_junction_object.Atomsobject = layer_down
            self.add_hetero_junction(Hetero_junction_object)
            try:
                parent = self.project_tree.selectedItems()[0].parent()
                stackchild = QtWidgets.QTreeWidgetItem(parent)
                stackchild.setText(0, "stack")
                stackchild.setText(1, self.dirkey)
                if self.project_tree.selectedItems()[0].text(2) != "":
                    text1 = self.project_tree.selectedItems()[0].text(2)
                else:
                    text1 = self.project_tree.selectedItems()[0].text(1)
                if self.project_tree.selectedItems()[1].text(2) != "":
                    text2 = self.project_tree.selectedItems()[1].text(2)
                else:
                    text2 = self.project_tree.selectedItems()[1].text(1)
                stackchild.setText(2, text1 + '-' + text2)
            except Exception as e:
                print(e)
                stackchild.setText(2, self.project_tree.selectedItems()[0].text(2))
        except Exception as e:
            print(e)

    def solve_twist_little_angle(self, vector1a, vector1b, vector2a, vector2b, minimum_gama,
                                 normal_strain_limit, shear_strain_limit, area_error_limit, len_max, maximul_gama, ab_relative_error):
        len_error = np.sqrt(normal_strain_limit ** 2 + shear_strain_limit ** 2)
        up_lst = []
        number = int(max(len_max / np.linalg.norm(vector1a), len_max / np.linalg.norm(vector1b),
                         len_max / np.linalg.norm(vector2a),
                         len_max / np.linalg.norm(vector2b)))
        for i in range(-1, 2 * number):
            for j in range(-1, 2 * number):
                up_pos = vector2a * i + vector2b * j
                if np.linalg.norm(up_pos) < len_max:
                    up_lst.append(up_pos)
        peidui_lst = []
        matrix = np.array([vector1a, vector1b])
        for up_pos in up_lst:
            min_error = 10000
            r = up_pos @ np.linalg.inv(matrix)
            r1 = np.floor(r[0])
            r2 = np.floor(r[1])
            close_down_pos = None
            down_pos = r1 * vector1a + r2 * vector1b
            if np.linalg.norm(down_pos - up_pos) < len_error * np.linalg.norm(down_pos)\
                    and np.linalg.norm(down_pos - up_pos) < min_error * np.linalg.norm(down_pos):
                if abs(down_pos @ down_pos - down_pos @ up_pos) < normal_strain_limit * (down_pos @ down_pos):
                    close_down_pos = down_pos
                    min_error = np.linalg.norm(down_pos - up_pos) / np.linalg.norm(down_pos)
            down_pos = (r1 + 1) * vector1a + r2 * vector1b
            if np.linalg.norm(down_pos - up_pos) < len_error * np.linalg.norm(down_pos) and \
                    np.linalg.norm(down_pos - up_pos) < min_error * np.linalg.norm(down_pos):
                if abs(down_pos @ down_pos - down_pos @ up_pos) < normal_strain_limit * (down_pos @ down_pos):
                    close_down_pos = down_pos
                    min_error = np.linalg.norm(down_pos - up_pos) / np.linalg.norm(down_pos)
            down_pos = (r1 + 1) * vector1a + (r2 + 1) * vector1b
            if np.linalg.norm(down_pos - up_pos) < len_error * np.linalg.norm(down_pos) and \
                    np.linalg.norm(down_pos - up_pos) < min_error * np.linalg.norm(down_pos):
                if abs(down_pos @ down_pos - down_pos @ up_pos) < normal_strain_limit * (down_pos @ down_pos):
                    close_down_pos = down_pos
                    min_error = np.linalg.norm(down_pos - up_pos) / np.linalg.norm(down_pos)
            down_pos = r1 * vector1a + (r2 + 1) * vector1b
            if np.linalg.norm(down_pos - up_pos) < len_error * np.linalg.norm(down_pos) and \
                    np.linalg.norm(down_pos - up_pos) < min_error * np.linalg.norm(down_pos):
                if abs(down_pos @ down_pos - down_pos @ up_pos) < normal_strain_limit * (down_pos @ down_pos):
                    close_down_pos = down_pos
            if close_down_pos is not None:
                peidui_lst.append([close_down_pos, up_pos])
        if len(peidui_lst) < 2:
            return None, None, None, None
        else:
            min_squrare = 100000000000000000000
            truevector1a = ...
            truevector1b = ...
            truevector2a = ...
            truevector2b = ...
            for indexi in range(len(peidui_lst)):
                for indexj in range(indexi + 1, len(peidui_lst)):
                    if abs(np.linalg.norm(peidui_lst[indexi][0]) - np.linalg.norm(peidui_lst[indexj][0])) \
                            < ab_relative_error * np.linalg.norm(peidui_lst[indexj][0]):
                        n_vector1 = np.cross(peidui_lst[indexi][0], peidui_lst[indexj][0])
                        square1 = abs(n_vector1)
                        square2 = abs(np.cross(peidui_lst[indexi][1], peidui_lst[indexj][1]))
                        if abs(square1 - square2) < area_error_limit * square1:

                            if maximul_gama * np.linalg.norm(peidui_lst[indexi][0]) * \
                                np.linalg.norm(peidui_lst[indexj][0]) > square1 > \
                                minimum_gama * np.linalg.norm(peidui_lst[indexi][0]) * np.linalg.norm(peidui_lst[indexj][0]):

                                shear_strain = abs(np.arccos((peidui_lst[indexi][0] @ peidui_lst[indexj][0]) /
                                                             np.linalg.norm(peidui_lst[indexi][0]) / np.linalg.norm(peidui_lst[indexj][0])) -
                                                   np.arccos((peidui_lst[indexi][1] @ peidui_lst[indexj][1]) /
                                                             np.linalg.norm(peidui_lst[indexi][1]) / np.linalg.norm(peidui_lst[indexj][1])))
                                if shear_strain < shear_strain_limit:
                                    if square1 < min_squrare:
                                        min_squrare = square1
                                        if n_vector1 > 0:
                                            truevector1a = peidui_lst[indexi][0]
                                            truevector1b = peidui_lst[indexj][0]
                                            truevector2a = peidui_lst[indexi][1]
                                            truevector2b = peidui_lst[indexj][1]

                                        else:
                                            truevector1a = peidui_lst[indexj][0]
                                            truevector1b = peidui_lst[indexi][0]
                                            truevector2a = peidui_lst[indexj][1]
                                            truevector2b = peidui_lst[indexi][1]
            if min_squrare != 100000000000000000000:
                return truevector1a, truevector1b, truevector2a, truevector2b
            else:
                return None, None, None, None

    def classification_main(self):
        """ To judge whether the crys is 2D materials or not."""
        try:
            self.classify_window.close()
        except:
             pass
        try:
            self.classify_window = Classification_window()
            self.classify_window.signal_emit_dirs.connect(self.after_classify_window)
        except Exception as e:
            print(e)

    def after_classify_window(self, cifdir, outputdir):
        try:
            objs = FileOperations(cifdir)  # 创建FileOperations类
            objs = objs.read_files()  # 打开dir_path， 返回创建的 atoms对象列表
            zongnum = len(objs)
            num = 0
            for cry_obj in objs:
                try:
                    ClassificationCrystalOperationsNum1(cry_obj, outputdir)  # 创一个对象
                    num += 1
                    self.classify_window.pbar.setValue(int(num / zongnum * 100))
                except Exception as e:
                    print(e)
            self.classify_window.close()
        except Exception as e:
            print(e)

class Application(QtWidgets.QApplication):
    def __init__(self, argv):
        QtWidgets.QApplication.__init__(self, argv)

    def _slot_setStyle(self):
        self.setStyleSheet('')
        tmp = self.sender().objectName()
        print('tmp = ', tmp)
        if tmp == 'Day':
            app.setStyle('Fusion')
        else:
            app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())


if __name__ == "__main__":
    # app = QtWidgets.QApplication(argv)
    app = Application(argv)
    app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
    matLego = mywindow()  # 主窗口实例化
    matLego.actionNightmode.triggered.connect(app._slot_setStyle)
    # matLego.actionDaytimemode.triggered.connect(app._slot_setStyle)
    exit(app.exec_())





