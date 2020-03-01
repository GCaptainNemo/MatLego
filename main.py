"""
Main Function of Matlego 1.5, depend mainly on ase and pyqtgraph.
MatLego:1.5
Finish time: 2019/3/5
Main function: to build up model for hetero-junction, especially 2D electron Devices.

"""

from src.GUI import *
from src.database import *
from PyQt5.QtWidgets import QFileDialog
from ase.build import cut
from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifWriter
import pyqtgraph.opengl as gl
from pymatgen.io.ase import AseAtomsAdaptor
from ase.spacegroup import get_spacegroup
import re
import os
from PyQt5.QtCore import QEvent
from ase.io import read
from ase.atom import Atom
# import src.machine_learning_knn_classification as smlkc
from OpenGL.GLU import gluProject
import numpy as np
from fractions import Fraction
import qdarkstyle
import sys

class mywindow(QtWidgets.QMainWindow, Ui_MainWindow):  # 实现代码与界面的分离
    """  The fuction of the Mainwindow. """
    datanumber = 1 # database
    row_lst = []  # database
    true_id_lst = [] # database
    dic_Formula_Atoms = {}  # key = formula value = Atoms_object
    cell_view_num = 0        # whether draw cell or not
    coordinatesystem_view_num = 2  # Draw cartestian coordinate first
    grid_itemnum = 0         # 是否画网格
    judge2layer_num = 0      # 判断两层器件
    judge3layer_num = 0      # 判断三层器件
    judgecustomize_num = 0   # 判断人工判定器件
    plot_num = 0             # Ball-stick(1), Ball-filling(0)
    up_down = False
    def __init__(self):
        # self.view_widget.size
        self.keylock = True
        self.key_double_click = True
        super(mywindow, self).__init__()
        self.setupUi(self)
        self.view_widget.installEventFilter(self)
        self.actionNew.triggered.connect(self.newfile)
        self.actionNewToolBar.triggered.connect(self.newfile)
        self.actionOpen.triggered.connect(self.openfile)
        self.actionOpenToolBar.triggered.connect(self.openfile)
        self.actionExit.triggered.connect(self.close)
        # self.actionJudgecustomize.setIconVisibleInMenu()
        self.actionCut.triggered.connect(self.cut_customize)
        self.actionRotate.triggered.connect(self.Rotate)
        self.actionRotate_cutomize.triggered.connect(self.rotate_customize)
        # Judge electronic devices
        self.actionJudge3layer.triggered.connect(self.Judge3layer)
        self.actionJudge2layer.triggered.connect(self.Judge2layer)
        self.actionJudgecustomize.triggered.connect(self.Judgecustomize)
        # 批量操作
        self.actionTraversal_rotate.triggered.connect(self.traversalrotate)
        self.actionTraversal_Stack1.triggered.connect(self.traversalstack1)
        self.actionTraversal_Stack2.triggered.connect(self.traversalstack2)
        self.actionTraversal_Cut.triggered.connect(self.traversalcut)
        #
        self.actionCutToolBar.triggered.connect(self.judgetool)
        self.actionExport.triggered.connect(self.Export)
        self.actionImport.triggered.connect(self.Import)
        self.actionStack.triggered.connect(self.Stack_main_one)     # one stack direction is known
        self.actionStackToolbar.triggered.connect(self.Stack_main_two)  # two stack direction is known
        # self.actionClassification.triggered.connect(self.classification_main)
        # 两种箭头
        self.actiondrag_rectangle.triggered.connect(self.dragrectangle)
        self.actionNormalchoose.triggered.connect(self.normalchoose)
        self.actionSetlayer.triggered.connect(self.setlayer)         # tool setlayer
        self.actionMovelayer.triggered.connect(self.movelayer)       # tool movelayer
        self.actionleft_screen.triggered.connect(self.left_screen)
        self.actionright_screen.triggered.connect(self.right_screen)
        self.actionup_screen.triggered.connect(self.up_screen)
        self.actiondown_screen.triggered.connect(self.down_screen)
        self.actionmiddle_screen.triggered.connect(self.middle_screen)
        #
        self.actiona.triggered.connect(self.viewfroma)
        self.actionb.triggered.connect(self.viewfromb)
        self.actionc.triggered.connect(self.viewfromc)
        self.actionSave.triggered.connect(self.export_cif)
        self.actionViewCell.triggered.connect(self.view_cell)
        self.actionViewCoordinate_System.triggered.connect(self.viewCoordinate_system)
        self.actionViewCoordinate_cell.triggered.connect(self.viewCoordinatcell_system)
        self.actionViewGrid.triggered.connect(self.viewgrid)
        # 鼠标右击，right_menubar
        self.view_widget.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.view_widget.customContextMenuRequested['QPoint'].connect(self.rightMenuShow)
        self.actionViewDatabase.triggered.connect(self.viewdatabase)
        self.actionViewProject.triggered.connect(self.viewproject)
        self.actionViewText.triggered.connect(self.viewtext)
        #
        self.tree.clicked.connect(self.onTreeClicked)
        self.tree.doubleClicked.connect(self.doubleclickedontree)
        self.datatable.itemDoubleClicked.connect(self.double_clicked_on_datatable)
        self.datatable.itemClicked.connect(self.clicked_on_datatable)
        self.tab3_table_widget2.itemClicked.connect(self.clicked_on_tab3wg2)

        # self.view_widget.mousePressEvent()
        self.tablewidget.itemClicked.connect(self.handleItemClick)  # object box click
        self.calculate_distance.clicked.connect(self.calculatedistance)
        self.calculate_degree.clicked.connect(self.calculatedegree)
        self.calculate_vector.clicked.connect(self.calculatevector)
        self.filterbutton.clicked.connect(self.filterfrom_to)
        self.ballfilling_button.clicked.connect(self.ballfillingtype)
        self.ball_stick_button.clicked.connect(self.ball_and_stick)
        self.searchLineEdit.textChanged.connect(self.search)
        self.searchdatabase_LineEdit.textChanged.connect(self.searchdatabase)
        self.btn1.toggled.connect(lambda: self.btnstate(self.btn1))
        self.btn2.toggled.connect(lambda: self.btnstate(self.btn2))
        self.btn3.toggled.connect(lambda: self.btnstate(self.btn3))
        self.timer = QtCore.QTimer(self)
        self.timer.start(10)
        self.timer.timeout.connect(self.timeout_slot)
        # Two viewports
        self.actionPerspective.triggered.connect(self.viewport_perspective)
        self.actionOrthogonal.triggered.connect(self.viewport_orthogonal)
        # database
        self.action_db_connect.triggered.connect(self.connect_db)
        self.action_cif_to_db.triggered.connect(self.database_cif_to_db_main)    # create a database
        # self.action_db_to_cifs.triggered.connect(self.database_to_cif_files)    #  composite a cif file
        self.action_add_data.triggered.connect(self.use_cif_add_data)                 # Add data to database
        self.action_setting_atom.triggered.connect(self.edit_atom_para)
        self.Atoms_color = sf.crys_data.atom_color
        self.Atom_radii = sf.crys_data.atom_radii_lst
        self.datatable.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)  ######允许右键产生子菜单
        self.datatable.customContextMenuRequested['QPoint'].connect(self.Right_database_menushow)  ####右键菜单

    def edit_atom_para(self):
        try:
            self.edit_atom_para_window = edit_atom_para_window(self.Atoms_color, self.Atom_radii)
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

    def Right_database_menushow(self):
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
        self.copyline_Action.setText("Copy a line")
        self.copyline_Action.triggered.connect(self.copy_linedatabase)

        self.paste_Action = QtWidgets.QAction(self)
        self.paste_Action.setText("Paste")
        self.paste_Action.triggered.connect(self.paste_database)

        self.append_Action = QtWidgets.QAction(self)
        self.append_Action.setText("Append data")
        self.append_Action.triggered.connect(self.append_data)
        self.add_line_Action = QtWidgets.QAction(self)
        self.add_line_Action.setText("Add line")
        self.add_line_Action.triggered.connect(self.add_line)
        # self.to_cif_action = QtWidgets.QAction(self)
        # self.to_cif_action.setText("To CIF file")
        # self.to_cif_action.triggered.connect(self.data_to_cif)
        rightMenu.addAction(self.Edit_Action)
        rightMenu.addAction(self.copy_Action)
        rightMenu.addAction(self.copyline_Action)
        rightMenu.addAction(self.paste_Action)
        rightMenu.addAction(self.append_Action)
        rightMenu.addAction(self.add_line_Action)
        # rightMenu.addAction(self.to_cif_action)
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


    # def data_to_cif(self):
    #     try:
    #         # item = self.datatable.item(self.row_num, self.column_num)
    #         dir = QtGui.QFileDialog.getExistingDirectory(self, 'Select directory')
    #         if not self.row_lst:
    #             row = self.row_num
    #         else:
    #             row = self.row_lst[self.row_num]
    #         Formula = self.database_class.info_data_base[row][1]
    #         cif_dir = dir + '/' + Formula + '.cif'
    #         num = 1
    #         while os.path.exists(cif_dir):
    #             cif_dir = dir + '/' + Formula + '({}).cif'.format(num)
    #             num += 1
    #         print(Formula)
    #         self.database_set_cif(self.database_dir, cif_dir, Formula)
    #     except Exception as e:
    #         QtWidgets.QMessageBox.warning(self, 'error', str(e))
    #         print(e)


    def beforeEdit(self):
        try:
            item = self.datatable.item(self.row_num, self.column_num)
            self.double_clicked_on_datatable(item)
        except Exception as e:
            print(e)

    def append_data(self):
        try:
            self.database_class.exhibit_all_info(fenkuai=False)
            ID_lst = list(set([self.database_class.info_data_base[i][1]
                               for i in range(len(self.database_class.info_data_base))]))
            self.append_data_window = Append_data_window(ID_lst)
            self.append_data_window.signal_emit.connect(self.after_append_data_window)
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error',
                                          'Connect database please.')
    def after_append_data_window(self, ID, formula):
        try:
            info = self.database_class.info_data_base
            TrueID = max([info[i][0] for i in range(len(info))]) + 1
            self.database_class.insert_row_into_test(TrueID, ID, formula)
            self.data_base_init()
        except Exception as e:
            print(e)

    def clicked_on_tab3wg2(self, item):
        try:
            print("text = ", item.text())
            if item.column() == 4:
                if item.text() != "":
                    address_lst = []
                    address = ""
                    for letter in item.text():
                        if letter != "+":
                            address += letter
                        else:
                            address_lst.append(address)
                            address = ""
                    address_lst.append(address)
                    print("address_lst = ", address_lst)
                    self.databaseshow(3)
                    self.show_image_table.setColumnCount(1)
                    self.show_image_table.setRowCount(len(address_lst))
                    self.show_image_table.setIconSize(QtCore.QSize(300, 200))
                    self.show_image_table.setColumnWidth(0, 300)
                    for i in range(len(address_lst)):
                        newitem = QtWidgets.QTableWidgetItem()
                        try:
                            Icon = QtGui.QIcon(address_lst[i])
                            newitem.setIcon(QtGui.QIcon(Icon))
                        except Exception as e:
                            print(e)
                        newitem.setFlags(QtCore.Qt.ItemIsEnabled)
                        self.show_image_table.setItem(i, 0, newitem)
                        self.show_image_table.setRowHeight(i, 200)
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
                            print("address =", address)
                            address_lst.append(address)
                    print("address_lst = ", address_lst)
                    self.show_image_table.setColumnCount(1)
                    self.show_image_table.setRowCount(len(address_lst))
                    self.show_image_table.setIconSize(QtCore.QSize(300, 200))
                    self.show_image_table.setColumnWidth(0, 300)
                    for i in range(len(address_lst)):
                        newitem = QtWidgets.QTableWidgetItem()
                        try:
                            Icon = QtGui.QIcon(address_lst[i])
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
            ID_lst = list(set([self.database_class.info_data_base[i][1]
                               for i in range(len(self.database_class.info_data_base))]))
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
            # ID_repeatlist = []
            # for c in self.database_class.info_data_base:
            #     if c[1] == data[1]:
            #         ID_repeatlist.append
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
                    information_lst = copy.deepcopy(self.informatino_lst)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error',
                                                  'Different shape.')

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
                l = []
                for c in datafenkuai_false:
                    if c[1] == data[1]:
                        l.append(c[0])
                for True_ID in l:
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
                l = []
                for c in datafenkuai_false:
                    if c[1] == data[1]:
                        l.append(c[0])
                for true_ID in l:
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
            self.database_dir, filetype = QFileDialog.getOpenFileName(self, '选择文件', '', 'files(*.db)')
            self.data_base_init()
            if self.database_dir != "":              # Not empty
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
            # self.searchdatabase()
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
        cif_dir, filetype = QFileDialog.getOpenFileName(self, '选择文件', '', 'files(*.cif)')
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
            re.findall('_symmetry_space_group_name_H-M(.*?)\n', content)[0].strip()
            a = float(re.findall('_cell_length_a(.*?)\n', content)[0].strip())
            b = float(re.findall('_cell_length_b(.*?)\n', content)[0].strip())
            c = float(re.findall('_cell_length_c(.*?)\n', content)[0].strip())
            angle_alpha = float(re.findall('_cell_angle_alpha(.*?)\n', content)[0].strip())
            angle_beta = float(re.findall('_cell_angle_beta(.*?)\n', content)[0].strip())
            angle_gamma = float(re.findall('_cell_angle_gamma(.*?)\n', content)[0].strip())
            volume = float(re.findall('_cell_volume(.*?)\n', content)[0].strip())
            conn = sqlite3.connect(db_dir)
            cn = conn.cursor()
            # cn.execute(
            #     "INSERT INTO test (ID, Formula, A, B, C, Angle_a, Angle_b, Angle_C, Volume) VALUES (%d, %s, %g, %g, %g, %g, %g, %g, %g)" % (
            #         i, formula, a, b, c, angle_alpha, angle_beta, angle_gamma, volume))
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
        # cn.execute('''CREATE TABLE test
        #                                          (ID       INT PRIMARY KEY    NOT NULL,
        #                                           Formula             TEXT    NOT NULL,
        #                                           A                   TEXT    ,
        #                                           B                   TEXT    ,
        #                                           C                   TEXT    ,
        #                                           Angle_a             TEXT    ,
        #                                           Angle_b             TEXT    ,
        #                                           Angle_c             TEXT    ,
        #                                           Volume              TEXT    ,
        #
        #                                           Hetero_junction     TEXT,
        #                                           Optimal_Match       TEXT,
        #                                           Layer               TEXT,
        #                                           Binding_energy      TEXT,
        #                                           Schottky_barrier    TEXT,
        #
        #                                           Device              Text,
        #                                           Optimal_MatchD      Text,
        #                                           LayerD              Text,
        #                                           Schottky_barrierD   TEXT,
        #                                           Image_dir           Text)''')
        conn.commit()
        conn.close()

    # table_site
    # def database_cif_to_db(self, db_dir, cif_dir):
    #     try:
    #         f = open(cif_dir, 'r')  # cif_dir = '/Users/mac/Desktop/1.cif'
    #         content = f.read()
    #         Atomsobject = read(cif_dir)
    #         r_formula = Atomsobject.get_chemical_formula(mode='hill')
    #         print(r_formula)
    #         formula = "'" + r_formula + "'"
    #         r_sites = re.findall('_atom_site_occupancy(.*?)loop_', content, re.S)[0]
    #         sites = r_sites.split('\n')
    #         sites.remove('')
    #
    #         conn = sqlite3.connect(db_dir)  # '/Users/mac/Desktop/crystal.db'
    #         cnn = conn.cursor()
    #         cnn.execute('''CREATE TABLE %s
    #                                       (Element             TEXT    NOT NULL,
    #                                        Num                 TEXT    NOT NULL,
    #                                        U1                  INT     NOT NULL,
    #                                        X                   REAL    NOT NULL,
    #                                        Y                   REAL    NOT NULL,
    #                                        Z                   REAL    NOT NULL,
    #                                        U2                  INT     NOT NULL)''' % formula)
    #         conn.commit()
    #         for x in sites:
    #             if x is '':
    #                 sites.remove(x)
    #         for i in range(0, len(sites)):
    #             row_lst = sites[i].split()
    #             print(row_lst)
    #             element = "'" + row_lst[0] + "'"
    #             order = row_lst[0] + str(i + 1)
    #             num = "'" + order + "'"
    #             cnn.execute("INSERT INTO %s (Element,Num,U1,X,Y,Z,U2) VALUES (%s, %s, 1, %g, %g, %g, 1)" % (
    #                 formula, element, num, float(row_lst[3]), float(row_lst[4]), float(row_lst[5])))
    #             conn.commit()
    #         conn.close()
    #     except Exception as e:
    #         print(e)

    # def database_to_cif_files(self):
    #     try:
    #         # database_dir, filetype = QFileDialog.getOpenFileName(self, '选择文件', '', 'files(*.db)')
    #         self.data_to_cif_window = database_to_cif(self.database_dir)
    #         self.data_to_cif_window.signal_emit_formula.connect(self.database_set_cif)
    #     except Exception as e:
    #         print(e)
    #         QtWidgets.QMessageBox.warning(self, 'error',
    #                                       'Connect database please.')

    # def database_set_cif(self, db_dir, cif_dir, r_formula):  # db_dir:database.py中创立的数据库的路径，cif_dir：新建立的cif文件的路径
    #     try:
    #         num = 0
    #         formula = "'" + r_formula + "'"
    #         conn = sqlite3.connect(db_dir)  # dir = "/Users/mac/Desktop/crystal.db"
    #         debuglogger("1")
    #         cnn = conn.cursor()
    #         # prepare
    #         cnn.execute("SELECT * FROM test WHERE Formula={}".format(formula))
    #         row = list(cnn.fetchone())
    #         num = 1
    #         print(row)
    #         debuglogger("1.5")
    #         cnn.execute("SELECT * FROM {0}".format(formula))
    #         sites = [list(x) for x in cnn.fetchall()]
    #         num = 2
    #         # write
    #         if not cif_dir:
    #             f = open(cif_dir, 'a')  # cif_dir = "/Users/mac/Desktop/foo.cif"
    #         num = 3
    #         f.write("# generated using pymatgen\n")
    #         debuglogger("2")
    #         f.write("data_{}\n".format(formula.strip("'")))
    #         f.write("_cell_length_a\t{0}\n".format(row[3]))
    #         f.write("_cell_length_b\t{0}\n".format(row[4]))
    #         f.write("_cell_length_c\t{0}\n".format(row[5]))
    #         f.write("_cell_angle_alpha\t{0}\n".format(row[6]))
    #         f.write("_cell_angle_beta\t{0}\n".format(row[7]))
    #         f.write("_cell_angle_gamma\t{0}\n".format(row[8]))
    #         f.write("_symmetry_Int_Tables_number   1\n")
    #         f.write("_cell_volume\t{0}\n".format(row[9]))
    #         f.write("loop_\n")
    #         f.write(" _symmetry_equiv_pos_site_id\n")
    #         f.write(" _symmetry_equiv_pos_as_xyz\n")
    #         f.write("  1  'x, y, z'\n")
    #         f.write("loop_\n")
    #         f.write(" _atom_site_type_symbol\n")
    #         f.write(" _atom_site_label\n")
    #         f.write(" _atom_site_symmetry_multiplicity\n")
    #         f.write(" _atom_site_fract_x\n")
    #         f.write(" _atom_site_fract_y\n")
    #         f.write(" _atom_site_fract_z\n")
    #         f.write(" _atom_site_occupancy\n")
    #         debuglogger("3")
    #         print(sites)
    #         for site in sites:
    #             f.write("%s %s %d %g %g %g %g\n" % (site[0], site[1], site[2], site[3], site[4], site[5], site[6]))
    #         f.write("loop_\n")
    #         f.write(" _atom_site_moment_label\n")
    #         f.write(" _atom_site_moment_crystalaxis_x\n")
    #         f.write(" _atom_site_moment_crystalaxis_y\n")
    #         f.write(" _atom_site_moment_crystalaxis_z\n")
    #         f.close()
    #         QtWidgets.QMessageBox.information(self, 'Information:', 'File address:' + cif_dir)
    #     except Exception as e:
    #         print(e, type(e))
    #         if num == 0:
    #             QtWidgets.QMessageBox.warning(self, 'error', "There is no crystaline information.")
    #         elif num == 1:
    #             QtWidgets.QMessageBox.warning(self, 'error', "There is no atoms' information.")



    def stdout(self, dir):  # database.py中创建的数据库的地址
        '''
        需要在数据库面板中显示的性质
        :param dir: 数据库的路径
        :return: None
        '''
        try:
            conn = sqlite3.connect(dir)  # dir:"/Users/mac/Desktop/crystal.db"
            cnn = conn.cursor()

            r_formula = input("你想要查询哪种材料的性质，请输入它的化学式: ")  # 想要查询材料的化学式，eg.Ag
            formula = "'" + r_formula + "'"
            cnn.execute("SELECT * FROM test WHERE Formula={}".format(formula))
            row = list(cnn.fetchone())
            print("Formula: %s" % row[1])
            print("A: %g" % row[2])
            print("B: %g" % row[3])
            print("C: %g" % row[4])
            print("Angle_alpha: %g" % row[5])
            print("Angle_beta: %g" % row[6])
            print("Angle_gamma: %g" % row[7])
            print("Volumn: %s" % row[8])
            print("Hetero_junction: %s" % row[9])
            print("Optimal_Match: %s" % row[10])
            print("Layer: %s" % row[11])
            print("Binding_energy: %g" % row[12])
            print("Schottky_barrier: %g" % row[13])
        except Exception as e:
            print(e)

    def viewport_perspective(self):
        try:
            self.view_widget.opts['distance'] = 10
            self.view_widget.opts['fov'] = 60
            self.plot(self.Atomsobject)
        except Exception as e:
            print(e)

    def viewport_orthogonal(self):
        try:
            self.view_widget.opts['distance'] = 2000
            self.view_widget.opts['fov'] = 2
            self.plot(self.Atomsobject)
        except Exception as e:
            print(e)

    def rotate_customize(self):
        try:
            formula_lst = []
            itemlst = []
            for item in self.tree.selectedItems():
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
                self.scatter_plot = scatter_plot_customize(self.rotate_obj1, self.rotate_obj2)
                self.scatter_plot.emit_information_signal.connect(self.plot_after_scatter)
        except Exception as e:
            print(e)

    def Rotate(self):
        """To stack along different orientations."""
        try:
            formula_lst = []
            itemlst = []
            for item in self.tree.selectedItems():
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
                self.scatter_plot = scatter_plot(self.rotate_obj1, self.rotate_obj2)
                self.scatter_plot.emit_information_signal.connect(self.plot_after_scatter)
        except Exception as e:
            print(e)

    def plot_after_scatter(self, degree, a_vector, b_vector, layer_distance, vaccum_distance):
        try:
            # layer_down
            cell1 = self.rotate_obj1.get_cell()
            cell1_a, cell1_b = [cell1[0][0], cell1[0][1]], [cell1[1][0], cell1[1][1]]
            A = np.array([cell1_a, cell1_b]).T
            r1a = np.linalg.solve(A, np.array(a_vector))  # 求解线性方程组，简单,直角坐标系下----用晶胞坐标系表示
            r1a = np.append(r1a, 0)
            r1b = np.linalg.solve(A, np.array(b_vector))
            r1b = np.append(r1b, 0)
            layer_down = self.deal_with_rotate(cut(self.rotate_obj1, r1a, r1b, [0, 0, 1], origo=[0, 0, 0]))
            # self.layer_down = self.deal_with_rotate(layer_down)
            # layer_up
            cell2_before_rotate = self.rotate_obj2.get_cell()
            cell2_a, cell2_b = [cell2_before_rotate[0][0], cell2_before_rotate[0][1]], \
                               [cell2_before_rotate[1][0], cell2_before_rotate[1][1]]
            # degree = degree / 180 * np.pi
            rotate_matrix = np.array(
                [[np.cos(degree), np.sin(degree)], [-np.sin(degree), np.cos(degree)]])
            new_cell2_a = cell2_a @ rotate_matrix
            new_cell2_b = cell2_b @ rotate_matrix
            A = np.array([new_cell2_a, new_cell2_b]).T
            r2a = np.linalg.solve(A, np.array(a_vector))
            r2a = np.append(r2a, 0)
            r2b = np.linalg.solve(A, np.array(b_vector))  # 求解线性方程组，简单,直角坐标系下----用晶胞坐标系表示
            r2b = np.append(r2b, 0)
            layer_up = cut(self.rotate_obj2,  r2a, r2b, [0, 0, 1], origo=[0, 0, 0])
            self.layer_up = self.deal_with_rotate(layer_up)
            # start stacking
            layer_down_cell = layer_down.get_cell_lengths_and_angles()
            self.layer_up.translate(np.array([0, 0, layer_distance + layer_down.get_cell()[2][2]]))
            layer_down_cell[2] += self.layer_up.get_cell_lengths_and_angles()[2] + layer_distance + 0.01
            layer_down.extend(self.layer_up)
            layer_down.set_cell(layer_down_cell)
            # add_vacuum _layer
            c_length = layer_down.get_cell_lengths_and_angles()[2]
            cell_par = layer_down.get_cell()
            cell_par[2] *= (c_length + vaccum_distance) / c_length
            layer_down.set_cell(cell_par)
            layer_down.translate(cell_par[2] / (c_length + vaccum_distance) * vaccum_distance / 2)
            # plot
            self.plot(layer_down, clear=True, globalAtomsobject=False, dictionary=True)
            try:
                parent = self.tree.selectedItems()[0].parent()
                stackchild = QtWidgets.QTreeWidgetItem(parent)
                stackchild.setText(0, "stack")
                stackchild.setText(1, self.dirkey)
                if self.tree.selectedItems()[0].text(2) != "":
                    text1 = self.tree.selectedItems()[0].text(2)
                else:
                    text1 = self.tree.selectedItems()[0].text(1)
                if self.tree.selectedItems()[1].text(2) != "":
                    text2 = self.tree.selectedItems()[1].text(2)
                else:
                    text2 = self.tree.selectedItems()[1].text(1)
                stackchild.setText(2, text1 + '-' + text2)
            except Exception as e:
                print(e)
                stackchild.setText(2, self.tree.selectedItems()[0].text(2))
            finally:
                self.scatter_plot.close()
        except Exception as e:
            print(e)

    def middle_screen(self):
        """ Put the model in the middle of the screen."""
        try:
            print("before_viewport = ", self.view_widget.getViewport())
            height = self.view_widget.height()
            width = self.view_widget.width()
            viewport = list(self.view_widget.getViewport())
            viewport[0] = 0
            viewport[1] = 0
            viewport[2] = width
            viewport[3] = height
            viewport = tuple(viewport)
            self.view_widget.opts['viewport'] = viewport
        except Exception as e:
            print(e)


    def down_screen(self):
        """ Move the model to the bottom of the screen."""
        try:
            height = self.view_widget.height()
            width = self.view_widget.width()
            viewport = list(self.view_widget.getViewport())
            viewport[1] -= 30
            viewport[2] = width - viewport[0]
            viewport[3] = height - viewport[1]
            viewport = tuple(viewport)
            self.view_widget.opts['viewport'] = viewport
        except Exception as e:
            print(e)


    def up_screen(self):
        """ Move the model to the top of the screen."""
        try:
            height = self.view_widget.height()
            width = self.view_widget.width()
            viewport = list(self.view_widget.getViewport())
            viewport[1] += 30
            viewport[2] = width - viewport[0]
            viewport[3] = height - viewport[1]
            viewport = tuple(viewport)
            self.view_widget.opts['viewport'] = viewport
        except Exception as e:
            print(e)

    def right_screen(self):
        """ Move the model to the right of the screen."""
        try:
            height = self.view_widget.height()
            width = self.view_widget.width()
            viewport = list(self.view_widget.getViewport())
            viewport[0] += 30
            viewport[2] = width - viewport[0]
            viewport[3] = height - viewport[1]
            viewport = tuple(viewport)
            self.view_widget.opts['viewport'] = viewport
        except Exception as e:
            print(e)


    def left_screen(self):
        """ Move the model to the left of the screen."""
        try:
            height = self.view_widget.height()
            width = self.view_widget.width()
            viewport = list(self.view_widget.getViewport())
            viewport[0] -= 30
            viewport[2] = width - viewport[0]
            viewport[3] = height - viewport[1]
            viewport = tuple(viewport)
            self.view_widget.opts['viewport'] = viewport
        except Exception as e:
            print(e)

    def ball_and_stick(self):
        """ Set the atomic model to the ball-stick model."""
        try:
            current = self.tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.tree.selectedItems()) == 1:  # 用户只选择一个的情况
                    self.plot_num = 1
                    self.plot(self.Atomsobject)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def ballfillingtype(self):
        """ Set the atomic model to the ball-filling model."""
        try:
            current = self.tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.tree.selectedItems()) == 1:  # 用户只选择一个的情况
                    self.plot_num = 0
                    self.plot(self.Atomsobject)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)



    def Judge3layer(self):
        """ To judge whether the model is 3 layer device(conductor-semiconductor-conductor) or not"""
        try:
            global judge3conductivity_lst
            try:
                for item in judge3conductivity_lst:
                    for i in range(3):
                        item.setBackground(i, QtGui.QBrush(QtGui.QColor(35, 38, 41)))
            except Exception as e:
                print(e)
            if self.judge3layer_num % 2 == 0:
                it = QtWidgets.QTreeWidgetItemIterator(self.tree)
                judge3conductivity_lst = []
                while it.value():
                    if it.value().text(2) != "":         # 用project栏heterostructure一栏判断
                        Heterostructure_text = it.value().text(2)
                        l = []
                        layerletter = []
                        for letter_index in range(len(Heterostructure_text)):
                            letter = Heterostructure_text[letter_index]
                            if letter == "-":
                                l.append(layerletter)
                                layerletter = []
                            elif letter.isupper() and letter_index < len(Heterostructure_text) - 1 and \
                                    Heterostructure_text[letter_index + 1].islower():
                                layerletter.append(letter + Heterostructure_text[letter_index + 1])
                            elif letter.isupper() and letter_index < len(Heterostructure_text) - 1 and \
                                    (not Heterostructure_text[letter_index + 1].islower()):
                                layerletter.append(letter)
                            elif letter.isupper() and letter_index == len(Heterostructure_text) - 1:
                                layerletter.append(letter)
                        l.append(layerletter)
                        if len(l) == 3:
                            print(l)
                            conductlist = [sf.crys_data.conductor_Elemental_composition,
                                           sf.crys_data.semiconductor_Element_composition,
                                           sf.crys_data.Insulator_Element_composition]
                            str_conductivity = ""
                            for Elementlist in l:
                                for i in range(3):
                                    for conduct_element in conductlist[i]:
                                        if len(conduct_element) == len(Elementlist):
                                            for element in Elementlist:  # Atomsobject 的元素
                                                if element not in conduct_element:
                                                    break
                                                elif element == Elementlist[-1]:
                                                    if i == 0:
                                                        num = 1
                                                    if i == 1:
                                                        num = 2
                                                    if i == 2:
                                                        num = 3
                                if num == 1:
                                    str_conductivity += "Conductor"
                                elif num == 2:
                                    str_conductivity += "Semiconductor"
                                elif num == 3:
                                    str_conductivity += "Insulator"
                                elif num == 0:
                                    str_conductivity += ""
                            if str_conductivity == "ConductorSemiconductorConductor":
                                for i in range(3):
                                    it.value().setBackground(i, QtGui.QBrush(QtGui.QColor(50, 0, 50)))
                                judge3conductivity_lst.append(it.value())
                    it += 1
            if len(judge3conductivity_lst) == 0:
                QtWidgets.QMessageBox.warning(self, 'error', "There's no Conductor-semiconductor-Conductor heterostructure (3layer electonic devices).")
            self.judge3layer_num += 1
        except Exception as e:
            print(e)

    def Judgecustomize(self):
        """ To judge the model according to the clients."""
        try:
            if self.judgecustomize_num % 2 == 1:
                global judgecustomize_lst
                try:
                    for item in judgecustomize_lst:
                        for i in range(3):
                            item.setBackground(i, QtGui.QBrush(QtGui.QColor(35, 38, 41)))
                except Exception as e:
                    print(e)
            else:
                self.Judge_customize_window = customize_judge()
                self.Judge_customize_window.signal_emit_conductivity.connect(self.afterJudgecustomize_window)
            self.judgecustomize_num += 1
        except Exception as e:
            print(e)

    def afterJudgecustomize_window(self, string_lst):
        try:
            global judgecustomize_lst
            try:
                for item in judgecustomize_lst:
                    for i in range(3):
                        item.setBackground(i, QtGui.QBrush(QtGui.QColor(35, 38, 41)))
            except Exception as e:
                print(e)
            string = ""
            for i in string_lst:
                string += i
            it = QtWidgets.QTreeWidgetItemIterator(self.tree)
            judgecustomize_lst = []
            while it.value():
                if it.value().text(2) != "":         # 用project栏heterostructure一栏判断
                    Heterostructure_text = it.value().text(2)
                    l = []
                    layerletter = []
                    for letter_index in range(len(Heterostructure_text)):
                        letter = Heterostructure_text[letter_index]
                        if letter == "-":
                            l.append(layerletter)
                            layerletter = []
                        elif letter.isupper() and letter_index < len(Heterostructure_text) - 1 and Heterostructure_text[letter_index + 1].islower():
                            layerletter.append(letter + Heterostructure_text[letter_index + 1])
                        elif letter.isupper() and letter_index < len(Heterostructure_text) - 1 and (not Heterostructure_text[letter_index + 1].islower()):
                            layerletter.append(letter)
                        elif letter.isupper() and letter_index == len(Heterostructure_text) - 1:
                            layerletter.append(letter)
                    l.append(layerletter)
                    if len(l) == len(string_lst):
                        print(l)
                        conductlist = [sf.crys_data.conductor_Elemental_composition,
                                       sf.crys_data.semiconductor_Element_composition,
                                       sf.crys_data.Insulator_Element_composition]
                        str_conductivity = ""
                        for Elementlist in l:
                            for i in range(3):
                                for conduct_element in conductlist[i]:
                                    if len(conduct_element) == len(Elementlist):
                                        for element in Elementlist:  # Atomsobject 的元素
                                            if element not in conduct_element:
                                                break
                                            elif element == Elementlist[-1]:
                                                if i == 0:
                                                    num = 1
                                                if i == 1:
                                                    num = 2
                                                if i == 2:
                                                    num = 3
                            if num == 1:
                                str_conductivity += "Conductor"
                            elif num == 2:
                                str_conductivity += "Semiconductor"
                            elif num == 3:
                                str_conductivity += "Insulator"
                            elif num == 0:
                                str_conductivity += ""
                        if str_conductivity == string:
                            for i in range(3):
                                it.value().setBackground(i, QtGui.QBrush(QtGui.QColor(0, 50, 50)))
                            judgecustomize_lst.append(it.value())
                it += 1
            if len(judgecustomize_lst) == 0:
                QtWidgets.QMessageBox.warning(self, 'error', "There's no customize heterostructure ({}layer electonic devices).".format(str(len(string_lst))))
        except Exception as e:
            print(e)


    def Judge2layer(self):
        """ To judge whether the model is 2 layer device(semiconductor-semiconductor) or not"""
        try:
            global judge2conductivity_lst
            try:
                for item in judge2conductivity_lst:
                    for i in range(3):
                        item.setBackground(i, QtGui.QBrush(QtGui.QColor(35, 38, 41)))
            except Exception as e:
                print(e)
            if self.judge2layer_num % 2 == 0:
                it = QtWidgets.QTreeWidgetItemIterator(self.tree)
                judge2conductivity_lst = []
                while it.value():
                    if it.value().text(2) != "":         # 用project栏heterostructure一栏判断
                        Heterostructure_text = it.value().text(2)
                        l = []
                        layerletter = []
                        for letter_index in range(len(Heterostructure_text)):
                            letter = Heterostructure_text[letter_index]
                            if letter == "-":
                                l.append(layerletter)
                                layerletter = []
                            elif letter.isupper() and letter_index < len(Heterostructure_text) - 1 and Heterostructure_text[letter_index + 1].islower():
                                layerletter.append(letter + Heterostructure_text[letter_index + 1])
                            elif letter.isupper() and letter_index < len(Heterostructure_text) - 1 and (not Heterostructure_text[letter_index + 1].islower()):
                                layerletter.append(letter)
                            elif letter.isupper() and letter_index == len(Heterostructure_text) - 1:
                                layerletter.append(letter)
                        l.append(layerletter)
                        if len(l) == 2:
                            print(l)
                            conductlist = [sf.crys_data.conductor_Elemental_composition,
                                           sf.crys_data.semiconductor_Element_composition,
                                           sf.crys_data.Insulator_Element_composition]
                            str_conductivity = ""
                            for Elementlist in l:
                                for i in range(3):
                                    for conduct_element in conductlist[i]:
                                        if len(conduct_element) == len(Elementlist):
                                            for element in Elementlist:  # Atomsobject 的元素
                                                if element not in conduct_element:
                                                    break
                                                elif element == Elementlist[-1]:
                                                    if i == 0:
                                                        num = 1
                                                    if i == 1:
                                                        num = 2
                                                    if i == 2:
                                                        num = 3
                                if num == 1:
                                    str_conductivity += "Conductor"
                                elif num == 2:
                                    str_conductivity += "Semiconductor"
                                elif num == 3:
                                    str_conductivity += "Insulator"
                                elif num == 0:
                                    str_conductivity += ""
                            if str_conductivity == "Semiconductor"*2:
                                for i in range(3):
                                    it.value().setBackground(i, QtGui.QBrush(QtGui.QColor(0, 50, 50)))
                                judge2conductivity_lst.append(it.value())
                    it += 1
            if len(judge2conductivity_lst) == 0:
                QtWidgets.QMessageBox.warning(self, 'error', "There's no semiconductor-semiconductor heterostructure (2layer electonic devices).")
            self.judge2layer_num += 1
        except Exception as e:
            print(e)

    def viewtext(self):
        """ To (not) view the text widget."""
        try:
            global Text_widget_height
            if self.tab.isHidden():
                self.tab.setVisible(True)
                self.view_widget.resize(self.view_widget.width(), self.view_widget.height() - Text_widget_height)
            else:
                self.tab.setVisible(False)
                self.view_widget.resize(self.view_widget.width(), self.view_widget.height() + Text_widget_height)
        except Exception as e:
            print(e)

    def viewproject(self):
        """ To (not) view the project widget."""
        try:
            global project_widget_width
            if self.tab_tree.isHidden():
                self.tab_tree.setVisible(True)
                self.view_widget.resize(self.view_widget.width() - project_widget_width, self.view_widget.height())
            else:
                self.tab_tree.setVisible(False)
                self.view_widget.resize(self.view_widget.width() + project_widget_width, self.view_widget.height())
        except Exception as e:
            print(e)

    def viewdatabase(self):
        """ To (not) view the database widget."""
        try:
            global view_databas_widget_width
            if self.vertical_widget.isHidden():
                self.vertical_widget.setVisible(True)
                self.view_widget.resize(self.view_widget.width() - view_databas_widget_width, self.view_widget.height())
            else:
                self.vertical_widget.setVisible(False)
                self.view_widget.resize(self.view_widget.width() + view_databas_widget_width, self.view_widget.height())
        except Exception as e:
            print(e)

    def timeout_slot(self):
        """ To adjust the size of screen. """
        try:
            global view_databas_widget_width, project_widget_width, Text_widget_height
            if self.vertical_widget.width() != 0:
                view_databas_widget_width = self.vertical_widget.width()
            if self.tab_tree.width() != 0:
                project_widget_width = self.tab_tree.width()
            if self.tab.height() != 0:
                Text_widget_height = self.tab.height()
        except Exception as e:
            print(e)
        try:
            viewport = list(self.view_widget.getViewport())
            width = self.view_widget.width()
            height = self.view_widget.height()
            viewport[2] = width - viewport[0]
            viewport[3] = height - viewport[1]
            self.view_widget.opts['viewport'] = tuple(viewport)
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
            # self.searchdatabase()
            if number == 1:        # Crystal
                self.datanumber = 1
                self.datatable.setColumnCount(9)
                self.datatable.setHorizontalHeaderLabels(['ID', 'Formula', 'a(Å)', 'b(Å)', 'c(Å)',
                                                          'α(°)', 'β(°)', 'γ(°)', 'Volume'])
                self.database_class.exhibit_all_info(fenkuai=True)
                data = self.database_class.info_data_base
                if not self.true_id_lst:    # search里没有信息,true_id_lst为空
                    data_info = [data[i][0] for i in range(len(data))]
                    self.datatable.setRowCount(len(data_info))
                    for i in range(len(data_info)):
                        for k in range(1, 10):
                            newitem = QtWidgets.QTableWidgetItem(str(list(data_info[i])[k]))
                            newitem.setTextAlignment(5 | 5)
                            self.datatable.setItem(i, k - 1, newitem)
                else:
                    self.database_class.exhibit_all_info(fenkuai=False)
                    data_fenkuaiFalse = self.database_class.info_data_base
                    true_id_lst = [data[i][0][0] for i in range(len(data))]
                    # 交集
                    self.Reflect_id_lst = list(set(true_id_lst).intersection(set(self.true_id_lst)))
                    self.datatable.setRowCount(len(self.Reflect_id_lst))
                    id_lst = [data_fenkuaiFalse[i][0] for i in range(len(data_fenkuaiFalse))]
                    print("id_lst = ", id_lst)
                    for i in range(len(self.Reflect_id_lst)):
                        for k in range(1, 10):
                            newitem = QtWidgets.QTableWidgetItem(str(data_fenkuaiFalse[id_lst.index(self.Reflect_id_lst[i])][k]))
                            newitem.setTextAlignment(5 | 5)
                            self.datatable.setItem(i, k - 1, newitem)
                self.datatable.setShowGrid(True)
            elif number == 2:      # Hetero-junction
                self.datanumber = 2
                try:
                    self.database_class.exhibit_all_info(fenkuai=False)
                    data_info = self.database_class.info_data_base
                    self.Reflect_id_lst = []
                    self.datatable.setColumnCount(5)
                    self.datatable.setHorizontalHeaderLabels(['Hetero-junction', 'Optimal Match', '#Layer',
                                                              'Binding energy(eV)', 'Schottky_barrier(eV)'])
                    if not self.true_id_lst:
                        self.datatable.setRowCount(len(data_info))  # 添加信息
                        for i in range(len(data_info)):
                            for k in range(10, 15):
                                newitem = QtWidgets.QTableWidgetItem(str(list(data_info[i])[k]))
                                newitem.setTextAlignment(5 | 5)
                                self.datatable.setItem(i, k - 10, newitem)
                    else:
                        self.datatable.setRowCount(len(self.true_id_lst))
                        id_lst = [data_info[i][0] for i in range(len(data_info))]

                        for i in range(len(self.true_id_lst)):
                            for k in range(10, 15):
                                newitem = QtWidgets.QTableWidgetItem(str(data_info[id_lst.index(self.true_id_lst[i])][k]))
                                newitem.setTextAlignment(5 | 5)
                                self.datatable.setItem(i, k - 10, newitem)
                    self.datatable.setShowGrid(True)
                except Exception as e:
                    print(e)
            elif number == 3:
                self.Reflect_id_lst = []
                self.datanumber = 3
                try:
                    self.database_class.exhibit_all_info(fenkuai=False)
                    data_info = self.database_class.info_data_base
                    self.datatable.setColumnCount(5)
                    self.datatable.setHorizontalHeaderLabels(['Device', 'Optimal Match', '#Layer',
                                                              'Schottky barrier (eV)', 'I-V curve (Theory)'])
                    if not self.true_id_lst:
                        self.datatable.setRowCount(len(data_info))  # 添加信息
                        for i in range(len(data_info)):
                            for k in range(15, 20):
                                newitem = QtWidgets.QTableWidgetItem(str(list(data_info[i])[k]))
                                newitem.setTextAlignment(5 | 5)
                                self.datatable.setItem(i, k - 15, newitem)
                    else:
                        id_lst = [data_info[i][0] for i in range(len(data_info))]
                        self.datatable.setRowCount(len(self.true_id_lst))
                        for i in range(len(self.true_id_lst)):
                            for k in range(15, 20):
                                newitem = QtWidgets.QTableWidgetItem(str(data_info[id_lst.index(self.true_id_lst[i])][k]))
                                newitem.setTextAlignment(5 | 5)
                                self.datatable.setItem(i, k - 15, newitem)
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
            if len(text) != 0:
                self.true_id_lst = []
                self.row_lst = []
                textlength = len(text)
                self.database_class.exhibit_all_info(fenkuai=False)
                data_info = self.database_class.info_data_base
                for row in range(len(data_info)):
                    num = 0
                    itemtext = data_info[row][2]
                    for letter in text:
                        if letter not in itemtext:
                            break
                        else:
                            num += 1
                    if num == textlength:
                        self.true_id_lst.append(data_info[row][0])
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
            if len(text) != 0:
                it = QtWidgets.QTreeWidgetItemIterator(self.tree)
                searchlist = []
                while it.value():
                    textlength = len(text)
                    num = 0
                    if it.value().text(0) != 'bulk' and it.value().text(0) != 'layer' and it.value().text(0) != 'stack':
                        for letter in text:
                            if (letter not in it.value().text(1)) and (letter not in it.value().text(0)):
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
            it = QtWidgets.QTreeWidgetItemIterator(self.tree)
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
                for child in l:
                    self.clear_painter()
                    self.clear_text()
                    self.dic_Formula_Atoms.pop(child.text(1))  # 对象字典去除该对象
                    if child.parent().childCount() > 1:
                        child.parent().removeChild(child)
                    else:
                        self.tree.takeTopLevelItem(self.tree.indexOfTopLevelItem(child.parent()))
            except Exception as e:
                print(e)
        except Exception as e:
            print(e)

    def normalchoose(self):
        """ To change self.view_widget.num1 to 0, and thus change the PaintEvent of self.view_widget."""
        try:
            global x_rectangle
            try:  # 消除之前画的图
                for obj in x_rectangle.values():
                    self.view_widget.removeItem(obj)
            except Exception as e:
                print(e)
            self.view_widget.num1 = 0  # 改变原来的view_widget 函数,从而重写glviewwidget
        except Exception as e:
            print(e)

    def dragrectangle(self):
        try:
            current = self.tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.tree.selectedItems()) == 1:  # 用户只选择一个的情况
                    self.view_widget.signal_release.connect(self.release)
                    self.view_widget.num1 = 1  # 改变原来的view_widget 函数,从而重写glviewwidget
                else:
                    self.view_widget.num1 = 0  # 改变原来的view_widget 函数,从而重写glviewwidget
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def release(self, initpos, currentpos):
        try:
            global x_rectangle
            try:  # 消除之前画的图
                for obj in x_rectangle.values():
                    self.view_widget.removeItem(obj)
            except Exception as e:
                print(e)
            self.initpos = initpos
            self.currentpos = currentpos
            x_rectangle = {}  # 字典
            for positions in self.positions_lst:  # 获得位置在屏幕上的投影的屏幕坐标
                screentuple = gluProject(positions[0] + self.translate_pos[0], positions[1] + self.translate_pos[1],
                                         positions[2] + self.translate_pos[2])
                if (self.currentpos[0] < screentuple[0] < self.initpos[0] or self.currentpos[0] > screentuple[0] >
                    self.initpos[0]) and \
                        (self.currentpos[1] < self.view_widget.height() - screentuple[1] < self.initpos[1] or
                         self.currentpos[1] > self.view_widget.height() - screentuple[1] > self.initpos[1]):
                    # 两个坐标系的y轴反向
                    # global x_rectangle
                    index = self.positions_lst.index(positions)
                    Atomic_number = self.Atomsobject.get_atomic_numbers()[index]

                    try:
                        md = gl.MeshData.sphere(rows=10, cols=20)
                        x = gl.GLMeshItem(meshdata=md, smooth=True,
                                          color=self.Atoms_color[self.Atoms_num_lst[index]], shader='shaded',
                                          drawEdges=True)
                        x.translate(self.positions_lst[index][0], self.positions_lst[index][1],
                                    self.positions_lst[index][2])
                        x.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])
                        if self.plot_num == 0:
                            x.scale(self.Atom_radii[Atomic_number] + 0.01,
                                    self.Atom_radii[Atomic_number] + 0.01,
                                    self.Atom_radii[Atomic_number] + 0.01)
                        elif self.plot_num == 1:
                            x.scale(self.Atom_radii[Atomic_number]/2 + 0.005,
                                    self.Atom_radii[Atomic_number]/2 + 0.005,
                                    self.Atom_radii[Atomic_number]/2 + 0.005)
                        self.view_widget.addItem(x)
                        x_rectangle[index] = x
                    except Exception as e:
                        print(e)
        except Exception as e:
            print(e)

    def setlayer(self):
        """ To set layer chosen by drag rectangle mode."""
        global x_rectangle
        try:
            if len(x_rectangle) == 0:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please use drag rectangle to select atoms')
            else:
                qm = QtWidgets.QMessageBox.question(self, 'Question',
                                                    'Would yo like to set {} atoms you choose as a layer?'.format(
                                                        str(len(x_rectangle))))
                if qm == QtWidgets.QMessageBox.Yes:
                    max_height = 0
                    min_height = 100
                    origo = ...
                    self.view_widget.num1 = 0
                    for index in x_rectangle.keys():
                        if self.positions_lst[index][2] > max_height:
                            max_height = self.positions_lst[index][2]
                        if self.positions_lst[index][2] < min_height:       # 取最低的原子作为origo
                            min_height = self.positions_lst[index][2]
                            origo = index
                    height = max_height - min_height
                    num = 0
                    if height > 0.01:
                        while True:  # 防止由于误差漏掉一层原子
                            layer = cut(self.Atomsobject, [1, 0, 0], [0, 1, 0], origo=origo,
                                        clength=height + 0.01 * num)
                            if abs(len(layer.get_positions()) - len(x_rectangle)) / len(x_rectangle) < 0.05:
                                break
                            else:
                                num += 1
                    else:           # c=cross(a,b)
                        layer = cut(self.Atomsobject, [1, 0, 0], [0, 1, 0], origo=origo, clength=height + 0.01)

                    self.plot(layer, clear=True, dictionary=True, globalAtomsobject=True, object=True)
                    text3column = self.judgeconductivity(layer)
                    childx = QtWidgets.QTreeWidgetItem(self.tree.currentItem().parent())
                    childx.setText(1, self.dirkey)
                    childx.setText(0, 'layer')
                    childx.setText(3, text3column)
        except Exception as e:
            print(e)

    def movelayer(self):
        """ To move layer."""
        try:
            self.up_down = True
            self.supercell1_cell = self.Atomsobject.get_cell_lengths_and_angles()
            self.view_widget.num1 = 0
            global x_rectangle
            if len(x_rectangle) == 0:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please use drag rectangle to select atoms')
            else:
                atom_num = len(self.positions_lst)
                index_lst = [i for i in range(atom_num)]
                client_lst = list(x_rectangle.keys())
                client_lst.sort()
                orijin1 = copy.deepcopy(self.Atomsobject)  # 用户选的原子
                orijin2 = copy.deepcopy(self.Atomsobject)
                num1 = 0
                num2 = 0
                for index in index_lst:  # 所有的索引
                    if index not in client_lst:
                        orijin1.pop(index - num1)
                        num1 += 1
                    elif index in client_lst:
                        orijin2.pop(index - num2)
                        num2 += 1
                self.crys1_supercell = copy.deepcopy(orijin2)
                self.crys2_supercell = copy.deepcopy(orijin1)  # 用户选的在上层
                self.superplot(self.crys1_supercell, clear=True, layer='bottom')
                self.superplot(self.crys2_supercell, clear=False, layer='top')
                self.keylock = False  # 释放键盘
                xmin = 100
                ymin = 100
                for i in range(len(self.superpositions_lst)):
                    if self.superpositions_lst[i][0] < xmin and self.superpositions_lst[i][1] < ymin:
                        xmin = self.superpositions_lst[i][0]  # superpositions_lst是堆垛在上面的超胞的位置列表, 取出[0,0,h]
                        ymin = self.superpositions_lst[i][1]  # superpositions_lst是堆垛在上面的超胞的位置列表, 取出[0,0,h]
                        self.index0 = i  # 取距离（0， 0）位置最近的点用于定位
        except Exception as e:
            print(e)

    def traversalcut(self):
        try:
            it = QtWidgets.QTreeWidgetItemIterator(self.tree)
            l = []
            s = []
            while it.value():
                v = it.value().text(0)
                if v == "bulk":
                    l.append(self.dic_Formula_Atoms[it.value().text(1)])
                    s.append(it.value())
                it += 1
            for i in range(len(l)):
                try:
                    self.Atomsobject = l[i]
                    self.cut_exf_layer(s[i].parent())
                except Exception as e:
                    print(e)
        except Exception as e:
            print(e)

    def traversalrotate(self):
        try:
            self.traversalrotate_window = traversal_rotate()
            self.traversalrotate_window.signaltraversalstack.connect(self.aftertraversalrotate_window)
        except Exception as e:
            print(e)

    def aftertraversalrotate_window(self, number):
        """Two sides don't know."""
        try:
            self.len_max = self.traversalrotate_window.max_length_var
            self.layerlength = self.traversalrotate_window.layerdistance_var
            self.vacuum_distance = self.traversalrotate_window.vacuumdistance_var
            self.traversalstack_tolerance = self.traversalrotate_window.tolerance_var / 100
            l = []  # 用来储存用户要求堆垛的晶胞,layer
            s = []  # stack
            objl = []  # Project栏的it.value,layer
            objs = []  # Project栏的it.value,stack
            it = QtWidgets.QTreeWidgetItemIterator(self.tree)
            while it.value():
                v = it.value().text(0)
                if number == 1:
                    if v == "layer":
                        l.append(self.dic_Formula_Atoms[it.value().text(1)])
                        objl.append(it.value())
                elif number == 2:  # layer - stack
                    if v == "layer":
                        l.append(self.dic_Formula_Atoms[it.value().text(1)])
                        objl.append(it.value())
                    elif v == 'stack':
                        s.append(self.dic_Formula_Atoms[it.value().text(1)])
                        objs.append(it.value())
                elif number == 3:  # all
                    if v == "layer" or v == 'stack':
                        l.append(self.dic_Formula_Atoms[it.value().text(1)])
                        objl.append(it.value())
                it += 1
            if number == 1 or number == 3:  # layer-layer or All
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
                            self.obj1 = copy.deepcopy(l[i])
                            self.orijincrys1 = copy.deepcopy(self.obj1)  # Atoms1，对象
                            self.obj2 = copy.deepcopy(l[j])
                            fail = self.rotate_all_angles(self.obj1, self.obj2)
                            print("fail = ", fail)
                            if fail != "fail":
                                num += 1
                                self.traversalrotate_window.signal_start_pbar.emit((num + failnum) / zongnum * 100)
                                successmessage = successmessage + l[i].get_chemical_formula(mode='hill') + "-" + l[
                                    j].get_chemical_formula(mode='hill') + "\t" + "Rotate {}°".format(fail / np.pi * 180) + "\n"
                            else:
                                failnum += 1
                                failmessage = failmessage + l[i].get_chemical_formula(mode='hill') \
                                              + "-" + l[j].get_chemical_formula(mode='hill') + "\n"
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
            elif number == 2:  # layer - stack
                zongnum = len(l) * len(s) * 2
                failnum = 0
                num = 0
                failmessage = ""
                successmessage = ""
                for i in range(len(l)):  # layer上叠stack
                    for j in range(len(s)):
                        try:
                            self.it_value1 = objl[i]  # l[i],l[j]的原子对象对应的it.value()
                            self.it_value2 = objs[j]
                            layer = l[i]
                            stack = s[j]
                            self.obj1 = copy.deepcopy(layer)  # 没有转换之前的
                            self.obj2 = copy.deepcopy(stack)  # 没有转换之前的
                            self.orijincrys1 = copy.deepcopy(self.obj1)
                            fail = self.rotate_all_angles(self.obj1, self.obj2)
                            if fail != "fail":
                                num += 1
                                self.traversalrotate_window.signal_start_pbar.emit((num + failnum) / zongnum * 100)
                                successmessage = successmessage + layer.get_chemical_formula(
                                    mode='hill') + "(layer)" \
                                                 + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" + "\t" +\
                                                 "Rotate {}°".format(fail / np.pi * 180) + "\n"
                            else:
                                failnum += 1
                                failmessage = failmessage + layer.get_chemical_formula(mode='hill') + "(layer)" \
                                              + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" + "\n"

                        except Exception as e:
                            failnum += 1
                            failmessage = failmessage + layer.get_chemical_formula(mode='hill') + "(layer)" \
                                          + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" + "\n"
                            print(e)
                for i in range(len(s)):  # stack 上叠layer
                    for j in range(len(l)):
                        try:
                            layer = l[j]
                            stack = s[i]
                            self.it_value1 = objs[i]
                            self.it_value2 = objl[j]  # l[i],l[j]的原子对象对应的it.value()
                            self.obj1 = copy.deepcopy(stack)
                            self.obj2 = copy.deepcopy(layer)
                            self.orijincrys1 = copy.deepcopy(stack)  # 没有转换之前的
                            fail = self.rotate_all_angles(self.obj1, self.obj2)

                            if fail != "fail":
                                num += 1
                                self.traversalrotate_window.signal_start_pbar.emit((num + failnum) / zongnum * 100)
                                successmessage = successmessage + layer.get_chemical_formula(mode='hill') + "(layer)" \
                                                 + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" + "\t" +\
                                                 "Rotate {}°".format(fail / np.pi * 180)+ "\n"
                            else:
                                failnum += 1
                                failmessage = failmessage + layer.get_chemical_formula(mode='hill') + "(layer)" \
                                              + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" + "\n"
                        except Exception as e:
                            failnum += 1
                            failmessage = failmessage + layer.get_chemical_formula(mode='hill') + "(layer)" \
                                          + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" + "\n"
                            print(e)

        except Exception as e:
            print(e)

    def rotate_all_angles(self, obj1, obj2):
        try:
            rotate_obj1 = copy.deepcopy(obj1)
            rotate_obj2 = copy.deepcopy(obj2)
            cell1 = list(rotate_obj1.get_cell())
            vectora1 = np.array([cell1[0][0], cell1[0][1]])
            vectorb1 = np.array([cell1[1][0], cell1[1][1]])
            cell2 = list(rotate_obj2.get_cell())
            vectora2 = np.array([cell2[0][0], cell2[0][1]])
            vectorb2 = np.array([cell2[1][0], cell2[1][1]])
            self.pos_lst_orijin = []
            self.pos_lst_orijin_up = []
            for i in range(- 10, 10):
                for j in range(-10,  10):
                    pos = vectora2 * i + vectorb2 * j
                    if np.linalg.norm(pos) <= self.len_max:
                        self.pos_lst_orijin_up.append(pos)
            for i in range(-10, 10):
                for j in range(0, 10):
                    pos = vectora1 * i + vectorb1 * j
                    if np.linalg.norm(pos) <= self.len_max:
                        self.pos_lst_orijin.append(pos)
            theta = 0
            unit = np.pi / 180
            compare_square = 10000
            a_down = ...
            b_down = ...
            theta_true = ...
            cell_2_par = rotate_obj2.get_cell_lengths_and_angles().tolist()
            if abs(cell_2_par[0] - cell_2_par[1]) < 10e-6:
                if abs(cell_2_par[5] - 60) < 1e-6 or abs(cell_2_par[5] - 120) < 1e-6:
                    division = 6
                elif abs(cell_2_par[5] - 90) < 1e-6:
                    division = 4
            else:
                division = 1
            while theta < np.pi/division:
                rotate_matrix = np.array(
                    [[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
                self.new_spot_up_lst = []
                for i in range(len(self.pos_lst_orijin_up)):
                    new_pos_up = self.pos_lst_orijin_up[i] @ rotate_matrix
                    self.new_spot_up_lst.append(new_pos_up)
                a, b, min_square = self.find_min()
                if a is not None and b is not None and min_square is not None:
                    if min_square < compare_square:
                        compare_square = min_square
                        a_down = a
                        b_down = b
                        theta_true = theta
                theta += unit
            if compare_square == 10000:
                return "fail"
            else:
                cell1 = rotate_obj1.get_cell()
                cell1_a = [cell1[0][0], cell1[0][1]]
                cell1_b = [cell1[1][0], cell1[1][1]]
                A = np.array([cell1_a, cell1_b]).T
                r1a = np.linalg.solve(A, np.array(a_down))  # 求解线性方程组，简单,直角坐标系下----用晶胞坐标系表示
                r1a = np.append(r1a, 0)
                r1b = np.linalg.solve(A, np.array(b_down))
                r1b = np.append(r1b, 0)
                layer_down = cut(rotate_obj1, r1a, r1b, [0, 0, 1], origo=[0, 0, 0])
                self.layer_down = self.deal_with_rotate(layer_down)
                # layer_up
                cell2_before_rotate = rotate_obj2.get_cell()
                cell2_a = [cell2_before_rotate[0][0], cell2_before_rotate[0][1]]
                cell2_b = [cell2_before_rotate[1][0], cell2_before_rotate[1][1]]
                # degree = degree / 180 * np.pi
                rotate_matrix = np.array(
                    [[np.cos(theta_true), np.sin(theta_true)], [-np.sin(theta_true), np.cos(theta_true)]])
                new_cell2_a = cell2_a @ rotate_matrix
                new_cell2_b = cell2_b @ rotate_matrix
                A = np.array([new_cell2_a, new_cell2_b]).T
                r2a = np.linalg.solve(A, np.array(a_down))
                r2a = np.append(r2a, 0)
                r2b = np.linalg.solve(A, np.array(b_down))  # 求解线性方程组，简单,直角坐标系下----用晶胞坐标系表示
                r2b = np.append(r2b, 0)
                layer_up = cut(rotate_obj2, r2a, r2b, [0, 0, 1], origo=[0, 0, 0])
                self.layer_up = self.deal_with_rotate(layer_up)
                # start stacking
                layer_down_cell = self.layer_down.get_cell_lengths_and_angles()
                self.layer_up.translate(np.array([0, 0, self.layerlength + self.layer_down.get_cell()[2][2]]))
                layer_down_cell[2] += self.layer_up.get_cell_lengths_and_angles()[2] + self.layerlength + 0.01
                self.layer_down.extend(self.layer_up)
                self.layer_down.set_cell(layer_down_cell)
                # add_vacuum _layer
                c_length = self.layer_down.get_cell_lengths_and_angles()[2]
                cell_par = self.layer_down.get_cell()
                cell_par[2] *= (c_length + self.vacuum_distance) / c_length
                self.layer_down.set_cell(cell_par)
                self.layer_down.translate(cell_par[2] / (c_length + self.vacuum_distance) * self.vacuum_distance / 2)
                # plot
                self.plot(self.layer_down, plot=False, object=False, tab3=False, clear=True, globalAtomsobject=False, dictionary=True)
                key = list(self.dic_Formula_Atoms.keys())[
                    list(self.dic_Formula_Atoms.values()).index(self.orijincrys1)]
                it = QtWidgets.QTreeWidgetItemIterator(self.tree)
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


    def find_min(self):
        try:
            peidui_lst = []
            for i in range(len(self.new_spot_up_lst)):
                for j in range(len(self.pos_lst_orijin)):
                    if (self.pos_lst_orijin[j][1] > self.traversalstack_tolerance * self.len_max and self.new_spot_up_lst[i][1] <= 0) \
                            or self.new_spot_up_lst[i][1] < -self.traversalstack_tolerance * self.len_max:
                        break
                    else:
                        pos1 = self.pos_lst_orijin[j]
                        pos2 = self.new_spot_up_lst[i]
                        dis = np.linalg.norm(pos1 - pos2)
                        if dis / np.linalg.norm(pos1) < self.traversalstack_tolerance:
                            peidui_lst.append(pos1)
            if len(peidui_lst) < 2:
                return None, None, None
            else:
                min = 1000000
                a = ...
                b = ...
                for i in range(len(peidui_lst)):
                    for j in range(i + 1, len(peidui_lst)):
                        n_vector = np.cross(peidui_lst[i], peidui_lst[j])
                        square = abs(n_vector)
                        if 1e-6 < square < min:
                            if n_vector > 0:     # 保证右手系
                                a = peidui_lst[i]
                                b = peidui_lst[j]
                            else:
                                b = peidui_lst[i]
                                a = peidui_lst[j]
                            min = square
                if min == 1000000:
                    return None, None, None
                else:
                    return a, b, min
        except Exception as e:
            print(e)



    def traversalstack1(self):
        try:
            self.traversalstack_window = traversal_stack()
            self.traversalstack_window.side_known = 1
            self.traversalstack_window.signaltraversalstack.connect(self.aftertraversalstack_window)
        except Exception as e:
            print(e)

    def traversalstack2(self):
        """ To traversal stack multiple crys. """
        try:
            self.traversalstack_window = traversal_stack()
            self.traversalstack_window.side_known = 2
            self.traversalstack_window.signaltraversalstack.connect(self.aftertraversalstack_window)
        except Exception as e:
            print(e)


    def aftertraversalstack_window(self, number):
        """ To traversal stack multiple crys. """
        try:
            self.layerlength = self.traversalstack_window.layerdistance_var
            self.vacuum_distance = self.traversalstack_window.vacuumdistance_var
            self.traversalstack_tolerance = self.traversalstack_window.tolerance_var / 100
            l = []  # 用来储存用户要求堆垛的晶胞,layer
            s = []  # stack
            objl = []  # Project栏的it.value,layer
            objs = []  # Project栏的it.value,stack
            it = QtWidgets.QTreeWidgetItemIterator(self.tree)
            while it.value():
                v = it.value().text(0)
                if number == 1:
                    if v == "layer":
                        l.append(self.dic_Formula_Atoms[it.value().text(1)])
                        objl.append(it.value())
                elif number == 2:  # layer - stack
                    if v == "layer":
                        l.append(self.dic_Formula_Atoms[it.value().text(1)])
                        objl.append(it.value())
                    elif v == 'stack':
                        s.append(self.dic_Formula_Atoms[it.value().text(1)])
                        objs.append(it.value())
                elif number == 3:  # all
                    if v == "layer" or v == 'stack':
                        l.append(self.dic_Formula_Atoms[it.value().text(1)])
                        objl.append(it.value())
                it += 1
            if self.traversalstack_window.side_known == 2:      # two-side know
                if number == 1 or number == 3:  # layer-layer or All
                    zongnum = len(l) * (len(l) - 1) / 2
                    failnum = 0
                    num = 0
                    failmessage = ""
                    successmessage = ""
                    for i in range(len(l)):
                        for j in range(i + 1, len(l)):
                            try:
                                self.compare_cell(l[i], l[j], traversalstack=True)
                                self.orijincrys1 = copy.deepcopy(l[i])         # 没有转换之前的
                                self.orijincrys2 = copy.deepcopy(l[j])         # 没有转换之前的
                                self.it_value1 = objl[i]         # l[i],l[j]的原子对象对应的it.value()
                                self.it_value2 = objl[j]
                                if self.crys1 is not None and self.crys2 is not None:
                                    self.layer_transform(traversalstack=True)
                                    num += 1
                                    self.traversalstack_window.signal_start_pbar.emit((num+failnum) / zongnum * 100)
                                    successmessage = successmessage + l[i].get_chemical_formula(mode='hill') + "-" + l[j].get_chemical_formula(mode='hill') + "\n"
                                else:
                                    failnum += 1
                                    failmessage = failmessage + l[i].get_chemical_formula(mode='hill')  \
                                                  + "-" + l[j].get_chemical_formula(mode='hill')  + "\n"
                            except Exception as e:
                                failnum += 1
                                failmessage = failmessage + l[i].get_chemical_formula(mode='hill') + "-" + l[j].get_chemical_formula(mode='hill') + "\n"
                                print(e)
                    message1 = "There are {} stack situations in total".format(int(zongnum)) \
                               + "\n" + "{} succeeded".format(num) + "\n" + "{} failed".format(failnum)
                    message2 = "\n" + "-" * 20 + "\n" + "failed situations:" + "\n" + failmessage
                    message3 = "\n" + "-" * 20 + "\n" + "succeeded situations:" + "\n" + successmessage
                    message = message1 + message2 + message3
                    self.traversalstack_window.Text.setText(message)
                elif number == 2:  # layer - stack
                    zongnum = len(l) * len(s) * 2
                    failnum = 0
                    num = 0
                    failmessage = ""
                    successmessage = ""
                    for i in range(len(l)):  # layer上叠stack
                        for j in range(len(s)):
                            try:
                                layer = l[i]
                                stack = s[j]
                                self.orijincrys1 = copy.deepcopy(layer)  # 没有转换之前的
                                self.orijincrys2 = copy.deepcopy(stack)  # 没有转换之前的
                                self.it_value1 = objl[i]  # l[i],l[j]的原子对象对应的it.value()
                                self.it_value2 = objs[j]
                                self.compare_cell(layer, stack, traversalstack=True)
                                if self.crys1 is not None and self.crys2 is not None:
                                    self.layer_transform(traversalstack=True)
                                    num += 1
                                    self.traversalstack_window.signal_start_pbar.emit((num+failnum) / zongnum * 100)
                                    successmessage = successmessage + layer.get_chemical_formula(mode='hill') + "(layer)" \
                                                  + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" + "\n"
                                else:
                                    failnum += 1
                                    failmessage = failmessage + layer.get_chemical_formula(mode='hill') + "(layer)" \
                                                  + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" + "\n"

                            except Exception as e:
                                failnum += 1
                                failmessage = failmessage + layer.get_chemical_formula(mode='hill') + "(layer)" \
                                              + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" +"\n"
                                print(e)
                    for i in range(len(s)):  # stack 上叠layer
                        for j in range(len(l)):
                            try:
                                layer = l[j]
                                stack = s[i]
                                self.orijincrys1 = copy.deepcopy(stack)  # 没有转换之前的
                                self.orijincrys2 = copy.deepcopy(layer)  # 没有转换之前的
                                self.it_value1 = objs[i]
                                self.it_value2 = objl[j]  # l[i],l[j]的原子对象对应的it.value()
                                self.compare_cell(stack, layer, traversalstack=True)
                                if self.crys1 is not None and self.crys2 is not None:
                                    self.layer_transform(traversalstack=True)
                                    num += 1
                                    self.traversalstack_window.signal_start_pbar.emit((num + failnum) / zongnum * 100)
                                    successmessage = successmessage + layer.get_chemical_formula(mode='hill') + "(layer)" \
                                                 + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" + "\n"
                                else:
                                    failnum += 1
                                    failmessage = failmessage + layer.get_chemical_formula(mode='hill') + "(layer)" \
                                                  + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" + "\n"
                            except Exception as e:
                                failnum += 1
                                failmessage = failmessage + layer.get_chemical_formula(mode='hill') + "(layer)" \
                                              + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" + "\n"
                                print(e)
                    message1 = "There are {} stack situations in total".format(int(zongnum)) \
                               + "\n" + "{} succeeded".format(num) + "\n" + "{} failed".format(failnum)
                    message2 = "\n" + "-"*20 + "\n" + "failed situations:" + "\n" + failmessage
                    message3 = "\n" + "-"*20 + "\n" + "succeeded situations:" + "\n" + successmessage
                    message = message1 + message2 + message3
                    self.traversalstack_window.Text.setText(message)
            else:    # one-sid-known
                if number == 1 or number == 3:  # layer-layer or All
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
                                self.obj1 = copy.deepcopy(l[i])
                                self.orijincrys1 = copy.deepcopy(self.obj1)  # Atoms1，对象
                                self.obj2 = copy.deepcopy(l[j])
                                fail = self.stack_one(self.traversalstack_tolerance, self.layerlength,
                                                           self.vacuum_distance, traversal=True)
                                print("fail = ", fail)
                                if fail != "fail":
                                    num += 1
                                    self.traversalstack_window.signal_start_pbar.emit((num + failnum) / zongnum * 100)
                                    successmessage = successmessage + l[i].get_chemical_formula(mode='hill') + "-" + l[
                                        j].get_chemical_formula(mode='hill') + "\n"
                                else:
                                    failnum += 1
                                    failmessage = failmessage + l[i].get_chemical_formula(mode='hill') \
                                                  + "-" + l[j].get_chemical_formula(mode='hill') + "\n"
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
                    self.traversalstack_window.Text.setText(message)
                elif number == 2:  # layer - stack
                    zongnum = len(l) * len(s) * 2
                    failnum = 0
                    num = 0
                    failmessage = ""
                    successmessage = ""
                    for i in range(len(l)):  # layer上叠stack
                        for j in range(len(s)):
                            try:
                                self.it_value1 = objl[i]  # l[i],l[j]的原子对象对应的it.value()
                                self.it_value2 = objs[j]
                                layer = l[i]
                                stack = s[j]
                                self.obj1 = copy.deepcopy(layer)  # 没有转换之前的
                                self.obj2 = copy.deepcopy(stack)  # 没有转换之前的
                                self.orijincrys1 = copy.deepcopy(self.obj1)
                                fail = self.stack_one(self.traversalstack_tolerance, self.layerlength,
                                                           self.vacuum_distance, traversal=True)
                                if fail != "fail":
                                    num += 1
                                    self.traversalstack_window.signal_start_pbar.emit((num + failnum) / zongnum * 100)
                                    successmessage = successmessage + layer.get_chemical_formula(
                                        mode='hill') + "(layer)" \
                                                     + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" + "\n"
                                else:
                                    failnum += 1
                                    failmessage = failmessage + layer.get_chemical_formula(mode='hill') + "(layer)" \
                                                  + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" + "\n"

                            except Exception as e:
                                failnum += 1
                                failmessage = failmessage + layer.get_chemical_formula(mode='hill') + "(layer)" \
                                              + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" + "\n"
                                print(e)
                    for i in range(len(s)):  # stack 上叠layer
                        for j in range(len(l)):
                            try:
                                layer = l[j]
                                stack = s[i]
                                self.it_value1 = objs[i]
                                self.it_value2 = objl[j]  # l[i],l[j]的原子对象对应的it.value()
                                self.obj1 = copy.deepcopy(stack)
                                self.obj2 = copy.deepcopy(layer)
                                self.orijincrys1 = copy.deepcopy(stack)  # 没有转换之前的
                                fail = self.stack_one(self.traversalstack_tolerance, self.layerlength,
                                                           self.vacuum_distance, traversal=True)

                                if fail != "fail":
                                    self.layer_transform(traversalstack=True)
                                    num += 1
                                    self.traversalstack_window.signal_start_pbar.emit((num + failnum) / zongnum * 100)
                                    successmessage = successmessage + layer.get_chemical_formula(mode='hill') + "(layer)" \
                                                 + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" + "\n"
                                else:
                                    failnum += 1
                                    failmessage = failmessage + layer.get_chemical_formula(mode='hill') + "(layer)" \
                                                  + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" + "\n"
                            except Exception as e:
                                failnum += 1
                                failmessage = failmessage + layer.get_chemical_formula(mode='hill') + "(layer)" \
                                              + "-" + stack.get_chemical_formula(mode='hill') + "(stack)" + "\n"
                                print(e)
                    message1 = "There are {} stack situations in total".format(int(zongnum)) \
                               + "\n" + "{} succeeded".format(num) + "\n" + "{} failed".format(failnum)
                    message2 = "\n" + "-"*20 + "\n" + "failed situations:" + "\n" + failmessage
                    message3 = "\n" + "-"*20 + "\n" + "succeeded situations:" + "\n" + successmessage
                    message = message1 + message2 + message3
                    self.traversalstack_window.Text.setText(message)
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

    def handleItemClick(self, item):
        """  To handle Atom clicked in the objectbox. """
        current = self.tree.currentItem()
        try:
            if current is not None and current.parent() is not None:
                if len(self.tree.selectedItems()) == 1:  # 用户只选择一个的情况
                    try:
                        global x
                        self.view_widget.removeItem(x)  # 移走双击选中的原子
                    except Exception as e:
                        print(e)
                    global index_lst
                    index_lst = []
                    for item in self.tablewidget.selectedItems():
                        index_lst.append(item.row())
                    index_lst = [index_lst[i] for i in range(len(index_lst)) if i % 3 == 1]
                    self.tab.setCurrentWidget(self.tab2)
                    try:
                        global obj_lst_table
                        for x in obj_lst_table:
                            self.view_widget.removeItem(x)
                    except Exception as e:
                        print(e)
                    try:
                        obj_lst_table = []
                        text = ''
                        for i in range(len(index_lst)):
                            index = index_lst[i]
                            Atomic_number = self.Atomsobject.get_atomic_numbers()[index]
                            md = gl.MeshData.sphere(rows=10, cols=20)
                            x = gl.GLMeshItem(meshdata=md, smooth=True,
                                              color=self.Atoms_color[self.Atoms_num_lst[index]],
                                              shader='shaded', drawEdges=True)
                            x.translate(self.positions_lst[index][0], self.positions_lst[index][1],
                                        self.positions_lst[index][2])
                            x.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])
                            if self.plot_num == 0:
                                x.scale(self.Atom_radii[Atomic_number] + 0.01,
                                        self.Atom_radii[Atomic_number] + 0.01,
                                        self.Atom_radii[Atomic_number] + 0.01)
                            elif self.plot_num == 1:
                                x.scale(self.Atom_radii[Atomic_number]/2 + 0.005,
                                        self.Atom_radii[Atomic_number]/2 + 0.005,
                                        self.Atom_radii[Atomic_number]/2 + 0.005)

                            self.view_widget.addItem(x)
                            j = sf.crys_data.atom_formula[Atomic_number] + str(index + 1) + "\n" + "Atomic number:" + str(
                                Atomic_number) + "\n" + "position:" + str(self.positions_lst[index]) + '\n' + "-" * 32 + '\n'
                            text += j
                            obj_lst_table.append(x)
                        self.text_widgt2.setText(text)
                    except Exception as e:
                        print(e)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:  # parent()是self.tree的情况
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
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

    def viewCoordinate_system(self):  # 笛卡尔坐标（直角坐标）
        """ Whether plot coordinate system or not."""
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
        except Exception as e:
            print(e)

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
        except Exception as e:
            print(e)

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
                if tantheta < 0.0 and abs(theta - 0) < 0.0001:
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
        except Exception as e:
            print(e)

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
                if tantheta < 0.0 and abs(theta - 0) < 0.0001:
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
        except Exception as e:
            print(e)

    # 创建新文件
    def newfile(self):
        """ To creat new cif file"""
        self.new_window = New_file()
        self.new_window.signal_determine.connect(self.opennewfile)
        self.new_window.show_win()

    # 打开新建文件
    def opennewfile(self, str):
        """ To open file."""
        self.fileName_choose = str
        self.openfile(new=True)

    def doubleclickedontree(self):
        """ To react to the MouseDbCLICKED on tree."""
        global current
        try:
            if current is not None and current.parent() is not None:
                self.setTexttree_window = tree_set_text()
                self.setTexttree_window.signal_edit_text.connect(self.edit_tree_text)
                self.setTexttree_window.signal_edit_text1.connect(self.edit_tree_text_conduct)

            elif current == None:  # parent()是self.tree的情况
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def edit_tree_text_conduct(self, text):
        """ To edit the text of tree. (Conductor, Semiconductor, Insulator)"""
        current = self.tree.currentItem()
        try:
            current.setText(3, text)
        except Exception as e:
            print(e)

    def edit_tree_text(self, text):
        """ To edit the text of tree. (bulk, layer, stack)"""
        current = self.tree.currentItem()
        try:
            current.setText(0, text)
            if text == 'bulk' or text == 'stack':      # 如果是bulk或stack没有conductivity
                current.setText(3, "")
            elif text == 'layer':
                if current.text(3) == '':
                    Text3column = self.judgeconductivity(self.dic_Formula_Atoms[current.text(1)])
                    current.setText(3, Text3column)
        except Exception as e:
            print(e)

    # 点击树
    def onTreeClicked(self):
        """ Reaction of click on self.treewidget """
        try:
            if self.view_widget.num1 == 0:
                try:  # 移去drag rectangle选中的Item
                    global x_rectangle, current
                    for obj in x_rectangle.values():
                        self.view_widget.removeItem(obj)
                    x_rectangle = {}
                except Exception as e:
                    print(e)
                current = self.tree.currentItem()
                self.plot(self.dic_Formula_Atoms[current.text(1)], dictionary=False, globalAtomsobject=True)  # 不加对象

            elif self.view_widget.num1 == 1:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose drag mode')
        except Exception as e:
            print(e)

    # 右击出现的菜单
    def rightMenuShow(self):
        """ Mouse right button click on self.view_widget to show a menubar."""
        try:
            if self.view_widget.num1 == 0:
                rightMenu = QtWidgets.QMenu(self.menuBar)
                self.add_cell_Action = QtWidgets.QAction(self)
                self.add_cell_Action.setText("Create supercell")
                self.add_cell_Action.triggered.connect(self.plot_add_cell)
                rightMenu.addAction(self.add_cell_Action)
                self.export_cif_Action = QtWidgets.QAction(self)
                self.export_cif_Action.setText("Export cif")
                self.export_cif_Action.triggered.connect(self.export_cif)
                rightMenu.addAction(self.export_cif_Action)
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
                self.delete = QtWidgets.QAction(self)
                self.delete.setText("Delete")
                self.delete.triggered.connect(self.deleteobject)
                self.delete.setEnabled(True)
                rightMenu.addAction(self.delete)
                self.add_vacuum_layer = QtWidgets.QAction(self)
                self.add_vacuum_layer.setText("Add vacuum layer")
                self.add_vacuum_layer.triggered.connect(self.addvacuumlayer)
                self.add_vacuum_layer.setEnabled(True)
                rightMenu.addAction(self.add_vacuum_layer)
                self.changegama = QtWidgets.QAction(self)
                self.changegama.setText("Change γ degree")
                self.changegama.triggered.connect(self.change_gama)
                self.changegama.setEnabled(True)
                rightMenu.addAction(self.changegama)

                rightMenu.exec_(QtGui.QCursor.pos())
        except Exception as e:
            print(e)

    def change_gama(self):
        """ To change gama according to the clients."""
        try:
            self.change_gama_window = Change_gama_window()
            self.change_gama_window.signal_gama_error.connect(self.trans_gama)
        except Exception as e:
            print(e)

    def trans_gama(self, gama_error):
        """ To change gamma according to object-gama and gama_error."""
        current = self.tree.currentItem()
        try:
            if current is not None and current.parent() is not None:
                if len(self.tree.selectedItems()) == 1:  # 用户只选择一个的情况
                    cell_vector = self.Atomsobject.get_cell()
                    a_vector = copy.deepcopy(cell_vector[0])
                    b_vector = copy.deepcopy(cell_vector[1])
                    newvector = self.transcrys(a_vector, b_vector, gama_error[0], gama_error[1], 0)
                    A = cell_vector.T
                    new_cut_vector = []
                    for vector in newvector:
                        b = vector.T
                        r = np.linalg.solve(A, b)  # 求解线性方程组，简单,直角坐标系下----用晶胞坐标系表示
                        new_cut_vector.append(r)
                    Atomsobject = copy.deepcopy(self.Atomsobject)
                    Atomsobject = cut(Atomsobject, a=new_cut_vector[0], b=new_cut_vector[1], c=[0, 0, 1])
                    cell_vector[0] = newvector[0]
                    cell_vector[1] = newvector[1]
                    Atomsobject.set_cell(cell_vector)
                    Atomsobject = self.set_rotate(Atomsobject)
                    self.plot(Atomsobject, clear=True, globalAtomsobject=False, dictionary=True)
                    childx = QtWidgets.QTreeWidgetItem(self.tree.currentItem().parent())
                    childx.setText(1, self.dirkey)
                    childx.setText(0, current.text(0))
                    childx.setText(3, current.text(3))
                    childx.setText(2, current.text(2))
                    self.change_gama_window.close()
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:  # parent()是self.tree的情况
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def set_rotate(self, Atoms):
        """ To set the layer lay on the x-y plane."""
        current = self.tree.currentItem()
        try:
            if current is not None and current.parent() is not None:
                if len(self.tree.selectedItems()) == 1:  # 用户只选择一个的情况
                    Atomsobject = copy.deepcopy(Atoms)
                    return self.deal_with_rotate(Atomsobject)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:  # parent()是self.tree的情况
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def deal_with_rotate(self, Atomsobject):
        try:
            cell_par = Atomsobject.get_cell_lengths_and_angles()
            pos_lst = Atomsobject.get_positions()
            A = Atomsobject.get_cell().T
            cell_coordinat_system_par = []  # 晶胞坐标系下的原子坐标
            for pos in pos_lst:
                b = pos.T
                r = np.linalg.solve(A, b)  # 求解线性方程组，简单,直角坐标系下----用晶胞坐标系表示
                cell_coordinat_system_par.append(r.T)
            a = [cell_par[0], 0, 0]
            b = [cell_par[1] * np.cos(cell_par[5] / 180 * np.pi),
                 cell_par[1] * np.sin(cell_par[5] / 180 * np.pi), 0]
            c = [cell_par[2] * np.cos(cell_par[4] / 180 * np.pi),
                 cell_par[2] * np.cos(cell_par[3] / 180 * np.pi),
                 cell_par[2] * np.sqrt(
                     1 - np.cos(cell_par[4] / 180 * np.pi) ** 2 - np.cos(cell_par[3] / 180 * np.pi) ** 2)]
            Atomsobject.set_cell([a, b, c])
            s = []
            for pos in cell_coordinat_system_par:
                new_pos = (np.array(pos[0]) * a + np.array(pos[1]) * b + np.array(pos[2]) * c)
                s.append(new_pos)
            Atomsobject.set_positions(s)
            return Atomsobject
        except Exception as e:
            print(e)


    def addvacuumlayer(self):
        """ To add vacuum layer."""
        current = self.tree.currentItem()
        try:
            if current is not None and current.parent() is not None:
                if len(self.tree.selectedItems()) == 1:  # 用户只选择一个的情况
                    self.vacuum_window = add_vacuum_layer_window()
                    self.vacuum_window.sinal_vacuum.connect(self.set_vacuum_layer)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:  # parent()是self.tree的情况
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
            Atomsobject = copy.deepcopy(self.Atomsobject)
            c_length = Atomsobject.get_cell_lengths_and_angles()[2]
            cell_par = Atomsobject.get_cell()
            cell_par[2] *= (c_length + vacuum_dis) / c_length
            Atomsobject.set_cell(cell_par)
            Atomsobject.translate(cell_par[2] / (c_length + vacuum_dis) * vacuum_dis / 2)
            if addproject == True:
                self.plot(Atomsobject, clear=True, globalAtomsobject=False, dictionary=True)
                childx = QtWidgets.QTreeWidgetItem(self.tree.currentItem().parent())
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
        global current
        try:
            if len(self.tree.selectedItems()) == 1:  # 用户只选择一个的情况
                if current is not None and current.parent() is not None:  # parent()不是self.tree的情况
                    self.clear_painter()
                    self.clear_text()
                    self.dic_Formula_Atoms.pop(current.text(1))  # 对象字典去除该对象
                    print("after delete = ", self.dic_Formula_Atoms)
                    if current.parent().childCount() > 1:
                        current.parent().removeChild(current)
                        current = None
                    else:
                        self.tree.takeTopLevelItem(self.tree.indexOfTopLevelItem(current.parent()))
                        current = None
                elif current == None:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
                elif current.parent() == None:  # parent()是self.tree的情况
                    # if not current.text(0):
                    for i in range(current.childCount()):
                        self.dic_Formula_Atoms.pop(current.child(i).text(1))
                    self.tree.takeTopLevelItem(self.tree.indexOfTopLevelItem(current))
                    self.clear_text()
                    self.clear_painter()
                    current = None

            elif len(self.tree.selectedItems()) > 1:
                for item in self.tree.selectedItems():
                    if item is not None and item.parent() is not None:  # parent()不是self.tree的情况
                        self.clear_painter()
                        self.clear_text()
                        self.dic_Formula_Atoms.pop(item.text(1))  # 对象字典去除该对象
                        if item.parent().childCount() > 1:
                            item.parent().removeChild(item)
                        else:
                            self.tree.takeTopLevelItem(self.tree.indexOfTopLevelItem(item.parent()))
                    elif item.parent() == None:  # parent()是self.tree的情况
                        for i in range(item.childCount()):
                            self.dic_Formula_Atoms.pop(item.child(i).text(1))
                        self.tree.takeTopLevelItem(self.tree.indexOfTopLevelItem(item))
                        self.clear_text()
                        self.clear_painter()
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)




    def remove_atom(self):
        """To remove the Atoms which is chosen by the client."""
        try:
            current = self.tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.tree.selectedItems()) == 1:  # 用户只选择一个的情况
                    Atomsobject = copy.deepcopy(self.Atomsobject)
                    Atomsobject.pop(self.remove_atom_index_num)
                    self.remove_atom_index_num = None
                    self.plot(Atomsobject, clear=True, dictionary=True)
                    child1 = QtWidgets.QTreeWidgetItem(current.parent())
                    child1.setText(0, current.text(0))
                    child1.setText(1, self.dirkey)
                    child1.setText(2, current.text(2))
                    child1.setText(3, current.text(3))
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def replace_atom(self):
        """To replace the Atoms which is chosen by the client."""
        try:
            current = self.tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.tree.selectedItems()) == 1:  # 用户只选择一个的情况
                    self.replace_window = Replace_Atom()
                    self.replace_window.signal_atomic_number.connect(self.afterreplace_window)
                else:
                    QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose in Projectbox')
        except Exception as e:
            print(e)

    def afterreplace_window(self, atomicnumber):
        """ To replace an Atom."""
        try:
            current = self.tree.currentItem()
            Atomsobject = copy.deepcopy(self.Atomsobject)
            Atomsobject.numbers[self.remove_atom_index_num] = atomicnumber
            self.remove_atom_index_num = None
            self.plot(Atomsobject, clear=True, dictionary=True)
            child1 = QtWidgets.QTreeWidgetItem(current.parent())
            child1.setText(0, current.text(0))
            child1.setText(1, self.dirkey)

        except Exception as e:
            print(e)

    # 导出cif文件

    def export_cif(self):  # 右击菜单显示的export cif
        """To export one cif file which client chooses in the self.tree widget."""
        try:
            current = self.tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.tree.selectedItems()) == 1:  # 用户只选择一个的情况
                    output_file = QFileDialog.getSaveFileName(self, "Save as cif", "{}.cif".format(
                        self.Atomsobject.get_chemical_formula(mode='hill')))
                    f = AseAtomsAdaptor.get_structure(self.Atomsobject)
                    print(AseAtomsAdaptor.get_structure(self.Atomsobject))
                    f1 = str(f) + '\n'
                    j_pv_lst = re.findall('abc\s\s\s:(.*?)\n', f1)[0]  # abc   :  19.257300  19.569178  21.133988
                    j1_pv_lst = j_pv_lst.split(' ')  # abc   :   6.419100   6.523059   7.044663
                    while '' in j1_pv_lst:
                        j1_pv_lst.remove('')
                    par_lst_matrix = self.Atomsobject.get_cell()

                    # 物质种类（比如：Lu2 Al4）
                    y1 = re.findall('Full\sFormula\s(.*?)\n', f1)[0]
                    y1_lst = y1.lstrip('(').rstrip(')')
                    material = re.findall('Full\sFormula\s(.*?)\n', f1)[0].lstrip('(').rstrip(')')
                    elements = material.split(' ')
                    print(elements)  # (Re4 S8)\(Re108 S216)
                    zmb_lst = [chr(i) for i in range(97, 123)] + [chr(i) for i in range(65, 91)]
                    szb_lst = [str(i) for i in range(0, 10)]
                    element_lst = []
                    number_lst = []
                    c = []
                    for element in elements:  # 这个循环之后，得到原子与对应的原子数两个列表
                        letter_lst = list(element)
                        symbol_lst = []
                        num_lst = []
                        element_lst1 = []
                        number_lst1 = []
                        # print(letter_lst)
                        for i in range(len(letter_lst)):
                            if letter_lst[i] in szb_lst:
                                num_lst.append(letter_lst[i])
                                # print(num_lst)
                            if letter_lst[i] in zmb_lst:
                                symbol_lst.append(letter_lst[i])
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
                    print(element_lst)
                    print(number_lst)
                    par_lst_species = []  # 用于Cifwrite参数(species)的
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
            current = self.tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.tree.selectedItems()) == 1:  # 用户只选择一个的情况
                    self.plot_cell = False
                    self.c = add_cell_plot()
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
        """To create crys after adding cell."""
        try:
            self.root = QtWidgets.QTreeWidgetItem(self.tree)
            self.plot(self.Atomsobject, dictionary=True)
            self.root.setText(0, self.dirkey)
            child = QtWidgets.QTreeWidgetItem(self.root)
            child.setText(0, "bulk")
            child.setText(1, self.dirkey)
            self.c.close()
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
                after_repeat = self.Atomsobject.repeat((self.repeat_a, self.repeat_b, self.repeat_c))
                self.plot(after_repeat, dictionary=False, globalAtomsobject=True, object=False)  # 中间变量，不用加入object
                self.plot_cell = False
                self.c.numberSpinBoxa.setValue(1)
                self.c.numberSpinBoxb.setValue(1)
                self.c.numberSpinBoxa.setValue(1)
        except Exception as e:
            print(e)

    def eventFilter(self, watched, event):
        """ Two Fuctions:
        1. MouseDblClick to choose an Atom.
        2. KeyPress to remove the superplot top layer.('a', 's', 'd', 'w', 'h', 'l')
        :param watched: the widget where the event happens.
        :param event: ...
        """
        if watched == self.view_widget:
            if event.type() == QEvent.MouseButtonDblClick and self.key_double_click == True:
                try:
                    mousex = event.pos().x()
                    mousey = event.pos().y()
                    try:
                        length_min = 25
                        choose_position = ...
                        for positions in self.positions_lst:  # 获得位置在屏幕上的投影的屏幕坐标
                            screentuple = gluProject(positions[0] + self.translate_pos[0], positions[1] +
                                                     self.translate_pos[1], positions[2] + self.translate_pos[2])
                            screenlen = np.sqrt((screentuple[0] - mousex) ** 2 +
                                                ((self.view_widget.height() - screentuple[1]) - mousey) ** 2)
                            if screenlen < length_min:
                                length_min = screenlen
                                choose_position = positions
                        if length_min != 25:
                            choose_position_index = self.positions_lst.index(choose_position)
                            self.tab.setCurrentWidget(self.tab2)
                            self.showmessage(choose_position_index)
                        else:
                            global x
                            try:
                                self.view_widget.removeItem(x)
                                self.tab.setCurrentWidget(self.tab1)
                                self.text_widgt2.setText("")
                                self.objectTextwidget.setText("")
                            except Exception as e:
                                print(e)

                            try:
                                global obj_lst_table
                                for obj in obj_lst_table:
                                    self.view_widget.removeItem(obj)
                                    self.tab.setCurrentWidget(self.tab1)
                                    self.text_widgt2.setText("")
                                    self.objectTextwidget.setText("")
                            except Exception as e:
                                print(e)

                            try:
                                global x_rectangle
                                for obj in x_rectangle.values():
                                    self.view_widget.removeItem(obj)
                                    self.tab.setCurrentWidget(self.tab1)
                                    self.text_widgt2.setText("")
                                    self.objectTextwidget.setText("")
                                x_rectangle = {}
                            except Exception as e:
                                print(e)
                    except Exception as e:
                        print(e)
                except Exception as e:
                    print(e)
            # 按键操作
            elif event.type() == QEvent.KeyPress and self.keylock == False:
                if event.key() == 87:  # 按下w                   # 向y轴正方向移动
                    try:
                        if self.superpositions_lst[self.index0][1] < 0:
                            for i in range(len(self.superobj_lst)):
                                self.superobj_lst[i].translate(0, 0.1, 0)
                                self.superpositions_lst[i][1] += 0.1

                    except Exception as e:
                        print(e)
                elif event.key() == 83:  # 按下s,向y轴负方向移动
                    try:
                        if self.superpositions_lst[self.index0][1] > -self.supercell1_cell[1] * abs(
                                np.sin(self.supercell1_cell[5] / 180 * np.pi)):

                            for i in range(len(self.superobj_lst)):
                                self.superobj_lst[i].translate(0, -0.1, 0)  # 向y轴负方向移动
                                self.superpositions_lst[i][1] += -0.1

                    except Exception as e:
                        print(e)
                elif event.key() == 65:  # 按下a      向x轴负方向移动
                    try:
                        if self.superpositions_lst[self.index0][0] > -self.supercell1_cell[0]:
                            for i in range(len(self.superobj_lst)):
                                self.superobj_lst[i].translate(-0.1, 0, 0)  # 向x轴负方向移动
                                self.superpositions_lst[i][0] += -0.1

                    except Exception as e:
                        print(e)

                elif event.key() == 68:  # 按下d
                    try:
                        if self.superpositions_lst[self.index0][0] < 0:
                            for i in range(len(self.superobj_lst)):
                                self.superobj_lst[i].translate(0.1, 0, 0)  # 向x轴正方向移动
                                self.superpositions_lst[i][0] += 0.1

                    except Exception as e:
                        print(e)
                elif event.key() == 72 and self.up_down == True:    # 按下h（high）
                    try:
                        for i in range(len(self.superobj_lst)):
                            self.superobj_lst[i].translate(0, 0, 0.1)  # 向x轴正方向移动
                            self.superpositions_lst[i][2] += 0.1
                    except Exception as e:
                        print(e)
                elif event.key() == 76 and self.up_down == True:    # 按下l（low）
                    try:
                        for i in range(len(self.superobj_lst)):
                            self.superobj_lst[i].translate(0, 0, -0.1)  # 向x轴正方向移动
                            self.superpositions_lst[i][2] -= 0.1
                    except Exception as e:
                        print(e)
                elif event.key() == 16777220:  # 按下回车
                    try:
                        self.afterenter()
                        self.keylock = True  # 锁住键盘
                        self.up_down = False # 锁住上下键
                        self.text_widgt3.setText("")
                    except:
                        pass
        return gl.GLViewWidget.eventFilter(self, watched, event)

    # 双击原子后显示的原子信息
    def showmessage(self, index, clear=True):
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
            global x
            try:
                self.view_widget.removeItem(x)
            except Exception as e:
                print(e)
        try:
            md = gl.MeshData.sphere(rows=10, cols=20)
            x = gl.GLMeshItem(meshdata=md, smooth=True, color=self.Atoms_color[self.Atoms_num_lst[index]],
                              shader='shaded', drawEdges=True)
            x.translate(self.positions_lst[index][0], self.positions_lst[index][1], self.positions_lst[index][2])
            x.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])
            if self.plot_num == 0:
                x.scale(self.Atom_radii[Atomic_number] + 0.01,
                        self.Atom_radii[Atomic_number] + 0.01,
                        self.Atom_radii[Atomic_number] + 0.01)
            elif self.plot_num == 1:
                x.scale(self.Atom_radii[Atomic_number]/2 + 0.005,
                        self.Atom_radii[Atomic_number]/2 + 0.005,
                        self.Atom_radii[Atomic_number]/2 + 0.005)
            self.view_widget.addItem(x)
            zhijiao_system = self.Atomsobject.get_cell()
            A = zhijiao_system.T
            b = np.array(self.positions_lst[index]).T
            r = np.linalg.solve(A, b)  # 求解线性方程组，直角坐标系下----用晶胞坐标系表示
            text = sf.crys_data.atom_formula[Atomic_number] + str(index + 1) + "\n" + "Atomic number:" + str(
                Atomic_number) + "\n" + "position:" + str(self.positions_lst[index]) + \
                   "\n" + "[uvw]:" + str(r)
            self.text_widgt2.setText(text)
            self.remove_atom_index_num = index
        except Exception as e:
            print(e)

    def openfile(self, new=False):
        """
        To open a file, draw it and put it in the self.tree widget
        :param new:create two forms, open file and new file
        """
        if new == False:
            self.fileName_choose, filetype = QFileDialog.getOpenFileName(self, '选择文件', '', 'files(*.cif , *.vasp)')
        if self.fileName_choose == "":
            print("\n取消选择")
        else:
            try:
                self.Atomsobject = copy.deepcopy(read(self.fileName_choose, index=None, format=None, parallel=True))
                # Text3column = self.judgeconductivity(self.Atomsobject)
                self.root = QtWidgets.QTreeWidgetItem(self.tree)
                self.plot(self.Atomsobject, dictionary=True)
                self.root.setText(0, self.dirkey)
                child = QtWidgets.QTreeWidgetItem(self.root)
                child.setText(0, "bulk")
                child.setText(1, self.dirkey)
                self.tree.expandAll()
            except Exception as e:
                print(e)


    def superplot(self, Atomsobject, clear=True, layer='top'):  # 对超胞plot 不会创造对象self.Atomsobject,
        """
        To plot a supercell.
        :param Atomsobject: the Atoms object which is going to be painted
        :param clear: clear = True to clear the painter
        :param layer: to judge whether the Atomsobject is in the top layer or not
        :return:
        """
        # clear==True时清空画布
        try:
            if clear == True:
                self.clear_painter()
            if layer == 'top':
                self.g = gl.GLGridItem()  # 画网格
                self.g.scale(1, 1, 1)
                self.view_widget.addItem(self.g)
            Atoms_num_lst = Atomsobject.get_atomic_numbers()
            self.superpositions_lst = Atomsobject.get_positions().tolist()
            # 画原子和位置
            md = gl.MeshData.sphere(rows=10, cols=20)
            if layer == 'top':  # 上面的超胞
                self.superobj_lst = []  # 存原子对象的列表
                for i in range(len(self.superpositions_lst)):
                    self.x = gl.GLMeshItem(meshdata=md, smooth=True, color=self.Atoms_color[Atoms_num_lst[i]],
                                           shader='shaded', drawFaces=True)
                    self.x.translate(self.superpositions_lst[i][0], self.superpositions_lst[i][1], self.superpositions_lst[i][2])
                    self.x.scale(self.Atom_radii[Atoms_num_lst[i]],
                                 self.Atom_radii[Atoms_num_lst[i]],
                                 self.Atom_radii[Atoms_num_lst[i]])
                    self.view_widget.addItem(self.x)
                    self.superobj_lst.append(self.x)
            elif layer == 'bottom':  # 底下的超胞
                self.superobj_lst1 = []  # 存原子对象的列表
                for i in range(len(self.superpositions_lst)):
                    self.x = gl.GLMeshItem(meshdata=md, smooth=True, color=self.Atoms_color[Atoms_num_lst[i]],
                                           shader='shaded',
                                           drawFaces=True)
                    self.x.translate(self.superpositions_lst[i][0], self.superpositions_lst[i][1], self.superpositions_lst[i][2])
                    self.x.scale(self.Atom_radii[Atoms_num_lst[i]],
                                 self.Atom_radii[Atoms_num_lst[i]],
                                 self.Atom_radii[Atoms_num_lst[i]])
                    self.view_widget.addItem(self.x)
                    self.superobj_lst1.append(self.x)
        except Exception as e:
            print(e)

    def clear_text(self):
        """ To clear the textwidget."""
        try:
            self.objectTextwidget.setText("")
        except Exception as e:
            print(e)
        try:
            self.text_widgt.setText("")
        except Exception as e:
            print(e)
        try:
            self.text_widgt2.setText("")
        except Exception as e:
            print(e)
        try:
            self.text_widgt3.setText("")
        except Exception as e:
            print(e)

    def clear_painter(self):
        """ To clear the painter."""
        try:
            global x
            self.view_widget.removeItem(x)
        except:
            pass
        try:
            self.view_widget.removeItem(self.g)  # 清除网格线
        except:
            pass
        try:
            self.tablewidget.setRowCount(0)  # 添加信息
        except:
            pass
        try:
            for item in self.cell_line_obj:  # 清空画布
                self.view_widget.removeItem(item)
        except:
            pass
        try:
            for item in self.axis_obj:
                self.view_widget.removeItem(item)
        except:
            pass

        try:
            for item in self.stick_lst:
                self.view_widget.removeItem(item)
        except:
            pass
        try:
            for item in self.obj_lst:
                self.view_widget.removeItem(item)
        except:
            pass
        try:
            for item in self.superobj_lst:  # 上面的超胞
                self.view_widget.removeItem(item)
        except:
            pass
        try:
            for item in self.superobj_lst1:  # 底下的超胞
                self.view_widget.removeItem(item)
        except:
            pass
        try:
            global obj_lst_table
            for item in obj_lst_table:
                self.view_widget.removeItem(item)
        except:
            pass
        try:
            global x_rectangle
            for obj in x_rectangle.values():
                self.view_widget.removeItem(obj)
            x_rectangle = {}
        except:
            pass

    def plot(self, Atomsobject, plot=True, clear=True, dictionary=True, globalAtomsobject=True,
             object=True, tab3=True):  # plot 会创造对象self.Atomsobject,目前有bug
        """
        Normally, to plot an ase.Atoms object.
        :param Atomsobject:the ase.Atoms object to be plotted
        :param clear: to clear the painter if clear == True
        :param dictionary: add the Atomsobject to self.dic_Formula_Atoms(dictionary)
        :param globalAtomsobject: whether create self.Atomsobject or not.
        :param object: whether add atoms to self.tablewidget or not
        """

        # clear==True时清空画布
        Number = len(Atomsobject.get_positions())
        if clear == True:
            self.clear_painter()
            self.clear_text()
        if globalAtomsobject == True:  # 若Atomsobject =True创建全局变量
            self.Atomsobject = Atomsobject
        if dictionary == True:  # 在字典中加入该对象
            self.dirkey = Atomsobject.get_chemical_formula(mode='hill')
            orijin_dirkey = self.dirkey
            num = 1
            while True:
                if self.dirkey in self.dic_Formula_Atoms.keys():
                    self.dirkey = orijin_dirkey + '({})'.format(num)
                    num += 1
                else:
                    break
            self.dic_Formula_Atoms[self.dirkey] = Atomsobject  # 在self.dic_formula_Atoms里添加Atomsobject
        if object == True:  # 在objectTab中录入晶体信息
            if Number < 500:
                Atomsobject_atomic_numbers_lst = Atomsobject.get_atomic_numbers()
                f = len(Atomsobject_atomic_numbers_lst)
                self.tablewidget.setRowCount(f)  # 添加信息
                for i in range(f):
                    newitem = QtWidgets.QTableWidgetItem(sf.crys_data.atom_formula[Atomsobject_atomic_numbers_lst[i]])
                    newitem.setTextAlignment(5 | 5)
                    self.tablewidget.setItem(i, 0, newitem)
                for i in range(f):
                    newitem = QtWidgets.QTableWidgetItem("")
                    color = self.Atoms_color[Atomsobject_atomic_numbers_lst[i]]
                    colortrans = (color[0] * 255, color[1] * 255, color[2] * 255)
                    self.tablewidget.setItem(i, 1, newitem)
                    self.tablewidget.item(i, 1).setBackground(QtGui.QBrush(
                        QtGui.QColor(colortrans[0], colortrans[1], colortrans[2])))

                for i in range(f):
                    newitem = QtWidgets.QTableWidgetItem(
                        str(self.Atom_radii[Atomsobject_atomic_numbers_lst[i]]))
                    newitem.setTextAlignment(5 | 5)
                    self.tablewidget.setItem(i, 2, newitem)
                self.tablewidget.setShowGrid(True)
        if tab3 == True:
            try:
                formula = Atomsobject.get_chemical_formula(mode='hill')
                self.database_class.exhibit_all_info(fenkuai=False)
                data_info = self.database_class.info_data_base
                tab3_info = []
                for i in range(len(data_info)):
                    if data_info[i][2] == formula:
                        tab3_info.append(data_info[i])
                self.tab3_table_widget1.setRowCount(len(tab3_info))
                self.tab3_table_widget2.setRowCount(len(tab3_info))
                print("tab3_info = ", tab3_info)
                for i in range(len(tab3_info)):
                    for k in range(10, 15):
                        newitem = QtWidgets.QTableWidgetItem(str(list(tab3_info[i])[k]))
                        newitem.setTextAlignment(5 | 5)
                        self.tab3_table_widget1.setItem(i, k - 10, newitem)
                    for k in range(15, 20):
                        newitem = QtWidgets.QTableWidgetItem(str(list(tab3_info[i])[k]))
                        newitem.setTextAlignment(5 | 5)
                        self.tab3_table_widget2.setItem(i, k - 15, newitem)
                self.tab3_table_widget2.setShowGrid(True)
                self.tab3_table_widget1.setShowGrid(True)
            except Exception as e:
                print(e)

        if plot == True:
            cell = Atomsobject.get_cell()
            print("cell = ", cell)
            self.translate_pos = -(cell[0] + cell[1] + cell[2]) / 2
            positions = Atomsobject.get_positions()
            self.Atoms_num_lst = Atomsobject.get_atomic_numbers()
            self.positions_lst = positions.tolist()
            self.cell_array = Atomsobject.get_cell()
            if Number < 900:

                self.key_double_click = True
                # 画原子和位置
                if self.plot_num == 0:
                    md = gl.MeshData.sphere(rows=10, cols=20)
                    self.obj_lst = []
                    for i in range(len(positions)):
                        self.x = gl.GLMeshItem(meshdata=md, smooth=True, color=self.Atoms_color[self.Atoms_num_lst[i]],
                                               shader='shaded',
                                               drawFaces=True)
                        self.x.translate(positions[i][0], positions[i][1], positions[i][2])
                        self.x.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])


                        self.x.scale(self.Atom_radii[self.Atoms_num_lst[i]],
                                     self.Atom_radii[self.Atoms_num_lst[i]],
                                     self.Atom_radii[self.Atoms_num_lst[i]])
                        self.view_widget.addItem(self.x)
                        self.obj_lst.append(self.x)
                elif self.plot_num == 1:
                    md = gl.MeshData.sphere(rows=10, cols=20)
                    self.obj_lst = []
                    self.stick_lst = []
                    for i in range(len(positions)):
                        self.x = gl.GLMeshItem(meshdata=md, smooth=True, color=self.Atoms_color[self.Atoms_num_lst[i]],
                                               shader='shaded',
                                               drawFaces=True)
                        self.x.translate(positions[i][0], positions[i][1], positions[i][2])
                        self.x.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])

                        self.x.scale(self.Atom_radii[self.Atoms_num_lst[i]]/2,
                                     self.Atom_radii[self.Atoms_num_lst[i]]/2,
                                     self.Atom_radii[self.Atoms_num_lst[i]]/2)
                        self.view_widget.addItem(self.x)
                        self.obj_lst.append(self.x)
                    for i in range(len(positions)):
                        for j in range(len(positions)):
                            if i == j:
                                continue
                            else:
                                vander1 = sf.crys_data.vander_wals_radii[self.Atoms_num_lst[i]]
                                vander2 = sf.crys_data.vander_wals_radii[self.Atoms_num_lst[j]]
                                dis = np.linalg.norm(positions[i] - positions[j])
                                if vander1 + vander2 - 1.3 > dis:
                                    # 不画金属键
                                    if sf.crys_data.metal_or_not[self.Atoms_num_lst[i]] == 0 or sf.crys_data.metal_or_not[self.Atoms_num_lst[j]] ==0:
                                        cylinder = gl.MeshData.cylinder(rows=10, cols=20, radius=[0.1, 0.1], length=dis/2, offset=False)
                                        self.y = gl.GLMeshItem(meshdata=cylinder, smooth=True, color=self.Atoms_color[self.Atoms_num_lst[i]],
                                                       shader='shaded',
                                                       drawFaces=True)
                                        pos1 = np.array(positions[i])
                                        pos2 = np.array(positions[j])
                                        pos_vec = pos2 - pos1
                                        if pos_vec[0] != 0 or pos_vec[1] != 0:
                                            if pos_vec[1] != 0:
                                                # vec_vertical = np.array([1, -pos_vec[0] / pos_vec[1], 0])
                                                if pos_vec[2] != 0:
                                                    fa_vec = np.array([0, 0, pos_vec[2]])
                                                    vec_vertical = np.cross(fa_vec, pos_vec)
                                                    project_vec = np.array([pos_vec[0], pos_vec[1], 0])
                                                    theta = np.pi/2 - np.arccos(pos_vec @ project_vec /
                                                                                np.linalg.norm(pos_vec) / np.linalg.norm(project_vec))
                                                    if pos_vec[2] > 0:
                                                        self.y.rotate(theta/np.pi*180, vec_vertical[0], vec_vertical[1], vec_vertical[2])
                                                    else:
                                                        self.y.rotate(theta/np.pi*180 + 180, vec_vertical[0], vec_vertical[1], vec_vertical[2])
                                                else:
                                                    fa_vec = np.array([0, 0, 1])
                                                    vec_vertical = np.cross(fa_vec, pos_vec)
                                                    self.y.rotate(90, vec_vertical[0], vec_vertical[1], vec_vertical[2])

                                            elif pos_vec[0] != 0:
                                                if pos_vec[2] != 0:
                                                    fa_vec = np.array([0, 0, pos_vec[2]])
                                                    vec_vertical = np.cross(fa_vec, pos_vec)
                                                    project_vec = np.array([pos_vec[0], pos_vec[1], 0])
                                                    theta = np.pi/2 - np.arccos(pos_vec @ project_vec / np.linalg.norm(pos_vec) / np.linalg.norm(project_vec))
                                                    if pos_vec[2] > 0:
                                                        self.y.rotate(theta/np.pi*180, vec_vertical[0], vec_vertical[1], vec_vertical[2])
                                                    else:
                                                        self.y.rotate(theta/np.pi*180 + 180, vec_vertical[0], vec_vertical[1], vec_vertical[2])
                                                else:
                                                    fa_vec = np.array([0, 0, 1])
                                                    vec_vertical = np.cross(fa_vec, pos_vec)
                                                    self.y.rotate(90, vec_vertical[0], vec_vertical[1], vec_vertical[2])

                                                # self.y.rotate(theta/np.pi*180, vec_vertical[0], vec_vertical[1], vec_vertical[2])
                                        # self.y.rotate()
                                        self.y.translate(positions[i][0], positions[i][1], positions[i][2])
                                        self.y.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])
                                        self.view_widget.addItem(self.y)
                                        self.stick_lst.append(self.y)
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Too many Atoms to show!')
                self.key_double_click = False
            # 画坐标轴
            if self.coordinatesystem_view_num == 2:  # 直角坐标系
                self.axis_obj = []
                vector_trans = np.array([-5, -5, 0])
                vector1_x = np.linspace(0, 5, 100)
                vector1_z = vector1_x * 0
                vector1_y = vector1_x * 0
                pts = np.vstack([vector1_x, vector1_y, vector1_z]).transpose()
                pts += vector_trans
                plt = gl.GLLinePlotItem(pos=pts, color=(10, 0, 0, 1), width=8 / 10., antialias=True)
                plt.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])

                self.view_widget.addItem(plt)
                self.axis_obj.append(plt)
                vector1_y = np.linspace(0, 5, 100)
                vector1_z = vector1_x * 0
                vector1_x = vector1_x * 0
                pts = np.vstack([vector1_x, vector1_y, vector1_z]).transpose()
                pts += vector_trans
                plt = gl.GLLinePlotItem(pos=pts, color=(0, 10, 0, 1), width=8 / 10., antialias=True)
                plt.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])
                self.view_widget.addItem(plt)
                self.axis_obj.append(plt)
                vector1_z = np.linspace(0, 5, 100)
                vector1_x = vector1_x * 0
                vector1_y = vector1_x * 0
                pts = np.vstack([vector1_x, vector1_y, vector1_z]).transpose()
                pts += vector_trans
                plt = gl.GLLinePlotItem(pos=pts, color=(0, 0, 10, 1), width=8 / 10., antialias=True)
                plt.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])

                self.view_widget.addItem(plt)
                self.axis_obj.append(plt)
            elif self.coordinatesystem_view_num == 1:  # 1的时候画晶体坐标系
                self.axis_obj = []
                vector_trans = np.array([-5, -5, 0])
                for i in range(3):
                    if self.cell_array[i][0] != 0:
                        vector1_x = np.linspace(0, self.cell_array[i][0], 100)
                        vector1_y = vector1_x / self.cell_array[i][0] * self.cell_array[i][1]
                        vector1_z = vector1_x / self.cell_array[i][0] * self.cell_array[i][2]
                        pts = np.vstack([vector1_x, vector1_y, vector1_z]).transpose()
                        pts += vector_trans
                        if i == 0:
                            plt = gl.GLLinePlotItem(pos=pts, color=(1, 0, 0, 1), width=5 / 10., antialias=True)
                        elif i == 1:
                            plt = gl.GLLinePlotItem(pos=pts, color=(0, 1, 0, 1), width=5 / 10., antialias=True)
                        elif i == 2:
                            plt = gl.GLLinePlotItem(pos=pts, color=(0, 0, 1, 1), width=5 / 10., antialias=True)

                    elif self.cell_array[i][2] != 0:
                        vector1_z = np.linspace(0, self.cell_array[i][2], 100)
                        vector1_y = vector1_z / self.cell_array[i][2] * self.cell_array[i][1]
                        vector1_x = vector1_z / self.cell_array[i][2] * self.cell_array[i][0]
                        pts = np.vstack([vector1_x, vector1_y, vector1_z]).transpose()
                        pts += vector_trans
                        if i == 0:
                            plt = gl.GLLinePlotItem(pos=pts, color=(1, 0, 0, 1), width=5 / 10., antialias=True)
                        elif i == 1:
                            plt = gl.GLLinePlotItem(pos=pts, color=(0, 1, 0, 1), width=5 / 10., antialias=True)
                        elif i == 2:
                            plt = gl.GLLinePlotItem(pos=pts, color=(0, 0, 1, 1), width=5 / 10., antialias=True)
                    elif self.cell_array[i][1] != 0:
                        vector1_y = np.linspace(0, self.cell_array[i][1], 100)
                        vector1_z = vector1_y / self.cell_array[i][1] * self.cell_array[i][2]
                        vector1_x = vector1_y / self.cell_array[i][1] * self.cell_array[i][0]
                        pts = np.vstack([vector1_x, vector1_y, vector1_z]).transpose()
                        pts += vector_trans
                        if i == 0:
                            plt = gl.GLLinePlotItem(pos=pts, color=(1, 0, 0, 1), width=5 / 10., antialias=True)
                        elif i == 1:
                            plt = gl.GLLinePlotItem(pos=pts, color=(0, 1, 0, 1), width=5 / 10., antialias=True)
                        elif i == 2:
                            plt = gl.GLLinePlotItem(pos=pts, color=(0, 0, 1, 1), width=5 / 10., antialias=True)
                    plt.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])
                    self.view_widget.addItem(plt)
                    self.axis_obj.append(plt)
                    self.view_widget.addItem(plt)
                    self.axis_obj.append(plt)

            # 画晶胞,画线
            if self.cell_view_num % 2 == 0:
                cell_array_line = []
                self.cell_line_obj = []
                for i in range(3):
                    if self.cell_array[i][0] != 0:
                        vector1_x = np.linspace(0, self.cell_array[i][0], 100)
                        vector1_y = vector1_x / self.cell_array[i][0] * self.cell_array[i][1]
                        vector1_z = vector1_x / self.cell_array[i][0] * self.cell_array[i][2]
                        pts = np.vstack([vector1_x, vector1_y, vector1_z]).transpose()
                        cell_array_line.append(pts)
                        plt = gl.GLLinePlotItem(pos=pts, color=pg.glColor((1, 2)), width=5 / 10., antialias=True)
                        plt.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])

                        self.view_widget.addItem(plt)
                        self.cell_line_obj.append(plt)
                    elif self.cell_array[i][2] != 0:
                        vector1_z = np.linspace(0, self.cell_array[i][2], 100)
                        vector1_y = vector1_z / self.cell_array[i][2] * self.cell_array[i][1]
                        vector1_x = vector1_z / self.cell_array[i][2] * self.cell_array[i][0]
                        pts = np.vstack([vector1_x, vector1_y, vector1_z]).transpose()
                        plt = gl.GLLinePlotItem(pos=pts, color=pg.glColor((1, 2)), width=5 / 10., antialias=True)
                        plt.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])

                        self.view_widget.addItem(plt)
                        cell_array_line.append(pts)
                        self.cell_line_obj.append(plt)
                    elif self.cell_array[i][1] != 0:
                        vector1_y = np.linspace(0, self.cell_array[i][1], 100)
                        vector1_z = vector1_y / self.cell_array[i][1] * self.cell_array[i][2]
                        vector1_x = vector1_y / self.cell_array[i][1] * self.cell_array[i][0]
                        pts = np.vstack([vector1_x, vector1_y, vector1_z]).transpose()
                        plt = gl.GLLinePlotItem(pos=pts, color=pg.glColor((1, 2)), width=5 / 10., antialias=True)
                        plt.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])

                        self.view_widget.addItem(plt)
                        cell_array_line.append(pts)
                        self.cell_line_obj.append(plt)
                # 画线

                for i in range(3):
                    lst = [0, 1, 2]
                    lst.remove(i)
                    x = cell_array_line[i] + self.cell_array[lst[0]]
                    plt = gl.GLLinePlotItem(pos=x, color=pg.glColor((1, 2)), width=5 / 10., antialias=True)
                    plt.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])

                    self.cell_line_obj.append(plt)
                    self.view_widget.addItem(plt)
                    y = cell_array_line[i] + self.cell_array[lst[1]]
                    plt = gl.GLLinePlotItem(pos=y, color=pg.glColor((1, 2)), width=5 / 10., antialias=True)
                    plt.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])

                    self.cell_line_obj.append(plt)
                    self.view_widget.addItem(plt)
                    z = cell_array_line[i] + self.cell_array[lst[0]] + self.cell_array[lst[1]]
                    plt = gl.GLLinePlotItem(pos=z, color=pg.glColor((1, 2)), width=5 / 10., antialias=True)
                    plt.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])

                    self.cell_line_obj.append(plt)
                    self.view_widget.addItem(plt)
            else:
                pass

            if self.grid_itemnum % 2 == 0:
                self.g = gl.GLGridItem()  # 画网格
                self.g.scale(1, 1, 1)
                self.g.translate(self.translate_pos[0], self.translate_pos[1], self.translate_pos[2])
                self.view_widget.addItem(self.g)
            else:
                ...
            text = self.all_info(Atomsobject)
            self.tab.setCurrentWidget(self.tab1)
            self.text_widgt.setText(text)

    def all_info(self, Atomsobject):
        """
        To get the information of ase.Atoms object.
        :param Atomsobject: ase.Atoms object
        :return: crysinfo
        """
        crys = Atomsobject
        f = str(AseAtomsAdaptor.get_structure(crys))
        num = re.findall('Sites (.*?)\n', f)[0]
        lst = f.split('Sites {}'.format(num))
        f1 = lst[0]
        f2 = lst[1]
        spacegroup = str(get_spacegroup(crys))
        space_group = 'Space group: ' + re.findall('\s\s(.*)\n', spacegroup)[0]
        crys_info = f1 + '\n' + space_group
        volume = str(crys.get_volume())
        unit_volume = 'Unit cell volume: ' + str(volume)
        crys_info = crys_info + '\n' + unit_volume + ' Å^3'
        f2 = 'Cell params:' + f2
        crys_info = crys_info + '\n\n' + f2
        return crys_info

    def judgetool(self):
        """ To cut self.Atomsobject.(defalut form)"""
        try:
            if self.Atomsobject == None:
                print("Please choose a cif file")
            else:
                current = self.tree.currentItem()
                if current is not None and current.parent() is not None:
                    if len(self.tree.selectedItems()) == 1:  # 用户只选择一个的情况
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
            current = self.tree.currentItem()
            if current is not None and current.parent() is not None:
                if len(self.tree.selectedItems()) == 1:  # 用户只选择一个的情况
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
            customize_cut_object = copy.deepcopy(self.Atomsobject)
            if number == 1:         # "xyz"
                A = customize_cut_object.get_cell().T
                cell_coordinat_system_par = []
                for i in cell_par:
                    b = np.array(i).T
                    r = np.linalg.solve(A, b)  # 求解线性方程组，直角坐标系下----用晶胞坐标系表示
                    cell_coordinat_system_par.append(r.T)
                layer = cut(customize_cut_object, a=cell_coordinat_system_par[0],
                                b=cell_coordinat_system_par[1],
                                c=cell_coordinat_system_par[2], origo=(0, 0, 0))
            else:         # "uvw"
                layer = cut(customize_cut_object, a=cell_par[0],
                                            b=cell_par[1],
                                            c=cell_par[2], origo=(0, 0, 0))
            Atoms = self.set_rotate(layer)
            self.plot(Atoms)
            childx = QtWidgets.QTreeWidgetItem(current.parent())
            childx.setText(1, self.dirkey)
            childx.setText(0, 'layer')
        except Exception as e:
            print(e)

    def cut_exf_layer(self, parent):
        """ Default cut, remove the exfoliable layer distance. (vanderwals1 + vanderwals2 -1.3 < dis())"""
        crystype = sf.CrystalOperationsNum1(self.Atomsobject)
        self.Atomsobject = copy.deepcopy(crystype.gamaToAcuteangle())   # 在切的时候统一转成锐角
        after_add_cell_self_ATO = copy.deepcopy(self.Atomsobject.repeat((1, 1, 2)))
        orijin = copy.deepcopy(after_add_cell_self_ATO)
        atomic_numbers_array = after_add_cell_self_ATO.get_atomic_numbers()
        f = after_add_cell_self_ATO.get_positions().tolist()
        gouzaolist = [f[i] for i in range(len(f))]
        for k in range(len(gouzaolist)):
            gouzaolist[k].append(k)
        gouzaolist.sort(key=operator.itemgetter(2))         # 根据z轴由小到大排序

        # 起点序号,z最小
        origo = gouzaolist[0][3]
        layer_zong_lst = []
        layer_lst = []
        for i in range(len(gouzaolist) - 1):
            if abs(gouzaolist[i+1][2] - gouzaolist[i][2]) < 0.1:       # z方向小于0.1埃就认为是一层的
                if i+1 == len(gouzaolist)-1:
                    layer_lst.append(Atom(atomic_numbers_array[gouzaolist[i][3]], (gouzaolist[i][0], gouzaolist[i][1],
                                                                                   gouzaolist[i][2])))
                    layer_lst.append(Atom(atomic_numbers_array[gouzaolist[i+1][3]], (gouzaolist[i+1][0],
                                                                                     gouzaolist[i+1][1], gouzaolist[i+1][2])))
                    layer_zong_lst.append(layer_lst)
                    break
                else:
                    layer_lst.append(Atom(atomic_numbers_array[gouzaolist[i][3]], (gouzaolist[i][0], gouzaolist[i][1],
                                                                                   gouzaolist[i][2])))
            else:
                if i+1 == len(gouzaolist)-1:
                    layer_lst.append(Atom(atomic_numbers_array[gouzaolist[i][3]], (gouzaolist[i][0], gouzaolist[i][1],
                                                                                   gouzaolist[i][2])))
                    layer_zong_lst.append(layer_lst)
                    layer_lst = []
                    layer_lst.append(Atom(atomic_numbers_array[gouzaolist[i+1][3]], (gouzaolist[i+1][0],
                                                                                     gouzaolist[i+1][1], gouzaolist[i+1][2])))
                    layer_zong_lst.append(layer_lst)
                    break
                else:
                    layer_lst.append(Atom(atomic_numbers_array[gouzaolist[i][3]], (gouzaolist[i][0],
                                                                                   gouzaolist[i][1], gouzaolist[i][2])))
                    layer_zong_lst.append(layer_lst)
                    layer_lst = []
        # layer_zong_lst: 每一层一个元素
        layer_layer_zong_lst = []
        for i in range(len(layer_zong_lst)-1):
            layer_Atom_lst1 = layer_zong_lst[i]
            layer_Atom_lst2 = layer_zong_lst[i+1]
            dismin = 1000000
            for layer_Atom1 in layer_Atom_lst1:
                for layer_Atom2 in layer_Atom_lst2:
                    pos_vector = layer_Atom2.position - layer_Atom1.position
                    dis = np.linalg.norm(pos_vector)
                    if dis < dismin:      # 找每两层之间的最小距离
                        dismin = dis
                        layer_layer = [layer_Atom1.number, layer_Atom2.number, dis,pos_vector]
            layer_layer_zong_lst.append(layer_layer)
        zhouqidis = 0
        exfo_dis = 0
        for i in range(len(layer_layer_zong_lst)):
            dis = layer_layer_zong_lst[i][2]
            vector = layer_layer_zong_lst[i][3]
            vander_wals1 = sf.crys_data.vander_wals_radii[layer_layer_zong_lst[i][0]]
            vander_wals2 = sf.crys_data.vander_wals_radii[layer_layer_zong_lst[i][1]]
            if dis > vander_wals1 + vander_wals2 - 1.3 and vector[2]/dis > 0.5:     # z轴占比在0.8
                exfo_dis += layer_layer_zong_lst[i][3] @ np.array([0, 0, 1])
                break
            else:
                zhouqidis += layer_layer_zong_lst[i][3] @ np.array([0, 0, 1])
        if exfo_dis == 0:
            QtWidgets.QMessageBox.warning(self, 'error', "Can't exfoliate!")
        else:
            if zhouqidis == 0:
                layer = cut(orijin, [1, 0, 0], [0, 1, 0], origo=origo, clength=0.01)
                self.layer = layer
            else:
                self.layer = cut(orijin, [1, 0, 0], [0, 1, 0], origo=origo, clength=zhouqidis + 0.01)
            self.plot(self.layer, dictionary=True, clear=True, globalAtomsobject=False)
            Text3column = self.judgeconductivity(self.layer)
            childx = QtWidgets.QTreeWidgetItem(parent)
            childx.setText(1, self.dirkey)
            childx.setText(0, 'layer')
            childx.setText(3, Text3column)

    def Stack_main_two(self):
        """ To stack two layer according to the client."""
        try:
            self.point1_index = None
            self.point2_index = None
            formula_lst = []
            itemlst = []
            for item in self.tree.selectedItems():
                formula_lst.append(item.text(1))
                itemlst.append(item)
            if len(formula_lst) != 2:
                QtWidgets.QMessageBox.warning(self, 'error',
                                              'Please choose 2 crystal in Projectbox' + '\n' + '(Press Ctrl or Shift).')
            elif itemlst[0].text(0) == 'bulk' or itemlst[1].text(0) == 'bulk':
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose stack or layer(doubleclick to change)')
            else:
                self.point2_index, self.point1_index = None, None
                self.zhexianplot = zhexian_plot()
                self.zhexianplot.signal_determine.connect(self.layer_transform)  # 按确定以后，执行layer_transform
                self.zhexianplot.signal_point1_index.connect(self.findnum1)
                self.zhexianplot.signal_point2_index.connect(self.findnum2)
                self.zhexianplot.signal_layer_length.connect(self.get_layerlength)
                self.obj1 = self.dic_Formula_Atoms[formula_lst[0]]  # Atoms1，对象
                self.obj2 = self.dic_Formula_Atoms[formula_lst[1]]  # Atoms2，对象
                self.orijincrys1 = copy.deepcopy(self.obj1)
                self.orijincrys2 = copy.deepcopy(self.obj2)
                try:
                    self.compare_cell(self.obj1, self.obj2)
                    if self.crys1 == None or self.crys2 == None:
                        QtWidgets.QMessageBox.warning(self, 'error', "Can't stack, for the γ difference is too big. ")
                    else:
                        Atomslst = [self.crys1, self.crys2]
                        self.zhexianplot.set_text_init(Atomslst)
                        self.zhexianplot.plot(self.axis_lst)  # 画折线plot
                except Exception as e:
                    print(e)
        except Exception as e:
            print(e)

    def Stack_main_one(self):
        try:
            formula_lst = []
            itemlst = []
            for item in self.tree.selectedItems():
                formula_lst.append(item.text(1))
                itemlst.append(item)
            if len(formula_lst) != 2:
                QtWidgets.QMessageBox.warning(self, 'error',
                                              'Please choose 2 crystal in Projectbox' + '\n' + '(Press Ctrl or Shift).')
            elif itemlst[0].text(0) == 'bulk' or itemlst[1].text(0) == 'bulk':
                QtWidgets.QMessageBox.warning(self, 'error', 'Please choose stack or layer(doubleclick to change)')
            else:
                self.stack_one_window = Stack_main_one_window()
                self.stack_one_window.signal_emit.connect(self.stack_one)
                self.orijincrys1 = copy.deepcopy(self.dic_Formula_Atoms[formula_lst[0]])  # Atoms1，对象
                obj2 = self.dic_Formula_Atoms[formula_lst[1]]  # Atoms2，对象
                crystype1 = sf.CrystalOperationsNum1(self.orijincrys1)
                self.obj1 = copy.deepcopy(crystype1.gamaToAcuteangle())  # 在切的时候统一转成锐角
                crystype2 = sf.CrystalOperationsNum1(obj2)
                self.obj2 = copy.deepcopy(crystype2.gamaToAcuteangle())  # 在切的时候统一转成锐角
        except Exception as e:
            print(e)

    def stack_one(self, strain, layer_dis, vacuum_dis, traversal=False):
        try:
            cell1 = self.obj1.get_cell()
            cell2 = self.obj2.get_cell()
            cell1_par = self.obj1.get_cell_lengths_and_angles()
            cell2_par = self.obj2.get_cell_lengths_and_angles()
            self.cell1 = np.array([[cell1[0][0], cell1[0][1]],
                                   [cell1[1][0], cell1[1][1]]])
            self.cell2 = np.array([[cell2[0][0], cell2[0][1]],
                                   [cell2[1][0], cell2[1][1]]])
            a1 = cell1_par[0]
            a2 = cell2_par[0]
            a1_times, a2_times = self.chase(a1, a2, strain)
            b1 = cell1_par[1]
            b2 = cell2_par[1]
            b_vertical1 = b1 * np.sin(cell1_par[5] / 180 * np.pi)
            b_vertical2 = b2 * np.sin(cell2_par[5] / 180 * np.pi)
            b1_times, b2_times = self.chase(b_vertical1, b_vertical2, strain)
            # print("b1_times = ", b1_times)
            # print("b2_times = ", b2_times)
            # print("b1_vlength =", b1_times * b_vertical1, "b2_vlength = ", b2_times * b_vertical2)
            init_b1 = copy.deepcopy(self.cell1[1] * b1_times)
            init_b2 = copy.deepcopy(self.cell2[1] * b2_times)
            init_a1 = copy.deepcopy(self.cell1[0])
            init_a2 = copy.deepcopy(self.cell2[0])
            number1 = - int(b1 * np.cos(cell1_par[5] / 180 * np.pi) * b1_times / a1) - 1
            number2 = - int(b2 * np.cos(cell2_par[5] / 180 * np.pi) * b2_times / a2) - 1
            number1, number2 = self.find_a_var(init_b1,  init_b2,
                                                       init_a1, init_a2, strain, number1, number2)
            if (number1 is None or number2 is None) and traversal_stack == False:
                QtWidgets.QMessageBox.warning(self, 'error',
                                              "Don't find.")
            elif (number1 is None or number2 is None) and traversal_stack == True:
                return "fail"
            else:
                a1 = (a1_times, 0, 0)
                b1 = (number1, b1_times, 0)
                a2 = (a2_times, 0, 0)
                b2 = (number2, b2_times, 0)
                self.Atomsobject = cut(self.obj1, a=a1, b=b1, c=[0, 0, 1], origo=(0, 0, 0))
                crys1_cell = self.Atomsobject.get_cell_lengths_and_angles()
                crys2_new = cut(self.obj2, a=a2, b=b2, c=[0, 0, 1], origo=(0, 0, 0))
                crys2_cell = crys2_new.get_cell_lengths_and_angles()
                cell1_new = self.Atomsobject.get_cell()
                crys2_new.translate(np.array([0, 0, layer_dis + cell1_new[2][2]]))
                crys1_cell[2] = crys1_cell[2] + crys2_cell[2] + layer_dis + 0.01
                self.Atomsobject.extend(crys2_new)
                self.Atomsobject.set_cell(crys1_cell)
                Atomsobject = self.set_vacuum_layer(vacuum_dis, addproject=False)
                self.plot(Atomsobject, plot=False, object=False, clear=False, globalAtomsobject=False, dictionary=True,
                          tab3=False)
                key = list(self.dic_Formula_Atoms.keys())[
                    list(self.dic_Formula_Atoms.values()).index(self.orijincrys1)]
                it = QtWidgets.QTreeWidgetItemIterator(self.tree)

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
                    parent = self.tree.selectedItems()[0].parent()
                    stackchild = QtWidgets.QTreeWidgetItem(parent)
                    stackchild.setText(0, "stack")
                    stackchild.setText(1, self.dirkey)
                    if self.tree.selectedItems()[0].text(2) != "":
                        text1 = self.tree.selectedItems()[0].text(2)
                    else:
                        text1 = self.tree.selectedItems()[0].text(1)
                    if self.tree.selectedItems()[1].text(2) != "":
                        text2 = self.tree.selectedItems()[1].text(2)
                    else:
                        text2 = self.tree.selectedItems()[1].text(1)
                    stackchild.setText(2, text1 + '-' + text2)
                except Exception as e:
                    print(e)
                    stackchild.setText(2, self.tree.selectedItems()[0].text(2))
        except Exception as e:
            print(e)


    def find_a_var(self, b1, b2, a1, a2, strain, number1, number2):
        try:
            print("b1 = ", b1)
            print("b2 = ", b2)
            while number1 < 20 and number2 < 20:
                b1_new = copy.deepcopy(b1 + number1 * a1)
                b2_new = copy.deepcopy(b2 + number2 * a2)
                relative_error = np.linalg.norm(b1_new - b2_new) / np.linalg.norm(b1_new)
                if relative_error < strain:
                    return number1, number2
                elif b1_new[0] > b2_new[0]:
                    number2 += 1
                else:
                    number1 += 1
        except Exception as e:
            print(e)


    def get_layerlength(self, length):
        """ globalize layerlength"""
        self.layerlength = length

    def chase_lst(self, float1, float2):
        lst = [Fraction(1, 1)]
        orijin1 = float1
        orijin2 = float2
        strain = abs(float1 - float2) / float1
        num1 = 1
        num2 = 1
        num = 0
        zongnum = int((strain - 0.005) / 0.001)
        while num < zongnum and strain > 0.005 and num1 < 50 and num2 < 50:
            if float1 > float2:
                float2 += orijin2
                num2 += 1
                strain_new = abs(float1 - float2) / float1
                if strain - strain_new > 0.001:
                    lst.append(Fraction(num2, num1))
                    strain = strain_new
                    num += 1
                    print(Fraction(num2, num1))
            else:
                float1 += orijin1
                num1 += 1
                strain_new = abs(float1 - float2) / float1
                if strain - strain_new > 0.001:
                    lst.append(Fraction(num2, num1))
                    strain = strain_new
                    num += 1
                    print(Fraction(num2, num1))

        return lst

    def chase(self, float1, float2, strain_object):
        # lst = [Fraction(1, 1)]
        orijin1 = float1
        orijin2 = float2
        strain = abs(float1 - float2) / float1
        num1 = 1
        num2 = 1
        while strain > strain_object:
            if float1 > float2:
                float2 += orijin2
                num2 += 1
                strain = abs(float1 - float2) / float1
            else:
                float1 += orijin1
                num1 += 1
                strain = abs(float1 - float2) / float1
        return num1, num2


    def compare_cell(self, crys1, crys2, traversalstack=False):  # 绘图    # 计算两原胞a，b轴,返回应变-超胞a,b长的列表
        """
        To compare two crys.
        :param crys1: the bottom crys
        :param crys2: the top crys
        :param traversalstack: whether is traversalstack form or not
        :return: lista(b) ,which have strain1, strain2, strain_sum, lcm,
        int1/100(crys1'a-axis), int2/100(cry2'a-axis) information
        """
        self.crys1 = copy.deepcopy(crys1)         # globalize crys1,crys2
        self.crys2 = copy.deepcopy(crys2)
        crys1_cell_lst = crys1.get_cell_lengths_and_angles().tolist()
        crys2_cell_lst = crys2.get_cell_lengths_and_angles().tolist()
        shear_strain = abs(crys1_cell_lst[5] - crys2_cell_lst[5]) / 180 * np.pi
        if shear_strain < 0.001745:  # 若两晶体的gama角之差小于0.1°，则认为属于一种晶系
            crys1_a_axis_len,  crys2_a_axis_len, crys1_b_axis_len, crys2_b_axis_len = \
                crys1_cell_lst[0], crys2_cell_lst[0], crys1_cell_lst[1], crys2_cell_lst[1]
            # lista:strain1, strain2, strain_sum, lcm, int1(晶体1的a-axis), int2(晶体2的a-axis)
            if traversalstack == False:
                a_ratio_of_elong = self.chase_lst(crys1_a_axis_len, crys2_a_axis_len)
                b_ratio_of_elong = self.chase_lst(crys1_b_axis_len, crys2_b_axis_len)
                # 用追逐法
                self.lista = self.make_listab(crys1_a_axis_len, crys2_a_axis_len, a_ratio_of_elong)
                self.listb = self.make_listab(crys1_b_axis_len, crys2_b_axis_len, b_ratio_of_elong)
                self.lista.sort(key=operator.itemgetter(2))  # lista以strainsum排序
                self.listb.sort(key=operator.itemgetter(2))
                functiona = [[self.lista[i][2], self.lista[i][3]] for i in range(len(self.lista))]
                functionb = [[self.listb[i][2], self.listb[i][3]] for i in range(len(self.listb))]
                # functiona = [strainsum, lcma]
                # functionb = [strainsum, lcmb]
                self.axis_lst = [functiona, functionb]
            else:     # traversal_stack
                self.a1_times, self.a2_times = self.chase(crys1_a_axis_len, crys2_a_axis_len, self.traversalstack_tolerance)
                self.b1_times, self.b2_times = self.chase(crys1_b_axis_len, crys2_b_axis_len, self.traversalstack_tolerance)

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

    def findnum1(self, number1):
        """ client chooses the strainsum - lcma dot.(GLOBAL THE NUMBER)"""
        self.zhexianplot1_info = ""
        self.point1_index = number1
        crys1_a = str(self.lista[self.point1_index][4])  # 晶体1的a-axis
        crys2_a = str(self.lista[self.point1_index][5])  # 晶体2的a-axis
        strain_sum = str(self.lista[self.point1_index][2])
        strain1 = str(self.lista[self.point1_index][0])
        strain2 = str(self.lista[self.point1_index][1])
        lcma = str(self.lista[self.point1_index][3])
        self.zhexianplot1_info = 'StrainSum-a = ' + strain_sum + '    lcm-a = ' + lcma + \
                                 '\n' + '-' * 32 + '\n' + 'crys1-a: ' + crys1_a + '    strain1-a: ' + \
                                 strain1 + '\n' + 'crys2-a: ' + crys2_a + '    strain2-a: ' + \
                                 strain2
        try:
            if self.point2_index is not None:
                self.zhexianplot.text_widgt.setText(self.zhexianplot1_info + '\n' + '\n' + self.zhexianplot2_info)
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
                at_num_text = "Atom's number = " + str(atom_num_sum)
                self.zhexianplot.text_widgt.setText(
                    at_num_text + '\n' + self.zhexianplot1_info + '\n' + '\n' + self.zhexianplot2_info)
        except Exception as e:
            self.zhexianplot.text_widgt.setText(self.zhexianplot1_info)
            print(e)

    def findnum2(self, number2):
        """ client chooses the strainsum - lcmb dot.(GLOBAL THE INDEX)"""
        self.zhexianplot2_info = ""
        self.point2_index = number2
        crys1_b = str(self.listb[self.point2_index][4])  # 晶体1的b-axis
        crys2_b = str(self.listb[self.point2_index][5])  # 晶体2的b-axis
        strain_sum = str(self.listb[self.point2_index][2])
        strain1 = str(self.listb[self.point2_index][0])
        strain2 = str(self.listb[self.point2_index][1])
        lcmb = str(self.listb[self.point2_index][3])
        self.zhexianplot2_info = 'StrainSum-b = ' + strain_sum + '    lcm-b = ' + lcmb + \
                                 '\n' + '-' * 32 + '\n' + 'crys1-b: ' + crys1_b + '    strain1-b: ' + \
                                 strain1 + '\n' + 'crys2-b: ' + crys2_b + '    strain2-b: ' + \
                                 strain2
        try:
            if self.point1_index is not None:
                self.zhexianplot.text_widgt.setText(
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
                at_num_text = "Atom's number = " + str(atom_num_sum)
                self.zhexianplot.text_widgt.setText(
                    at_num_text + '\n' + self.zhexianplot1_info + '\n' + '\n' + self.zhexianplot2_info)
        except Exception as e:
            self.zhexianplot.text_widgt.setText(self.zhexianplot2_info)
            print(e)

    def layer_transform(self, traversalstack=False):
        """ To change two crys' a,b axis."""
        try:
            if self.crys1 is not None and self.crys2 is not None:
                if traversalstack == False:
                    orijin1 = copy.deepcopy(self.crys1)  # 不改变crys1的晶胞状态
                    orijin2 = copy.deepcopy(self.crys2)
                    crys1_cell_lst = orijin1.get_cell_lengths_and_angles().tolist()  # 返回存有[应变， 最小公倍数]的列表
                    crys2_cell_lst = orijin2.get_cell_lengths_and_angles().tolist()
                    if self.point1_index is None or self.point2_index is None:
                        print("请选择", "Num1 = ", self.point1_index, "Num2=", self.point2_index)
                    else:
                        obj_strain_lcma = self.axis_lst[0][self.point1_index]  # client选择的strain-lcma列表
                        obj_strain_lcmb = self.axis_lst[1][self.point2_index]  # client选择的strain-lcmb列表
                        crys1_cell_lst[0] = self.lista[self.point1_index][4]  # 晶体1的a-axis
                        crys1_cell_lst[1] = self.listb[self.point2_index][4]  # 晶体1的b-axis
                        crys2_cell_lst[0] = self.lista[self.point1_index][5]  # 晶体2的a-axis
                        crys2_cell_lst[1] = self.listb[self.point2_index][5]  # 晶体2的b-axis
                        orijin1.set_cell(crys1_cell_lst, scale_atoms=True)
                        orijin2.set_cell(crys2_cell_lst, scale_atoms=True)  # 按照client需求改变晶胞大小,里面的原子做仿射变换
                        multi_a_crys1 = round(obj_strain_lcma[1] / crys1_cell_lst[0])
                        multi_a_crys2 = round(obj_strain_lcma[1] / crys2_cell_lst[0])
                        multi_b_crys1 = round(obj_strain_lcmb[1] / crys1_cell_lst[1])
                        multi_b_crys2 = round(obj_strain_lcmb[1] / crys2_cell_lst[1])
                        self.crys1_supercell = orijin1.repeat((multi_a_crys1, multi_b_crys1, 1))
                        self.crys2_supercell = orijin2.repeat((multi_a_crys2, multi_b_crys2, 1))
                        self.zhexianplot2_info = ""  # 恢复初始情况
                        self.zhexianplot1_info = ""
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
            self.supercell1_cell = self.crys1_supercell.get_cell_lengths_and_angles()
            self.supercell2_cell = self.crys2_supercell.get_cell_lengths_and_angles()
            supercell1_cell_ndarray = self.crys1.get_cell()
            self.crys2_supercell.translate(np.array([0, 0, layerlength + supercell1_cell_ndarray[2][2]]))
            self.supercell1_cell[2] = self.supercell1_cell[2] + self.supercell2_cell[2] + layerlength + 0.01
            if traversalstack == 0:  # 非traversalstack
                self.superplot(self.crys1_supercell, clear=True, layer='bottom')
                self.superplot(self.crys2_supercell, clear=False, layer='top')
                self.keylock = False  # 释放键盘
                xmin = 100
                ymin = 100
                for i in range(len(self.superpositions_lst)):
                    if self.superpositions_lst[i][0] < xmin and self.superpositions_lst[i][1] < ymin:
                        xmin = self.superpositions_lst[i][0]  # superpositions_lst是堆垛在上面的超胞的位置列表, 取出[0,0,h]
                        ymin = self.superpositions_lst[i][1]  # superpositions_lst是堆垛在上面的超胞的位置列表, 取出[0,0,h]
                        self.index0 = i  # 取距离（0， 0）位置最近的点用于定位
            else:  # 批量处理，用户没法自定义
                self.crys1_supercell.extend(self.crys2_supercell)
                self.crys1_supercell.set_cell(self.supercell1_cell)
                self.Atomsobject = copy.deepcopy(self.crys1_supercell)
                self.crys1_supercell = self.set_vacuum_layer(self.vacuum_distance, addproject=False)
                self.plot(self.crys1_supercell, plot=False, object=False, clear=False, globalAtomsobject=False,
                          dictionary=True, tab3=False)
                key = list(self.dic_Formula_Atoms.keys())[list(self.dic_Formula_Atoms.values()).index(self.orijincrys1)]
                it = QtWidgets.QTreeWidgetItemIterator(self.tree)
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

    def afterenter(self):  # 按下回车键后触发
        """ After client pressing down 'enter', create a supercell."""
        try:  # 阔胞后切割
            self.crys2_supercell.set_positions(self.superpositions_lst)
            self.crys1_supercell.extend(self.crys2_supercell)  # 把crys2加在crys1上
            self.crys1_supercell.set_cell(self.supercell1_cell)
            self.crys1_supercell = copy.deepcopy(cut(self.crys1_supercell, a=(1, 0, 0), b=(0, 1, 0), c=(0, 0, 1),
                                                     origo=(0, 0, 0)))
        except Exception as e:
            print(e)
        # 画出切割后的超胞
        try:
            self.plot(self.crys1_supercell, clear=True, globalAtomsobject=True, dictionary=True)
            parent = self.tree.selectedItems()[0].parent()
            stackchild = QtWidgets.QTreeWidgetItem(parent)
            stackchild.setText(0, "stack")
            stackchild.setText(1, self.dirkey)
            if self.tree.selectedItems()[0].text(2) != "":
                text1 = self.tree.selectedItems()[0].text(2)
            else:
                text1 = self.tree.selectedItems()[0].text(1)
            if self.tree.selectedItems()[1].text(2) != "":
                text2 = self.tree.selectedItems()[1].text(2)
            else:
                text2 = self.tree.selectedItems()[1].text(1)
            stackchild.setText(2, text1 + '-' + text2)

        except Exception as e:
            print(e)
            stackchild.setText(2, self.tree.selectedItems()[0].text(2))

    def Export(self):
        """ To export cif files as a dictionary."""
        try:
            output_file = QFileDialog.getSaveFileName(self, "Save as cif",
                                                      "untitled.dic")
            if output_file[0] == "":
                print("请重新选择")
            else:
                folder = os.path.exists(output_file[0])
                if not folder:  # 判断是否存在文件夹如果不存在则创建为文件夹	,建文件时如果路径不存在会创建这个路径
                    os.makedirs(output_file[0])
                it = QtWidgets.QTreeWidgetItemIterator(self.tree)
                while it.value():
                    if it.value().text(1) == '':  # 根目录
                        try:
                            root = output_file[0] + '/{}'.format(it.value().text(0))  # 根目录地址
                            os.makedirs(root)
                        except Exception as e:
                            print(e)
                    else:
                        self.export_to_cif(self.dic_Formula_Atoms[it.value().text(1)], root, it.value().text(1))
                    it += 1
        except Exception as e:
            print(e)

    def export_to_cif(self, crys, dir, name):  # 输出crys.cif文件
        """ To export cif files as a dictionary."""
        f = AseAtomsAdaptor.get_structure(crys)
        f1 = str(f) + '\n'
        j_pv_lst = re.findall('abc\s\s\s:(.*?)\n', f1)[0]  # abc   :  19.257300  19.569178  21.133988
        j1_pv_lst = j_pv_lst.split(' ')  # abc   :   6.419100   6.523059   7.044663
        while '' in j1_pv_lst:
            j1_pv_lst.remove('')
        par_lst_matrix = crys.get_cell()
        # 物质种类（比如：Lu2 Al4）
        y1 = re.findall('Full\sFormula\s(.*?)\n', f1)[0]
        y1.lstrip('(').rstrip(')')
        material = re.findall('Full\sFormula\s(.*?)\n', f1)[0].lstrip('(').rstrip(')')
        mat_name = name
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
            for i in range(len(letter_lst)):
                if letter_lst[i] in szb_lst:
                    num_lst.append(letter_lst[i])
                if letter_lst[i] in zmb_lst:
                    symbol_lst.append(letter_lst[i])
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
        for i in range(len(element_lst)):
            num = int(number_lst[i])
            for j in range(num):
                par_lst_species.append(element_lst[i])
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
        slab.write_file(dir + "\{}.cif".format(mat_name))

    def Import(self):  # 两层或一层文件夹
        """ To import dictionary with cif files."""
        try:
            directory = QtGui.QFileDialog.getExistingDirectory(self, 'Select directory')
            names = []
            dir_lst = []
            for dirpath, dirs, files in os.walk(directory):  # 递归遍历当前目录和所有子目录的文件和目录
                for name in files:  # files保存的是所有的文件名
                    if os.path.splitext(name)[1] == '.cif' or os.path.splitext(name)[1] == '.vasp':
                        file_path = os.path.join(dirpath, name)  # 加上路径，dirpath是遍历时文件对应的路径
                        names.append(name)
                        dir_lst.append(file_path)
            for i in range(len(names)):
                try:
                    self.Atomsobject = copy.deepcopy(read(dir_lst[i]))
                    if i != len(names) - 1:
                        self.plot(self.Atomsobject, plot=False, object=False, clear=False, dictionary=True,
                                  globalAtomsobject=False, tab3=False)
                    else:
                        self.plot(self.Atomsobject, dictionary=True, globalAtomsobject=False)
                    self.root = QtWidgets.QTreeWidgetItem(self.tree)
                    self.root.setText(0, names[i])
                    # Textcolumn3 = self.judgeconductivity(self.Atomsobject)   # 判断导电性
                    child = QtWidgets.QTreeWidgetItem(self.root)
                    child.setText(0, "bulk")
                    child.setText(1, self.dirkey)
                except Exception as e:
                    print(e)
            self.tree.expandAll()
        except Exception as e:
            print(e)

    def judgeconductivity(self, Atomsobject):
        """ To judge the conductivity of an Atoms object."""
        Numberlist = Atomsobject.get_atomic_numbers()
        Numberset_lst = list(set(Numberlist))
        Elementlist = []
        for number in Numberset_lst:
            Elementlist.append(sf.crys_data.atom_formula[number])

        num = 0
        conductlist = [sf.crys_data.conductor_Elemental_composition, sf.crys_data.semiconductor_Element_composition,
                       sf.crys_data.Insulator_Element_composition]
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

    # def classification_main(self):
    #     """ To judge whether the crys is 2D materials or not."""
    #     try:
    #         self.classify_window = Classification_window()
    #         self.classify_window.signal_emit_dirs.connect(self.after_classify_window)
    #     except Exception as e:
    #         print(e)
    #
    # def after_classify_window(self, cifdir, outputdir):
    #     # dir_path = cifdir
    #     try:
    #         objs = sf.FileOperations(cifdir)  # 创建FileOperations类
    #         objs = objs.read_files()  # 打开dir_path， 返回创建的 atoms对象列表
    #         zongnum = len(objs)
    #         num = 0
    #         for cry_obj in objs:
    #             try:
    #                 smlkc.CrystalOperationsNum1(cry_obj, outputdir)  # 创一个对象
    #                 num += 1
    #                 self.classify_window.pbar.setValue(int(num / zongnum * 100))
    #             except Exception as e:
    #                 print(e)
    #         self.classify_window.close()
    #     except Exception as e:
    #         print(e)

    def transcrys(self, a_vector, b_vector, object_angle, episilon, num):
        """ To change crys' γ degree. To refind a,b axis-vector."""
        try:
            global orijina, orijinb
            if abs(object_angle - 90) <= 71:
                if num == 0:
                    orijina = copy.deepcopy(a_vector)
                    orijinb = copy.deepcopy(b_vector)
                gama = sf.MathOperation.get_theat(a_vector, b_vector)
                objectgama = object_angle / 180 * np.pi
                episilon = episilon / 180 * np.pi

                if abs(objectgama - gama) <= episilon:
                    print(a_vector, b_vector, "gamma = ", gama)
                    return a_vector, b_vector
                else:
                    if gama < objectgama:
                        if num % 2 == 0:            # b轴变化
                            while gama < objectgama:
                                orijin = copy.deepcopy(b_vector)
                                b_vector -= orijina
                                gama = sf.MathOperation.get_theat(a_vector, b_vector)
                                if abs(objectgama - gama) <= episilon:
                                    return a_vector, b_vector
                                    break
                            episilon = episilon * 180 / np.pi
                            gama_orijin = sf.MathOperation.get_theat(a_vector, orijin)
                            gama_b_vector = sf.MathOperation.get_theat(a_vector, b_vector)
                            if abs(gama_orijin - objectgama) < abs(gama_b_vector - objectgama):
                                return self.transcrys(a_vector, orijin + orijinb, object_angle, episilon, num + 1)
                            else:
                                return self.transcrys(a_vector, b_vector + orijinb, object_angle, episilon, num + 1)
                        else:            # a轴变化
                            while gama < objectgama:
                                orijin = copy.deepcopy(a_vector)
                                a_vector -= orijinb
                                gama = sf.MathOperation.get_theat(a_vector, b_vector)
                                if abs(objectgama - gama) <= episilon:
                                    return a_vector, b_vector
                                    break
                            episilon = episilon * 180 / np.pi
                            gama_orijin = sf.MathOperation.get_theat(orijin, b_vector)
                            gama_a_vector = sf.MathOperation.get_theat(a_vector, b_vector)

                            if abs(gama_orijin - objectgama) < abs(gama_a_vector - objectgama):
                                return self.transcrys(orijin + orijina, b_vector, object_angle, episilon, num + 1)
                            else:
                                return self.transcrys(a_vector + orijina, b_vector, object_angle, episilon, num + 1)
                    elif gama > objectgama:
                        if num % 2 == 0:      # b轴变化
                            while gama > objectgama:
                                orijin = copy.deepcopy(b_vector)
                                b_vector += orijina
                                print("orijin = ", orijin)
                                gama = sf.MathOperation.get_theat(a_vector, b_vector)
                                if abs(objectgama - gama) <= episilon:
                                    return a_vector, b_vector
                                    break
                            episilon = episilon * 180 / np.pi
                            gama_orijin = sf.MathOperation.get_theat(a_vector, orijin)
                            gama_b_vector = sf.MathOperation.get_theat(a_vector, b_vector)
                            if abs(gama_orijin - objectgama) < abs(gama_b_vector - objectgama):
                                return self.transcrys(a_vector, orijin + orijinb, object_angle, episilon, num + 1)
                            else:
                                return self.transcrys(a_vector, b_vector + orijinb, object_angle, episilon, num + 1)
                        else:    # a轴变化
                            while gama > objectgama:
                                orijin = copy.deepcopy(a_vector)
                                a_vector += orijinb
                                print("orijin = ", orijin)
                                gama = sf.MathOperation.get_theat(a_vector, b_vector)
                                if abs(objectgama - gama) <= episilon:
                                    return a_vector, b_vector
                                    break
                            episilon = episilon * 180 / np.pi
                            gama_orijin = sf.MathOperation.get_theat(orijin, b_vector)
                            gama_a_vector = sf.MathOperation.get_theat(a_vector, b_vector)
                            if abs(gama_orijin - objectgama) < abs(gama_a_vector - objectgama):
                                return self.transcrys(orijin + orijina, b_vector, object_angle, episilon, num + 1)
                            else:
                                return self.transcrys(a_vector + orijina, b_vector, object_angle, episilon, num + 1)
            else:
                try:
                    if num == 0:
                        # global orijina, orijinb
                        orijina = copy.deepcopy(a_vector)
                        orijinb = copy.deepcopy(b_vector)
                    gama = sf.MathOperation.get_theat(a_vector, b_vector)
                    objectgama = object_angle / 180 * np.pi
                    episilon = episilon / 180 * np.pi
                    if abs(objectgama - gama) <= episilon:
                        return a_vector, b_vector
                    else:
                        if gama < objectgama:
                            if num % 2 == 0:
                                while gama < objectgama:
                                    orijin = copy.deepcopy(b_vector)
                                    b_vector -= orijina
                                    gama = sf.MathOperation.get_theat(a_vector, b_vector)
                                episilon = episilon * 180 / np.pi
                                return self.transcrys(a_vector, b_vector + orijin, object_angle, episilon, num + 1)
                            else:
                                while gama < objectgama:
                                    orijin = copy.deepcopy(a_vector)
                                    a_vector -= orijinb
                                    gama = sf.MathOperation.get_theat(a_vector, b_vector)
                                episilon = episilon * 180 / np.pi
                                return self.transcrys(a_vector + orijin, b_vector, object_angle, episilon, num + 1)

                        elif gama > objectgama:
                            if num % 2 == 0:
                                while True:
                                    orijin = copy.deepcopy(b_vector)
                                    b_vector += orijina
                                    print("orijin = ", orijin)
                                    gama = sf.MathOperation.get_theat(a_vector, b_vector)
                                    if gama <= objectgama:
                                        break
                                episilon = episilon * 180 / np.pi
                                return self.transcrys(a_vector, b_vector + orijin, object_angle, episilon, num + 1)
                            else:
                                while True:
                                    orijin = copy.deepcopy(a_vector)
                                    a_vector += orijinb
                                    print("orijin = ", orijin)
                                    gama = sf.MathOperation.get_theat(a_vector, b_vector)
                                    if gama <= objectgama:
                                        break
                                episilon = episilon * 180 / np.pi
                                return self.transcrys(a_vector + orijin, b_vector, object_angle, episilon, num + 1)
                except Exception as e:
                    print(e)
        except Exception as e:
            print(e)

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
    myshow = mywindow()  # 主窗口实例化
    sys.exit(app.exec_())

