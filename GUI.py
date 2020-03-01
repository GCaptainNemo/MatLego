"""
The GUI of matlego 1.5.
"""


from PyQt5 import QtCore, QtGui, QtWidgets
import pyqtgraph.opengl as gl
from PyQt5.QtCore import pyqtSignal
import numpy as np
import pyqtgraph as pg
from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifWriter
import src.function as sf
import copy
from pyqtgraph.parametertree import Parameter, ParameterTree


class mytreewidget(QtWidgets.QTreeWidget):
    """ To rewrite the Qtwidgets.QTreeWidget. """
    def __init__(self):
        super(mytreewidget, self).__init__()
        self.setAnimated(True)
        self.setDragDropMode(QtGui.QAbstractItemView.InternalMove)
        self.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)


class Myviewwidget(gl.GLViewWidget):      # 重写glviewwidget
    """ To rewrite the gl.GLViewwidget to achieve the drawrectangle function of self.view_widget. """
    num1 = 0
    signal_release = pyqtSignal(tuple, tuple)
    def __init__(self):
        super(Myviewwidget, self).__init__()
        self.current_pos = None
        self.init_pos = None

    def wheelEvent(self, event):
        if self.num1 == 1:
            ...
        elif self.num1 == 0:
            return super().wheelEvent(event)

    def mousePressEvent(self, event):
        try:
            if self.num1 == 1:
                if event.button() == 1:     # 鼠标左键按下
                    self.init_pos = (event.pos().x(), event.pos().y())
                elif event.button() != 1:  # 不是鼠标左键按下
                    self.init_pos = None
            elif self.num1 == 0:
                return super().mousePressEvent(event)
        except Exception as e:
            print(e)

    def mouseReleaseEvent(self, event):
        try:
            if self.num1 == 1:
                self.signal_release.emit(self.init_pos, self.current_pos)
                self.current_pos = copy.deepcopy(self.init_pos)
                self.update()
                self.num1 = 0
                self.repaint()
                self.num1 = 1
                self.init_pos == None
                self.current_pos == None
            elif self.num1 == 0:
                return super().mousePressEvent(event)
        except Exception as e:
            print(e)

    def focusInEvent(self, event):
        if self.num1 == 1:
            print('focusInEvent')
        elif self.num1 == 0:
            return super().focusInEvent(event)

    def focusOutEvent(self, event):
        if self.num1 == 1:
            print('focusOutEvent')
        elif self.num1 == 0:
            return super().focusOutEvent(event)

    def mouseMoveEvent(self, event):
        try:
            if self.num1 == 1:
                print('mousemoveEvent')
                self.current_pos = (event.pos().x(), event.pos().y())
                self.update()   # 调用paintevent
                self.num1 = 0
                self.repaint()
                self.num1 = 1
            elif self.num1 == 0:
                return super().mouseMoveEvent(event)
        except Exception as e:
            print(e)

    def paintEvent(self, event):
        try:
            if self.num1 == 1:
                self.q = QtGui.QPainter(self)
                self.q.begin(self)
                self.drawRect(self.q)
                self.q.end()
            elif self.num1 == 0:
                return super().paintEvent(event)
        except Exception as e:
            print(e)

    def drawRect(self, qp):
        try:
            qp.setPen(QtGui.QColor(255, 255, 255))
            self.recwidth = self.current_pos[0] - self.init_pos[0]
            self.recheight = self.current_pos[1] - self.init_pos[1]
            self.num1 = 0
            self.repaint()
            qp.drawRect(self.init_pos[0], self.init_pos[1], self.recwidth, self.recheight)
            self.num1 = 1
        except Exception as e:
            print(e)

    def enterEvent(self, event):   # 鼠标移入label
        if self.num1 == 1:
            print('enterEvent')
        elif self.num1 == 0:
            return super().enterEvent(event)


class Ui_MainWindow(object):  # 主窗口
    """ The layout of the mainwindow. """
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.showFullScreen()
        font = QtGui.QFont()
        font.setPointSize(9)
        MainWindow.setFont(font)
        # toolbar
        self.toolbarBox1 = QtGui.QToolBar(self)
        self.toolbarBox2 = QtGui.QToolBar(self)
        self.addToolBar(self.toolbarBox1)
        self.addToolBar(self.toolbarBox2)
        self.sonwindow = QtWidgets.QWidget()
        MainWindow.setCentralWidget(self.sonwindow)
        layout = QtWidgets.QHBoxLayout()
        self.sonwindow.setLayout(layout)
        splitter1 = QtWidgets.QSplitter(QtCore.Qt.Horizontal)    # 水平
        splitter2 = QtWidgets.QSplitter(QtCore.Qt.Vertical)     # 垂直
        # parametertree
        self.tree = mytreewidget()
        self.tree.setColumnCount(4)
        self.tree.setHeaderLabels(['Files', 'Cell Formula', 'Heterostructure', 'Conductivity'])
        self.tree.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.tree.setStyleSheet("""
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
        # self.tab_tree.setMaximumSize(500, 1000)
        self.tab_tree_project = QtWidgets.QWidget()
        self.tab_tree.addTab(self.tab_tree_project, "Project")

        tab_tree_widgetlayout = QtWidgets.QVBoxLayout()
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
        tab_tree_widgetlayout.addWidget(self.tree)
        # 给projectTab设置一个搜索框
        searchlayout = QtWidgets.QHBoxLayout()
        self.searchLineEdit = QtWidgets.QLineEdit()
        self.searchLineEdit.setPlaceholderText("Search...")
        searchlayout.addWidget(self.searchLineEdit)
        tab_tree_widgetlayout.addLayout(searchlayout)
        self.tab_tree_project.setLayout(tab_tree_widgetlayout)
        # 左边控件 - object
        Verticallayout = QtWidgets.QVBoxLayout()
        self.calculate_distance = QtWidgets.QPushButton('Distance')
        self.calculate_degree = QtWidgets.QPushButton('Degree')
        self.calculate_vector = QtWidgets.QPushButton('Vector')
        Horizonlayout = QtWidgets.QHBoxLayout()
        Horizonlayout.addWidget(self.calculate_distance)
        Horizonlayout.addWidget(self.calculate_vector)
        Horizonlayout.addWidget(self.calculate_degree)
        Verticallayout.addLayout(Horizonlayout)
        Splitter_object = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        self.tablewidget = QtWidgets.QTableWidget()
        self.tablewidget.setStyleSheet("color: rgb(205, 205, 205);background-color: rgb(35, 38, 41);")
        self.tablewidget.setColumnCount(3)
        self.tablewidget.setHorizontalHeaderLabels(["Atom", 'Color', 'Atom radii(Å)'])
        self.tablewidget.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.tablewidget.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.tablewidget.setColumnWidth(1, 70)
        Splitter_object.addWidget(self.tablewidget)
        self.objectTextwidget = QtWidgets.QTextEdit()
        self.objectTextwidget.setStyleSheet("color: rgb(205, 205, 205);background-color: rgb(35, 38, 41);")
        Splitter_object.addWidget(self.objectTextwidget)
        Verticallayout.addWidget(Splitter_object)
        self.ballfilling_button = QtWidgets.QPushButton('Space filling')
        self.ball_stick_button = QtWidgets.QPushButton('Ball-and-stick')
        Horizonlayout1 = QtWidgets.QHBoxLayout()
        Horizonlayout1.addWidget(self.ballfilling_button)
        Horizonlayout1.addWidget(self.ball_stick_button)
        Verticallayout.addLayout(Horizonlayout1)
        self.tab_tree_Object = QtWidgets.QWidget()
        self.tab_tree.addTab(self.tab_tree_Object, "Object")
        self.tab_tree_Object.setLayout(Verticallayout)
        splitter1.addWidget(self.tab_tree)
        ## 右边的控件
        self.view_widget = Myviewwidget()
        self.view_widget.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        self.view_widget.setBackgroundColor((35, 38, 41, 0))
        self.widget_view_widget = QtWidgets.QWidget()
        view_widget_vertical_layout = QtWidgets.QVBoxLayout()
        view_widget_vertical_layout.addWidget(self.view_widget)
        self.widget_view_widget.setLayout(view_widget_vertical_layout)
        self.splitterview_widget = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        # self.splitterview_widget.addWidget(self.view_widget)
        self.splitterview_widget.addWidget(self.widget_view_widget)
        self.vertical_widget = QtWidgets.QWidget()
        database_vertical_layout = QtWidgets.QVBoxLayout()
        self.vertical_widget.setLayout(database_vertical_layout)
        self.datatable = QtWidgets.QTableWidget()
        self.datatable.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        # self.datatable.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)  ######允许右键产生子菜单
        # self.datatable.customContextMenuRequested.connect(self.generateMenu)  ####右键菜单
        self.datatable.setColumnCount(9)
        # self.datatable.setHorizontalHeaderLabels(['Formula', 'a(Å)', 'b(Å)', 'c(Å)', 'α(°)', 'β(°)', 'γ(°)', 'Volume'])
        self.datatable.setHorizontalHeaderLabels(['ID', 'Formula', 'a(Å)', 'b(Å)', 'c(Å)',
                                                  'α(°)', 'β(°)', 'γ(°)', 'Volume'])

        self.show_image_table = QtWidgets.QTableWidget()
        self.show_image_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)


        self.tab_view_widget = QtWidgets.QTabWidget()
        self.tab_view_widget_init = QtWidgets.QWidget()
        self.tab_view_widget.addTab(self.tab_view_widget_init, "MatLego")
        tab_view_widgetlayout = QtWidgets.QVBoxLayout()
        # database 里的搜索框
        searchdatabase_layout = QtWidgets.QHBoxLayout()
        self.searchdatabase_LineEdit = QtWidgets.QLineEdit()
        self.searchdatabase_LineEdit.setPlaceholderText("Search...")
        searchdatabase_layout.addWidget(self.searchdatabase_LineEdit)
        # databaseradiobutton
        Hlayout0 = QtWidgets.QHBoxLayout()
        self.btn1 = QtWidgets.QRadioButton("Crystal")
        self.btn1.setChecked(True)
        Hlayout0.addWidget(self.btn1)
        self.btn2 = QtWidgets.QRadioButton("Hetero-junction")
        self.btn2.setChecked(False)
        Hlayout0.addWidget(self.btn2)
        self.btn3 = QtWidgets.QRadioButton("Device")
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
        self.tab_view_widget_init.setLayout(tab_view_widgetlayout)
        splitter2.addWidget(self.tab_view_widget)
        # 右下
        self.tab = QtWidgets.QTabWidget()
        self.tab1 = QtWidgets.QWidget()
        self.tab2 = QtWidgets.QWidget()
        self.tab3 = QtWidgets.QWidget()
        self.tab.addTab(self.tab1, "Crystal information")
        self.tab.addTab(self.tab2, "Atom information")
        self.tab.addTab(self.tab3, "Hetero-junction/Device information")
        # tab1的布局

        self.text_widgt = QtWidgets.QTextEdit()
        self.text_widgt.setReadOnly(True)
        self.text_widgt.resize(24, 20)
        self.text_widgt.setStyleSheet("color: rgb(205, 205, 205);background-color: rgb(35, 38, 41);")

        # self.text_widgt.setStyleSheet("color: rgb(205, 205, 205);")

        tab1layout = QtWidgets.QVBoxLayout()
        tab1layout.addWidget(self.text_widgt)
        self.tab1.setLayout(tab1layout)
        # tab2的布局
        self.text_widgt2 = QtWidgets.QTextEdit()
        self.text_widgt2.setReadOnly(True)
        self.text_widgt2.resize(24, 20)
        self.text_widgt2.setStyleSheet("color: rgb(205, 205, 205);background-color: rgb(35, 38, 41);")
        # self.text_widgt2.setStyleSheet("color: rgb(205, 205, 205);background-color: rgb(0, 0, 0);")

        tab2layout = QtWidgets.QVBoxLayout()
        tab2layout.addWidget(self.text_widgt2)
        self.tab2.setLayout(tab2layout)
        # tab3的布局
        self.tab3_table_widget1 = QtWidgets.QTableWidget()
        self.tab3_table_widget1.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.tab3_table_widget1.setStyleSheet("color: rgb(255, 255, 255);background-color: rgb(35, 38, 41);")
        self.tab3_table_widget1.setColumnCount(5)
        self.tab3_table_widget1.setHorizontalHeaderLabels(['Hetero-junction', 'Optimal Match', '#Layer',
                                                           'Binding  energy (eV/MoS2)', 'Schottky barrier (eV)'])
        self.tab3_table_widget1.resizeColumnsToContents()
        self.tab3_table_widget2 = QtWidgets.QTableWidget()
        self.tab3_table_widget2.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.tab3_table_widget2.setStyleSheet("color: rgb(255, 255, 255);background-color: rgb(35, 38, 41);")
        self.tab3_table_widget2.setColumnCount(5)
        self.tab3_table_widget2.setHorizontalHeaderLabels(['Device', 'Optimal Match', '#Layer',
                                                           'Schottky barrier (eV)', 'I-V curve (Theory)'])
        self.tab3_table_widget2.resizeColumnsToContents()
        tab3layout = QtWidgets.QVBoxLayout()
        tab3splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)    # 水平
        tab3splitter.addWidget(self.tab3_table_widget1)
        tab3splitter.addWidget(self.tab3_table_widget2)
        tab3layout.addWidget(tab3splitter)
        self.tab3.setLayout(tab3layout)

        # splitter2.addWidget(self.text_widgt)
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
        # self.menuBar.setGeometry(QtCore.QRect(0, 0, 1274, 30))
        self.document_menu = QtWidgets.QMenu(self.menuBar)
        # self.document_menu.setGeometry(QtCore.QRect(0, 0, 800, 300))
        self.function_menu = QtWidgets.QMenu(self.menuBar)
        self.view_menu = QtWidgets.QMenu(self.menuBar)
        self.database_menu = QtWidgets.QMenu(self.menuBar)
        self.Setting_menu = QtWidgets.QMenu(self.menuBar)
        self.help_menu = QtWidgets.QMenu(self.menuBar)
        MainWindow.setMenuBar(self.menuBar)
        self.mainToolBar = QtWidgets.QToolBar(MainWindow)
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.mainToolBar)
        # New
        self.actionNew = QtWidgets.QAction(MainWindow)
        self.actionNewToolBar = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\New.png"), "New", MainWindow)

        self.toolbarBox1.addAction(self.actionNewToolBar)
        # Open
        self.actionOpen = QtWidgets.QAction(MainWindow)
        self.actionOpenToolBar = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\open1.png"), "Open", MainWindow)
        self.toolbarBox1.addAction(self.actionOpenToolBar)
        # Save
        self.actionSave = QtWidgets.QAction(MainWindow)
        # exit
        self.actionExit = QtWidgets.QAction(MainWindow)
        # Cut
        self.actionCut = QtWidgets.QAction(MainWindow)
        self.actionTraversal_Cut = QtWidgets.QAction(MainWindow)
        # Cut ToolBar
        self.actionCutToolBar = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\Cut7.png"), "Exfoliate", MainWindow)
        self.toolbarBox1.addAction(self.actionCutToolBar)
        # Judge3layer
        self.actionJudge2layer = QtWidgets.QAction(MainWindow)
        self.actionJudge3layer = QtWidgets.QAction(MainWindow)
        self.actionJudgecustomize = QtWidgets.QAction(MainWindow)
        #
        self.actionGuide = QtWidgets.QAction(MainWindow)
        self.action_setting_atom = QtWidgets.QAction(MainWindow)
        # Stack
        self.actionStack = QtWidgets.QAction(MainWindow)
        self.actionStackToolbar = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\Stack2.png"), "Stack", MainWindow)
        self.toolbarBox1.addAction(self.actionStackToolbar)
        # Rotate
        self.actionRotate = QtWidgets.QAction(MainWindow)
        self.actionRotate_cutomize = QtWidgets.QAction(MainWindow)
        # drag_rectangle
        self.actiondrag_rectangle = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\dragrectangle0.png"), "Drag rectangle", MainWindow)
        self.actionNormalchoose = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\drag6.png"), "Drag", MainWindow)
        self.actionSetlayer = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\setlayer.png"), "Set layer", MainWindow)
        self.actionMovelayer = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\movelayer.png"), "Movelayer", MainWindow)
        self.actionleft_screen = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\left.png"), "left screen", MainWindow)
        self.actionright_screen = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\right.png"), "right screen", MainWindow)
        self.actionup_screen = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\up.png"), "up screen", MainWindow)
        self.actiondown_screen = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\down.png"), "down screen", MainWindow)
        self.actionmiddle_screen = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\middle.png"), "middle screen", MainWindow)
        self.toolbarBox1.addAction(self.actiondrag_rectangle)
        self.toolbarBox1.addAction(self.actionNormalchoose)
        self.toolbarBox1.addAction(self.actionSetlayer)
        self.toolbarBox1.addAction(self.actionMovelayer)
        self.toolbarBox2.addAction(self.actionleft_screen)
        self.toolbarBox2.addAction(self.actionright_screen)
        self.toolbarBox2.addAction(self.actionup_screen)
        self.toolbarBox2.addAction(self.actiondown_screen)
        self.toolbarBox2.addAction(self.actionmiddle_screen)
        self.actionTraversal_rotate = QtWidgets.QAction(MainWindow)
        self.actionTraversal_Stack1 = QtWidgets.QAction(MainWindow)
        self.actionTraversal_Stack2 = QtWidgets.QAction(MainWindow)
        self.actionClassification = QtWidgets.QAction(MainWindow)
        self.actionImport = QtWidgets.QAction(MainWindow)
        self.actionExport = QtWidgets.QAction(MainWindow)
        self.actiona = QtWidgets.QAction(MainWindow)
        self.actionb = QtWidgets.QAction(MainWindow)
        self.actionc = QtWidgets.QAction(MainWindow)
        self.actionViewCell = QtWidgets.QAction(QtGui.QIcon(r"D:\coding_python\2DMaterial\src\Photo\dot1.png"), "View cell", MainWindow)
        self.actionViewCell.setIconVisibleInMenu(True)
        self.actionViewCoordinate_System = QtWidgets.QAction(MainWindow)
        self.actionViewCoordinate_cell = QtWidgets.QAction(MainWindow)
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
        # self.function_menu.addSeparator()
        self.function_menu.addAction(self.actionStack)
        self.function_menu.addAction(self.actionRotate)
        self.function_menu.addAction(self.actionRotate_cutomize)

        self.function_menu.addSeparator()
        self.function_menu.addAction(self.actionTraversal_rotate)
        self.function_menu.addAction(self.actionTraversal_Stack1)
        self.function_menu.addAction(self.actionTraversal_Stack2)
        self.function_menu.addAction(self.actionTraversal_Cut)
        self.Judge = self.function_menu.addMenu("Judge")
        self.Judge.addAction(self.actionJudge2layer)
        self.Judge.addAction(self.actionJudge3layer)
        self.Judge.addAction(self.actionJudgecustomize)
        self.function_menu.addSeparator()
        self.function_menu.addAction(self.actionClassification)
        self.function_menu.addSeparator()
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
        self.view_from.addAction(self.actiona)
        self.view_from.addAction(self.actionb)
        self.view_from.addAction(self.actionc)
        self.view_menu.addAction(self.actionViewCell)
        self.coordinate_system_menu = self.view_menu.addMenu("Coordinate System")
        self.coordinate_system_menu.addAction(self.actionViewCoordinate_System)
        self.coordinate_system_menu.addAction(self.actionViewCoordinate_cell)
        self.view_menu.addAction(self.actionViewGrid)
        # help menu
        self.help_menu.addAction(self.actionGuide)
        # Setting menu
        self.Setting_menu.addAction(self.action_setting_atom)
        # database_menu
        self.action_db_connect = QtWidgets.QAction(MainWindow)
        self.action_add_data =  QtWidgets.QAction(MainWindow)
        self.action_cif_to_db = QtWidgets.QAction(MainWindow)
        self.action_db_to_cifs = QtWidgets.QAction(MainWindow)
        self.database_menu.addAction(self.action_db_connect)
        self.database_menu.addAction(self.action_add_data)
        self.database_menu.addAction(self.action_cif_to_db)
        self.database_menu.addAction(self.action_db_to_cifs)
        #
        self.menuBar.addAction(self.document_menu.menuAction())
        self.menuBar.addAction(self.function_menu.menuAction())
        self.menuBar.addAction(self.view_menu.menuAction())
        self.menuBar.addAction(self.help_menu.menuAction())
        self.menuBar.addAction(self.Setting_menu.menuAction())
        self.menuBar.addAction(self.database_menu.menuAction())
        self.init_database()
        # calculate menu
        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)


    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.document_menu.setTitle(_translate("MainWindow", "File"))
        self.function_menu.setTitle(_translate("MainWindow", "Function"))
        self.view_menu.setTitle(_translate("MainWindow", "View"))
        self.database_menu.setTitle(_translate("MainWindow", "Database"))
        self.help_menu.setTitle(_translate("MainWindow", "Help"))
        self.Setting_menu.setTitle(_translate("MainWindow", "Settings"))
        # Judge2layer
        self.actionJudge2layer.setText(_translate("MainWindow", "Two - layer electronic device"))
        self.actionJudge3layer.setText(_translate("MainWindow", "Three - layer electronic device"))
        self.actionJudgecustomize.setText(_translate("MainWindow", "Customize"))
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
        self.actionRotate_cutomize.setText(_translate("MainWindow", "Rotate(customize)"))

        # actionCut
        self.actionCut.setText(_translate("MainWindow", "Exfoliate"))
        # actionGuide
        self.actionGuide.setText(_translate("MainWindow", "Guide"))
        # Setting
        self.action_setting_atom.setText(_translate("MainWindow", "Atom parameter"))

        # actionStack
        self.actionStack.setText(_translate("MainWindow", "Stack"))
        # actionTraversal
        self.actionTraversal_rotate.setText(_translate("MainWindow", "Traversal-Rotate"))

        self.actionTraversal_Stack1.setText(_translate("MainWindow", "Traversal-Stack(1)"))
        self.actionTraversal_Stack2.setText(_translate("MainWindow", "Traversal-Stack(2)"))

        # TraversalCut
        self.actionTraversal_Cut.setText(_translate("MainWindow", "Traversal-Cut"))
        # actionClassification
        self.actionClassification.setText(_translate("MainWindow", "Classification"))
        # actionImport
        self.actionImport.setText(_translate("MainWindow", "Import"))
        # actionExport
        self.actionExport.setText(_translate("MainWindow", "Export"))

        self.actionViewProject.setText(_translate("MainWindow", "Project"))
        self.actiona.setText(_translate("MainWindow", "a-axis"))
        self.actionb.setText(_translate("MainWindow", "b-axis"))
        self.actionc.setText(_translate("MainWindow", "c-axis"))
        self.actionViewCell.setText(_translate("MainWindow", "Cell"))
        self.actionViewCoordinate_System.setText(_translate("MainWindow", "Cartesian"))
        self.actionViewCoordinate_cell.setText(_translate("MainWindow", "Cell Coordinate"))
        self.actionViewGrid.setText(_translate("MainWindow", "Grid"))
        self.actionViewDatabase.setText(_translate("MainWindow", "Database"))
        self.actionViewText.setText(_translate("MainWindow", "Text"))
        self.actionOrthogonal.setText(_translate("MainWindow", "Orthogonal viewport"))
        self.actionPerspective.setText(_translate("MainWindow", "Perspective viewport"))
        self.action_db_connect.setText(_translate("MainWindow", "Connect"))
        self.action_add_data.setText(_translate("MainWindow", "Add CIF"))
        self.action_cif_to_db.setText(_translate("MainWindow", "Create a database"))
        self.action_db_to_cifs.setText(_translate("MainWindow", "Composite CIF file"))

    def init_database(self):
        self.datatable.setRowCount(10)  # 添加信息
        self.datatable.setShowGrid(True)
        self.vertical_widget.setMaximumWidth(1200)
        # self.vertical_widget.(1000, 550)


class create_database_window(QtWidgets.QWidget):
    Signal_emit_number = pyqtSignal(float)
    def __init__(self, parent=None):
        super(create_database_window, self).__init__(parent)
        layoutH = QtWidgets.QHBoxLayout(self)
        label = QtWidgets.QLabel("Progress:")
        self.pbar = QtWidgets.QProgressBar()
        self.pbar.setValue(15)
        layoutH.addWidget(label)
        layoutH.addWidget(self.pbar)
        self.setLayout(layoutH)
        self.show()

    def pbar_number(self, value):
        self.pbar.setValue(value)

class Append_data_window(QtWidgets.QWidget):
    signal_emit = pyqtSignal(int, str)
    def __init__(self, ID_lst, parent=None):
        super(Append_data_window, self).__init__(parent)
        self.ID_lst = ID_lst
        Hlayout = QtWidgets.QHBoxLayout()
        label1 = QtWidgets.QLabel("ID:")
        self.linedit1 = QtWidgets.QLineEdit()
        label2 = QtWidgets.QLabel("Formula:")
        self.linedit2 = QtWidgets.QLineEdit()
        Hlayout.addWidget(label1)
        Hlayout.addWidget(self.linedit1)
        Hlayout.addWidget(label2)
        Hlayout.addWidget(self.linedit2)
        self.button = QtWidgets.QPushButton("OK")
        self.button.clicked.connect(self.determine)
        Hlayout.addWidget(self.button)
        self.setLayout(Hlayout)
        self.show()

    def determine(self):
        try:
            self.ID = int(self.linedit1.text())
            if self.ID not in self.ID_lst:
                self.Formula = self.linedit2.text()
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
    signal_emit_information = pyqtSignal(int, list, list)
    signal_emit_delete = pyqtSignal(int)
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
            for i in range(len(self.information_lst)):
                newitem = QtWidgets.QTableWidgetItem(self.information_lst[i])
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
        self.setLayout(layout1)
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

class edit_atom_para_window(QtWidgets.QWidget):
    signal_emit_atom_par = pyqtSignal(list, list)
    def __init__(self, color_lst, radii_lst, parent=None):
        super(edit_atom_para_window, self).__init__(parent)
        self.color_lst = copy.deepcopy(color_lst)
        self.radii_lst = copy.deepcopy(radii_lst)
        self.color_lst.pop(0)
        self.radii_lst.pop(0)
        Hlayout = QtWidgets.QHBoxLayout()
        self.button_cancel = QtWidgets.QPushButton("Cancel")
        self.button_cancel.clicked.connect(self.close)
        self.button_ok = QtWidgets.QPushButton("OK")
        self.button_ok.clicked.connect(self.determine)
        Hlayout.addWidget(self.button_cancel)
        Hlayout.addWidget(self.button_ok)
        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(Hlayout)
        self.tablewidget = QtWidgets.QTableWidget()
        self.tablewidget.setColumnCount(2)
        self.tablewidget.setHorizontalHeaderLabels(['Color', 'Atom radii'])
        self.tablewidget.resizeColumnsToContents()
        self.formula_lst = copy.deepcopy(list(sf.crys_data.atom_formula))
        self.formula_lst.pop(0)
        self.atom_number = len(self.formula_lst)
        self.tablewidget.setRowCount(self.atom_number)
        self.tablewidget.setVerticalHeaderLabels(self.formula_lst)
        layout.addWidget(self.tablewidget)
        self.setLayout(layout)
        self.tablewidget.itemDoubleClicked.connect(self.double_clicked_on_table)
        self.init()
        self.show()

    def init(self):
        for i in range(0, self.atom_number):
            try:
                newitem_color = QtWidgets.QTableWidgetItem("")
                newitem_radii = QtWidgets.QTableWidgetItem(str(self.radii_lst[i]))
                color = self.color_lst[i]
                colortrans = (color[0] * 255, color[1] * 255, color[2] * 255)
                self.tablewidget.setItem(i, 0, newitem_color)
                self.tablewidget.setItem(i, 1, newitem_radii)
                self.tablewidget.item(i, 0).setBackground(QtGui.QBrush(
                    QtGui.QColor(colortrans[0], colortrans[1], colortrans[2])))
            except Exception as e:
                print(e)

    def double_clicked_on_table(self, item):
        try:
            self.row = item.row()
            color = self.tablewidget.item(self.row, 0).background().color().getRgb()
            print("color = ", color)
            # color = hex(color)
            print("color = ", color)
            radii = float(self.tablewidget.item(item.row(), 1).text())
            self.edit_information_window = edit_atom_para(radii, color)
            self.edit_information_window.signal_emit_change.connect(self.change)
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error',
                                          'Connect database please.')

    def change(self, lst):
        try:
            self.tablewidget.item(self.row, 0).setBackground(lst[1])
            Rgb = list(lst[1].getRgb())
            self.color_lst[self.row] = (Rgb[0]/255, Rgb[1] / 255, Rgb[2]/255, 1)
        except Exception as e:
            print(e)
        try:
            self.tablewidget.item(self.row, 1).setText(str(lst[0]))
            self.radii_lst[self.row] = lst[0]
        except Exception as e:
            print(e)

    def determine(self):
        try:
            self.color_lst.insert(0, (0, 0, 0, 0))
            self.radii_lst.insert(0, np.nan)
            self.signal_emit_atom_par.emit(self.color_lst, self.radii_lst)
            self.close()
        except Exception as e:
            print(e)




class edit_atom_para(QtWidgets.QWidget):
    signal_emit_change = pyqtSignal(list)
    def __init__(self, radii, color, parent=None):
        super(edit_atom_para, self).__init__(parent)
        self.radii = radii
        self.color = color
        print("color = ", self.color)
        self.params = [
            {'name': 'Basic parameter data types', 'type': 'group', 'children': [
                {'name': 'Atom radii', 'type': 'float', 'value': self.radii, 'step': 0.1},
                {'name': 'Color', 'type': 'color', 'value': self.color, 'tip': "This is a color button"}
            ]
             }]
        self.p = Parameter.create(name='params', type='group', children=self.params)
        self.p.sigTreeStateChanged.connect(self.change)
        self.t = ParameterTree()
        self.t.setParameters(self.p, showTop=False)
        layout = QtWidgets.QVBoxLayout()
        self.setLayout(layout)
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
                # self.color = data.getRgb()
                self.color = data

                print('self.color =', self.color)
                print('  parameter: %s' % childName)
                print('  change:    %s' % change)
                print('  data:      %s' % str(data))
                print('  ----------')
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
    signal_emit_ID = pyqtSignal(str, str, int)
    def __init__(self, cif_dir, database_dir, parent=None):
        super(add_data_window, self).__init__(parent)
        self.cif_dir = cif_dir
        self.database_dir = database_dir
        layout1 = QtWidgets.QHBoxLayout(self)
        label = QtWidgets.QLabel("ID：")
        self.ID_widget = QtWidgets.QLineEdit()
        self.button1 = QtWidgets.QPushButton('Ok')
        self.button1.clicked.connect(self.determine)
        layout1.addWidget(label)
        layout1.addWidget(self.ID_widget)
        layout1.addWidget(self.button1)
        self.setLayout(layout1)
        self.show()

    def determine(self):
        try:
            ID = int(self.ID_widget.text())
            self.signal_emit_ID.emit(self.cif_dir, self.database_dir, ID)
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, 'error', 'ID input error!')


class database_to_cif(QtWidgets.QWidget):
    signal_emit_formula = pyqtSignal(str, str, str)
    def __init__(self, database_dir, parent=None):
        super(database_to_cif, self).__init__(parent)
        self.database_dir = database_dir
        layout = QtWidgets.QHBoxLayout(self)
        label = QtWidgets.QLabel("Formula：")
        self.button = QtWidgets.QPushButton('Ok')
        self.button.clicked.connect(self.determine)
        self.formula_widget = QtWidgets.QLineEdit()
        layout.addWidget(label)
        layout.addWidget(self.formula_widget)
        layout.addWidget(self.button)
        self.setLayout(layout)
        self.show()

    def determine(self):
        cif_dir = QtGui.QFileDialog.getExistingDirectory(self, 'Select cif files directory')
        self.signal_emit_formula.emit(self.database_dir, cif_dir + '/' + str(self.formula_widget.text()) + '.cif',
                                      str(self.formula_widget.text()))
        self.close()


class customize_judge(QtWidgets.QWidget):
    signal_emit_conductivity = pyqtSignal(list)
    def __init__(self, parent=None):
        super(customize_judge, self).__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)  # 垂直布局
        layout0 = QtWidgets.QHBoxLayout(self)
        label = QtWidgets.QLabel("Number of layers：")
        self.layernumber = QtWidgets.QLineEdit()
        self.button = QtWidgets.QPushButton('Ok')
        layout0.addWidget(label)
        layout0.addWidget(self.layernumber)
        layout0.addWidget(self.button)
        self.tablewidget = QtWidgets.QTableWidget()
        self.tablewidget.setColumnCount(1)
        self.tablewidget.setHorizontalHeaderLabels(['Conductivity (from bottom to top)'])
        self.tablewidget.resizeColumnsToContents()
        layout.addLayout(layout0)
        layout.addWidget(self.tablewidget)
        self.setLayout(layout)
        self.layernumber.textChanged.connect(self.tablechange)
        self.button.clicked.connect(self.determine)
        self.show()

    def tablechange(self):
        try:
            if self.layernumber.text() == "":
                self.tablewidget.setRowCount(0)
            else:
                num = int(self.layernumber.text())
                self.tablewidget.setRowCount(num)
                self.l = []
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
            l = []
            for i in range(self.tablewidget.rowCount()):
                l.append(self.l[i].currentText())
            try:
                self.signal_emit_conductivity.emit(l)
            except Exception as e:
                print(e)
        except Exception as e:
            print(e)


class create_2D_devices(QtWidgets.QWidget):
    """ The layout of create 2D devices"""
    signal_emit_conductivity = pyqtSignal(list)
    signal_emit_vacuum_layer = pyqtSignal(float)
    signal_devices = pyqtSignal(float)
    def __init__(self, parent=None):
        super(create_2D_devices, self).__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)  # 垂直布局
        layout0 = QtWidgets.QHBoxLayout(self)
        layout1 = QtWidgets.QHBoxLayout(self)
        layout2 = QtWidgets.QHBoxLayout(self)
        self.pbar = QtWidgets.QProgressBar()
        self.button = QtWidgets.QPushButton('Start')
        layout2.addWidget(self.pbar)
        layout2.addWidget(self.button)
        label = QtWidgets.QLabel("Number of layers：")
        self.layernumber = QtWidgets.QLineEdit()
        layout0.addWidget(label)
        layout0.addWidget(self.layernumber)
        # layout0.addWidget(self.button)
        label1 = QtWidgets.QLabel("Vacuum distance(Å) ：")
        self.vacuum_distance = QtWidgets.QLineEdit()
        layout1.addWidget(label1)
        layout1.addWidget(self.vacuum_distance)
        self.tablewidget = QtWidgets.QTableWidget()
        self.resize(530, 400)
        self.tablewidget.setColumnCount(2)
        self.tablewidget.setHorizontalHeaderLabels(['Conductivity (from bottom to top)', 'layer distance(Å)'])
        self.tablewidget.resizeColumnsToContents()
        layout.addLayout(layout0)
        layout.addLayout(layout1)
        layout.addLayout(layout2)
        layout.addWidget(self.tablewidget)
        self.setLayout(layout)
        self.layernumber.textChanged.connect(self.tablechange)
        self.button.clicked.connect(self.determine)
        self.signal_devices.connect(self.pbar_set)
        self.show()

    def pbar_set(self, value):
        self.pbar.setValue(value)

    def tablechange(self):
        try:
            if self.layernumber.text() == "":
                self.tablewidget.setRowCount(0)
            else:
                num = int(self.layernumber.text())
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
            try:
                vacuum_dis = float(self.vacuum_distance.text())
                self.signal_emit_vacuum_layer.emit(vacuum_dis)
            except Exception as e:
                print(e)
            l = []
            layer_dis_lst = []
            for i in range(self.tablewidget.rowCount()):
                l.append(self.l[i].currentText())
            try:
                for i in range(1, self.tablewidget.rowCount()):
                    layer_dis = float(self.tablewidget.item(i, 1).text())
                    layer_dis_lst.append(layer_dis)
                self.signal_emit_conductivity.emit([l, layer_dis_lst])
            except Exception as e:
                QtWidgets.QMessageBox.warning(self, 'error', 'layer distance Input error!')
        except Exception as e:
            print(e)


class tree_set_text(QtWidgets.QWidget):
    """ The layout of the self.tree set text window. """
    signal_edit_text = pyqtSignal(str)
    signal_edit_text1 = pyqtSignal(str)
    def __init__(self, parent=None):
        super(tree_set_text, self).__init__(parent)
        # layout = QtWidgets.QVBoxLayout(self)
        layout0 = QtWidgets.QHBoxLayout(self)  # 水平布局
        # layout1 = QtWidgets.QHBoxLayout(self)
        label = QtWidgets.QLabel("Set Text ：")
        self.treetext = QtWidgets.QComboBox()
        self.treetext.addItem("bulk")
        self.treetext.addItem("layer")
        self.treetext.addItem("stack")
        label1 = QtWidgets.QLabel("Set Text ：")
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
        layout0.addWidget(label1)
        layout0.addWidget(self.treetext1)
        layout0.addWidget(self.determine_button1)
        self.determine_button.clicked.connect(self.determine)
        self.determine_button1.clicked.connect(self.determine1)
        self.setLayout(layout0)
        self.show()

    def determine1(self):
        try:
            text = self.treetext1.currentText()
            self.signal_edit_text1.emit(text)
            # self.close()
        except Exception as e:
            print(e)

    def determine(self):
        try:
            text = self.treetext.currentText()
            self.signal_edit_text.emit(text)
            # self.close()
        except Exception as e:
            print(e)

class Stack_main_one_window(QtWidgets.QWidget):
    signal_emit = pyqtSignal(float, float, float)
    def __init__(self, parent=None):
        super(Stack_main_one_window, self).__init__(parent)
        Vlayout = QtWidgets.QVBoxLayout(self)
        self.setWindowTitle("A stacking direction is known(a1, a2).")
        self.button = QtWidgets.QPushButton('Ok')
        self.button.clicked.connect(self.determine)
        Hlayout1  = QtWidgets.QHBoxLayout(self)
        label1 = QtWidgets.QLabel("Error tolerance(%)：")
        self.strain_var_widget = QtWidgets.QLineEdit()
        Hlayout1.addWidget(label1)
        Hlayout1.addWidget(self.strain_var_widget)
        Hlayout1.addWidget(self.button)
        #
        Hlayout2 = QtWidgets.QHBoxLayout(self)
        label2 = QtWidgets.QLabel("layer distance(Å)：")
        self.layer_distance_widget = QtWidgets.QLineEdit()
        Hlayout2.addWidget(label2)
        Hlayout2.addWidget(self.layer_distance_widget)
        Hlayout3 = QtWidgets.QHBoxLayout(self)
        label3 = QtWidgets.QLabel("Vacuum distance(Å)：")
        self.vacuum_distance_widget = QtWidgets.QLineEdit()
        Hlayout3.addWidget(label3)
        Hlayout3.addWidget(self.vacuum_distance_widget)
        Vlayout.addLayout(Hlayout2)
        Vlayout.addLayout(Hlayout3)
        Vlayout.addLayout(Hlayout1)
        self.setLayout(Vlayout)
        self.show()

    def determine(self):
        try:
            self.strain_var = float(self.strain_var_widget.text()) / 100
            self.vacuum_distance_var = float(self.vacuum_distance_widget.text())
            self.layer_distance_var = float(self.layer_distance_widget.text())
            self.signal_emit.emit(self.strain_var, self.layer_distance_var, self.vacuum_distance_var)
            self.close()
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')
            print(e)

class traversal_stack(QtWidgets.QWidget):
    """ The layout of the traversal_stack window. """
    signaltraversalstack = pyqtSignal(int)
    signal_start_pbar = pyqtSignal(float)
    side_known = ...
    def __init__(self, parent=None):
        super(traversal_stack, self).__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)
        Hlayout0 = QtWidgets.QHBoxLayout(self)  # 水平布局
        Hlayout1 = QtWidgets.QHBoxLayout(self)  # 水平布局
        Hlayout2 = QtWidgets.QHBoxLayout(self)  # 水平布局
        Hlayout3 = QtWidgets.QHBoxLayout(self)  # 水平布局
        Hlayout4 = QtWidgets.QHBoxLayout(self)  # 水平布局
        Hlayout5 = QtWidgets.QHBoxLayout(self)  # 水平布局
        self.btn1 = QtWidgets.QRadioButton("layer-layer")
        self.btn1.setChecked(False)
        self.btn1.toggled.connect(lambda: self.btnstate(self.btn1))
        Hlayout0.addWidget(self.btn1)
        self.btn2 = QtWidgets.QRadioButton("layer-stack")
        self.btn2.setChecked(False)
        self.btn2.toggled.connect(lambda: self.btnstate(self.btn2))
        Hlayout0.addWidget(self.btn2)
        self.btn3 = QtWidgets.QRadioButton("All")
        self.btn3.setChecked(False)
        self.btn3.toggled.connect(lambda: self.btnstate(self.btn3))
        Hlayout0.addWidget(self.btn3)
        label1 = QtWidgets.QLabel("Set vacuum distance(Å) ：")
        self.vacuum_distance = QtWidgets.QLineEdit()
        Hlayout1.addWidget(label1)
        Hlayout1.addWidget(self.vacuum_distance)
        label2 = QtWidgets.QLabel("Set layer distance(Å)  ：")
        self.layer_distance = QtWidgets.QLineEdit()
        label3 = QtWidgets.QLabel("Set error tolerance(%)  ：")
        self.tolerance = QtWidgets.QLineEdit()
        self.determine_button = QtWidgets.QPushButton('Start')
        self.determine_button.clicked.connect(self.determine)
        self.pbar = QtWidgets.QProgressBar()
        self.Text = QtWidgets.QTextEdit()
        Hlayout2.addWidget(label2)
        Hlayout2.addWidget(self.layer_distance)
        Hlayout3.addWidget(label3)
        Hlayout3.addWidget(self.tolerance)
        Hlayout4.addWidget(self.pbar)
        Hlayout4.addWidget(self.determine_button)
        Hlayout5.addWidget(self.Text)
        layout.addLayout(Hlayout0)
        layout.addLayout(Hlayout1)
        layout.addLayout(Hlayout2)
        layout.addLayout(Hlayout3)
        layout.addLayout(Hlayout4)
        layout.addLayout(Hlayout5)
        self.signal_start_pbar.connect(self.start_pbar)
        self.show()

    def start_pbar(self, num):
        self.pbar.setValue(num)

    def btnstate(self, btn):
        try:
            self.choosebutton = btn
        except Exception as e:
            print(e)

    def determine(self):
        try:
            test = 0
            self.layerdistance_var = float(self.layer_distance.text())
            test = 1
            self.vacuumdistance_var = float(self.vacuum_distance.text())
            test = 2
            self.tolerance_var = float(self.tolerance.text())
            test = 3
            if self.choosebutton.text() == "layer-layer":
                self.signaltraversalstack.emit(1)
            elif self.choosebutton.text() == 'layer-stack':
                self.signaltraversalstack.emit(2)
            elif self.choosebutton.text() == 'All':
                self.signaltraversalstack.emit(3)
            # self.close()

        except Exception as e:
            print(e)
            if test == 0:
                QtWidgets.QMessageBox.warning(self, 'error', 'layer distance Input error!')
            elif test == 1:
                QtWidgets.QMessageBox.warning(self, 'error', 'vacuum distance Input error!')
            elif test == 2:
                QtWidgets.QMessageBox.warning(self, 'error', 'error tolerance Input error!')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please Choose in RadioButtons')

class traversal_rotate(QtWidgets.QWidget):
    """ The layout of the traversal_stack window. """
    signaltraversalstack = pyqtSignal(int)
    signal_start_pbar = pyqtSignal(float)
    side_known = ...
    def __init__(self, parent=None):
        super(traversal_rotate, self).__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)
        Hlayout0 = QtWidgets.QHBoxLayout(self)  # 水平布局
        Hlayout1 = QtWidgets.QHBoxLayout(self)  # 水平布局
        Hlayout2 = QtWidgets.QHBoxLayout(self)  # 水平布局
        Hlayout3 = QtWidgets.QHBoxLayout(self)  # 水平布局
        Hlayout4 = QtWidgets.QHBoxLayout(self)  # 水平布局
        Hlayout5 = QtWidgets.QHBoxLayout(self)  # 水平布局
        Hlayout3_4 = QtWidgets.QHBoxLayout(self)
        self.btn1 = QtWidgets.QRadioButton("layer-layer")
        self.btn1.setChecked(False)
        self.btn1.toggled.connect(lambda: self.btnstate(self.btn1))
        Hlayout0.addWidget(self.btn1)
        self.btn2 = QtWidgets.QRadioButton("layer-stack")
        self.btn2.setChecked(False)
        self.btn2.toggled.connect(lambda: self.btnstate(self.btn2))
        Hlayout0.addWidget(self.btn2)
        self.btn3 = QtWidgets.QRadioButton("All")
        self.btn3.setChecked(False)
        self.btn3.toggled.connect(lambda: self.btnstate(self.btn3))
        Hlayout0.addWidget(self.btn3)
        label1 = QtWidgets.QLabel("Set vacuum distance(Å) ：")
        self.vacuum_distance = QtWidgets.QLineEdit()
        Hlayout1.addWidget(label1)
        Hlayout1.addWidget(self.vacuum_distance)
        label2 = QtWidgets.QLabel("Set layer distance(Å)  ：")
        self.layer_distance = QtWidgets.QLineEdit()
        label3 = QtWidgets.QLabel("Set error tolerance(%)  ：")
        self.tolerance = QtWidgets.QLineEdit()
        label4 = QtWidgets.QLabel("Set a,b max length(Å)  ：")
        self.max_length = QtWidgets.QLineEdit()

        self.layer_distance = QtWidgets.QLineEdit()
        self.determine_button = QtWidgets.QPushButton('Start')
        self.determine_button.clicked.connect(self.determine)
        self.pbar = QtWidgets.QProgressBar()
        self.Text = QtWidgets.QTextEdit()
        Hlayout2.addWidget(label2)
        Hlayout2.addWidget(self.layer_distance)
        Hlayout3.addWidget(label3)
        Hlayout3.addWidget(self.tolerance)
        Hlayout3_4.addWidget(label4)
        Hlayout3_4.addWidget(self.max_length)
        Hlayout4.addWidget(self.pbar)
        Hlayout4.addWidget(self.determine_button)
        Hlayout5.addWidget(self.Text)
        layout.addLayout(Hlayout0)
        layout.addLayout(Hlayout1)
        layout.addLayout(Hlayout2)
        layout.addLayout(Hlayout3)
        layout.addLayout(Hlayout3_4)
        layout.addLayout(Hlayout4)
        layout.addLayout(Hlayout5)
        self.signal_start_pbar.connect(self.start_pbar)
        self.show()

    def start_pbar(self, num):
        self.pbar.setValue(num)

    def btnstate(self, btn):
        try:
            self.choosebutton = btn
        except Exception as e:
            print(e)

    def determine(self):
        try:
            test = 0
            self.layerdistance_var = float(self.layer_distance.text())
            test = 1
            self.vacuumdistance_var = float(self.vacuum_distance.text())
            test = 2
            self.tolerance_var = float(self.tolerance.text())
            test = 3
            self.max_length_var = float(self.max_length.text())
            if self.choosebutton.text() == "layer-layer":
                self.signaltraversalstack.emit(1)
            elif self.choosebutton.text() == 'layer-stack':
                self.signaltraversalstack.emit(2)
            elif self.choosebutton.text() == 'All':
                self.signaltraversalstack.emit(3)
            # self.close()

        except Exception as e:
            print(e)
            if test == 0:
                QtWidgets.QMessageBox.warning(self, 'error', 'layer distance Input error!')
            elif test == 1:
                QtWidgets.QMessageBox.warning(self, 'error', 'vacuum distance Input error!')
            elif test == 2:
                QtWidgets.QMessageBox.warning(self, 'error', 'error tolerance Input error!')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Please Choose in RadioButtons')


class scatter_plot_customize(QtWidgets.QWidget):
    len_max = 15
    error_limit = 0.02
    emit_information_signal = QtCore.pyqtSignal(float, list, list, float, float)
    def __init__(self, crys1, crys2, parent=None):
        super(scatter_plot_customize, self).__init__(parent)
        pg.setConfigOption('background', (35, 38, 41))
        self.resize(1000, 2000)
        self.crys1 = crys1  # down
        self.crys2 = crys2  # up
        layout = QtWidgets.QVBoxLayout(self)  # 垂直布局
        Hlayout = QtWidgets.QHBoxLayout(self)
        splitter2 = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        self.determine_button = QtWidgets.QPushButton('OK')
        self.determine_button2 = QtWidgets.QPushButton('OK')
        self.error_limit_widget = QtWidgets.QLineEdit()
        self.error_limit_widget.setPlaceholderText(str(self.error_limit * 100))
        error_label = QtWidgets.QLabel("Error tolerance(%)：：")
        error_label.setBuddy(self.error_limit_widget)
        self.max_length_widget = QtWidgets.QLineEdit()
        self.max_length_widget.setPlaceholderText(str(self.len_max))
        max_length_label = QtWidgets.QLabel("Max length(Å): ")
        max_length_label.setBuddy(self.max_length_widget)

        self.view = pg.GraphicsLayoutWidget()  ## GraphicsView with GraphicsLayout inserted by default
        self.p1 = self.view.addPlot()
        self.vb = self.p1.vb
        self.s1 = pg.ScatterPlotItem(size=10, pen=pg.mkPen(None))
        # self.s1.sigClicked.connect(self.clicked)
        self.text_widgt = QtWidgets.QTextEdit()
        self.text_widgt.setReadOnly(True)
        self.text_widgt.resize(24, 20)
        self.text_widgt.setStyleSheet("color: rgb(255, 255, 255);background-color: rgb(35, 38, 41);")
        splitter2.addWidget(self.view)
        splitter2.addWidget(self.text_widgt)
        self.slider = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.slider.setFixedSize(700, 50)  # 2
        self.slider.setRange(0, 1800)  # 3
        self.slider.valueChanged.connect(self.on_change_func)  # 5
        self.label_slider = QtWidgets.QLabel('0', self)
        self.label_slider.setFont(QtGui.QFont('Arial Black', 20))
        self.h_layout = QtWidgets.QHBoxLayout()
        self.h_layout.addWidget(self.slider)
        self.h_layout.addWidget(self.label_slider)
        self.h_layout.addWidget(self.determine_button2)
        labela = QtWidgets.QLabel("a:")
        labelb = QtWidgets.QLabel("b:")
        self.atextedit = QtWidgets.QLineEdit()
        self.atextedit.setPlaceholderText('[u, v]')
        self.btextedit = QtWidgets.QLineEdit()
        self.btextedit.setPlaceholderText('[u, v]')

        H2Layout = QtWidgets.QHBoxLayout()
        H2Layout.addWidget(labela)
        H2Layout.addWidget(self.atextedit)
        H2Layout.addWidget(labelb)
        H2Layout.addWidget(self.btextedit)
        self.atextedit.textChanged.connect(self.a_change)
        self.btextedit.textChanged.connect(self.b_change)
        layout.addWidget(splitter2)
        Hlayout.addWidget(error_label)
        Hlayout.addWidget(self.error_limit_widget)
        Hlayout.addWidget(max_length_label)
        Hlayout.addWidget(self.max_length_widget)
        Hlayout.addWidget(self.determine_button)
        layout.addLayout(Hlayout)
        layout.addLayout(self.h_layout)
        layout.addLayout(H2Layout)
        self.determine_button.clicked.connect(self.determine)
        self.determine_button2.clicked.connect(self.determine2)
        self.setLayout(layout)
        self.plot()

    def a_change(self):
        try:
            text = self.atextedit.text()
            if text[0] == '[' and text[-1] == ']':
                c = ''
                for letter in text:
                    if sf.is_number(letter) or letter == ',' or letter == '-':
                        c += letter
                l = c.split(',')
                if len(l) == 2:
                    # for i in range(len(l)):
                    l[0] = int(l[0])
                    l[1] = int(l[1])
                    self.a = self.vectora1 * l[0] + self.vectorb1 * l[1]
                    if self.b is not None:
                        self.x_axisb = [0, self.a[0], self.a[0] + self.b[0], self.b[0], 0]
                        self.y_axisb = [0, self.a[1], self.a[1] + self.b[1], self.b[1], 0]
                        array_bx = np.array(self.x_axisb)
                        array_by = np.array(self.y_axisb)
                        self.p1.clear()
                        self.p1.addItem(self.s1)
                        self.p1.plot(array_bx, array_by)
                        self.s1.setPoints(self.new_spot)
                else:
                    self.a_clent = None

            else:
                self.a_clent = None
        except Exception as e:
            print(e)

    def b_change(self):
        try:
            text = self.btextedit.text()
            if text[0] == '[' and text[-1] == ']':
                c = ''
                for letter in text:
                    if sf.is_number(letter) or letter == ',' or letter == '-':
                        c += letter
                l = c.split(',')
                if len(l) == 2:
                    self.b = self.vectora1 * int(l[0]) + self.vectorb1 * int(l[1])
                    if self.a is not None:
                        self.x_axisb = [0, self.a[0], self.a[0] + self.b[0], self.b[0], 0]
                        self.y_axisb = [0, self.a[1], self.a[1] + self.b[1], self.b[1], 0]
                        array_bx = np.array(self.x_axisb)
                        array_by = np.array(self.y_axisb)
                        self.p1.clear()
                        self.p1.addItem(self.s1)
                        self.p1.plot(array_bx, array_by)
                        self.s1.setPoints(self.new_spot)
                else:
                    self.b = None
            else:
                self.b = None
        except Exception as e:
            print(e)

    def on_change_func(self):
        self.theta = float(self.slider.value() / 10) / 180 * np.pi
        self.label_slider.setText(str(self.slider.value() / 10) + "°")
        rotate_matrix = np.array([[np.cos(self.theta), np.sin(self.theta)], [-np.sin(self.theta), np.cos(self.theta)]])
        self.new_x_up = []
        self.new_y_up = []
        for i in range(len(self.xu)):
            spot = np.array([self.xu[i], self.yu[i]])
            new_pos = list(spot @ rotate_matrix)
            self.new_x_up.append(new_pos[0])
            self.new_y_up.append(new_pos[1])

        self.find_in_tolerance()
        if self.near_index:
            x_near = []
            y_near = []
            for index in self.near_index:
                x_near.append(self.xd[index])
                y_near.append(self.yd[index])
            self.new_x = np.append(self.xd, self.new_x_up)  # 所有x_down点 + x_up点
            self.new_y = np.append(self.yd, self.new_y_up)
            self.new_new_x = np.append(self.new_x, x_near)
            self.new_new_y = np.append(self.new_y, y_near)
            pos = np.array([self.new_new_x, self.new_new_y])
            self.new_spot = [{'pos': pos[:, i], 'data': 1, 'brush': pg.mkBrush(255, 250, 250, 120)} for i in
                             range(len(self.xd))] + \
                            [{'pos': pos[:, i], 'data': 1, 'brush': pg.mkBrush(255, 20, 147, 120)} for i in
                             range(len(self.xd), len(self.new_x))] + \
                            [{'pos': pos[:, i], 'data': 1, 'brush': pg.mkBrush(0, 255, 0, 120)} for i in
                             range(len(self.new_x), len(self.new_new_x))]
        else:
            self.new_x = np.append(self.xd, self.new_x_up)  # 所有x_down点 + x_up点
            self.new_y = np.append(self.yd, self.new_y_up)
            pos = np.array([self.new_x, self.new_y])
            self.new_spot = [{'pos': pos[:, i], 'data': 1, 'brush': pg.mkBrush(255, 250, 250, 120)} for i in
                             range(len(self.xd))] + \
                            [{'pos': pos[:, i], 'data': 1, 'brush': pg.mkBrush(255, 20, 147, 120)} for i in
                             range(len(self.xd), len(self.new_x))]
        self.plot_cell()

    def plot(self):               #  orijin
        self.p1.clear()
        self.s1 = pg.ScatterPlotItem(size=10, pen=pg.mkPen(None))
        self.lattice_enlarge()
        self.xd = [self.pos_lst_orijin[i][0] for i in range(len(self.pos_lst_orijin))]
        self.yd = [self.pos_lst_orijin[i][1] for i in range(len(self.pos_lst_orijin))]
        self.xu = [self.pos_lst_orijin_up[i][0] for i in range(len(self.pos_lst_orijin_up))]
        self.yu = [self.pos_lst_orijin_up[i][1] for i in range(len(self.pos_lst_orijin_up))]
        self.lattice_plot()
        self.p1.addItem(self.s1)
        self.show()

    def lattice_enlarge(self):
        # self.crys1
        cell1 = list(self.crys1.get_cell())
        vectora1 = np.array([cell1[0][0], cell1[0][1]])
        self.vectora1 = vectora1
        vectorb1 = np.array([cell1[1][0], cell1[1][1]])
        self.vectorb1 = vectorb1

        a_len1 = np.linalg.norm(vectora1)
        b_len1 = np.linalg.norm(vectorb1)
        if a_len1 < b_len1:
            n1 = int(self.len_max / a_len1)
        else:
            n1 = int(self.len_max / b_len1)
        # self.crys2
        cell2 = list(self.crys2.get_cell())
        vectora2 = np.array([cell2[0][0], cell2[0][1]])
        vectorb2 = np.array([cell2[1][0], cell2[1][1]])
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
            for j in range(0, 2 * n1):
                pos = vectora1 * i + vectorb1 * j
                if np.linalg.norm(pos) <= self.len_max:
                    self.pos_lst_orijin.append(pos)

    def lattice_plot(self):
        try:
            x = np.append(self.xd, self.xu)
            y = np.append(self.yd, self.yu)
            pos = np.array([x, y])
            spots = [{'pos': pos[:, i], 'data': 1, 'brush': pg.mkBrush(255, 250, 250, 120)} for i in
                     range(len(self.xd))] + \
                    [{'pos': pos[:, i], 'data': 1, 'brush': pg.mkBrush(255, 20, 147, 120)} for i in
                     range(len(self.xd), len(x))]
            self.s1.addPoints(spots)
        except Exception as e:
            print(e)

    def determine(self):
        try:
            self.len_max = float(self.max_length_widget.text())
            self.error_limit = float(self.error_limit_widget.text()) / 100
            self.error_limit_widget.setPlaceholderText(str(self.error_limit * 100))
            self.max_length_widget.setPlaceholderText(str(self.len_max))
            self.label_slider.setText("0")
            self.slider.setValue(0)
            self.theta = 0
            self.plot()
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, 'error', "Input error!")
            print(e)

    def determine2(self):
        try:
            self.after_scatter = after_scatter_plot()
            self.after_scatter.signal_distances.connect(self.distances)
        except Exception as e:
            print(e)

    def distances(self, layer_distance, vacuum_distance):
        try:
            self.emit_information_signal.emit(self.theta, list(self.a), list(self.b), layer_distance, vacuum_distance)
        except Exception as e:
            print(e)

    def plot_cell(self):
        try:
            self.p1.clear()
            self.p1.addItem(self.s1)
            self.s1.setPoints(self.new_spot)
            self.a_change()
        except Exception as e:
            print(e)

    def find_in_tolerance(self):
        try:
            self.near_index = []
            for i in range(len(self.new_x_up)):
                for j in range(len(self.xd)):
                    if (self.yd[j] > self.error_limit * self.len_max and self.new_y_up[i] <= 0) or \
                            self.new_y_up[i] < - self.error_limit * self.len_max:
                        break
                    else:
                        pos1 = np.array([self.xd[j], self.yd[j]])
                        pos2 = np.array([self.new_x_up[i], self.new_y_up[i]])
                        dis = np.linalg.norm(pos1 - pos2)
                        # if dis < 0.1:
                        if dis < self.error_limit * np.linalg.norm(pos1):
                            self.near_index.append(j)
            orientation_txt = ""
            for index in self.near_index:
                mubiao_vec = np.array([self.xd[index], self.yd[index]])
                B = np.array([self.vectora1, self.vectorb1]).T
                r = list(np.dot(np.linalg.inv(B), mubiao_vec))
                r = '[{}, {}]'.format(str(round(r[0])), str(round(r[1])))
                orientation_txt = orientation_txt + r + '\n'
            self.text = "Rotate {}°".format(self.theta / np.pi * 180) + "\n" +\
                        "-" * 12 + "Near points orientation:" + '\n' + orientation_txt
            self.text_widgt.setText(self.text)
            self.plot_cell()
        except Exception as e:
            print(e)



class scatter_plot(QtWidgets.QWidget):
    len_max = 15
    error_limit = 0.02
    emit_information_signal = QtCore.pyqtSignal(float, list, list, float, float)
    def __init__(self, crys1, crys2, parent=None):
        super(scatter_plot, self).__init__(parent)
        pg.setConfigOption('background', (35, 38, 41))
        self.resize(1000, 2000)
        self.crys1 = crys1      # down
        self.crys2 = crys2      # up
        layout = QtWidgets.QVBoxLayout(self)  # 垂直布局
        Hlayout = QtWidgets.QHBoxLayout(self)
        splitter2 = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        self.determine_button = QtWidgets.QPushButton('OK')
        self.determine_button2 = QtWidgets.QPushButton('OK')
        self.error_limit_widget = QtWidgets.QLineEdit()
        self.error_limit_widget.setPlaceholderText(str(self.error_limit * 100))
        error_label = QtWidgets.QLabel("Error tolerance(%)：：")
        error_label.setBuddy(self.error_limit_widget)
        self.max_length_widget = QtWidgets.QLineEdit()
        self.max_length_widget.setPlaceholderText(str(self.len_max))

        max_length_label = QtWidgets.QLabel("Max length(Å): ")
        max_length_label.setBuddy(self.max_length_widget)
        self.view = pg.GraphicsLayoutWidget()  ## GraphicsView with GraphicsLayout inserted by default
        self.p1 = self.view.addPlot()
        self.vb = self.p1.vb
        self.s1 = pg.ScatterPlotItem(size=10, pen=pg.mkPen(None))
        pg.setConfigOption('background', (35, 38, 41))
        self.text_widgt = QtWidgets.QTextEdit()
        self.text_widgt.setReadOnly(True)
        self.text_widgt.resize(24, 20)
        self.text_widgt.setStyleSheet("color: rgb(255, 255, 255);background-color: rgb(35, 38, 41);")
        splitter2.addWidget(self.view)
        splitter2.addWidget(self.text_widgt)
        self.slider = QtWidgets.QSlider(QtCore.Qt.Horizontal,self)
        self.slider.setFixedSize(700, 50)  # 2
        self.slider.setRange(0, 1800)  # 3
        self.slider.valueChanged.connect(self.on_change_func)  # 5
        self.label_slider = QtWidgets.QLabel('0', self)
        self.label_slider.setFont(QtGui.QFont('Arial Black', 20))
        self.h_layout = QtWidgets.QHBoxLayout()
        self.h_layout.addWidget(self.slider)
        self.h_layout.addWidget(self.label_slider)
        self.h_layout.addWidget(self.determine_button2)
        layout.addWidget(splitter2)
        Hlayout.addWidget(error_label)
        Hlayout.addWidget(self.error_limit_widget)
        Hlayout.addWidget(max_length_label)
        Hlayout.addWidget(self.max_length_widget)
        Hlayout.addWidget(self.determine_button)
        layout.addLayout(Hlayout)
        layout.addLayout(self.h_layout)
        self.determine_button.clicked.connect(self.determine)
        self.determine_button2.clicked.connect(self.determine2)
        self.setLayout(layout)
        self.plot()

    def on_change_func(self):
        self.theta = float(self.slider.value()/10) / 180 * np.pi
        self.label_slider.setText(str(self.slider.value()/10) + "°")
        rotate_matrix = np.array([[np.cos(self.theta), np.sin(self.theta)], [-np.sin(self.theta), np.cos(self.theta)]])
        self.new_x_up = []
        self.new_y_up = []
        for i in range(len(self.xu)):
            spot = np.array([self.xu[i], self.yu[i]])
            new_pos = list(spot @ rotate_matrix)
            self.new_x_up.append(new_pos[0])
            self.new_y_up.append(new_pos[1])
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
        self.xd = [self.pos_lst_orijin[i][0] for i in range(len(self.pos_lst_orijin))]
        self.yd = [self.pos_lst_orijin[i][1] for i in range(len(self.pos_lst_orijin))]
        self.xu = [self.pos_lst_orijin_up[i][0] for i in range(len(self.pos_lst_orijin_up))]
        self.yu = [self.pos_lst_orijin_up[i][1] for i in range(len(self.pos_lst_orijin_up))]
        self.lattice_plot()
        self.p1.addItem(self.s1)
        self.show()

    def lattice_enlarge(self):
        # self.crys1
        cell1 = list(self.crys1.get_cell())
        vectora1 = np.array([cell1[0][0], cell1[0][1]])
        vectorb1 = np.array([cell1[1][0], cell1[1][1]])
        a_len1 = np.linalg.norm(vectora1)
        b_len1 = np.linalg.norm(vectorb1)
        if a_len1 < b_len1:
            n1 = int(self.len_max / a_len1)
        else:
            n1 = int(self.len_max / b_len1)
        # self.crys2
        cell2 = list(self.crys2.get_cell())
        vectora2 = np.array([cell2[0][0], cell2[0][1]])
        vectorb2 = np.array([cell2[1][0], cell2[1][1]])
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
            for j in range(0, 2 * n1):
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

    def determine(self):
        try:
            self.len_max = float(self.max_length_widget.text())
            self.error_limit = float(self.error_limit_widget.text()) / 100
            self.error_limit_widget.setPlaceholderText(str(self.error_limit * 100))
            self.max_length_widget.setPlaceholderText(str(self.len_max))
            self.label_slider.setText("0")
            self.slider.setValue(0)
            self.theta = 0
            self.plot()
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, 'error', "Input error!")
            print(e)

    def determine2(self):
        try:
            self.after_scatter = after_scatter_plot()
            self.after_scatter.signal_distances.connect(self.distances)
        except Exception as e:
            print(e)

    def distances(self, layer_distance, vacuum_distance):
        try:
            self.emit_information_signal.emit(self.theta, list(self.a), list(self.b), layer_distance, vacuum_distance)
        except Exception as e:
            print(e)


    def find_min(self):
        try:
            peidui_lst = []
            for i in range(len(self.new_x_up)):
                for j in range(len(self.xd)):
                    if (self.yd[j] > self.error_limit * self.len_max and self.new_y_up[i] <= 0) or \
                            self.new_y_up[i] < - self.error_limit * self.len_max:
                        break
                    else:
                        pos1 = np.array([self.xd[j], self.yd[j]])
                        pos2 = np.array([self.new_x_up[i], self.new_y_up[i]])
                        dis = np.linalg.norm(pos1 - pos2)
                        if dis < self.error_limit *  np.linalg.norm(pos1):
                            peidui_lst.append(pos1)
            if len(peidui_lst) < 2:
                self.text = ""
                self.p1.clear()
                self.p1.addItem(self.s1)
                self.s1.setPoints(self.new_spot)
                self.text_widgt.setText(self.text)
            else:
                self.min = 1000000
                self.a = ...
                self.b = ...
                for i in range(len(peidui_lst)):
                    for j in range(i + 1, len(peidui_lst)):
                        n_vector = np.cross(peidui_lst[i], peidui_lst[j])
                        square = abs(n_vector)
                        if 1e-6 < square < self.min:
                            if n_vector > 0:     # 保证右手系
                                self.a = peidui_lst[i]
                                self.b = peidui_lst[j]
                            else:
                                self.b = peidui_lst[i]
                                self.a = peidui_lst[j]
                            self.min = square
                if self.min == 1000000:
                    self.text = ""
                    self.p1.clear()
                    self.p1.addItem(self.s1)
                    self.s1.setPoints(self.new_spot)
                    self.text_widgt.setText(self.text)
                else:
                    self.text = "Rotate {}°, a = {}, b = {}".format(self.theta / np.pi * 180, self.a, self.b) + "\n" +\
                                "Minimum square = {}".format(self.min)
                    self.text_widgt.setText(self.text)
                    self.plot_cell()
        except Exception as e:
            print(e)

    def plot_cell(self):
        try:
            self.x_axisb = [0, self.a[0], self.a[0] + self.b[0], self.b[0], 0]
            self.y_axisb = [0, self.a[1], self.a[1] + self.b[1], self.b[1], 0]
            array_bx = np.array(self.x_axisb)
            array_by = np.array(self.y_axisb)
            self.p1.clear()
            self.p1.addItem(self.s1)
            self.p1.plot(array_bx, array_by)
            self.s1.setPoints(self.new_spot)
        except Exception as e:
            print(e)

class after_scatter_plot(QtWidgets.QWidget):
    signal_distances = QtCore.pyqtSignal(float, float)
    def __init__(self, parent=None):
        super(after_scatter_plot, self).__init__(parent)
        hlayout1 = QtWidgets.QHBoxLayout(self)
        self.layer_distance = QtWidgets.QLineEdit()
        label_layer = QtWidgets.QLabel("layer distance(Å)：")
        self.vacuum_distance = QtWidgets.QLineEdit()
        label_vacuum = QtWidgets.QLabel("Vacuum layer distance(Å)：")
        self.determine_button = QtWidgets.QPushButton('OK')
        hlayout1.addWidget(label_layer)
        hlayout1.addWidget(self.layer_distance)
        hlayout1.addWidget(label_vacuum)
        hlayout1.addWidget(self.vacuum_distance)
        hlayout1.addWidget(self.determine_button)
        self.setLayout(hlayout1)
        self.determine_button.clicked.connect(self.determine)
        self.show()

    def determine(self):
        try:
            self.layer_distance_value = float(self.layer_distance.text())
            self.vaccum_distance_value = float(self.vacuum_distance.text())
            self.signal_distances.emit(self.layer_distance_value, self.vaccum_distance_value)
            self.close()
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')


class zhexian_plot(QtWidgets.QWidget):
    """ The layout and function of zhexian_plot window. """
    signal_judge_point1 = QtCore.pyqtSignal(float, float)
    signal_judge_point2 = QtCore.pyqtSignal(float, float)
    signal_point1_index = QtCore.pyqtSignal(int)
    signal_point2_index = QtCore.pyqtSignal(int)
    signal_determine = QtCore.pyqtSignal()
    signal_layer_length = QtCore.pyqtSignal(float)
    num1 = 0
    num2 = 0
    def __init__(self, parent=None):
        pg.setConfigOption('background', (35, 38, 41))
        super(zhexian_plot, self).__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)    # 垂直布局
        Hlayout = QtWidgets.QHBoxLayout(self)
        splitter2 = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        self.determine_button = QtWidgets.QPushButton('OK')
        self.layer_distance = QtWidgets.QLineEdit()
        label = QtWidgets.QLabel("layer distance(Å)：")
        label.setBuddy(self.layer_distance)
        self.lastClicked = []
        self.lastClicked2 = []
        self.signal_judge_point1.connect(self.judge_point1)
        self.signal_judge_point2.connect(self.judge_point2)
        self.view = pg.GraphicsLayoutWidget()  ## GraphicsView with GraphicsLayout inserted by default
        self.plot_window1 = self.view.addPlot(title="lcma - StrainaSum plot")
        self.plot_window2 = self.view.addPlot(title="lcmb - StrainbSum plot")
        self.plot_window1.setLabel('left', "a-length", units='Å')
        self.plot_window1.setLabel('bottom', "StrainaSum", units='')
        self.plot_window2.setLabel('left', "b-length", units='Å')
        self.plot_window2.setLabel('bottom', "StrainbSum", units='')
        # self.buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok|QtWidgets.QDialogButtonBox.Cancel, QtCore.Qt.Horizontal, self))
        self.text_widgt = QtWidgets.QTextEdit()
        self.text_widgt.setReadOnly(True)
        self.text_widgt.resize(24, 20)
        self.text_widgt.setStyleSheet("color: rgb(255, 255, 255);background-color: rgb(35, 38, 41);")
        splitter2.addWidget(self.view)
        splitter2.addWidget(self.text_widgt)
        layout.addWidget(splitter2)
        Hlayout.addWidget(label)
        Hlayout.addWidget(self.layer_distance)
        Hlayout.addWidget(self.determine_button)
        layout.addLayout(Hlayout)
        self.determine_button.clicked.connect(self.determine)
        self.setLayout(layout)

    def plot(self, plotlist):        # 画折线图
        self.list = plotlist
        self.x_axisa = [self.list[0][i][0] for i in range(len(self.list[0]))]
        self.y_axisa = [self.list[0][i][1] for i in range(len(self.list[0]))]
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
        self.x_axisb = [self.list[1][i][0] for i in range(len(self.list[1]))]
        self.y_axisb = [self.list[1][i][1] for i in range(len(self.list[1]))]
        array_bx = np.array(self.x_axisb)
        array_by = np.array(self.y_axisb)
        self.plot_window2.plot(array_bx, array_by)
        self.scene2.setData(x=np.array(self.x_axisb), y=np.array(self.y_axisb))
        self.plot_window2.addItem(self.scene2)
        self.scene2.sigClicked.connect(self.clicked2)
        self.show()         # 绘图

    def set_text_init(self, Atomslst):
        crys1 = Atomslst[0]
        crys2 = Atomslst[1]
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
        for i in range(len(self.x_axisa)):
            length = np.sqrt((self.x_axisa[i] * 10000 - float1 * 10000) ** 2 + (self.y_axisa[i] - float2) ** 2)
            if length < length_min:
                length_min = length
                point1_index = i
        self.signal_point1_index.emit(point1_index)
        self.num1 += 1

    def judge_point2(self, float1, float2):
        length_min = 1000         # 判断用于点击的点的index
        point2_index = ...
        for i in range(len(self.x_axisb)):
            length = np.sqrt((self.x_axisb[i] * 10000 - float1 * 10000) ** 2 + (self.y_axisb[i] - float2) ** 2)
            if length < length_min:
                length_min = length
                point2_index = i
        self.signal_point2_index.emit(point2_index)
        self.num2 += 1


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
                self.signal_judge_point1.emit(mousePoint1.x(), mousePoint1.y())
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
                self.signal_judge_point2.emit(mousePoint2.x(), mousePoint2.y())
            self.lastClicked2 = points
        except Exception as e:
            print(e)

    def determine(self):
        try:
            try:
                c = float(self.layer_distance.text())
                if self.num1 == 0 or self.num2 == 0:
                    QtWidgets.QMessageBox.warning(self, 'error', 'please choose a，b length')
                else:
                    self.signal_layer_length.emit(c)  # 把用户输入的层间距赋值给self.layerlength
                    self.signal_determine.emit()  # 选定,与self.layertransform连接
                    self.close()
            except:
                QtWidgets.QMessageBox.warning(self, 'error', 'please Input valid layer distance')
        except Exception as e:
            print(e)

class Classification_window(QtWidgets.QWidget):
    signal_emit_dirs = pyqtSignal(str, str)
    input_directory = ""
    output_directory = ""
    def __init__(self, parent=None):
        super(Classification_window, self).__init__(parent)
        Vlayout = QtWidgets.QVBoxLayout()
        Hlayout1 = QtWidgets.QHBoxLayout()
        Hlayout2 = QtWidgets.QHBoxLayout()
        Hlayout3 = QtWidgets.QHBoxLayout()
        label1 = QtWidgets.QLabel("CIF Directory:")
        self.text_widget_cif_dir = QtWidgets.QLineEdit()
        self.button1 = QtWidgets.QPushButton("Browse")
        self.button1.clicked.connect(self.CIF_dir)
        self.text_widget_cif_dir.setEnabled(False)
        self.text_widget_cif_dir.setStyleSheet("color: rgb(255, 255, 255);background-color: rgb(35, 38, 41);")

        Hlayout1.addWidget(label1)
        Hlayout1.addWidget(self.text_widget_cif_dir)
        Hlayout1.addWidget(self.button1)
        label2 = QtWidgets.QLabel("Export address:")
        self.text_widget_export = QtWidgets.QLineEdit()
        self.text_widget_export.setStyleSheet("color: rgb(255, 255, 255);background-color: rgb(35, 38, 41);")

        self.text_widget_export.setEnabled(False)
        self.button2 = QtWidgets.QPushButton("Browse")
        self.button2.clicked.connect(self.export_dir)
        Hlayout2.addWidget(label2)
        Hlayout2.addWidget(self.text_widget_export)
        Hlayout2.addWidget(self.button2)
        self.button = QtWidgets.QPushButton("Start")
        self.button.clicked.connect(self.determine)
        self.pbar = QtWidgets.QProgressBar()
        Hlayout3.addWidget(self.button)
        Hlayout3.addWidget(self.pbar)
        Vlayout.addLayout(Hlayout1)
        Vlayout.addLayout(Hlayout2)
        Vlayout.addLayout(Hlayout3)
        self.setLayout(Vlayout)
        self.show()

    def CIF_dir(self):
        self.input_directory = QtGui.QFileDialog.getExistingDirectory(self, 'Select Classify directory')
        self.text_widget_cif_dir.setText(self.input_directory)

    def export_dir(self):
        self.output_directory = QtGui.QFileDialog.getExistingDirectory(self, 'Select Export directory')
        self.text_widget_export.setText(self.output_directory)

    def determine(self):


        try:
            if self.input_directory != "" and self.output_directory  != "":
                self.signal_emit_dirs.emit(self.input_directory, self.output_directory)
                # self.close()
            elif self.input_directory == "":
                QtWidgets.QMessageBox.warning(self, 'error', 'Choose CIF directory!')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Choose export directory!')
        except Exception as e:
            print(e)


class Change_gama_window(QtWidgets.QWidget):
    """ To change crys' γ"""

    signal_gama_error = pyqtSignal(list)
    def __init__(self, parent=None):
        super(Change_gama_window, self).__init__(parent)
        layout = QtWidgets.QHBoxLayout()
        label1 = QtWidgets.QLabel("Change γ to (degree)：")
        self.text1 = QtWidgets.QLineEdit()
        label2 = QtWidgets.QLabel("Error limit(°): ")
        self.text2 = QtWidgets.QLineEdit()
        self.button = QtWidgets.QPushButton("Ok")
        layout.addWidget(label1)
        layout.addWidget(self.text1)
        layout.addWidget(label2)
        layout.addWidget(self.text2)
        layout.addWidget(self.button)
        self.setLayout(layout)
        self.button.clicked.connect(self.determine)
        self.show()

    def determine(self):
        try:
            gamadegree = float(self.text1.text())
            errorlimit = float(self.text2.text())
            if 0<gamadegree<180 and errorlimit> 0:
                self.signal_gama_error.emit([gamadegree, errorlimit])
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error', 'Input error!')



class add_vacuum_layer_window(QtWidgets.QWidget):
    """ The layout of the add vacuum layer window. """
    sinal_vacuum = pyqtSignal(float)
    def __init__(self, parent=None):
        super(add_vacuum_layer_window, self).__init__(parent)
        self.setWindowTitle("Add vacuum layer")
        layout = QtWidgets.QHBoxLayout()
        label = QtWidgets.QLabel("vacuum layer distance(Å)：")
        self.text = QtWidgets.QLineEdit()
        self.determinebutton = QtWidgets.QPushButton("OK")
        layout.addWidget(label)
        layout.addWidget(self.text)
        layout.addWidget(self.determinebutton)
        self.determinebutton.clicked.connect(self.determine)
        self.setLayout(layout)
        self.show()

    def determine(self):
        try:
            vacuum_distance = float(self.text.text())
            self.sinal_vacuum.emit(vacuum_distance)
            self.close()
        except Exception as e:
            print(e)


class cut_window(QtWidgets.QWidget):
    """ The layout of the cut window. """
    # signal_cut_exf = pyqtSignal()
    signal_custom_cut_exf = pyqtSignal(list, int)
    num = 0
    def __init__(self, cell, parent=None):
        super(cut_window, self).__init__(parent)
        self.cell = cell
        self.resize(550, 500)
        zonglayout = QtWidgets.QVBoxLayout()
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
        self.setLayout(zonglayout)
        self.setWindowTitle("Cut Style")
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
                self.tablewidget.setHorizontalHeaderLabels(['x(Å)', 'y(Å)', 'z(Å)'])
                self.tablewidget.setVerticalHeaderLabels(['a-axis', 'b-axis', 'c-axis'])
            elif self.choosebutton.text() == 'uvw':
                self.tablewidget.setHorizontalHeaderLabels(['u', 'v', 'w'])
                self.tablewidget.setVerticalHeaderLabels(['a', 'b', 'c'])
        except Exception as e:
            print(e)

    def determine(self):
        try:
            cell_par = []
            for i in range(3):
                l = []
                for j in range(3):
                    l.append(float(self.tablewidget.item(i, j).text()))
                cell_par.append(l)
            if self.choosebutton.text() == "xyz":
                self.signal_custom_cut_exf.emit(cell_par, 1)
            elif self.choosebutton.text() == 'uvw':
                self.signal_custom_cut_exf.emit(cell_par, 2)
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.warning(self, 'error', 'Input Error')


class New_file(QtWidgets.QWidget):
    """ The layout of the New file window. """

    signal_determine = QtCore.pyqtSignal(str)
    def __init__(self, parent=None):
        super(New_file, self).__init__(parent)
        self.setWindowTitle('New cifFile')
        self.resize(300, 500)
        self.setMaximumSize(300, 1000)
        layout = QtWidgets.QVBoxLayout(self)    # 垂直布局
        self.tablewidget = QtWidgets.QTableWidget()
        self.tablewidget.setRowCount(100)
        self.tablewidget.setColumnCount(4)
        self.tablewidget.setHorizontalHeaderLabels(["Atom", 'x(Å)', 'y(Å)', 'z(Å)'])
        Hlayout1 = QtWidgets.QHBoxLayout(self)
        labela = QtWidgets.QLabel("a-axis length(Å): ")
        self.a = QtWidgets.QLineEdit()
        labelb = QtWidgets.QLabel("b-axis length(Å): ")
        self.b = QtWidgets.QLineEdit()
        labelc = QtWidgets.QLabel("c-axis length(Å): ")
        self.c = QtWidgets.QLineEdit()
        Hlayout1.addWidget(labela)
        Hlayout1.addWidget(self.a)
        Hlayout1.addWidget(labelb)
        Hlayout1.addWidget(self.b)
        Hlayout1.addWidget(labelc)
        Hlayout1.addWidget(self.c)
        layout.addLayout(Hlayout1)
        Hlayout2 =QtWidgets.QHBoxLayout(self)
        labelalpha = QtWidgets.QLabel("alpha(degree): ")
        self.alpha = QtWidgets.QLineEdit()
        labelbetta = QtWidgets.QLabel("beta(degree): ")
        self.beta = QtWidgets.QLineEdit()
        labelgamma = QtWidgets.QLabel("gama(degree): ")
        self.gama = QtWidgets.QLineEdit()
        Hlayout2.addWidget(labelalpha)
        Hlayout2.addWidget(self.alpha)
        Hlayout2.addWidget(labelbetta)
        Hlayout2.addWidget(self.beta)
        Hlayout2.addWidget(labelgamma)
        Hlayout2.addWidget(self.gama)
        layout.addLayout(Hlayout2)
        self.determine_button = QtWidgets.QPushButton('OK')
        self.determine_button.clicked.connect(self.determine)
        Hlayout3 = QtWidgets.QHBoxLayout(self)
        self.name_label = QtWidgets.QLabel("File Name: ")
        self.name = QtWidgets.QLineEdit()
        Hlayout3.addWidget(self.name_label)
        Hlayout3.addWidget(self.name)
        Hlayout3.addWidget(self.determine_button)
        layout.addWidget(self.tablewidget)
        layout.addLayout(Hlayout3)
        self.setLayout(layout)

    def show_win(self):
        self.show()

    def determine(self):
        try:
            # cell parameter
            test = 0
            alpha = float(self.alpha.text())
            betta = float(self.beta.text())
            gama = float(self.gama.text())
            a = float(self.a.text())
            b = float(self.b.text())
            c = float(self.c.text())
            self.make_par_lst_matrix(alpha, betta, gama, a, b, c)  # 生成晶胞向量
            # table
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
            self.make_par_lst_species_par_lst_coords()

            self.cifname = self.name.text()
            if self.cifname == '':
                QtWidgets.QMessageBox.warning(self, 'error', 'Name Input error!')
            else:
                self.generate()
        except Exception as e:
            print(e)
            if test == 0:
                QtWidgets.QMessageBox.warning(self, 'error', 'Cell parameter Input error!')
            # elif test == 3:
            #     QtWidgets.QMessageBox.warning(self, 'error', 'Table Atom Input error!')
            else:
                QtWidgets.QMessageBox.warning(self, 'error', 'Table Pos Input error!')


    def coords_transform(self, A1, B1, C1):
        eps = 0.0001
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
        self.par_lst_matrix = [a_vec, b_vec, c_vec]

    def make_par_lst_species_par_lst_coords(self):
        self.par_lst_species = [self.atom_lst[i][0] for i in range(len(self.atom_lst))]
        self.par_lst_coords = [[self.atom_lst[i][1], self.atom_lst[i][2], self.atom_lst[i][3]] for i in range(len(self.atom_lst))]

    def generate(self):
        try:
            self.directory = QtGui.QFileDialog.getExistingDirectory(self, 'Select directory')
            structure = Structure(self.par_lst_matrix, self.par_lst_species, self.par_lst_coords)
            slab = CifWriter(structure, write_magmoms=True)
            slab.write_file(self.directory +'/' + self.cifname + '.cif')
            self.signal_determine.emit(self.directory +'/' + self.cifname + '.cif')
        except Exception as e:
            print(e)

class add_cell_plot(QtWidgets.QWidget):      # 扩胞
    """ The layout of the creat supercell window. """
    signala = pyqtSignal(int)
    signalb = pyqtSignal(int)
    signalc = pyqtSignal(int)
    def __init__(self, parent=None):
        super(add_cell_plot, self).__init__(parent)
        self.setWindowTitle('add cell')
        self.resize(300, 100)
        layout = QtWidgets.QVBoxLayout(self)    # 垂直布局
        Hlayout = QtWidgets.QHBoxLayout(self)
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
        self.setLayout(layout)
        self.show()

    def valuea(self):
        self.signala.emit(self.numberSpinBoxa.value())

    def valueb(self):
        self.signalb.emit(self.numberSpinBoxb.value())

    def valuec(self):
        self.signalc.emit(self.numberSpinBoxc.value())

class Replace_Atom(QtWidgets.QWidget):
    signal_atomic_number = pyqtSignal(int)
    def __init__(self):
        super(Replace_Atom, self).__init__()
        Hlayout0 = QtWidgets.QHBoxLayout(self)
        label_atomicnumber = QtWidgets.QLabel("atomic number: ")
        self.lineEdit = QtWidgets.QLineEdit()
        self.button = QtWidgets.QPushButton("Ok")
        Hlayout0.addWidget(label_atomicnumber)
        Hlayout0.addWidget(self.lineEdit)
        Hlayout0.addWidget(self.button)
        self.button.clicked.connect(self.determine)
        self.setLayout(Hlayout0)
        self.show()

    def determine(self):
        try:
            atomic_number = int(self.lineEdit.text())
            self.signal_atomic_number.emit(atomic_number)
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, 'error', 'Atomic number input error!')
            print(e)









