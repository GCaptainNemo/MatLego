import sqlite3
import operator
import copy

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
            self.info_data_base.sort(key=operator.itemgetter(1))     # 以ID排序，不是TrueID
        else:
            self.info_data_base.sort(key=operator.itemgetter(1))     # 以ID排序，不是TrueID
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
        self.info_data_base = copy.deepcopy(zong_lst)

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
            for i in range(len(self.info_data_base)):
                if self.info_data_base[i][0] == content:
                    formula = "'" + self.info_data_base[i][2] + "'"
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


