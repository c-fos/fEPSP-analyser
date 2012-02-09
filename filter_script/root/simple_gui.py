import sys
from PyQt4 import QtCore, QtGui
from simple import Ui_MainWindow
from graph import graphReconstruction
from fEPSPanalyser import fepspAnalyser
import pickle

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s
        
        
class MyForm(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.fsmodel = QtGui.QFileSystemModel()
        self.stmodel = QtGui.QStandardItemModel()
        self.root_index=self.fsmodel.setRootPath(self.ui.pathLine.text())
        QtCore.QObject.connect(self.ui.sourceList, QtCore.SIGNAL(_fromUtf8("doubleClicked(QModelIndex)")),\
         self.sourceListClicked)
        QtCore.QObject.connect(self.ui.processedList, QtCore.SIGNAL(_fromUtf8("doubleClicked(QModelIndex)")),\
         self.print_image_2)
        QtCore.QObject.connect(self.ui.pathLine, QtCore.SIGNAL(_fromUtf8("editingFinished()")),\
                 self.show_path)
        QtCore.QObject.connect(self.ui.createListButton, QtCore.SIGNAL(_fromUtf8("clicked()")),\
                 self.form_list)
        QtCore.QObject.connect(self.ui.graphButton, QtCore.SIGNAL(_fromUtf8("clicked()")),\
                 self.makeGraph)
        QtCore.QObject.connect(self.ui.fEPSP_button, QtCore.SIGNAL(_fromUtf8("clicked()")),\
                 self.fEPSP_start)
        QtCore.QObject.connect(self.ui.upButton, QtCore.SIGNAL(_fromUtf8("clicked()")),\
                 self.pathUp)
        QtCore.QObject.connect(self.ui.rmItemButton, QtCore.SIGNAL(_fromUtf8("clicked()")),\
                 self.rmItem)
        QtCore.QObject.connect(self.ui.clearDbButton, QtCore.SIGNAL(_fromUtf8("clicked()")),\
                 self.rmDbConfig)
        QtCore.QObject.connect(self.ui.saveDbButton, QtCore.SIGNAL(_fromUtf8("clicked()")),\
                 self.mkDbConfig)
        sourceList = self.ui.sourceList
        sourceList.setModel(self.fsmodel)
        sourceList.setRootIndex(self.root_index)
        sourceList.show()
        self.loadConfig()
    
    def rmDbConfig(self):
        print("this function has not implemented yet")
        
    def mkDbConfig(self):
        fd=open("dbConfig",'w')
        list=[str(self.ui.dbServerIpLine.text()),str(self.ui.dbNameLine.text()),str(self.ui.dbUserLine.text()),str(self.ui.dbPassLine.text())]
        pickle.dump(list,fd)
        fd.close()
    
    def loadConfig(self):
        try:
            fd=open("dbConfig",'r')
            dbAccessVars=pickle.load(fd)
            fd.close()
            self.ui.dbServerIpLine.setText(dbAccessVars[0])
            self.ui.dbNameLine.setText(dbAccessVars[1])
            self.ui.dbUserLine.setText(dbAccessVars[2])
            self.ui.dbPassLine.setText(dbAccessVars[3])
        except:
            pass
        
    def rmItem(self):
        self.ui.processedList.takeItem(self.ui.processedList.currentRow())

    def fEPSP_start(self):
        path=str(self.ui.pathLine.text())
        frequency=str(self.ui.frequency_line.text())
        substance=str(self.ui.substance_line.text())
        if self.ui.debugBox.isChecked():
            debug="1"
        else:
            debug="0"
        if self.ui.clusterizationBox.isChecked():
            cluster="0"
        else:
            cluster="1"
        if self.ui.database_checkBox.isChecked():
            write="1"
        else:
            write="0"
        if debug=="1":
            print((path,frequency,"data","1",substance,debug,write))
        analyserObject=fepspAnalyser([0,path,frequency,"data","1",substance,debug,write,cluster])
        
    def show_path(self):
        path = self.ui.pathLine.text()#QtGui.QDesktopServices.storageLocati$
        print(path)
        self.root_index = self.fsmodel.setRootPath(path)
        self.ui.sourceList.setRootIndex(self.root_index)
        self.ui.sourceList.update()
        
    def pathUp(self):
        path = self.ui.pathLine.text()#QtGui.QDesktopServices.storageLocati$
        print(path)
        index = self.fsmodel.setRootPath(path)
        newPathIndex=self.fsmodel.parent(index)
        newPath=self.fsmodel.filePath(newPathIndex)
        self.ui.sourceList.setRootIndex(newPathIndex)
        self.ui.pathLine.setText(newPath)
        
        
    def sourceListClicked(self):
        curIndex=self.ui.sourceList.currentIndex()
        if self.fsmodel.isDir(curIndex):
            self.ui.pathLine.setText(self.fsmodel.filePath(curIndex))
            self.ui.sourceList.setRootIndex(curIndex)
        else:
            self.print_image("source")
            
    def print_image_2(self):
            self.print_image("processed")  
        
    def print_image(self,ident):
        if ident=="source":
            imageIndex = self.ui.sourceList.currentIndex()
            imagePath=self.fsmodel.fileName(imageIndex)
        else:
            imageIndex = self.ui.processedList.currentItem()
            imagePath=imageIndex.text()
            print((imageIndex))
        imagePath.prepend(self.ui.pathLine.text())
        print((imagePath,"1"))
        pixmap = QtGui.QPixmap(imagePath)
        self.ui.imageLabel.setPixmap(pixmap)
        self.ui.imageLabel.resize(400, 300)
        self.show()   
        
    def form_list(self):
        fileIndexes = self.ui.sourceList.selectedIndexes()
        for i in range(len(fileIndexes)):
            self.ui.processedList.addItem(self.fsmodel.fileName(fileIndexes[i]))
            
    def makeGraph(self):
        start = int(self.ui.startLine.text())
        stop = int(self.ui.stopLine.text())
        rstride = int(self.ui.rstrideLine.text())
        cstride = int(self.ui.cstrideLine.text())
        if self.ui.debugBox2.isChecked():
            debug2="1"
        else:
            debug2="0"
            
        try:
            self.ui.processedList.selectAll()
            fileList = list(self.ui.processedList.selectedItems())
            for x in range(len(fileList)):
                fileList[x]=str(fileList[x].text().prepend(self.ui.pathLine.text()))
            fileList.sort()
            if debug2=="1":
                print(fileList)
            reconst=graphReconstruction(fileList,start=start,stop=stop,rstride=rstride,cstride=cstride,debug=debug2)
        except:
            pass   

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = MyForm()
    myapp.show()
    sys.exit(app.exec_())