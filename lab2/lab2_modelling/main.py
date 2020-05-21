import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMessageBox

from Modeller import Modeller
from Ui_mainwindow import Ui_MainWindow

class MyApp(QtWidgets.QMainWindow):
	def __init__(self):
		super(MyApp, self).__init__()
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)

		self.ui.run_button.clicked.connect(self.run)
		self.ui.set_default_btn.clicked.connect(self.set_defaults)

		self.defaults = {"R"      : 0.35,
						 "Le"     : 12,
						 "Lk"     : 187e-6,
						 "Ck"     : 268e-6,
						 "Rk"     : 0.25,
						 "Uc0"    : 1400,
						 "I0"     : 0.5,
						 "Tw"     : 2000,
						 "Tstart" : 0,
						 "Tend"   : 0.0006,
						 "Tstep"  : 1e-6}
		self.set_defaults()
		self.data = {}

	def set_defaults(self):
		self.ui.lineEdit_r.setText(str(self.defaults.get("R")))
		self.ui.lineEdit_le.setText(str(self.defaults.get("Le")))
		self.ui.lineEdit_lk.setText(str(self.defaults.get("Lk")))
		self.ui.lineEdit_ck.setText(str(self.defaults.get("Ck")))
		self.ui.lineEdit_rk.setText(str(self.defaults.get("Rk")))
		self.ui.lineEdit_uc0.setText(str(self.defaults.get("Uc0")))
		self.ui.lineEdit_i0.setText(str(self.defaults.get("I0")))
		self.ui.lineEdit_tw.setText(str(self.defaults.get("Tw")))
		self.ui.lineEdit_tstart.setText(str(self.defaults.get("Tstart")))
		self.ui.lineEdit_tend.setText(str(self.defaults.get("Tend")))
		self.ui.lineEdit_tstep.setText(str(self.defaults.get("Tstep")))

	def run(self):
		if self.check():
			mdlr = Modeller(self.data)
			isequal = self.ui.checkBox.isChecked()
			mdlr.compute(isequal)
		else:
			self.msg_box("Error!", "Error! Incorrect input!", QMessageBox.Critical)
	
	def msg_box(self, title, message, type):
		msg = QMessageBox(self)
		msg.setIcon(type)
		msg.setWindowTitle(title)
		msg.setText(message)
		msg.addButton('Ok', QMessageBox.AcceptRole)
		msg.exec()

	def check(self):
		try:
			self.data["R"] = float(self.ui.lineEdit_r.text())
			self.data["Le"] = float(self.ui.lineEdit_le.text())
			self.data["Lk"] = float(self.ui.lineEdit_lk.text())
			self.data["Ck"] = float(self.ui.lineEdit_ck.text())
			self.data["Rk"] = float(self.ui.lineEdit_rk.text())
			self.data["Uc0"] = float(self.ui.lineEdit_uc0.text())
			self.data["I0"] = float(self.ui.lineEdit_i0.text())
			self.data["Tw"] = float(self.ui.lineEdit_tw.text())
			self.data["Tstart"] = float(self.ui.lineEdit_tstart.text())
			self.data["Tend"] = float(self.ui.lineEdit_tend.text())
			self.data["Tstep"] = float(self.ui.lineEdit_tstep.text())
		except ValueError:
			return False
		return True


def main():
	app = QtWidgets.QApplication(sys.argv)
	window = MyApp()
	window.show()
	app.exec_()


if __name__ == '__main__':
	main()
