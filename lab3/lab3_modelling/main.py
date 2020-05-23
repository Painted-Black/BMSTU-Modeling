import sys

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMessageBox

from Ui_MainWindow import Ui_MainWindow

from Modeller import Modeller

class MyApp(QtWidgets.QMainWindow):
	def __init__(self):
		super(MyApp, self).__init__()
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)

		self.ui.set_def_button.clicked.connect(self.set_defaults)
		self.ui.run_button.clicked.connect(self.run)

		self.defaults = {
			"k0" : 0.4,
			"kN" : 0.1,
			"a0" : 0.05,
			"aN" : 0.01,
			"l"  : 10,
			"T0" : 300,
			"R"  : 0.5,
			"F0" : 50,
			"h"  : 0.1
		}

		self.data = {
			"k0" : None,
			"kN" : None,
			"a0" : None,
			"aN" : None,
			"l"  : None,
			"T0" : None,
			"R"  : None,
			"F0" : None,
			"h"  : None
		}

		self.set_defaults()
	
	def set_defaults(self):
		self.ui.lineEdit_k0.setText(str(self.defaults.get("k0")))
		self.ui.lineEdit_kN.setText(str(self.defaults.get("kN")))
		self.ui.lineEdit_alpha0.setText(str(self.defaults.get("a0")))
		self.ui.lineEdit_alphaN.setText(str(self.defaults.get("aN")))
		self.ui.lineEdit_l.setText(str(self.defaults.get("l")))
		self.ui.lineEdit_T0.setText(str(self.defaults.get("T0")))
		self.ui.lineEdit_R.setText(str(self.defaults.get("R")))
		self.ui.lineEdit_F0.setText(str(self.defaults.get("F0")))
		self.ui.lineEdit_h.setText(str(self.defaults.get("h")))

	def get_data(self):
		try:
			self.data["k0"] = float(self.ui.lineEdit_k0.text())
			self.data["kN"] = float(self.ui.lineEdit_kN.text())
			self.data["a0"] = float(self.ui.lineEdit_alpha0.text())
			self.data["aN"] = float(self.ui.lineEdit_alphaN.text())
			self.data["l"] = float(self.ui.lineEdit_l.text())
			self.data["T0"] = float(self.ui.lineEdit_T0.text())
			self.data["R"] = float(self.ui.lineEdit_R.text())
			self.data["F0"] = float(self.ui.lineEdit_F0.text())
			self.data["h"] = float(self.ui.lineEdit_h.text())
		except ValueError:
			return False
		return True

	def run(self):
		if self.get_data():
			print("Computing...")
			mdlr = Modeller(self.data)
			mdlr.compute()
			print("Finish.")
		else:
			self.msg_box("Error!", "Error! Incorrect input!", QMessageBox.Critical)

	def msg_box(self, title, message, type):
		msg = QMessageBox(self)
		msg.setIcon(type)
		msg.setWindowTitle(title)
		msg.setText(message)
		msg.addButton('Ok', QMessageBox.AcceptRole)
		msg.exec()


def main():
	app = QtWidgets.QApplication(sys.argv)
	window = MyApp()
	window.show()
	app.exec_()


if __name__ == '__main__':
	main()
