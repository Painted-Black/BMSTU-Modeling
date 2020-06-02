import sys

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMessageBox

from Ui_mainwindow import Ui_MainWindow

from Modeller import Modeller

class MyApp(QtWidgets.QMainWindow):
	def __init__(self):
		super(MyApp, self).__init__()
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)

		self.ui.set_def_button.clicked.connect(self.set_defaults)
		self.ui.run_button.clicked.connect(self.run)

		self.defaults = {
			"a1" 	 : 0.0134,
			"b1" 	 : 1,
			"c1" 	 : 4.35e-4,
			"m1" 	 : 1,
			"a2" 	 : 2.049,
			"b2" 	 : 0.563e-3,
			"c2" 	 : 0.528e5,
			"m2" 	 : 1,
			"alpha0" : 0.05,
			"alphaN" : 0.01,
			"l"	 	 : 10,
			"T0"	 : 300,
			"R"		 : 0.5,
			"F0"	 : 50,
			"h"		 : 0.001,
			"t"		 : 1
		}

		self.data = {
			"a1" 	 : None,
			"b1" 	 : None,
			"c1" 	 : None,
			"m1" 	 : None,
			"a2" 	 : None,
			"b2" 	 : None,
			"c2" 	 : None,
			"m2" 	 : None,
			"alpha0" : None,
			"alphaN" : None,
			"l"	 	 : None,
			"T0"	 : None,
			"R"		 : None,
			"F0"	 : None,
			"h"		 : None,
			"t"		 : None
		}

		self.set_defaults()
	
	def set_defaults(self):
		self.ui.lineEdit_a1.setText(str(self.defaults.get("a1")))
		self.ui.lineEdit_b1.setText(str(self.defaults.get("b1")))
		self.ui.lineEdit_c1.setText(str(self.defaults.get("c1")))
		self.ui.lineEdit_m1.setText(str(self.defaults.get("m1")))

		self.ui.lineEdit_a2.setText(str(self.defaults.get("a2")))
		self.ui.lineEdit_b2.setText(str(self.defaults.get("b2")))
		self.ui.lineEdit_c2.setText(str(self.defaults.get("c2")))
		self.ui.lineEdit_m2.setText(str(self.defaults.get("m2")))

		self.ui.lineEdit_alpha0.setText(str(self.defaults.get("alpha0")))
		self.ui.lineEdit_alphaN.setText(str(self.defaults.get("alphaN")))

		self.ui.lineEdit_l.setText(str(self.defaults.get("l")))
		self.ui.lineEdit_T0.setText(str(self.defaults.get("T0")))
		self.ui.lineEdit_R.setText(str(self.defaults.get("R")))
		self.ui.lineEdit_F0.setText(str(self.defaults.get("F0")))

		self.ui.lineEdit_h.setText(str(self.defaults.get("h")))
		self.ui.lineEdit_t.setText(str(self.defaults.get("t")))

	def get_data(self):
		try:
			self.data["a1"] = float(self.ui.lineEdit_a1.text())
			self.data["b1"] = float(self.ui.lineEdit_b1.text())
			self.data["c1"] = float(self.ui.lineEdit_c1.text())
			self.data["m1"] = float(self.ui.lineEdit_m1.text())

			self.data["a2"] = float(self.ui.lineEdit_a2.text())
			self.data["b2"] = float(self.ui.lineEdit_b2.text())
			self.data["c2"] = float(self.ui.lineEdit_c2.text())
			self.data["m2"] = float(self.ui.lineEdit_m2.text())

			self.data["alpha0"] = float(self.ui.lineEdit_alpha0.text())
			self.data["alphaN"] = float(self.ui.lineEdit_alphaN.text())

			self.data["l"] = float(self.ui.lineEdit_l.text())
			self.data["T0"] = float(self.ui.lineEdit_T0.text())
			self.data["R"] = float(self.ui.lineEdit_R.text())
			self.data["F0"] = float(self.ui.lineEdit_F0.text())

			self.data["h"] = float(self.ui.lineEdit_h.text())
			self.data["t"] = float(self.ui.lineEdit_t.text())

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
