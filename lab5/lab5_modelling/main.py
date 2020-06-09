from Modeller import Modeller

defaults = {
			"a1" 	 : 1.333,
			"b1" 	 : -3.333,
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
			"Fmax"	 : 0.1,
			"tmax"	 : 10,
			"h"		 : 0.01,
			"t"		 : 1,
			"period" : 1,
			"tu"	 : 0.9,
			"F0"	 : 50
		}

def main():
	print("Computing...")
	mdlr = Modeller(defaults)
	mdlr.compute()
	print("Finish.")


if __name__ == '__main__':
	main()
