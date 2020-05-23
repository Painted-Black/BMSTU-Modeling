import matplotlib.pyplot as plt

class Modeller():
	def __init__(self, data):
		self.data = data

		self.k0 = self.data.get("k0")
		self.kN = self.data.get("kN")
		self.a0 = self.data.get("a0")
		self.aN = self.data.get("aN")
		self.l = self.data.get("l")
		self.T0 = self.data.get("T0")
		self.R = self.data.get("R")
		self.F0 = self.data.get("F0")
		self.h = self.data.get("h")

		self.start_l = 0

		self.b = (self.kN * self.l) / (self.kN - self.k0)
		self.a = - self.k0 * self.b
		self.d = (self.aN * self.l) / (self.aN - self.a0)
		self.c = - self.a0 * self.d

		self.M0, self.K0, self.P0 = self.left_boundary_condition()
		self.MN, self.KN, self.PN = self.right_boundary_condition()

	# приближенно вычислим интеграл методом средних
	def chi_plus_half(self, func, x):
		return (func(x) + func(x + self.h)) / 2
	
	# приближенно вычислим интеграл методом средних
	def chi_minus_half(self, func, x):
		return (func(x) + func(x - self.h)) / 2

	def right_boundary_condition(self):
		chi = self.chi_minus_half(self.k, self.l)

		pn = self.p(self.l)
		fn = self.f(self.l)

		pn12 = (pn + self.p(self.start_l - self.h)) / 2
		fn12 = (fn + self.f(self.start_l - self.h)) / 2

		MN = chi / self.h - (self.h * pn12) / 8
		KN = -self.aN - chi / self.h - (self.h * pn) / 4 - (self.h * pn12) / 8
		PN = - self.aN * self.T0 - self.h / 4 * (fn + fn12)

		return MN, KN, PN
		
	def left_boundary_condition(self):
		chi = self.chi_minus_half(self.k, self.start_l)
		h_2 = self.h ** 2
		p0 = self.p(self.start_l)
		f0 = self.f(self.start_l)
		p12 = (p0 + self.p(self.start_l + self.h)) / 2
		f12 = (f0 + self.f(self.start_l + self.h)) / 2
		M0 = - chi + h_2 / 8 * p12
		K0 = chi + h_2 / 8 * p12 + h_2 / 4 * p0
		P0 = self.h * self.F0 + h_2 / 4 * (f12 + f0)

		return M0, K0, P0

	def compute(self):
		epsilon = [0, - self.M0 / self.K0]
		eta = [0, self.P0 / self.K0]
		xes = [0, self.h]

		# вычисляем значения прогоночных коэффицентов
		x = self.h
		n = 1
		while x + self.h < self.l:
			newEps = self.C(x) / (self.B(x) - self.A(x) * epsilon[n])
			epsilon.append(newEps)

			newEta = (self.D(x) + self.A(x) * eta[n]) / (self.B(x) - self.A(x) * epsilon[n])
			eta.append(newEta)

			x += self.h
			n += 1
			xes.append(x)

		# обратный ход

		T = [0] * (n + 1)
		T[n] = (self.PN - self.MN * eta[n]) / (self.KN + self.MN * epsilon[n])

		for i in range(n-1, -1, -1):
			T[i] = epsilon[i + 1] *T[i + 1] + eta[i+1]

		self.plot1(xes, T)

	def plot1(self, x, y):
		print(y)
		plt.plot(x, y)
		plt.xlabel("x, cm")
		plt.ylabel("t, K")
		plt.grid()
		plt.show()
	
	def p(self, x):
		return (2 * self.alpha(x)) / self.R

	def f(self, x):
		return (2 * self.T0 * self.alpha(x)) / self.R

	def k(self, x):
		return self.a / (x - self.b)

	def alpha(self, x):
		return self.c / (x - self.d)

	def A(self, n):
		return self.chi_plus_half(self.k, n) / self.h

	def C(self, n):
		return self.chi_minus_half(self.k, n) / self.h

	def B(self, n):
		return self.A(n) + self.C(n) + self.p(n) * self.h

	def D(self, n):
		return self.f(n) * self.h
