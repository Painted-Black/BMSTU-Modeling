import matplotlib.pyplot as plt
import numpy as np
from math import fabs, exp, fmod

class Modeller():
	def __init__(self, data):
		self.data = data

		self.a1 = self.data.get("a1")
		self.b1 = self.data.get("b1")
		self.c1 = self.data.get("c1")
		self.m1 = self.data.get("m1")

		self.a2 = self.data.get("a2")
		self.b2 = self.data.get("b2")
		self.c2 = self.data.get("c2")
		self.m2 = self.data.get("m2")
		
		self.alpha0 = self.data.get("alpha0")
		self.alphaN = self.data.get("alphaN")

		self.l = self.data.get("l")
		self.T0 = self.data.get("T0")
		self.R = self.data.get("R")

		self.F0_val = self.data.get("F0")
		self.Fmax = self.data.get("Fmax")
		self.tmax = self.data.get("tmax")

		self.period = self.data.get("period")
		self.tu = self.data.get("tu")

		self.h = self.data.get("h")
		self.t = self.data.get("t")

		self.d = (self.alphaN * self.l) / (self.alphaN - self.alpha0)
		self.c_koef = - self.alpha0 * self.d

		self.eps = 1e-4

	def F0(self, T):
		return self.F0_val
		#return 45
		#return (self.Fmax / self.tmax) * T * exp(-(T / self.tmax - 1))
		if T % self.period <= self.tu:
			return (self.Fmax / self.tmax) * T * exp(-(T / self.tmax - 1))
		else:
			return 0
		#nu = 1 / self.period
		#newT = round(fmod(T, nu) * 100) / 100.0
		#if (nu <= 0.01):
		#	newT = 10
		#	return (self.Fmax / self.tmax) * newT * exp(-(newT / self.tmax - 1))
		#return (self.Fmax / self.tmax) * T * exp(-(T / self.tmax - 1))

	def c(self, T):
		res = self.a2 + self.b2 * (T ** self.m2) - (self.c2 / (T ** 2))
		return res
		#return 0

	def f_plus_half(self, x, h, func):
		res = (func(x) + func(x + h)) / 2
		return res
	
	def f_minus_half(self, x, h, func):
		res = (func(x) + func(x - h)) / 2
		return res
	
	def old_k(self, x):
		res = self.a1 / (x - self.b1)
		return res

	def k(self, T):
		res = self.old_k(self.x)
		#res = self.a1 * (self.b1 + self.c1 * (T ** self.m1))
		return res

	def alpha(self, x):
		return self.c_koef / (x - self.d)

	def p(self, x):
		return 2 * self.alpha(x) / self.R

	def f(self, x):
		return 2 * self.alpha(x) * self.T0 / self.R

	def __left_boundary_condition(self, T, time):
		h8 = self.h / 8
		h4 = self.h / 4
		h2 = self.h / 2
		c_p12 = self.f_plus_half(T[0], self.t, self.c)
		chi_p12 = self.f_plus_half(T[0], self.t, self.k)
		t_over_h = self.t / self.h

		K0 = h8 * c_p12 + \
			 h4 * self.c(T[0]) + chi_p12 * t_over_h + \
			 self.t * h8 * self.p(h2) + (self.t * self.h) / 4 * self.p(0)

		M0 = h8 * c_p12 - chi_p12 * t_over_h + \
			 self.t * h8 * self.p(h2)

		P0 = h8 * c_p12 * (T[0] + T[1]) + \
			 h4 * self.c(T[0]) * T[0] + self.F0(time) * self.t + \
			 self.t * h8 * (3 * self.f(0) + self.f(self.h))
			 #self.t * h4 * (self.f(0) + self.f_plus_half(T[0], self.t, self.f)) # ?
		return K0, M0, P0

	def __right_boundary_condition(self, T):
		h8 = self.h / 8
		h4 = self.h / 4
		h2 = self.h / 2
		c_m12 = self.f_minus_half(T[-1], self.t, self.c)
		chi_m12 = self.f_minus_half(T[-1], self.t, self.k)

		h8_cm12 = h8 * c_m12 

		KN = h4 * self.c(T[-1]) + h8_cm12 + \
			self.t * self.alphaN + self.t * chi_m12 / self.h + \
			h4 * self.t * self.p(self.l) + h8 * self.t * self.p(self.l - h2)
		
		MN = h8_cm12 - \
	 		 self.t * chi_m12 / self.h + \
	 		 h8 * self.p(self.l - h2)
			
		PN = h4 * self.c(T[-1]) * T[-1] + h8_cm12 * T[-1] +\
	 		 h8_cm12 * T[-2] + \
	 		 self.t * self.alphaN * self.T0 + \
	 		 h4 * self.t * (self.f(self.l) + self.f(self.l - h2))
		
		return KN, MN, PN

	def A(self, T):
		return self.t / self.h * self.f_minus_half(T, self.t, self.k)


	def D(self, T):
		return self.t / self.h * self.f_plus_half(T, self.t, self.k)

	def B(self, x, T):
		return self.A(T) + self.D(T) + self.c(T) * self.h + \
			 self.p(x) * self.h * self.t

	def F(self, x, T):
		return self.f(x) * self.h * self.t + self.c(T) * T * self.h

	def progon(self, T, K0, M0, P0, KN, MN, PN):
		epsilon = [0, - M0 / K0]
		eta = [0, P0 / K0]

		x = self.h
		n = 1
		while x + self.h < self.l:
			An = self.A(T[n])
			Bn = self.B(x, T[n])

			newEps = self.D(T[n]) / (Bn - An * epsilon[n])
			epsilon.append(newEps)

			newEta = (self.F(x, T[n]) + An * eta[n]) / (Bn - An * epsilon[n])
			eta.append(newEta)

			x += self.h
			n += 1
		
		T = [0] * (n + 1)
		T[n] = (PN - MN * eta[n]) / (KN + MN * epsilon[n])

		for i in range(n-1, -1, -1):
			T[i] = epsilon[i + 1] * T[i + 1] + eta[i+1]

		return T

	def check_iter(self, T, T_next):
		# если нагрев замедлилися и изменение температуры за 
		# шаг по времени достаточно мало
		max = fabs((T[0] - T_next[0]) / T_next[0])
		for i, j in zip(T, T_next):
			d = fabs(i - j) / j
			if d > max:
				max = d
		return max < 1

	def check_epsilon(self, T, T_next):
		# если максимально изменившийся элемент
		# изменился меньше, чем на eps
		for i, j in zip(T, T_next):
			
			if fabs((i - j) / j) > self.eps:
				return True
		return False

	def __iterate(self):
		res = []
		n = int(self.l / self.h)
		print(n)
		T = [self.T0] * (n + 1)
		T_new = [0] * (n + 1)
		ti = 0
		res.append(T)
		while True:
			tmp = T
			self.x = self.h
			ti += self.t
			while True:
				K0, M0, P0 = self.__left_boundary_condition(tmp, ti)
				KN, MN, PN = self.__right_boundary_condition(tmp)
				T_new = self.progon(tmp, K0, M0, P0, KN, MN, PN)
				self.x += self.h

				if self.check_iter(tmp, T_new):
					break
				tmp = T_new

			res.append(T_new)
			#ti += self.t
			if not self.check_epsilon(T, T_new):
				break
			T = T_new
		return res, ti

	def compute(self):
		periods = [1]
		for curp in periods:
			result = []
			resulti = []
			self.period = curp
			res, ti = self.__iterate()
			print(ti)
			te = []

			i = 0
			while i < ti:
				te.append(i)
				i += self.t
			
			i = 0
			xfix = [temp[int(i / self.h)] for temp in res]
			result.append(xfix[:-1])
			resulti.append(te)
		i = 0
		for res, teres in zip(result, resulti):
			#l = 'nu = ' + str(1/periods[i])
			plt.plot(teres, res)
			i += 1
		
		plt.xlabel("t, sec")
		plt.ylabel("T, K")
		plt.grid()
		plt.legend()
		plt.show()

		'''
		Fmaxes = [10, 15, 25, 35, 45, 50, 75, 100, 125, 150]
		result = []
		resulti = []
		for curFmax in Fmaxes:
			self.Fmax = curFmax
			res, ti = self.__iterate()
			print(ti)
			te = []

			i = 0
			while i < ti:
				te.append(i)
				i += self.t
			
			i = 0
			xfix = [temp[int(i / self.h)] for temp in res]
			result.append(xfix[:-1])
			resulti.append(te)
		i = 0
		for res, teres in zip(result, resulti):
			l = 'Fmax = ' + str(Fmaxes[i])
			plt.plot(teres, res, label=l)
			i += 1
		
		plt.xlabel("t, sec")
		plt.ylabel("T, K")
		plt.grid()
		plt.legend()
		plt.show()
		'''
		'''
		a2s = [2.0493, 3, 4, 5]
		result = []
		resulti = []
		for cura in a2s:
			self.a2 = cura
			res, ti = self.__iterate()
			print(ti)
			te = []

			i = 0
			while i < ti:
				te.append(i)
				i += self.t
			
			i = 0
			xfix = [temp[int(i / self.h)] for temp in res]
			result.append(xfix[:-1])
			resulti.append(te)
		i = 0
		for res, teres in zip(result, resulti):
			l = 'a2 = ' + str(a2s[i])
			plt.plot(teres, res, label=l)
			i += 1

			#xes = [i for i in np.arange(0, self.l, self.h)]
			#ts = [i for i in range(0, int(ti), int(self.t))]
			#ts = [i for i in np.arange(0, ti, self.t)]
			
			#s = 0
			#p = [j[int(s / self.h)] for j in res]
			#plt.plot(ts, p)
			#s = 0
			#while s < self.l / 3:
			#	p = [j[int(s / self.h)] for j in res]
			#	plt.plot(ts, p[:-1])
			#	s += 0.1

		
		plt.xlabel("t, sec")
		plt.ylabel("T, K")
		plt.grid()
		plt.legend()
		plt.show()
		'''
		'''
		hes = [0.3, 0.2, 0.1, 0.01]
		#taus = [2, 1, 0.1, 0.05]
		

		result = []

		for hi in hes:
			self.h = hi
			res, ti = self.__iterate()
			print('h = ', self.h, 'tau = ', self.t)
			i = 0
			xfix = [temp[int(i / self.h)] for temp in res]
			result.append(xfix)

		#print(len(hes))
		#print(len(result))
		count = 0
		print('{:^10}|{:^10}|{:^10}|{:^10}|'.format(hes[0],hes[1],hes[2],hes[3]))
		for i in range(len(result[0])):
			st = ''
			for j in range(len(result)):
				if count % 20 == 0:
					st += '{:^10.3f}|'.format(result[j][i])

				#print('{:^10.3f}|{:^10.3f}|{:^10.3f}|{:^10.3f}'.format())
				#print('  %3.3f  |' % result[j][i], end = '')
			if count % 20 == 0:
				print(st)
			count += 1
		'''
		