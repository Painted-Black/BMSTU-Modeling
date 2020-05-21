from math import pi
from scipy.interpolate import InterpolatedUnivariateSpline, CubicSpline
from matplotlib import pyplot as plt
from numpy import arange
from scipy import integrate

class Modeller():
    def __init__(self, data):
        self.data = data

        self.itk = [[0.5, 6700, 0.5],
                    [1,	6790, 0.55],
                    [5,	7150, 1.7],
                    [10, 7270, 3],
                    [50, 8010, 11],
                    [200, 9185, 32],
                    [400, 10010, 40],
                    [800, 11140, 41],
                    [1200, 12010, 39]]
        
        self.tsigma = [
            [4000, 0.031],
            [5000, 0.27],
            [6000, 2.05],
            [7000, 6.06],
            [8000, 12.0],
            [9000, 19.9],
            [10000, 29.6],
            [11000, 41.1],
            [12000, 54.1],
            [13000, 67.7],
            [14000, 81.5]]

        self.m = None
        self.T0 = None
        #self.m4 = None
        #self.T04 = None
        self.simpsonN = 41

    def compute(self, isequal):
        tstart = self.data.get("Tstart")
        tend = self.data.get("Tend")
        tstep = self.data.get("Tstep")
        I0 = self.data.get("I0")
        Uc0 = self.data.get("Uc0")
        I04 = self.data.get("I0")
        Uc04 = self.data.get("Uc0")

        t_plot = []
        rp_plot = []
        rp_plot4 = []
        I_plot = []
        I_plot4 = []
        U_plot = []
        U_plot4 = []
        T0_plot = []
        IRp_plot = []
        IRp_plot4 = []

        for t in arange(tstart, tend, tstep):
            Rp = self.Rp(I0)
            if isequal:
                self.data['Rk'] = -Rp
            I0, Uc0 = self.runge_kutta2(I0, Uc0, tstep, Rp)

            t_plot.append(t)
            rp_plot.append(Rp)
            I_plot.append(I0)
            U_plot.append(Uc0)
            T0_plot.append(self.T0)
            IRp_plot.append(I0 * Rp)

            Rp4 = self.Rp(I04)
            if isequal:
                self.data['Rk'] = -Rp4
            I04, Uc04 = self.runge_kutta4(I04, Uc04, tstep, Rp4)
            rp_plot4.append(Rp4)
            I_plot4.append(I04)
            U_plot4.append(Uc04)
            IRp_plot4.append(I04 * Rp4)

        
        plt.figure(1)

        plt.suptitle("Методы Рунге-Кутта 2-го и 4-го порядка")
        plt.subplot(321)
        plt.plot(t_plot, rp_plot, label = '2-й порядок')
        plt.plot(t_plot, rp_plot4, label = '4-й порядок')
        plt.xlabel('t, с')
        plt.ylabel('Rp, Ом')
        plt.grid(True)
        plt.legend()

        plt.subplot(322)
        plt.plot(t_plot, I_plot, label = '2-й порядок')
        plt.plot(t_plot, I_plot4, label = '4-й порядок')
        plt.xlabel('t, с')
        plt.ylabel('I, А')
        plt.grid(True)
        plt.legend()

        plt.subplot(323)
        plt.plot(t_plot, U_plot, label = '2-й порядок')
        plt.plot(t_plot, U_plot4, label = '4-й порядок')
        plt.xlabel('t, с')
        plt.ylabel('U, В')
        plt.grid(True)
        plt.legend()

        plt.subplot(324)
        plt.plot(t_plot, T0_plot)
        plt.xlabel('t, с')
        plt.ylabel('T0, К')
        plt.grid(True)
        
        plt.subplot(325)
        plt.plot(t_plot, IRp_plot, label = '2-й порядок')
        plt.plot(t_plot, IRp_plot4, label = '4-й порядок')
        plt.xlabel('t, с')
        plt.ylabel('I * Rp, В')
        plt.grid(True)
        plt.legend()

        plt.show()

    def Rp(self, I):
        I_table = []
        itk_len = len(self.itk)
        for i in range(itk_len):
            I_table.append(self.itk[i][0])
        T0_table = []
        for i in range(itk_len):
            T0_table.append(self.itk[i][1])

        m_table = []
        for i in range(itk_len):
            m_table.append(self.itk[i][2])
        
        self.m = self.interpolate(I, I_table, m_table)
        self.T0 = self.interpolate(I, I_table, T0_table)

        intgl = self.simpson(self.integrand, 0, 1, self.simpsonN)
        #intgl = integrate.simps(self.integrand, 0, 1)
        #intgl = integrate.quad(self.integrand, 0, 1)
        Rp = self.data.get("Le") / (2 + pi * (self.data.get("R") ** 2) * intgl)
        return Rp

    def integrand(self, z):
        return self.sigma(self.Tz(z)) * z
        

    def Tz(self, z):
        return (self.T0 + (self.data.get("Tw") - self.T0) * (z ** self.m))
    
    def f(self, U, I, Rp):
        Rk = self.data.get("Rk")
        Lk = self.data.get("Lk")
        return ((U - (Rk + Rp) * I) / Lk)

    def phi(self, U, I):
        Ck = self.data.get("Ck")
        return - I / Ck
    
    def runge_kutta4(self, yn, zn, hn, Rp):
        hn2 = hn / 2

        k1 = hn * self.f(zn, yn, Rp)
        q1 = hn * self.phi(zn, yn)

        k2 = hn * self.f(zn + k1 / 2, yn + q1 / 2, Rp)
        q2 = hn * self.phi(zn + k1 / 2, yn + q1 / 2)
        
        k3 = hn * self.f(zn + k2 / 2, yn + q2 / 2, Rp)
        q3 = hn * self.phi(zn + k2 / 2, yn + q2 / 2)

        k4 = hn * self.f(zn + k3, yn + q3, Rp)
        q4 = hn * self.phi(zn + k3, yn + q3)

        y_next = yn + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        z_next = zn + (q1 + 2 * q2 + 2 * q3 + q4) / 6

        return y_next, z_next




    def runge_kutta2(self, yn, zn, hn, Rp): # I, U
        alpha = 0.5

        nh = hn / (2 * alpha)
        k1 = self.f(zn, yn, Rp)
        q1 = self.phi(zn, yn)
        k2 = self.f(zn + nh * q1, yn + nh * k1, Rp)
        q2 = self.phi(zn + nh * q1, yn + nh * k1)
        next_y = yn + hn * ((1 - alpha) * k1 + alpha * k2)
        next_z = zn + hn * ((1 - alpha) * q1 + alpha * q2)
        return next_y, next_z

    def sigma(self, T):
        T_table = []
        len_tsigma = len(self.tsigma)
        for i in range(len_tsigma):
            T_table.append(self.tsigma[i][0])
        sigma_table = []
        for i in range(len_tsigma):
            sigma_table.append(self.tsigma[i][1])
        
        return self.interpolate(T, T_table, sigma_table)
    
    def interpolate(self, I_known, kn_col1, kn_col2):
        #res = InterpolatedUnivariateSpline(kn_col1, kn_col2, k = 1)
        res = CubicSpline(kn_col1, kn_col2, extrapolate=True)
        return float(res(I_known))
    
    def simpson(self, func, a, b, N):
        h = float((b - a) / N)
        res_sum = 0
        for step in range(N):
            x1 = a + step * h
            x2 = a + (step + 1) * h
            res_sum += (x2 - x1) / 6.0 *(func(x1) + 4.0 * func(0.5 * (x1 + x2)) + func(x2))
        return res_sum
