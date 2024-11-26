# import numpy as np
# import matplotlib.pyplot as plt

# N = 50000

# a = np.zeros(N)
# b = np.zeros(N)
# c = np.zeros(N)
# p = np.zeros(N)
# d = np.zeros(N)

# a[0] = 0.98
# b[0] = 0.9
# c[0] = 1.1
# p[0] = 1.05
# d[0] = 0.95

# dt = 1e-1
# k1 = 1
# k3 = 1
# k5 = 1
# k7 = 0.001
# k9 = 0.001
# k = np.array([0,k1,k1,k3,k3,k5,k5,k7,k7,k9,k9])

# for i in range(1,N):
# 	p[i] = p[i-1] + dt*(k[1]*b[i-1]+k[8]*d[i-1]-k[2]*p[i-1]-k[7]*c[i-1]*p[i-1])
# 	c[i] = c[i-1] + dt*(k[3]*b[i-1]+k[8]*d[i-1]-k[4]*c[i-1]-k[7]*c[i-1]*p[i-1])
# 	b[i] = b[i-1] + dt*(k[5]*d[i-1]+k[10]*a[i-1]-k[6]*b[i-1]-k[9]*b[i-1]*d[i-1])
# 	d[i] = d[i-1] + dt*(k[7]*c[i-1]*p[i-1]-k[8]*d[i-1]+k[10]*a[i-1]-k[9]*b[i-1]*d[i-1])
# 	a[i] = a[i-1] + dt*(k[9]*b[i-1]*d[i-1]-k[10]*a[i-1])

# plt.plot(np.arange(N),b)
# plt.plot(np.arange(N),c)
# plt.plot(np.arange(N),p)
# plt.plot(np.arange(N),d)
# plt.plot(np.arange(N),a)
# plt.show()

# import numpy as np
# import matplotlib.pyplot as plt

# N = 25000

# b = np.zeros(N)
# c = np.zeros(N)
# p = np.zeros(N)
# d = np.zeros(N)

# X = 500

# b[0] = 0.98*X
# c[0] = 0.96*X
# p[0] = 1.02*X
# d[0] = 1.04*X

# dt = 1e-1
# alpha = 0.15789473684210525
# k = np.array([0.0005452861026334875,0.002,0.1,1.1,10.0,11.0])

# for i in range(1,N):
# 	b[i] = b[i-1] + dt*(1/(1+d[i-1]**alpha)-k[0]*b[i-1])
# 	c[i] = c[i-1] + dt*(k[2]*b[i-1]+d[i-1]-k[3]*c[i-1])
# 	p[i] = p[i-1] + dt*(k[4]*b[i-1]+d[i-1]-k[5]*p[i-1])
# 	d[i] = d[i-1] + dt*(k[1]*c[i-1]*p[i-1]-d[i-1])

# plt.plot(dt*np.arange(N),b,label="B")
# plt.plot(dt*np.arange(N),c,label="C")
# plt.plot(dt*np.arange(N),p,label="P")
# plt.plot(dt*np.arange(N),d,label="D")
# plt.legend()
# plt.show()

import numpy as np

# Forward Euler method to solve dimensionless ODE system
# X: stable value around which the species should fluctuate
# alpha, K2, and K4 are the free parameters
def dimensionless_model(X, alpha, K2, K4):
	N = 250 # Number of timesteps
	dt = 1e-1 # Timestep

	# K0 and K1 are fixed by the value around which the species fluctuate
	# K3 and K5 are fixed by K2 and K4, respectively
	k = [1/(X*(1+X**alpha)), 1/X, K2, K2, K4, K4]
	
	# Species will all start at 0
	b = np.zeros(N)
	c = np.zeros(N)
	p = np.zeros(N)
	d = np.zeros(N)

	for i in range(1,N):
		b[i] = b[i-1] + dt*(1/(1+(d[i-1]/1)**alpha)-k[0]*b[i-1])
		c[i] = c[i-1] + dt*(k[2]*b[i-1]+d[i-1]-k[3]*c[i-1]-k[1]*c[i-1]*p[i-1])
		p[i] = p[i-1] + dt*(k[4]*b[i-1]+d[i-1]-k[5]*p[i-1]-k[1]*c[i-1]*p[i-1])
		d[i] = d[i-1] + dt*(k[1]*c[i-1]*p[i-1]-d[i-1])

	return b, c, p, d


# from scipy.linalg import eig
# A = np.array([[-k[0], 0, 0, -alpha*X**(alpha-1)/((1+X**alpha)**2)],
# 						[k[2], -k[3], 0, 1],
# 						[k[4], 0, -k[5], 1],
# 						[0, k[1]*X, k[1]*X, -1]])
# evals, _ = eig(A)
# print(evals)

# e = (alpha*X**(alpha-1))*(2*K2*K4+K2+K4)/(1+X**alpha)**2+(K2*K4-1)/(X*(1+X**alpha))
# d = (alpha*X**(alpha-1))*(K2+K4)/(1+X**alpha)**2+(2*K2+2*K4+K2*K4+1)/(X*(1+X**alpha))+(K2*K4-1)
# c = (K2+K4+3)/(X*(1+X**alpha))+2*K2+2*K4+K2*K4+1
# b = (1)/(X*(1+X**alpha))+K2+K4+3
# print(f"Delta = {256*e**3-192*b*d*e**2-128*c**2*e**2+144*c*d**2*e-27*d**4+144*b**2*c*e**2-6*b**2*d**2*e-80*b*c**2*d*e+18*b*c*d**3+16*c**4*e-4*c**3*d**2-27*b**4*e**2+18*b**3*c*d*e-4*b**3*d**3-4*b**2*c**3*e+b**2*c**2*d**2}")
# print(f"D = {64*e-16*c**2+16*b**2*c-16*b*d-3*b**4}")
# print(f"P = {8*c-3*b**2}")


