import numpy as np
import matplotlib.pyplot as plt
# from scipy.optimize import fsolve
# from scipy.linalg import eig
# import warnings

N = 50000
dt = 1e-3
alpha = 10

K = np.array([0.5,1,1,2,1,2])
k0 = 1
k6 = 1
k7 = 1
k = np.array([k0, K[2]*k6**2*k0/k7, K[3]*k6, K[4]*k6**2*k0/k7, K[5]*k6, K[1]*k6/k0, k6, k7, K[0]*k6])

# K simulation
b = np.zeros(N)
c = np.zeros(N)
p = np.zeros(N)
d = np.zeros(N)

b[0] = 0.9
c[0] = 1.1
p[0] = 1.05
d[0] = 0.95

for i in range(1,N):
	b[i] = b[i-1] + dt*(1/(1+d[i-1]**alpha)-K[0]*b[i-1])
	c[i] = c[i-1] + dt*(K[2]*b[i-1]+d[i-1]-K[3]*c[i-1])
	p[i] = p[i-1] + dt*(K[4]*b[i-1]+d[i-1]-K[5]*p[i-1])
	d[i] = d[i-1] + dt*(K[1]*c[i-1]*p[i-1]-d[i-1])

plt.subplot(1,2,1)
plt.plot(np.arange(N),b)
plt.plot(np.arange(N),c)
plt.plot(np.arange(N),p)
plt.plot(np.arange(N),d)

# k simulation
b = np.zeros(N)
c = np.zeros(N)
p = np.zeros(N)
d = np.zeros(N)

b[0] = 0.9
c[0] = 1.1
p[0] = 1.05
d[0] = 0.95

for i in range(1,N):
	b[i] = b[i-1] + dt*(k[7]/(1+(d[i-1]/k[0])**alpha)-k[8]*b[i-1])
	c[i] = c[i-1] + dt*(k[1]*b[i-1]+k[6]*d[i-1]-k[2]*c[i-1])
	p[i] = p[i-1] + dt*(k[3]*b[i-1]+k[6]*d[i-1]-k[4]*p[i-1])
	d[i] = d[i-1] + dt*(k[5]*c[i-1]*p[i-1]-k[6]*d[i-1])

plt.subplot(1,2,2)
plt.plot(np.arange(N),b)
plt.plot(np.arange(N),c)
plt.plot(np.arange(N),p)
plt.plot(np.arange(N),d)
plt.show()

# warnings.simplefilter("ignore")

# alpha = 10

# def test_params(K0,K1,K2,K3,K4,K5):
# 	f_d = lambda d: K2*K4*(1/(K0*(1+d**alpha)))**2+(K2+K4)*(1/(K0*(1+d**alpha)))*d+d**2-K3*K5*d/K1
# 	d0 = fsolve(f_d, 1)[0]

# 	# Check that d0 converged
# 	if np.abs(f_d(d0)) < 1e-5:
# 		b0 = 1/(K0*(1+d0**alpha))
# 		c0 = (K2*b0+d0)/K3
# 		p0 = (K4*b0+d0)/K5

# 		A = np.array([[-K0, 0, 0, -alpha*d0**(alpha-1)/((1+d0**alpha)**2)],
# 						[K2, -K3, 0, 1],
# 						[K4, 0, -K5, 1],
# 						[0, K1*p0, K1*c0, -1]])
# 		evals, _ = eig(A)
# 		# Stable if all evals have 0/negative real part.
# 		# Stable orbital oscillation if all evals have 0 real part and +/- complex part.
# 		if (not np.all(evals.imag==0)) and np.all(np.abs(evals.real)<3):
# 			print(K0,K1,K2,K3,K4,K5)
# 			print(d0,b0,c0,p0)
# 			print(evals)
# 			print("\n")


# for k2 in np.logspace(-3,3,10):
# 	for k4 in np.logspace(-3,3,10):
# 		test_params(0.5, 1, k2, k2+1, k4, k4+1)



# import numpy as np
# from scipy.optimize import fsolve
# from scipy.linalg import eig
# import warnings

# warnings.simplefilter("ignore")

# alpha = 10


# def test_params(K0,K1,K2,K3,K4,K5):
# 	f_d = lambda d: K2*K4*(1/(K0*(1+d**alpha)))**2+(K2+K4)*(1/(K0*(1+d**alpha)))*d+d**2-K3*K5*d/K1
# 	d0 = fsolve(f_d, 1)[0]

# 	# Check that d0 converged
# 	if np.abs(f_d(d0)) < 1e-5:
# 		b0 = 1/(K0*(1+d0**alpha))
# 		c0 = (K2*b0+d0)/K3
# 		p0 = (K4*b0+d0)/K5

# 		A = np.array([[-K0, 0, 0, -alpha*d0**(alpha-1)/((1+d0**alpha)**2)],
# 						[K2, -K3, 0, 1],
# 						[K4, 0, -K5, 1],
# 						[0, K1*p0, K1*c0, -1]])
# 		evals, _ = eig(A)
# 		# Stable if all evals have 0/negative real part.
# 		# Stable orbital oscillation if all evals have 0 real part and +/- complex part.
# 		if (not np.all(evals.imag==0)) and np.all(np.abs(evals.real)<3):
# 			print(K0,K1,K2,K3,K4,K5)
# 			print(d0,b0,c0,p0)
# 			print(evals)
# 			print("\n")


# for k2 in np.logspace(-3,3,10):
# 	for k4 in np.logspace(-3,3,10):
# 		test_params(0.5, 1, k2, k2+1, k4, k4+1)





