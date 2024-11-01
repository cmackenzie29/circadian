# import numpy as np
# import matplotlib.pyplot as plt
# # from scipy.optimize import fsolve
# # from scipy.linalg import eig
# # import warnings

# N = 50000
# dt = 1e-3

# for alpha in range(1,20,100):
# 	for K2 in np.logspace(-3,3,100):
# 		for K4 in np.logspace(-3,3,100):
# 			K = np.array([0.5,1,K2,K2+1,K4,K4+1])

# 			b = np.zeros(N)
# 			c = np.zeros(N)
# 			p = np.zeros(N)
# 			d = np.zeros(N)

# 			b[0] = 0.9
# 			c[0] = 1.1
# 			p[0] = 1.05
# 			d[0] = 0.95

# 			for i in range(1,N):
# 				b[i] = b[i-1] + dt*(1/(1+d[i-1]**alpha)-K[0]*b[i-1])
# 				c[i] = c[i-1] + dt*(K[2]*b[i-1]+d[i-1]-K[3]*c[i-1])
# 				p[i] = p[i-1] + dt*(K[4]*b[i-1]+d[i-1]-K[5]*p[i-1])
# 				d[i] = d[i-1] + dt*(K[1]*c[i-1]*p[i-1]-d[i-1])

# 			# Test fit to real data

### ORIGINAL SYSTEM WITH FIVE SPECIES
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.linalg import eig
import warnings

N = 50000
dt = 1e-3

def test_params(k):
	# Assumes a0=b0=c0=d0=p0=1
	A = np.array([[-k[1]-k[6], -k[6], k[0], k[7], 0],
					[-k[6], -k[3]-k[6], k[2], k[7], 0],
					[0, 0, -k[5]-k[8], k[4]-k[8], k[9]],
					[k[6], k[6], -k[8], -k[7]-k[8], k[9]],
					[0, 0, k[8], k[8], -k[9]]])
	evals, _ = eig(A)
	# Stable if all evals have 0/negative real part.
	# Stable orbital oscillation if all evals have 0 real part and +/- complex part.
	if np.all(evals.real <= 0): #and (not np.all(evals.imag==0)) and np.all(np.abs(evals.real)<10):
		print(k)
		print(evals)
		print("\n")

for k1 in np.logspace(-3,3,5):
	for k3 in np.logspace(-3,3,5):
		for k5 in np.logspace(-3,3,5):
			for k7 in np.logspace(-3,3,5):
				for k9 in np.logspace(-3,3,5):
					test_params(np.array([k1,k1,k3,k3,k5,k5,k7,k7,k9,k9]))
		



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





