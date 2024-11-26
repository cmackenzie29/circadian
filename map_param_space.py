import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy.linalg import eig

warnings.simplefilter("ignore")

X = 500
K1 = 1/X
dt = 1e-1

def sim(K0,K1,K2,K3,K4,K5,alpha):
	N = 100000
	b = np.zeros(N)
	c = np.zeros(N)
	p = np.zeros(N)
	d = np.zeros(N)

	b[0] = 0
	c[0] = 0
	p[0] = 0
	d[0] = 0

	
	for i in range(1,N):
		b[i] = b[i-1] + dt*(1/(1+d[i-1]**alpha)-K0*b[i-1])
		c[i] = c[i-1] + dt*(K2*b[i-1]+d[i-1]-K3*c[i-1])
		p[i] = p[i-1] + dt*(K4*b[i-1]+d[i-1]-K5*p[i-1])
		d[i] = d[i-1] + dt*(K1*c[i-1]*p[i-1]-d[i-1])

	# Start at first per peak
	start = p.argmax()
	if start*2 < len(p): # The start is within the first half of the simulation
		b = b[start:]
		c = c[start:]
		p = p[start:]
		if np.all(periods(p) != 0): # The solution oscillates
			return [b, c, p]
	return np.zeros(0)

def periods(a):
	ps = np.zeros(0)
	N = 3
	for _ in range(N):
		a = a[np.argmin(a):]
		if len(a)<=1:
			return np.zeros(N)
		a = a[np.argmax(a):]
		if len(a)<=1:
			return np.zeros(N)
		ps = np.append(ps, dt*(np.argmax(a[np.argmin(a):])+np.argmin(a)))
	return ps

# What subspace of alpha, K2, K4 produces stable oscillations?
for alpha in [0.9]:#np.arange(2,50)/100+0.01:
	K0 = 1/(X*(1+X**alpha))
	K0_d0 = -alpha*X**(alpha-1)/((1+X**alpha)**2)
	for K2 in np.linspace(0.5,2,10000):
		# Can we find where K2=K4?
		A = np.array([[-K0, 0, 0, K0_d0],
						[K2, -K2-1, 0, 1],
						[K2, 0, -K2-1, 1],
						[0, 1, 1, -1]])
		evals, _ = eig(A)
		if (not np.all(evals.imag==0)) and np.all(evals.real<=0):
			print(K2)






