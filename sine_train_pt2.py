import numpy as np
import matplotlib.pyplot as plt
import warnings

warnings.simplefilter("ignore")

# Now train deterministic model to idealized oscillations:
T = 2*np.pi/22.9966
cs = [1.19380521, 1.69646003, 4.20973416]
def training_data(x, i):
	return 150*np.sin(x*T-cs[i])+500

N = 7500
X = 500

def sim(alpha, K2, K4, k6):
	k = [X, K2*k6, K2*k6, K4*k6, K4*k6, k6/X, k6, X*k6, k6/2]
	
	b = np.zeros(N)
	c = np.zeros(N)
	p = np.zeros(N)
	d = np.zeros(N)

	b[0] = 0
	c[0] = 0
	p[0] = 0
	d[0] = 0

	dt = 1e-1

	for i in range(1,N):
		b[i] = b[i-1] + dt*(k[7]/(1+(d[i-1]/k[0])**alpha)-k[8]*b[i-1])
		c[i] = c[i-1] + dt*(k[1]*b[i-1]+k[6]*d[i-1]-k[2]*c[i-1]-k[5]*c[i-1]*p[i-1])
		p[i] = p[i-1] + dt*(k[3]*b[i-1]+k[6]*d[i-1]-k[4]*p[i-1]-k[5]*c[i-1]*p[i-1])
		d[i] = d[i-1] + dt*(k[5]*c[i-1]*p[i-1]-k[6]*d[i-1])

	# Check the start (at first per peak) is within the first half of the simulation and the result oscillates
	start = p.argmax()
	if start*2 < len(p):
		b = b[start:]
		c = c[start:]
		p = p[start:]
		if np.all(periods(p) != 0): # The solution oscillates
			# Start on the third peak
			b = b[p.argmax():]
			c = c[p.argmax():]
			p = p[p.argmax():]
			b = b[p.argmax():]
			c = c[p.argmax():]
			p = p[p.argmax():]
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

# Initial parameters are, somewhat arbitrarily, chosen to be alpha = 10, K2 = 1.0, K4 = 1.5, k6 = 0.25
dt = 1e-1

best_mse = 1e15
best_params = None
for alpha in np.linspace(55,75,21):
	print(alpha)
	for K2 in np.linspace(8.5,10.5,9):
		for K4 in np.linspace(2.25,3.25,9):
			for k6 in np.linspace(0.06,0.14,9):
				simulation = sim(alpha, K2, K4, k6)
				if len(simulation) > 0:
					mse = 0
					for s in range(3):
						mse += np.sum((simulation[s]-training_data(dt*np.arange(len(simulation[s])),s))**2)
					if mse < best_mse:
						best_mse = mse
						best_params = [alpha, K2, K4, k6]
						print(mse)
						print(best_params)
						print("\n")

# Compare model to results
for s in range(3):
	simulation = sim(best_params[0], best_params[1], best_params[2], best_params[3])
	plt.plot(dt*np.arange(len(simulation[s])), training_data(dt*np.arange(len(simulation[s])),s))
	plt.plot(dt*np.arange(len(simulation[s])), simulation[s],color="black")
	plt.show()

# Results = [71.0, 9.25, 2.375, 0.09]

