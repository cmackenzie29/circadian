import numpy as np
import matplotlib.pyplot as plt
import pickle
import warnings

warnings.simplefilter("ignore")

def periods(a):
	ps = np.zeros(0)
	N = 3
	for _ in range(N):
		a = a[np.argmin(a):]
		if len(a)<=1:
			return np.zeros(3)
		a = a[np.argmax(a):]
		if len(a)<=1:
			return np.zeros(3)
		ps = np.append(ps, dt*(np.argmax(a[np.argmin(a):])+np.argmin(a)))
	return ps

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
	# If nothing returned by now, then simulation did not meet criteria with these parameters
	return np.zeros(0)

def evaluate(y, x, x_scale): # Takes a function discretized as numpy array and evaluates it at a continuous point, subject to scaling
	x /= x_scale
	left_x = int(np.floor(x))
	if left_x >= len(y)-1:
		return y[-1]
	left_y = y[left_x]
	right_y = y[left_x+1]
	return left_y + (x-left_x)*(right_y-left_y)

# Scale least squares by timing, so decay in importance as time goes on
dir = "/Users/cameron/Desktop/circadian/Data/"
train_data = [pickle.load(open(dir+'bmal.sav', 'rb')), pickle.load(open(dir+'cry.sav', 'rb')), pickle.load(open(dir+'per.sav', 'rb'))]

# Initial parameters are, somewhat arbitrarily, chosen to be alpha = 0.25, K2 = 0.1, K4 = 10
# Time scale needs to be stretched 
X = 500
dt = 1e-1
t_scale = dt/45

alpha = [1e-5,0.6]
K2 = [1e-3,100]
K4 = [1e-3,100]

best_wmse = 1e10

while best_wmse > 1:
	alpha_temp = alpha[0]+alpha[1]*np.random.rand()
	K2_temp = K2[0]+K2[1]*np.random.rand()
	K4_temp = K4[0]+K4[1]*np.random.rand()
	
	simulation = sim(1/(X*(1+X**alpha_temp)), 1/X, K2_temp, K2_temp+1, K4_temp, K4_temp+1, alpha_temp)
	if len(simulation) > 0:
		print("Oscillatory")
		wmse = 0
		for s, species in enumerate(train_data):
			sp_wmse = 0
			for train_series in species:
				for t in range(len(train_series)):
					update = np.exp(-t*t_scale/500)*(evaluate(simulation[s],train_series[t,0],t_scale)-train_series[t,1])**2
					sp_wmse += update
			wmse += sp_wmse**(1/len(species))
		if wmse < best_wmse:
			print(wmse)
			print(alpha_temp, K2_temp, K4_temp)
			print("\n")
			best_wmse = wmse


# Compare model to results
for s in range(len(train_data)):
	for train_series in train_data[s]:
		plt.plot(train_series[:,0], train_series[:,1])
	simulation = sim(1/(X*(1+X**alpha)), 1/X, K2, K2+1, K4, K4+1, alpha)
	plt.plot(t_scale*np.arange(len(simulation[s])),simulation[s],color="black")
	plt.show()

