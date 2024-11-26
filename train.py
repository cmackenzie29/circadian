import numpy as np
import matplotlib.pyplot as plt
import pickle
import warnings
from simulate import deterministic_model

warnings.simplefilter("ignore")

dt = 1e-1

# Starting parameters
alpha = 10
K2 = 5
K4 = 5
k6 = 0.1

# Function to find the first three periods of the series
def periods(a):
	pers = np.zeros(0)
	for _ in range(3):
		a = a[np.argmin(a):]
		if len(a)<=1:
			return np.zeros(3) # Return zeros if it doesn't oscillate 
		a = a[np.argmax(a):]
		if len(a)<=1:
			return np.zeros(3)
		pers = np.append(pers, dt*(np.argmax(a[np.argmin(a):])+np.argmin(a)))
	return pers

# Takes a function discretized as numpy array and evaluates it at a continuous point, subject to scaling
def evaluate(y, x, x_scale):
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

prev_wmse = 1e10
current_wmse = 1e9

delta = 0.01
while delta > 0.002:
	best_wmse = 1e10
	best_err_index = None
	for i in [0,-1,1]:
		for j in [0,-1,1]:
			for k in [0,-1,1]:
				for l in [0,-1,1]:
					alpha_temp = alpha * (1+delta*i)
					K2_temp = K2 * (1+delta*j)
					K4_temp = K4 * (1+delta*k)
					k6_temp = k6 * (1+delta*l)
					simulation = deterministic_model(alpha_temp, K2_temp, K4_temp, k6_temp)

					# Check the solution oscillates
					if np.all(periods(simulation[2]) != 0):
						wmse = 0
						for s, species in enumerate(train_data):
							sp_wmse = 0
							for train_series in species:
								# Loop through series for up to 7 days
								for t in range(round(min(np.floor(len(train_series)*dt), 7*24)/dt)):
									update = (evaluate(simulation[s],train_series[t,0],dt)-train_series[t,1])**2
									sp_wmse += update
							wmse += sp_wmse**(1/len(species))
						if wmse < best_wmse:
							best_wmse = wmse
							best_err_index = [i,j,k,l]

	if best_wmse >= current_wmse: # No improvements
		delta *= 0.5
		print(f"No improvement. Decrease delta to {delta}")
	else: # Improvements happened. Update parameters.
		alpha = alpha * (1+delta*best_err_index[0])
		K2 = K2 * (1+delta*best_err_index[1])
		K4 = K4 * (1+delta*best_err_index[2])
		k6 = k6 * (1+delta*best_err_index[3])
		prev_wmse = current_wmse
		current_wmse = best_wmse
		print(f"Updated with error {current_wmse}. Explore now with delta = {delta}.")
		print(f"Alpha = {alpha}, K2 = {K2}, K4 = {K4}, k6 = {k6}")
		print("\n")

# Compare model to results
for s in range(len(train_data)):
	for train_series in train_data[s]:
		plt.plot(train_series[:,0], train_series[:,1])
	simulation = sim(alpha, K2, K4, k6)
	plt.plot(dt*np.arange(len(simulation[s])),simulation[s],color="black")
	plt.show()


