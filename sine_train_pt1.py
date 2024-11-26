import numpy as np
import matplotlib.pyplot as plt
import pickle
import warnings

warnings.simplefilter("ignore")

# Scale least squares by timing, so decay in importance as time goes on
dir = "/Users/cameron/Desktop/circadian/Data/"
train_data = [pickle.load(open(dir+'bmal.sav', 'rb')), pickle.load(open(dir+'cry.sav', 'rb')), pickle.load(open(dir+'per.sav', 'rb'))]

# Estimate periods
Ts = np.zeros(0)
for species in train_data:
	for train_series in species:
		t_series = train_series[:,0]
		a = train_series[:,1]
		t = t_series[0]

		# Truncate noise during first two hours of series
		t = 2
		a = a[np.sum(t_series<t):]
		t_series = t_series[np.sum(t_series<t):]
		
		while (t_series[-1] - t) > 26: # More than 26 hours of the series left to analyze
			first_max = np.argmax(a[:np.sum(t_series<(t+15))])
			t_first_max = t_series[first_max]

			per_min = first_max + np.argmin(a[first_max:np.sum(t_series<(t_first_max+20))])
			t_min = t_series[per_min]
			
			next_max = per_min + np.argmax(a[per_min:np.sum(t_series<(t_min+15))])
			t_next_max = t_series[next_max]

			Ts = np.append(Ts, t_next_max-t_first_max)
			a = a[next_max:]
			t_series = t_series[next_max:]
			t = t_next_max
		
T = 2*np.pi/Ts.mean()
# Calculate the phase for each
cs = np.zeros(3)
for s in range(3):
	best_mse = 1e10
	best_c = 0
	for c in np.arange(100)*np.pi/50:
		for train_series in train_data[s]:
			mse = np.sum(((train_series[:,1]) - (150*np.sin(train_series[:,0]*T-c)+500))**2)
			if mse < best_mse:
				best_mse = mse
				best_c = c
	cs[s] = best_c

print(f"Results: sin(2pi*t/{Ts.mean()} + c), where c = {cs}")

# Compare model to results
for s in range(3):
	for train_series in train_data[s]:
		plt.plot(train_series[:,0], train_series[:,1])
		plt.plot(train_series[:,0], 150*np.sin(train_series[:,0]*T-cs[s])+500, color="black")
	plt.show()

