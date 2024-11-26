import numpy as np

# k6: parameter that scales time between dimensionless and full model
# scale: factor by which to multiply all rate constants for section 4.2
# scale_k: allows a single parameter to be multiplied by a scalar for section 4.1
def deterministic_model(alpha, K2, K4, k6, scale=1, scale_k=None):
	N = 4000 # Number of timesteps
	dt = 1e-1 # Timestep
	X = 500 # Stable value around which species should fluctuate

	k = np.array([X, K2*k6, K2*k6, K4*k6, K4*k6, k6/X, k6, X*k6, k6/2])*scale

	if scale_k != None: # Scale a rate constant kn if given some scale_k = [n, scale]
		k[scale_k[0]] *= scale_k[1]
		
	b = np.zeros(N)
	c = np.zeros(N)
	p = np.zeros(N)
	d = np.zeros(N)

	for i in range(1,N):
		# Forward Euler
		b[i] = b[i-1] + dt*(k[7]/(1+(d[i-1]/k[0])**alpha)-k[8]*b[i-1])
		c[i] = c[i-1] + dt*(k[1]*b[i-1]+k[6]*d[i-1]-k[2]*c[i-1]-k[5]*c[i-1]*p[i-1])
		p[i] = p[i-1] + dt*(k[3]*b[i-1]+k[6]*d[i-1]-k[4]*p[i-1]-k[5]*c[i-1]*p[i-1])
		d[i] = d[i-1] + dt*(k[5]*c[i-1]*p[i-1]-k[6]*d[i-1])

	# Return values starting from the first PER peak
	# Mutant simulations will not always oscillate, so turn this off in case of scaling
	if scale_k == None:
		b = b[np.argmax(p):]
		c = c[np.argmax(p):]
		p = p[np.argmax(p):]

	return b, c, p

def stochastic_model(alpha, K2, K4, k6, scale=1, scale_k=None):
	n_rxns = 600000 # Number of reactions to simulate
	X = 500 # Stable value around which species should fluctuate

	k_vals = np.array([X, K2*k6, K2*k6, K4*k6, K4*k6, k6/X, k6, X*k6, k6/2])*scale

	if scale_k != None: # Scale a rate constant kn if given some scale_k = [n, scale]
		k_vals[scale_k[0]] *= scale_k[1]

	A = np.zeros((n_rxns,5)) # 4 species, and keep track of the time at index 0
	A[0,1:] = np.array([0, 0, 0, 0]) #BMAL (B), CRY (C), PER (P), CRY-PER (D)

	# Draw random numbers to determine next reaction times and which reaction occurs at each time
	tau = -np.log(np.random.rand(n_rxns))
	next_rxn = np.random.rand(n_rxns)

	for i in range(0,n_rxns-1):
		# Relative probability/propensity of each possible reaction happening
		probs = np.array([k_vals[7]/(1+(A[i,4]/k_vals[0])**alpha),
							k_vals[1]*A[i,1],
							k_vals[2]*A[i,2],
							k_vals[3]*A[i,1],
							k_vals[4]*A[i,3],
							k_vals[5]*A[i,2]*A[i,3],
							k_vals[6]*A[i,4],
							k_vals[8]*A[i,1]])
		
		prob_sum = probs.cumsum()
		A[i+1, 0] = A[i, 0] + tau[i]/prob_sum[-1] # Fill in the time of the reaction
		A[i+1, 1:] = A[i, 1:] # Carry forward previous molecule counts at new time
		
		# Figure out which partition the random number falls in to determine which reaction happened
		match 8 - (prob_sum[-1]*next_rxn[i] < prob_sum).sum():
			case 0: # Reaction 1: +1 B
				A[i+1, 1] += 1
			case 1: # Reaction 2: +1 C
				A[i+1, 2] += 1
			case 2: # Reaction 3: -1 C
				A[i+1, 2] -= 1
			case 3: # Reaction 4: +1 P
				A[i+1, 3] += 1
			case 4: # Reaction 5: -1 P
				A[i+1, 3] -= 1
			case 5: # Reaction 6: C + P -> D
				A[i+1, 2] -= 1
				A[i+1, 3] -= 1
				A[i+1, 4] += 1
			case 6: # Reaction 7: D -> C + P
				A[i+1, 2] += 1
				A[i+1, 3] += 1
				A[i+1, 4] -= 1
			case 7: # Reaction 8: -1 B
				A[i+1, 1] -= 1

	return A


def get_periods(t_series, a):
	Ts = np.zeros(0)
	
	# Chop data to start at day 3
	t = 72
	a = a[np.sum(t_series<t):]
	t_series = t_series[np.sum(t_series<t):]

	while (t_series[-1] - t) > 26: # More than 26 hours of the series left to analyze
		# Find the maximum within the first 15 hours
		first_max = np.argmax(a[:np.sum(t_series<(t+15))])
		t_first_max = t_series[first_max]

		# Find the minimum within 20 hours after the maximum
		per_min = first_max + np.argmin(a[first_max:np.sum(t_series<(t_first_max+20))])
		t_min = t_series[per_min]
		
		# Find the next maximum within 15 hours of the minimum
		next_max = per_min + np.argmax(a[per_min:np.sum(t_series<(t_min+15))])
		t_next_max = t_series[next_max]

		# Period is the next max minus the first max
		Ts = np.append(Ts, t_next_max-t_first_max)
		
		# Now search for next period starting from the second max
		a = a[next_max:]
		t_series = t_series[next_max:]
		t = t_next_max
	return Ts







