import matplotlib.pyplot as plt
import numpy as np

n_rxns = 40000

alpha = 10
x = 1000 # Stability point
K2 = 1
K4 = 1.5

k_vals = np.array([1/(x*(1+x**alpha)), 1/x, K2, K2+1, K4, K4+1])

A = np.zeros((n_rxns,5)) # 4 species, and keep track of the time
A[0,1:] = np.array([950, 1100, 1050, 950]) #BMAL, CRY, PER, C-P

# Draw random numbers for next times and reactions
tau = -np.log(np.random.rand(n_rxns))
next_rxn = np.random.rand(n_rxns)

for i in range(0,n_rxns-1):
	probs = np.array([1/(1+A[i,4]**alpha),
						A[i,4],
						k_vals[0]*A[i,1],
						k_vals[1]*A[i,2]*A[i,3],
						k_vals[2]*A[i,1],
						k_vals[3]*A[i,2],
						k_vals[4]*A[i,1],
						k_vals[5]*A[i,3]])
	
	prob_sum = probs.cumsum()
	A[i+1, 0] = A[i, 0] + tau[i]/prob_sum[-1] # Fill the time of the reaction

	# Update molecule counts at new time
	A[i+1, 1:] = A[i, 1:]
	# Figure out which partition the random number falls in
	match 8 - (prob_sum[-1]*next_rxn[i] < prob_sum).sum():
		case 0:
			A[i+1, 1] += 1
		case 1:
			A[i+1, 2] += 1
			A[i+1, 3] += 1
			A[i+1, 4] -= 1
		case 2:
			A[i+1, 1] -= 1
		case 3:
			A[i+1, 4] += 1
		case 4:
			A[i+1, 2] += 1
		case 5:
			A[i+1, 2] -= 1
		case 6:
			A[i+1, 3] += 1
		case 7:
			A[i+1, 3] -= 1

plt.plot(A[:,0],A[:,1],label="B")
plt.plot(A[:,0],A[:,2],label="C")
plt.plot(A[:,0],A[:,3],label="P")
plt.plot(A[:,0],A[:,3],label="D")
plt.legend()
plt.show()










