import matplotlib.pyplot as plt
import numpy as np

n_rxns = 100000
k_vals = np.array([0.002, 1, 1, 1, 1, 1])

A = np.zeros((n_rxns,4)) # 3 species, and keep track of the time
A[0,1:] = np.array([500,500,500]) #X, Y

# Draw random numbers for next times and reactions
tau = -np.log(np.random.rand(n_rxns))
next_rxn = np.random.rand(n_rxns)

for i in range(0,n_rxns-1):
	probs = np.array([k_vals[0]*A[i,2]*A[i,3],
						k_vals[1]*A[i,1],
						k_vals[2]*A[i,3],
						k_vals[3]*A[i,2],
						k_vals[4]*A[i,2],
						k_vals[5]*A[i,3]])

	prob_sum = probs.cumsum()
	A[i+1, 0] = A[i, 0] + tau[i]/prob_sum[-1] # Fill the time of the reaction

	# Update molecule counts at new time
	A[i+1, 1:] = A[i, 1:]
	# Figure out which partition the random number falls in
	match len(k_vals) - (prob_sum[-1]*next_rxn[i] < prob_sum).sum():
		case 0:
			A[i+1, 1] += 1
		case 1:
			A[i+1, 1] -= 1
		case 2:
			A[i+1, 2] += 1
		case 3:
			A[i+1, 2] -= 1
		case 4:
			A[i+1, 3] += 1
		case 5:
			A[i+1, 3] -= 1

# plt.plot(A[:,1],A[:,2])
# plt.plot(A[:,1],A[:,3])
# plt.plot(A[:,2],A[:,3])
plt.plot(A[:,0],A[:,1],label="A")
plt.plot(A[:,0],A[:,2],label="B")
plt.plot(A[:,0],A[:,3],label="C")
# plt.plot(A[:,0],(A[:,1]+A[:,2]+A[:,3])/3,label="Total")
plt.legend()
plt.show()










