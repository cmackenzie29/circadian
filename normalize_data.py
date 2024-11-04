import numpy as np
import matplotlib.pyplot as plt
import os
import pickle

# Line up and normalize WT data from the various sources
# Data will be centered around 500 with range 300
# Data will be taken in 72 hour chunks, in 10 minute intervals.
# First chunk will start at the first PER1 peak for each WT result.
dir = "/Users/cameron/Desktop/circadian/Data/RawData/"
per = []
cry = []
bmal = []

# Ono
p1 = np.loadtxt(dir+"Ono_WT_PER1.dat")
p2 = np.loadtxt(dir+"Ono_WT_PER2.dat")
b = np.loadtxt(dir+"Ono_WT_BMAL1.dat")

T0 = p2[p2[:,1]==np.max(p2[:,1][p2[:,0]<20])][0][0] # First PER2 peak
p1[:,0] -= T0
p2[:,0] -= T0
b[:,0] -= T0

# 72 hour chunks
p1 = p1[(p1[:,0]>=0)&(p1[:,0]<72)]
p2 = p2[(p2[:,0]>=0)&(p2[:,0]<72)]
b = b[(b[:,0]>=0)&(b[:,0]<72)]

# Rescale
p1[:,1] += 5*p1[:,1].ptp()/3-p1[:,1].mean()
p1[:,1] *= 500/p1[:,1].mean()
p2[:,1] += 5*p2[:,1].ptp()/3-p2[:,1].mean()
p2[:,1] *= 500/p2[:,1].mean()
b[:,1] += 5*b[:,1].ptp()/3-b[:,1].mean()
b[:,1] *= 500/b[:,1].mean()

per.append(p1)
per.append(p2)
bmal.append(b)

for p in per:
	plt.plot(p[:,0],p[:,1])
plt.show()


# Maywood WT
p = np.loadtxt(dir+"Maywood_WT_PER2.dat")
c = np.loadtxt(dir+"Maywood_WT_CRY1.dat")

T0 = p[p[:,1]==np.max(p[:,1][p[:,0]>4])][0][0] # First real PER2 peak
p[:,0] -= T0
c[:,0] -= T0

# 72 hour chunks
for i in range(3):
	p_i = p[(p[:,0]>=72*i)&(p[:,0]<72*(i+1))]
	c_i = c[(c[:,0]>=72*i)&(c[:,0]<72*(i+1))]

	# Rescale
	p_i[:,1] += 5*p_i[:,1].ptp()/3-p_i[:,1].mean()
	p_i[:,1] *= 500/p_i[:,1].mean()
	c_i[:,1] += 5*c_i[:,1].ptp()/3-c_i[:,1].mean()
	c_i[:,1] *= 500/c_i[:,1].mean()

	per.append(p_i)
	cry.append(c_i)

for p in per:
	plt.plot(p[:,0],p[:,1])
plt.show()

# Nishide WT
p = np.loadtxt(dir+"Nishide_WT_PER2.dat")
b = np.loadtxt(dir+"Nishide_WT_BMAL1.dat")

T0 = p[p[:,1]==np.max(p[:,1][p[:,0]>12])][0][0] # First real PER2 peak
p[:,0] -= T0
b[:,0] -= T0

# 72 hour chunks
for i in range(5):
	p_i = p[(p[:,0]>=72*i)&(p[:,0]<72*(i+1))]
	b_i = b[(b[:,0]>=72*i)&(b[:,0]<72*(i+1))]

	# Rescale
	p_i[:,1] += 5*p_i[:,1].ptp()/3-p_i[:,1].mean()
	p_i[:,1] *= 500/p_i[:,1].mean()
	b_i[:,1] += 5*b_i[:,1].ptp()/3-b_i[:,1].mean()
	b_i[:,1] *= 500/b_i[:,1].mean()

	per.append(p_i)
	bmal.append(b_i)

for p in per:
	plt.plot(p[:,0],p[:,1])
plt.show()

# Noguchi WT
p = np.loadtxt(dir+"Noguchi_WT_PER2.dat")
b = np.loadtxt(dir+"Noguchi_WT_BMAL1.dat")

T0 = p[p[:,1]==np.max(p[:,1][p[:,0]<20])][0][0] # First PER2 peak
p[:,0] -= T0
b[:,0] -= T0

# 72 hour chunks
for i in range(2):
	p_i = p[(p[:,0]>=72*i)&(p[:,0]<72*(i+1))]
	b_i = b[(b[:,0]>=72*i)&(b[:,0]<72*(i+1))]

	# Rescale
	p_i[:,1] += 5*p_i[:,1].ptp()/3-p_i[:,1].mean()
	p_i[:,1] *= 500/p_i[:,1].mean()
	b_i[:,1] += 5*b_i[:,1].ptp()/3-b_i[:,1].mean()
	b_i[:,1] *= 500/b_i[:,1].mean()

	per.append(p_i)
	bmal.append(b_i)

#Binning
# t = np.linspace(0,71.9,720)
# y = np.zeros(720)
# for i in range(720):
# 	a = 0
# 	n = 0
# 	for p in per:
# 		pts = p[(p[:,0] >= t[i]) & (p[:,0] < t[i]+0.1)]

for p in per:
	plt.plot(p[:,0],p[:,1])
plt.show()






