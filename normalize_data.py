import numpy as np
import matplotlib.pyplot as plt
import os
import pickle

# Line up and normalize WT data from the various sources
# Data will be centered around 500 with range 300
# Data will be taken in 72 hour chunks, in 10 minute intervals.
# First chunk will start at the first PER1 peak for each WT result.
dir = ""
per = []
cry = []
bmal = []

# Ono
p1 = np.loadtxt(dir+"RawData/Ono_WT_PER1.dat")
p2 = np.loadtxt(dir+"RawData/Ono_WT_PER2.dat")
b = np.loadtxt(dir+"RawData/Ono_WT_BMAL1.dat")

T0 = p2[p2[:,1]==np.max(p2[:,1][p2[:,0]<20])][0][0] # First PER2 peak
p1[:,0] -= T0
p2[:,0] -= T0
b[:,0] -= T0

p1 = p1[p1[:,0]>=0]
p2 = p2[p2[:,0]>=0]
b = b[b[:,0]>=0]

# Detrend
p1_trend = np.polyfit(range(len(p1[:,1])),p1[:,1],1)
p1[:,1] -= (range(len(p1[:,1]))*p1_trend[0]+p1_trend[1])
p1[:,1] += 1
p2_trend = np.polyfit(range(len(p2[:,1])),p2[:,1],1)
p2[:,1] -= (range(len(p2[:,1]))*p2_trend[0]+p2_trend[1])
p2[:,1] += 1
b_trend = np.polyfit(range(len(b[:,1])),b[:,1],1)
b[:,1] -= (range(len(b[:,1]))*b_trend[0]+b_trend[1])
b[:,1] += 1

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


# Maywood WT
p = np.loadtxt(dir+"RawData/Maywood_WT_PER2.dat")
c = np.loadtxt(dir+"RawData/Maywood_WT_CRY1.dat")

T0 = p[p[:,1]==np.max(p[:,1][p[:,0]>25])][0][0] # Second real PER2 peak
p[:,0] -= T0
c[:,0] -= T0
p = p[p[:,0]>=0]
c = c[c[:,0]>=0]

# Detrend
p_trend = np.polyfit(range(len(p[:,1])),p[:,1],1)
p[:,1] -= (range(len(p[:,1]))*p_trend[0]+p_trend[1])
p[:,1] += 1
c_trend = np.polyfit(range(len(c[:,1])),c[:,1],1)
c[:,1] -= (range(len(c[:,1]))*c_trend[0]+c_trend[1])
c[:,1] += 1

# Rescale
p[:,1] += 5*p[:,1].ptp()/3-p[:,1].mean()
p[:,1] *= 500/p[:,1].mean()
c[:,1] += 5*c[:,1].ptp()/3-c[:,1].mean()
c[:,1] *= 500/c[:,1].mean()

per.append(p)
cry.append(c)


# Nishide WT
p = np.loadtxt(dir+"RawData/Nishide_WT_PER2.dat")
b = np.loadtxt(dir+"RawData/Nishide_WT_BMAL1.dat")

T0 = p[p[:,1]==np.max(p[:,1][p[:,0]>12])][0][0] # First real PER2 peak
p[:,0] -= T0
b[:,0] -= T0
p = p[p[:,0]>=0]
b = b[b[:,0]>=0]

# Detrend
p_trend = np.polyfit(range(len(p[:,1])),p[:,1],1)
p[:,1] -= (range(len(p[:,1]))*p_trend[0]+p_trend[1])
p[:,1] += 1
b_trend = np.polyfit(range(len(b[:,1])),b[:,1],1)
b[:,1] -= (range(len(b[:,1]))*b_trend[0]+b_trend[1])
b[:,1] += 1

# Rescale
p[:,1] += 5*p[:,1].ptp()/3-p[:,1].mean()
p[:,1] *= 500/p[:,1].mean()
b[:,1] += 5*b[:,1].ptp()/3-b[:,1].mean()
b[:,1] *= 500/b[:,1].mean()

per.append(p)
bmal.append(b)


# Noguchi WT
p = np.loadtxt(dir+"RawData/Noguchi_WT_PER2.dat")
b = np.loadtxt(dir+"RawData/Noguchi_WT_BMAL1.dat")

T0 = p[p[:,1]==np.max(p[:,1][p[:,0]<20])][0][0] # First PER2 peak
p[:,0] -= T0
b[:,0] -= T0
p = p[p[:,0]>=0]
b = b[b[:,0]>=0]

# Detrend
p_trend = np.polyfit(range(len(p[:,1])),p[:,1],1)
p[:,1] -= (range(len(p[:,1]))*p_trend[0]+p_trend[1])
p[:,1] += 1
b_trend = np.polyfit(range(len(b[:,1])),b[:,1],1)
b[:,1] -= (range(len(b[:,1]))*b_trend[0]+b_trend[1])
b[:,1] += 1

# Rescale
p[:,1] += 5*p[:,1].ptp()/3-p[:,1].mean()
p[:,1] *= 500/p[:,1].mean()
b[:,1] += 5*b[:,1].ptp()/3-b[:,1].mean()
b[:,1] *= 500/b[:,1].mean()

per.append(p)
bmal.append(b)


# for p in per:
# 	plt.plot(p[:,0],p[:,1])

# plt.show()

pickle.dump(per, open(dir+'per.sav', 'wb'))
pickle.dump(cry, open(dir+'cry.sav', 'wb'))
pickle.dump(bmal, open(dir+'bmal.sav', 'wb'))




