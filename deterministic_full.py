import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig

N = 2500
X = 500

alpha = 45.7
K2 = 10.2
K4 = 2.5
k6 = 0.1
scale = 1 # For testing how scaling up/down all factors affects the period
k = [scale*X, scale*K2*k6, scale*(K2+1)*k6, scale*K4*k6, scale*(K4+1)*k6, scale*k6/X, scale*k6, scale*X*k6, scale*k6/2]
	
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
	c[i] = c[i-1] + dt*(k[1]*b[i-1]+k[6]*d[i-1]-k[2]*c[i-1])
	p[i] = p[i-1] + dt*(k[3]*b[i-1]+k[6]*d[i-1]-k[4]*p[i-1])
	d[i] = d[i-1] + dt*(k[5]*c[i-1]*p[i-1]-k[6]*d[i-1])

# Start at first per peak
# b = b[np.argmax(p):]
# c = c[np.argmax(p):]
# p = p[np.argmax(p):]

plt.plot(dt*np.arange(len(b)),b,label="B")
plt.plot(dt*np.arange(len(b)),c,label="C")
plt.plot(dt*np.arange(len(b)),p,label="P")
#plt.plot(dt*np.arange(len(b)),d,label="D")
plt.show()


