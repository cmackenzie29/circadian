# import numpy as np
# import matplotlib.pyplot as plt

# N = 50000

# a = np.zeros(N)
# b = np.zeros(N)
# c = np.zeros(N)
# p = np.zeros(N)
# d = np.zeros(N)

# a[0] = 0.98
# b[0] = 0.9
# c[0] = 1.1
# p[0] = 1.05
# d[0] = 0.95

# dt = 1e-1
# k1 = 1
# k3 = 1
# k5 = 1
# k7 = 0.001
# k9 = 0.001
# k = np.array([0,k1,k1,k3,k3,k5,k5,k7,k7,k9,k9])

# for i in range(1,N):
# 	p[i] = p[i-1] + dt*(k[1]*b[i-1]+k[8]*d[i-1]-k[2]*p[i-1]-k[7]*c[i-1]*p[i-1])
# 	c[i] = c[i-1] + dt*(k[3]*b[i-1]+k[8]*d[i-1]-k[4]*c[i-1]-k[7]*c[i-1]*p[i-1])
# 	b[i] = b[i-1] + dt*(k[5]*d[i-1]+k[10]*a[i-1]-k[6]*b[i-1]-k[9]*b[i-1]*d[i-1])
# 	d[i] = d[i-1] + dt*(k[7]*c[i-1]*p[i-1]-k[8]*d[i-1]+k[10]*a[i-1]-k[9]*b[i-1]*d[i-1])
# 	a[i] = a[i-1] + dt*(k[9]*b[i-1]*d[i-1]-k[10]*a[i-1])

# plt.plot(np.arange(N),b)
# plt.plot(np.arange(N),c)
# plt.plot(np.arange(N),p)
# plt.plot(np.arange(N),d)
# plt.plot(np.arange(N),a)
# plt.show()

import numpy as np
import matplotlib.pyplot as plt

N = 25000

b = np.zeros(N)
c = np.zeros(N)
p = np.zeros(N)
d = np.zeros(N)

b[0] = 0.9
c[0] = 1.1
p[0] = 1.05
d[0] = 0.95

dt = 1e-3
alpha = 10
k2 = 1
k4 = 1.5
k = np.array([0.5, 1, k2, k2+1, k4, k4+1])

for i in range(1,N):
	b[i] = b[i-1] + dt*(1/(1+d[i-1]**alpha)-k[0]*b[i-1])
	c[i] = c[i-1] + dt*(k[2]*b[i-1]+d[i-1]-k[3]*c[i-1])
	p[i] = p[i-1] + dt*(k[4]*b[i-1]+d[i-1]-k[5]*p[i-1])
	d[i] = d[i-1] + dt*(k[1]*c[i-1]*p[i-1]-d[i-1])

plt.plot(np.arange(N),b,label="B")
plt.plot(np.arange(N),c,label="C")
plt.plot(np.arange(N),p,label="P")
plt.plot(np.arange(N),d,label="D")
plt.legend()
plt.show()


