import numpy as np
import matplotlib.pyplot as plt
import pickle
import seaborn as sns

sns.set(context='paper',style='white',font_scale=1.5,rc={"lines.linewidth":2.5})
dir = "/Users/cameron/Desktop/circadian/Data/RawData/"
dir_2 = "/Users/cameron/Desktop/circadian/Data/"
colors = sns.color_palette()

# Figure 2.2
from deterministic import dimensionless_model
b, c, p, d = dimensionless_model(X=1, alpha=20, K2=1.0, K4=2.5)
plt.plot(b,label="b")
plt.plot(c,label="c")
plt.plot(p,label="p")
plt.plot(d,label="d")
plt.xlabel(r"$\tau$")
plt.legend()
plt.show()


# Figure 2.3 - REMOVED
# from scipy.linalg import eigvals
# N = 200
# bound = 10
# K_range = np.linspace(bound/N, bound, N)

# #plt.subplot(1,2,1)
# alpha = 5
# K_grid = np.zeros((N,N))
# for i in range(N):
# 	for j in range(N):
# 		evals = eigvals(np.array([[-0.5,0,-alpha/4,0],[K_range[i], -K_range[i]-1, 1, -1],[0, 1, -1, 1],[K_range[j], -1, 1, -K_range[j]-1]]))
# 		if np.all(evals.real <= 0) and (not np.all(evals.imag == 0)):
# 			K_grid[i,j] = 1
# plt.imshow(K_grid, origin="lower", extent=(0,bound,0,bound), cmap="Greys", alpha=0.5)
# plt.title(fr"  $\alpha={alpha}$", loc="left")
# plt.xlabel(r"$K_2$")
# plt.ylabel(r"$K_4$")
# plt.show()

# plt.subplot(1,2,2)
# alpha = 10
# K_grid = np.zeros((N,N))
# for i in range(N):
# 	for j in range(N):
# 		evals = eigvals(np.array([[-0.5,0,-alpha/4,0],[K_range[i], -K_range[i]-1, 1, -1],[0, 1, -1, 1],[K_range[j], -1, 1, -K_range[j]-1]]))
# 		if np.all(evals.real <= 0) and (not np.all(evals.imag == 0)):
# 			K_grid[i,j] = 1
# plt.imshow(K_grid, origin="lower", extent=(0,bound,0,bound), cmap="Greys", alpha=0.5)
# plt.title(r"  $\alpha=10$", loc="left")
# plt.xlabel(r"$K_2$")
# plt.yticks([])
# plt.show()



# Figure 3.1
plt.subplot(2,2,2)
p1 = np.loadtxt(dir+"Ono_WT_PER1.dat")
plt.scatter(p1[:,0],p1[:,1])
plt.ylabel("Relative intensity")
plt.xlim(0,100)
plt.ylim(0,1)

plt.subplot(2,2,3)
p2 = np.loadtxt(dir+"Ono_WT_PER2.dat")
plt.scatter(p2[:,0],p2[:,1],c="red")
plt.xlabel("Time (hrs)")
plt.ylabel("Relative intensity")
plt.xlim(0,100)
plt.ylim(0,1)

plt.subplot(2,2,4)
b = np.loadtxt(dir+"Ono_WT_BMAL1.dat")
plt.scatter(b[:,0],b[:,1],c="green")
plt.xlabel("Time (hrs)")
plt.ylabel("Relative intensity")
plt.xlim(0,100)
plt.ylim(0,1)

plt.tight_layout()
plt.show()



# Figure 3.2
from simulate import deterministic_model, stochastic_model

alpha = 76.3
K2 = 3.45
K4 = 2.88
k6 = 0.110
dt = 1e-1 # Timestep

b, c, p = deterministic_model(alpha, K2, K4, k6)
A = np.loadtxt(dir+"stoch_sim.dat") # Load a previous stochastic simulation run

bmal_data = pickle.load(open(dir_2+'bmal.sav', 'rb'))[:3]
cry_data = pickle.load(open(dir_2+'cry.sav', 'rb'))[0]
per_data = pickle.load(open(dir_2+'per.sav', 'rb'))[3]

plt.subplot(2,3,1)
plt.title("A", loc="left")
plt.plot(dt*np.arange(0,1000),b[150:1150],color=colors[0])
tf_1 = (bmal_data[0][:,0] >= 15) & (bmal_data[0][:,0] <= 115)
tf_2 = (bmal_data[1][:,0] >= 15) & (bmal_data[1][:,0] <= 115)
tf_3 = (bmal_data[2][:,0] >= 15) & (bmal_data[2][:,0] <= 115)
plt.scatter(bmal_data[0][tf_1,0]-bmal_data[0][tf_1,0][0], bmal_data[0][tf_1,1], color=colors[3], s=6)
plt.scatter(bmal_data[1][tf_2,0]-bmal_data[1][tf_2,0][0], bmal_data[1][tf_2,1], color=colors[8], s=6)
plt.scatter(bmal_data[2][tf_3,0]-bmal_data[2][tf_3,0][0], bmal_data[2][tf_3,1], color=colors[6], s=6)
plt.xticks([])
plt.yticks([400,500,600])
plt.ylabel("Molec. Counts")

plt.subplot(2,3,4)
plt.title("D", loc="left")
timeframe = (A[:,0] >= 347) & (A[:,0] <= 447)
plt.plot(A[timeframe,0]-A[timeframe,0][0], A[timeframe,1],color=colors[0],linewidth=1.5)
plt.scatter(bmal_data[0][tf_1,0]-bmal_data[0][tf_1,0][0], bmal_data[0][tf_1,1], color=colors[3], s=6)
plt.scatter(bmal_data[1][tf_2,0]-bmal_data[1][tf_2,0][0], bmal_data[1][tf_2,1], color=colors[8], s=6)
plt.scatter(bmal_data[2][tf_3,0]-bmal_data[2][tf_3,0][0], bmal_data[2][tf_3,1], color=colors[6], s=6)
plt.yticks([400,500,600])
plt.ylabel("Molec. Counts")
plt.xlabel("Time (hrs)")

plt.subplot(2,3,2)
plt.title("B", loc="left")
plt.plot(dt*np.arange(0,1000),b[200:1200],color=colors[1])
tf = (cry_data[:,0] >= 20) & (cry_data[:,0] <= 120)
plt.scatter(cry_data[tf,0]-cry_data[tf,0][0], cry_data[tf,1], color=colors[9], s=6)
plt.xticks([])
plt.yticks([400,500,600])

plt.subplot(2,3,5)
plt.title("E", loc="left")
timeframe = (A[:,0] >= 76) & (A[:,0] <= 176)
plt.plot(A[timeframe,0]-A[timeframe,0][0], A[timeframe,2],color=colors[1],linewidth=1.5)
plt.scatter(cry_data[tf,0]-cry_data[tf,0][0], cry_data[tf,1], color=colors[9], s=6)
plt.yticks([400,500,600])
plt.xlabel("Time (hrs)")

plt.subplot(2,3,3)
plt.title("C", loc="left")
plt.plot(dt*np.arange(0,1000),p[1750:2750],label="P",color=colors[2])
tf = (per_data[:,0] >= 175) & (per_data[:,0] <= 275)
plt.scatter(per_data[tf,0]-per_data[tf,0][0], per_data[tf,1], color=colors[7], s=6)
plt.xticks([])
plt.yticks([400,500,600])

plt.subplot(2,3,6)
plt.title("F", loc="left")
timeframe = (A[:,0] >= 108) & (A[:,0] <= 208)
plt.plot(A[timeframe,0]-A[timeframe,0][0], A[timeframe,3],color=colors[2],linewidth=1.5)
plt.scatter(per_data[tf,0]-per_data[tf,0][0], per_data[tf,1], color=colors[7], s=6)
plt.yticks([400,500,600])
plt.xlabel("Time (hrs)")

plt.tight_layout()
plt.show()


# Figure 3.3
b, c, p = deterministic_model(alpha, K2, K4, k6)
A = np.loadtxt(dir+"stoch_sim.dat") # Load a previous stochastic simulation run

plt.subplot(2,2,1)
plt.title("A", loc="left")
plt.plot(dt*np.arange(0,1000),b[1200:2200],label="B")
plt.plot(dt*np.arange(0,1000),c[1200:2200],label="C")
plt.plot(dt*np.arange(0,1000),p[1200:2200],label="P")
plt.legend()
plt.yticks([400,600])
plt.xlabel("Time (hrs)")
plt.ylabel("Molec. Counts")

plt.subplot(2,2,2)
plt.title("B", loc="left")
timeframe = (A[:,0] >= 150) & (A[:,0] <= 250)
plt.plot(A[timeframe,0]-A[timeframe,0][0], A[timeframe,1],linewidth=1.5)
plt.plot(A[timeframe,0]-A[timeframe,0][0], A[timeframe,2],linewidth=1.5)
plt.plot(A[timeframe,0]-A[timeframe,0][0], A[timeframe,3],linewidth=1.5)
plt.yticks([400,600])
plt.xlabel("Time (hrs)")

plt.subplot(2,2,3, projection='3d')
plt.title("C", loc="left")
plt.plot(b, c, p, lw=1, c="black")
plt.xlabel("B")
plt.ylabel("C")

plt.subplot(2,2,4, projection= '3d')
plt.title("D", loc="left")
plt.plot(A[8000:,1], A[8000:,2], A[8000:,3], lw=0.25, c="black")
plt.xlabel("B")
plt.ylabel("C")

plt.tight_layout()
plt.show()



# Figure 3.4
from simulate import get_periods

np.random.seed(0)
Ts = np.zeros(0)
for _ in range(10): # Run 10 simulations to collect periods
	stoch = stochastic_model(alpha,K2,K4,k6)
	Ts = np.append(Ts, get_periods(stoch[:,0], stoch[:,3]))
print(Ts.mean())
print(Ts.std())
print(len(Ts))

fig, axs = plt.subplots(1, 2)
axs[0].set_title("A", loc="left")
axs[0].plot(stoch[80000:180000,0]-stoch[80000,0],stoch[80000:180000,1],label="B",linewidth=1.5)
axs[0].plot(stoch[80000:180000,0]-stoch[80000,0],stoch[80000:180000,2],label="C",linewidth=1.5)
axs[0].plot(stoch[80000:180000,0]-stoch[80000,0],stoch[80000:180000,3],label="P",linewidth=1.5)
axs[0].legend()
axs[0].set_ylabel("Molec. counts")
axs[0].set_xlabel("Time (hrs)")
axs[1].set_title("B", loc="left")
sns.stripplot([Ts], ax=axs[1], color="None", edgecolor="black", linewidth=1)
sns.barplot([Ts], errorbar="sd", color = "grey", edgecolor="black", capsize=0.4, ax=axs[1], width=0.6)
axs[1].set_xticks([])
axs[1].set_ylabel("Period (hrs)")
fig.tight_layout()
plt.show()



# Figure 4.1
plt.subplot(3,5,1)
plt.title("A      BMAL1 Het. KO", loc="left")
park_het_per = np.loadtxt(dir+"Park_BMAL1HetKO_PER2.dat")
plt.scatter(park_het_per[:,0], park_het_per[:,1], color = colors[2], s=6)
plt.xticks([])
plt.ylabel("Relative intensity")

plt.subplot(3,5,2)
plt.title("B      BMAL1 Hom. KO", loc="left")
park_hom_per = np.loadtxt(dir+"Park_BMAL1HomKO_PER2.dat")
plt.scatter(park_hom_per[:,0], park_hom_per[:,1], color = colors[2], s=6)
plt.xticks([])

plt.subplot(3,5,3)
plt.title("C          PER1 KO", loc="left")
per1_ko = np.loadtxt(dir+"Maywood_PER1KO_CRY2.dat")
plt.scatter(per1_ko[:,0], per1_ko[:,1], color = colors[1], s=6)
plt.xticks([])

plt.subplot(3,5,4)
plt.title("D         PER1/2 KO", loc="left")
per12_ko = np.loadtxt(dir+"Maywood_PER12KO_CRY2.dat")
plt.scatter(per12_ko[:,0], per12_ko[:,1], color = colors[1], s=6)
plt.xticks([])

plt.subplot(3,5,5)
plt.title("E         CRY1/2 KO", loc="left")
cry_ko = np.loadtxt(dir+"Honma_CRY12KO_PER2.dat")
plt.scatter(cry_ko[:,0], cry_ko[:,1], color = colors[2], s=6)
plt.xticks([])

plt.subplot(3,5,6)
plt.title("F", loc="left")
sim = deterministic_model(alpha,K2,K4,k6,scale_k=[7,0.5])
plt.plot(sim[0][900:900+round(park_het_per[-1,0]/dt)])
plt.plot(sim[1][900:900+round(park_het_per[-1,0]/dt)])
plt.plot(sim[2][900:900+round(park_het_per[-1,0]/dt)])
plt.xticks([])
plt.ylabel("Molec. counts")

plt.subplot(3,5,7)
plt.title("G", loc="left")
sim = deterministic_model(alpha,K2,K4,k6,scale_k=[7,0.2])
plt.plot(sim[0][:round(park_hom_per[-1,0]/dt)])
plt.plot(sim[1][:round(park_hom_per[-1,0]/dt)])
plt.plot(sim[2][:round(park_hom_per[-1,0]/dt)])
plt.xticks([])

plt.subplot(3,5,8)
plt.title("H", loc="left")
sim = deterministic_model(alpha,K2,K4,k6,scale_k=[3,0.5])
plt.plot(sim[0][:round(per1_ko[-1,0]/dt)])
plt.plot(sim[1][:round(per1_ko[-1,0]/dt)])
plt.plot(sim[2][:round(per1_ko[-1,0]/dt)])
plt.xticks([])
plt.yticks([0,500])

plt.subplot(3,5,9)
plt.title("I", loc="left")
sim = deterministic_model(alpha,K2,K4,k6,scale_k=[3,0.2])
plt.plot(sim[0][:round(per12_ko[-1,0]/dt)])
plt.plot(sim[1][:round(per12_ko[-1,0]/dt)])
plt.plot(sim[2][:round(per12_ko[-1,0]/dt)])
plt.xticks([])
plt.yticks([0,500])

plt.subplot(3,5,10)
plt.title("J", loc="left")
sim = deterministic_model(alpha,K2,K4,k6,scale_k=[1,0.2])
plt.plot(sim[0][:round(cry_ko[-1,0]/dt)])
plt.plot(sim[1][:round(cry_ko[-1,0]/dt)])
plt.plot(sim[2][:round(cry_ko[-1,0]/dt)])
plt.xticks([])
plt.yticks([0,500])

plt.subplot(3,5,11)
plt.title("K", loc="left")
np.random.seed(1)
sim = stochastic_model(alpha,K2,K4,k6,scale_k=[7,0.5])
tf = (sim[:,0] > 90) & (sim[:,0] < 90 + park_het_per[-1,0])
plt.plot(sim[tf,0]-sim[tf,0][0], sim[tf,1], linewidth=1.5)
plt.plot(sim[tf,0]-sim[tf,0][0], sim[tf,2], linewidth=1.5)
plt.plot(sim[tf,0]-sim[tf,0][0], sim[tf,3], linewidth=1.5)
plt.xlabel("Time (hrs)")
plt.ylabel("Molec. counts")

plt.subplot(3,5,12)
plt.title("L", loc="left")
sim = stochastic_model(alpha,K2,K4,k6,scale_k=[7,0.2])
tf = (sim[:,0] < park_hom_per[-1,0])
plt.plot(sim[tf,0]-sim[tf,0][0], sim[tf,1], linewidth=1.5)
plt.plot(sim[tf,0]-sim[tf,0][0], sim[tf,2], linewidth=1.5)
plt.plot(sim[tf,0]-sim[tf,0][0], sim[tf,3], linewidth=1.5)
plt.xlabel("Time (hrs)")

plt.subplot(3,5,13)
sim = stochastic_model(alpha,K2,K4,k6,scale_k=[3,0.5])
tf = (sim[:,0] < per1_ko[-1,0])
plt.plot(sim[tf,0]-sim[tf,0][0], sim[tf,1], linewidth=1.5)
plt.plot(sim[tf,0]-sim[tf,0][0], sim[tf,2], linewidth=1.5)
plt.plot(sim[tf,0]-sim[tf,0][0], sim[tf,3], linewidth=1.5)
plt.title("M", loc="left")
plt.xlabel("Time (hrs)")

plt.subplot(3,5,14)
sim = stochastic_model(alpha,K2,K4,k6,scale_k=[3,0.2])
tf = (sim[:,0] < per12_ko[-1,0])
plt.plot(sim[tf,0]-sim[tf,0][0], sim[tf,1], linewidth=1.5)
plt.plot(sim[tf,0]-sim[tf,0][0], sim[tf,2], linewidth=1.5)
plt.plot(sim[tf,0]-sim[tf,0][0], sim[tf,3], linewidth=1.5)
plt.title("N", loc="left")
plt.xlabel("Time (hrs)")
plt.yticks([0,500])

plt.subplot(3,5,15)
sim = stochastic_model(alpha,K2,K4,k6,scale_k=[1,0.2])
tf = (sim[:,0] < cry_ko[-1,0])
plt.plot(sim[tf,0]-sim[tf,0][0], sim[tf,1], linewidth=1.5)
plt.plot(sim[tf,0]-sim[tf,0][0], sim[tf,2], linewidth=1.5)
plt.plot(sim[tf,0]-sim[tf,0][0], sim[tf,3], linewidth=1.5)
plt.title("O", loc="left")
plt.xlabel("Time (hrs)")
plt.yticks([0,500])

plt.show()



# Figure 4.2
plt.subplot(2,2,1)
plt.title("A", loc="left")
scale = np.arange(0.5,1.5,0.02)
periods = np.zeros_like(scale)
for i, s in enumerate(scale):
	det = deterministic_model(alpha, K2, K4, k6, scale=s)
	periods[i] = get_periods(dt*np.arange(len(det[0])), det[0])[-1]
plt.scatter(scale,periods,color="black")
plt.xlabel("Scale Factor")
plt.ylabel("Period (hrs)")

np.random.seed(1)
scale = np.arange(0.5,1.5,0.05)
scales = np.zeros(0)
periods = np.zeros(0)
for s in scale:
	print(s)
	stoch = stochastic_model(alpha, K2, K4, k6, scale=s)
	pers = get_periods(stoch[:,0], stoch[:,1])
	for p in pers:
		scales = np.append(scales, s)
		periods = np.append(periods, p)
	# General sample plots
	if np.isclose(s, 0.8):
		plt.subplot(2,2,3)
		plt.title("C          scale: 0.8", loc="left")
		tf = (stoch[:,0] > 500) & (stoch[:,0] < 600)
		plt.plot(stoch[tf,0]-stoch[tf,0][0], stoch[tf,1], label="B", linewidth=1.5)
		plt.plot(stoch[tf,0]-stoch[tf,0][0], stoch[tf,2], label="C", linewidth=1.5)
		plt.plot(stoch[tf,0]-stoch[tf,0][0], stoch[tf,3], label="P", linewidth=1.5)
		plt.xlabel("Time (hrs)")
		plt.ylabel("Molec. counts")
	elif np.isclose(s, 1.2):
		plt.subplot(2,2,4)
		plt.title("D          scale: 1.2", loc="left")
		tf = (stoch[:,0] > 250) & (stoch[:,0] < 350)
		plt.plot(stoch[tf,0]-stoch[tf,0][0], stoch[tf,1], linewidth=1.5)
		plt.plot(stoch[tf,0]-stoch[tf,0][0], stoch[tf,2], linewidth=1.5)
		plt.plot(stoch[tf,0]-stoch[tf,0][0], stoch[tf,3], linewidth=1.5)
		plt.xlabel("Time (hrs)")
		plt.ylabel("Molec. counts")

plt.subplot(2,2,2)
plt.title("B", loc="left")
plt.scatter(scales,periods,color="black")
plt.xlabel("Scale Factor")
plt.ylabel("Period (hrs)")

plt.tight_layout()
plt.show()





