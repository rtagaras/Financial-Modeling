import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import lognorm

# configuring various figure elements
fig, (ax1, ax2) = plt.subplots(1,2)
#fig, ax2 = plt.subplots()
fig.suptitle("Jump-Diffusion Processes")

ax1.set_title("1000 Paths")
ax2.set_title("Distribution of Path Endpoints")

ax1.set(xlabel = "Days", ylabel = "Value")
ax2.set(xlabel = "End Value", ylabel = "Probability Density")

# plot the first 1000 paths only
for i in range(100000):
    if(i % 100 == 0):
        dataname = "./Data/data_" + str(i) + ".txt"
        timename = "./Data/times_" + str(i) + ".txt"
        data = np.loadtxt(dataname)
        times = np.loadtxt(timename)
        ax1.plot(times, data)

# Plot histogram of endpoints
data = np.loadtxt("./Data/endpoints.txt")

xmin = 0
xmax = 260
num_bins = 1000
bin_step = (xmax-xmin)/num_bins
bins = [xmin + bin_step*i for i in range(num_bins)]

ax2.hist(data, bins, weights=np.ones(len(data))/len(data), density = True)

# Plot exact distribution on histogram
s_0 = 100
mu = 0.03
sig = 0.4
T = 100/365
a = np.log(s_0)+mu*T-0.5*T*sig**2
b = sig*np.sqrt(T)

x = np.linspace(40,260, 100)
ax2.plot(x, lognorm.pdf(x, b, scale = np.exp(a)), lw=2)

plt.show()