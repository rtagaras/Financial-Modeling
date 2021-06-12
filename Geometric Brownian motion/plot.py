import numpy as np
import matplotlib.pyplot as plt

# configuring various figure elements
fig, (ax1, ax2) = plt.subplots(1,2)
fig.suptitle("Geometric Brownian Motion")

ax1.set_title("1000 GBM Paths")
ax2.set_title("Distribution of Path Endpoints")

ax1.set(xlabel = "Days", ylabel = "Value")
ax2.set(xlabel = "End Value", ylabel = "Probability Density")


# Adjusting units of x-axis ticks
x1_min = 0
x1_max = 60
x1_steps = 600
x1_ticks = [i/(x1_steps/x1_max) for i in range(x1_steps)]

# plot the first 10000 paths only
for i in range(10000):
    filename = "./Data/data" + str(i) + ".txt"
    data = np.loadtxt(filename)
    ax1.plot(x1_ticks, data)

# Plot histogram of endpoints
data = np.loadtxt("./Data/end_values.txt")

xmin = min(data)
xmax = max(data)
xmin = 0
xmax = 300
num_bins = 10000
bin_step = (xmax-xmin)/num_bins
bins = [xmin + bin_step*i for i in range(num_bins)]

ax2.hist(data, bins, weights=np.ones(len(data))/len(data))

plt.show()