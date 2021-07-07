import numpy as np
import matplotlib.pyplot as plt

# configuring various figure elements
fig, ax1 = plt.subplots()
fig.suptitle("Geometric Brownian Motion")

ax1.set_title("1000 GBM Paths")
ax1.set(xlabel = "Days", ylabel = "Value")

# Adjusting units of x-axis ticks
x1_min = 0
x1_max = 365
x1_steps = 3651
x1_ticks = [i/(x1_steps/x1_max) for i in range(x1_steps)]


data = np.loadtxt("data.txt")
ax1.plot(data)

plt.show()