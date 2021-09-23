import numpy as np
import matplotlib.pyplot as plt

#Plot Investigation 3
fig, ax = plt.subplots()
fig.suptitle("Ratio of Price to Volatility", fontsize=40)

# Adjusting x-axis ticks
x_steps = 10000
x_ticks = [i/10 for i in range(x_steps) if i%100==0]

axis_names = ["30", "50", "80", "100", "120", "150"]


# Config
ax.set_xlabel("$\sigma$", fontsize = 40)
ax.set_ylabel("S/$\sigma$", fontsize = 40)
ax.tick_params(axis='x', labelsize=30)
ax.tick_params(axis='y', labelsize=30)
ax.set_xticks(x_ticks)


# Load and plot each datafile
x = np.loadtxt("./Data/sigma.txt")

for i in range(len(axis_names)):
    data = np.loadtxt("./Data/ATM_ratios_" + axis_names[i] + ".txt")
    ax.plot(x, data, lw=4, label="$S_0$ = " + axis_names[i])

ax.legend(prop={'size': 30}, loc="upper right")

# Plot Investigation 4
fig2, ax2 = plt.subplots()

ax2.set_xlabel("S", fontsize = 40)
ax2.set_ylabel("Value", fontsize = 40)
ax2.tick_params(axis='x', labelsize=30)
ax2.tick_params(axis='y', labelsize=30)

data = np.loadtxt("./Data/put_prices.txt")
ax2.plot(data, lw=4, label="Option Value")

data = np.loadtxt("./Data/intrinsic_values.txt")
ax2.plot(data, lw=4, label="Intrinsic Value")

ax2.legend(prop={'size': 30}, loc="upper right")

plt.show()