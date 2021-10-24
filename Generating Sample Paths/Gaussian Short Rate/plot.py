import numpy as np
import matplotlib.pyplot as plt

# configuring various figure elements
fig, (ax1, ax2) = plt.subplots(1,2)
fig.suptitle("Short-Rate Models", fontsize = 50)

ax1.set_title("Ho-Lee Model", fontsize = 40)
ax2.set_title("Vasicek Model", fontsize = 40)

ax1.set_xlabel("Days", fontsize = 30)
ax1.set_ylabel("Value", fontsize = 30)
ax2.set_xlabel("Days", fontsize = 30)
ax2.set_ylabel("Value", fontsize = 30)

ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)

# Plot data
data = np.loadtxt("./Data/HL_data.txt")
ax1.plot(data, lw=3)

data = np.loadtxt("./Data/V_data.txt")
ax2.plot(data, lw=3)

plt.show()