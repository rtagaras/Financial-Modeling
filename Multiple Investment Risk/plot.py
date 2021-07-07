import numpy as np
import matplotlib.pyplot as plt

# configuring various figure elements
#fig, (ax1, ax2, ax3) = plt.subplots(1,3)
#fig.suptitle("Market Scenario Data")

fig, ax3 = plt.subplots()

# ax1.set_title("Interpolation")
# ax1.set(xlabel = "Days", ylabel = "Value")

# ax2.set_title("GRW")
# ax2.set(xlabel = "Days", ylabel = "Value")

ax3.set_title("Sum")
ax3.set(xlabel = "Days", ylabel = "Value")


# interpolation = np.loadtxt("./Data/interpolation.txt")
# GRW = np.loadtxt("./Data/GRW.txt")
for i in range(10):
    filename = "./Data/sum" + str(i) + ".txt"
    data = np.loadtxt(filename)
    ax3.plot(data)

# ax1.plot(interpolation)
# ax2.plot(GRW)
#ax3.plot(s)

plt.show()