import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.ticker as ticker

fig = plt.figure()
fig.suptitle("Dependence of Greeks on Various Parameters", fontsize=40)

# Adjusting x-axis ticks
x_steps = 20
x_ticks = [i for i in range(x_steps) if i%4==0]

# Files to load
filenames = ["Delta_spot", "Gamma_spot", "Kappa_spot", "Delta_ITM", "Delta_ATM", "Delta_OTM", "Kappa_time", "Kappa_vol"]

# x-values to load
x_values = ["S", "S", "S", "t", "t", "t", "t", "sigma"]

x_axis_names = ["S", "S", "S", "t", "t", "t", "t", "$\sigma$"]
y_axis_names = ["$\Delta$", "$\Gamma$", "K", "$\Delta$", "$\Delta$", "$\Delta$", "K", "K"]

# Create and configure each subplot, then plot data
fig.subplots_adjust(hspace=0.4, wspace=0.4)
for i in range(8):

    # Config
    ax = fig.add_subplot(4, 2, i+1)
    ax.set_xlabel(x_axis_names[i], fontsize = 40)
    ax.set_ylabel(y_axis_names[i], fontsize = 40)
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=30)
    #ax.set_xticks(x_ticks)

    # Load and plot each datafile
    data = np.loadtxt("./Data/" + filenames[i] + ".txt")
    x = np.loadtxt("./Data/" + x_values[i] + ".txt")
    ax.plot(x, data, lw=4)

fig2 = plt.figure()
fig2.suptitle("Convergence of Delta when Calculated Using Monte Carlo", fontsize=40)

filenames = ["Delta_1_vals", "Delta_2_vals", "likelihood_vals", "pathwise_vals"]
y_axis_names = ["$\Delta_1$", "$\Delta_2$", "$\Delta_3$", "$\Delta_4$"]

fig2.subplots_adjust(hspace=0.3, wspace=0.3)
for i in range(4):

    # Config
    ax = fig2.add_subplot(2, 2, i+1)
    ax.set_xlabel("Log$_2(N)$", fontsize=40, labelpad=30)
    ax.set_ylabel(y_axis_names[i], fontsize = 40, labelpad=50)
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=30)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=2))


    # Load and plot each datafile
    data = np.loadtxt("./Data/" + filenames[i] + ".txt")
    ax.plot(data, marker="o", markersize=10, lw=4)

    # Draw ellipse and lines on first subplot
    if(i==0):
        x1, y1 = [19.45, 18.42], [28, 335]
        x2, y2 = [10.51, 5.72], [29, 335]
        ax.plot(x1, y1, x2, y2, lw=2, color="red")
        ellipse = Ellipse(xy=(15, 10), width=9, height=150, edgecolor='r', fc='None', lw=2)
        ax.add_patch(ellipse)


# Add sub-sub-plot to subplot 1
# config
x = np.linspace(11,20,9)
ax = fig2.add_axes([0.23, 0.67, 0.2, 0.2])
ax.tick_params(axis='x', labelsize=25)
ax.tick_params(axis='y', labelsize=25)

# We only want the last few data points, so we load the file and truncate before plotting
data = np.loadtxt("./Data/Delta_1_vals.txt")
data = [data[i] for i in range(len(data)-9,len(data))]
ax.plot(x, data, marker="o", markersize=8, lw=3)

plt.show()