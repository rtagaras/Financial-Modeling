import numpy as np
import matplotlib.pyplot as plt

# configuring various figure elements
fig, ax = plt.subplots()
ax.set_title("Stock Price Correlated to Market Scenario")
ax.set(xlabel = "Days", ylabel = "Value")


market = np.loadtxt("./Data/market.txt")
stock = np.loadtxt("./Data/data.txt")

ax.plot(market)
ax.plot(stock)
plt.show()