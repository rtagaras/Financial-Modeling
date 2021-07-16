import numpy as np
import matplotlib.pyplot as plt

# configuring various figure elements
fig, ax = plt.subplots()
ax.set_title("Stock Price Correlated to Market Scenario")
ax.set(xlabel = "Days", ylabel = "Value")


market = np.loadtxt("./Data/market_scenario.txt")
stock1 = np.loadtxt("./Data/security_1.txt")
stock2 = np.loadtxt("./Data/security_2.txt")


ax.plot(market, label= "Market Scenario")
ax.plot(stock1, label= "Security 1")
ax.plot(stock2, label= "Security 2")
ax.legend()

plt.show()