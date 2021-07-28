from scipy.stats import t
import numpy as np

def batch_samples(samples, num_batches):
    batch_size = int(len(samples)/num_batches)
    batches = []

    if(len(samples) % num_batches != 0):
        print("Batch size does not divide size of sample set")

    for i in range(num_batches):
        b = []
        for j in range(batch_size*i, batch_size*(i+1)):
            b.append(samples[j])

        batches.append(b)

    return batches

def sample_variance(samples, num_batches):
    """
    return s_B^2 and theta
    """

    batches = batch_samples(samples, num_batches)
    averages = [sum(i)/len(i) for i in batches]
    
    theta = sum(averages)/num_batches
    s = 0

    for a in averages:
        s += (a-theta)**2 / (num_batches - 1)

    return s, theta

def confidence_interval(samples, percentage, num_batches):
    
    s, theta = sample_variance(samples, num_batches)
    t_a = t.ppf((1+percentage)/2, len(samples)-1)

    max = theta + t_a*np.sqrt(s/num_batches)
    min = theta - t_a*np.sqrt(s/num_batches)

    return min, max

option_types = ["European_GRW", "European_JD", "European_lattice", "Asian_GRW", "Barrier_GRW", "Basket_GRW", "Exchange_GRW"]
for n in option_types:
    filename = "./Data/" + n + ".txt"
    data = np.loadtxt(filename)

    mean = sum(data)/len(data)
    percentage = 0.95

    min, max = confidence_interval(data, percentage, len(data))
    print(f"Mean {n} value: {mean} \nWith {percentage*100}% certainty, the option value is in the interval ({min}, {max}). \n")