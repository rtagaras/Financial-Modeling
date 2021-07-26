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

    print(t_a)

    return min, max

x = [45, 55, 67, 45, 68, 79, 98, 87, 84, 82]
percentage = 0.999
min, max = confidence_interval(x, percentage, len(x))
print(f"With {percentage*100}% certainty, the value is in the interval ({min}, {max})")