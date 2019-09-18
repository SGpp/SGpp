import numpy as np
import scipy.integrate as integrate

mu = 0.1
sigma = 0.0161912
lower_limit = 0.05
upper_limit = 0.15
mean = integrate.quad(lambda x: x * (np.exp(-(x - mu) ** 2 / (2 * sigma ** 2)) / (np.sqrt(2 * np.pi) * sigma)), lower_limit, upper_limit)
print("mean: {}".format(mean))

mean2 = integrate.quad(lambda x: x ** 2 * (np.exp(-(x - mu) ** 2 / (2 * sigma ** 2)) / (np.sqrt(2 * np.pi) * sigma)), lower_limit, upper_limit)
var = mean2[0] - mean[0] ** 2
stdv = np.sqrt(var)
  
print("var: {}".format(var))
print("stdv: {}".format(stdv))
