import numpy as np
import scipy.integrate as integrate

mu = 7.71
sigma = 1.0056
lower_limit = 100
upper_limit = 50000
mean = integrate.quad(lambda x: x * (1. / (x * sigma * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu) ** 2 / (2 * sigma ** 2)), lower_limit, upper_limit)
print("mean: {:.16f}".format(mean[0]))

mean2 = integrate.quad(lambda x: x ** 2 * (1. / (x * sigma * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu) ** 2 / (2 * sigma ** 2)), lower_limit, upper_limit)
var = mean2[0] - mean[0] ** 2
stdv = np.sqrt(var)
  
print("var: {:.16f}".format(var))
print("stdv: {:.16f}".format(stdv))
