import numpy as np
import matplotlib.pyplot as plt

# Dakota
gridSizes = [18, 177, 1260 , 7169]  # , 34290]
means = [ 6.5632828919103816e-03, 6.5768859155912159e-03, 6.5770271305317547e-03, 6.5770283342825135e-03]  # , 6.5770283458260678e-03]
stdvs = [2.6322104837547540e-03, 2.8567031011856765e-03 , 2.8672908560281080e-03, 2.8675734257209411e-03]  # , 2.8675785527979727e-03]

scalingFactor = 1.0 / 0.9960008520179663

# SGpp regular
sgppsizes_regular = [ 17, 161, 1121, 6401]  # , 31745] 
SGppMeans_regular = [0.0065259130029695, 0.0065498295576200, 0.0065506962800845, 0.0065507250599070]  # , 0.0065507258125752]
scaledSGppMeans_regular = [mean * scalingFactor  for mean in SGppMeans_regular]
SGppStDvs_regular = [0.0026317311114392, 0.0028764501011042, 0.0028911784736958, 0.0028917664643657]

# SGpp adaptive numRefine 100 
sgppsizes_adaptive100 = [17, 161, 118, 1001, 7052]
SGppMeans_adaptive100 = [0.0065259130029695, 0.0065259130029695, 0.0065498295576200, 0.0065506962800852, 0.0065507257143593]
scaledSGppMeans_adaptive100 = [mean * scalingFactor  for mean in SGppMeans_adaptive100]
SGppStDvs_adaptive100 = [0.0026317311114392, 0.0026317311114392, 0.0028764501011042, 0.0028911784736984, 0.0028917851648883]

# SGpp adaptive numRefine 20 
sgppsizes_adaptive20 = [17, 161, 401, 1425, 5022]
SGppMeans_adaptive20 = [0.0065259130029695, 0.0065498295576200, 0.0065506879300093, 0.0065506719093779, 0.0065507256190988]
scaledSGppMeans_adaptive20 = [mean * scalingFactor  for mean in SGppMeans_adaptive20]
SGppStDvs_adaptive20 = [0.0026317311114392, 0.0028764501011042, 0.0028909238741513, 0.0028917629363247, 0.0028917850033505]

# SGpp adaptive numRefine 10
sgppsizes_adaptive10 = [17, 137, 415, 1461, 5099 ]
SGppMeans_adaptive10 = [0.0065259130029695, 0.0065498294188831, 0.0065506172787326, 0.0065507291304433, 0.0065507254857405 ]
scaledSGppMeans_adaptive10 = [mean * scalingFactor   for mean in SGppMeans_adaptive10]
SGppMeanSquared_adaptive10 = [ 4.951354916524374e-05, 5.117423059331098e-05, 5.126801880560545e-05, 5.127449990322398e-05, 5.127442491310185e-05 ]
SGppStDvs_adaptive10 = [ 0.0026317311114392, 0.0028764500998357, 0.0028909223568223, 0.0028917897162113, 0.0028917850064575]
scaledSGppStDvs_adaptive10 = [10] * len(sgppsizes_adaptive10)
for i in range(len(sgppsizes_adaptive10)):
    mean = SGppMeans_adaptive10[i]
    meanSquare = SGppMeanSquared_adaptive10[i]
    scaledSGppStDvs_adaptive10[i] = np.sqrt(meanSquare * scalingFactor - (mean * scalingFactor) ** 2)

# SGpp adaptive numRefine 5
sgppsizes_adaptive5 = [17, 88, 159, 415, 1390, 5044 ]
SGppMeans_adaptive5 = [0.0065259130029695, 0.0065979974213449, 0.0065517428806454, 0.0065506172787326, 0.0065506700795172, 0.0065507254803619 ]
scaledSGppMeans_adaptive5 = [mean * scalingFactor   for mean in SGppMeans_adaptive5]
SGppStDvs_adaptive5 = [0.0026317311114392, 0.0028736438042499, 0.0028877233431313, 0.0028909223568223, 0.0028917574987656, 0.0028917850039092]

# dakota 34k
# trueMean = 6.5770283458260678e-03 
# SGpp 35k, difference is 1e-11, convergence results are the same for both
trueMean = 0.00655072585279 * scalingFactor
dakErrM = [abs(mean - trueMean) for mean in means]
sgregErrM = [abs(mean - trueMean) for mean in scaledSGppMeans_regular]
sgada100ErrM = [abs(mean - trueMean) for mean in scaledSGppMeans_adaptive100]
sgada20ErrM = [abs(mean - trueMean) for mean in scaledSGppMeans_adaptive20]
sgada10ErrM = [abs(mean - trueMean) for mean in scaledSGppMeans_adaptive10]
sgada5ErrM = [abs(mean - trueMean) for mean in scaledSGppMeans_adaptive5]

plt.title('mean')
plt.semilogy(gridSizes, dakErrM, '-*b', label='Dakota')
# plt.semilogy(sgppsizes_regular, sgregErrM, '-*g')
# plt.semilogy(sgppsizes_adaptive100, sgada100ErrM, '-*y')
# plt.semilogy(sgppsizes_adaptive20, sgada20ErrM, '-*m')
plt.semilogy(sgppsizes_adaptive10, sgada10ErrM, '-*r', label='SGpp')
# plt.semilogy(sgppsizes_adaptive5, sgada5ErrM, '-*k')
plt.legend()

# I have not recorded meanSqure for all sgpp stuff, on ly for the best in case of mean, adaptive10
meanSGpp = 0.00655072585279
meanSquareSGpp = 5.127443183386146e-05

# dakota with 34k
# trueStDv = 2.8675785527979727e-03
# SGpp with 35 k, difference is 1e-11, convergence results are the same for both  
trueStDv = np.sqrt(meanSquareSGpp * scalingFactor - (meanSGpp * scalingFactor) ** 2)

dakErrS = [abs(stdv - trueStDv) for stdv in stdvs]
sgppErrS = [abs(stdv - trueStDv) for stdv in scaledSGppStDvs_adaptive10]

plt.figure()
plt.title('stdv')
plt.semilogy(gridSizes, dakErrS, '-*b', label='Dakota')
plt.semilogy(sgppsizes_adaptive10, sgppErrS, '-*r', label='SGpp')
plt.legend()

plt.show()

