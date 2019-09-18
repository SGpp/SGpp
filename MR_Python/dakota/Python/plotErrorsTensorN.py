import numpy as np
import matplotlib.pyplot as plt

# tensorMonomialN 6D
realMean = 1.
realStDv = np.sqrt(15624)

# Dakota
gridSizes = [13, 109, 629, 2909, 10901]
means = [9.9999999999999778e-01, 9.9999999999999634e-01, 9.9999999999999911e-01, 1.0000000000000058e+00, 9.9999999999999500e-01  ]
stdvs = [4.8989794855663531e+00, 1.6248076809271907e+01, 3.9293765408776970e+01, 7.3375745311376463e+01, 1.0736852425175611e+02]
dakotaMeanErrs = [abs(mean - realMean) for mean in means]
dakotaStDvErrs = [abs(stdv - realStDv) for stdv in stdvs]

# SGpp numRefine 100
sgppGridSizes_a = [13, 97, 545, 1376, 6609 ]
sgppMeans_a = [0.9999999999999784, 0.9999999999999803, 0.9999999999998878, 1.0024386526528750, 1.0000000154475477] 
sgppStDvs_a = [4.8989794855663069, 16.2480768092717440, 39.2937654087769275, 71.5295263022036352, 124.9959999366083423]
sgppMeanErrs_a = [abs(mean - realMean) for mean in sgppMeans_a]
sgppStDvErrs_a = [abs(stdv - realStDv) for stdv in sgppStDvs_a]

# SGpp numRefine 50
sgppGridSizes_a50 = [13, 97, 453, 2010 , 6363]
sgppMeans_a50 = [0.9999999999999784, 0.9999999999999803, 1.0438957475995576, 1.0062321123726998, 1.0022131371564538 ] 
sgppStDvs_a50 = [4.8989794855663069, 16.2480768092717440, 38.8045922059172526, 90.3217835090174930, 130.4847434354594782]
sgppMeanErrs_a50 = [abs(mean - realMean) for mean in sgppMeans_a50]
sgppStDvErrs_a50 = [abs(stdv - realStDv) for stdv in sgppStDvs_a50]

# SGpp numRefine 200
sgppGridSizes_a200 = [13, 97, 545, 1934, 6030 ]
sgppMeans_a200 = [0.9999999999999784, 0.9999999999999805, 0.9999999999998878, 0.9999999999694443, 0.9999999968126423 ] 
sgppStDvs_a200 = [4.8989794855663069, 16.2480768092717440, 39.2937654087769275, 73.3757453113785374, 124.9959999360369807]
sgppMeanErrs_a200 = [abs(mean - realMean) for mean in sgppMeans_a200]
sgppStDvErrs_a200 = [abs(stdv - realStDv) for stdv in sgppStDvs_a200]

# SGpp regular
sgppGridSizes_r = [1, 13, 97, 545, 2561]  # , 10625 ]
sgppMeans_r = [0.9999999999999787, 0.9999999999999784, 0.9999999999999785, 0.9999999999999778, 1.0000000000007787]
sgppStDvs_r = [ 0.0000001504943511, 4.8989794855663069, 16.2480768092717440, 39.2937654087768351, 73.3757453113780258]
sgppMeanErrs_r = [abs(mean - realMean) for mean in sgppMeans_r]
sgppStDvErrs_r = [abs(stdv - realStDv) for stdv in sgppStDvs_r]

plt.title('mean')
plt.semilogy(gridSizes, dakotaMeanErrs, '-*b', label='dakota')
plt.semilogy(sgppGridSizes_r, sgppMeanErrs_r, '-*r', label='regular')
plt.semilogy(sgppGridSizes_a, sgppMeanErrs_a, '-*g', label='100')
plt.semilogy(sgppGridSizes_a50, sgppMeanErrs_a50, '--*g', label='50')
plt.semilogy(sgppGridSizes_a50, sgppMeanErrs_a200, '-.*g', label='200')
plt.legend()

plt.figure()
plt.title('stdv')
plt.semilogy(gridSizes, dakotaStDvErrs, '-*b', label='dakota')
plt.semilogy(sgppGridSizes_r, sgppStDvErrs_r, '-*r', label='regular')
plt.semilogy(sgppGridSizes_a, sgppStDvErrs_a, '-*g', label='100')
plt.semilogy(sgppGridSizes_a50, sgppStDvErrs_a50, '--*g', label='50')
plt.semilogy(sgppGridSizes_a200, sgppStDvErrs_a200, '-.*g', label='200')
plt.legend()

plt.show()

