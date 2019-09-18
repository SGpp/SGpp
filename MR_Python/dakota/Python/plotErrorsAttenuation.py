import numpy as np
import matplotlib.pyplot as plt

# # attenuationU 3D
# realMean = 0.6150055735340820
# realStDv = 0.1028798173353860
# # Dakota attenuationU 3D
# gridSizes = [7, 31, 111, 303]
# means = [6.1496640320635632e-01, 6.1500551308585627e-01, 6.1500557353405572e-01, 6.1500557353408203e-01 ]
# stdvs = [ 1.0122898046748717e-01, 1.0286862415567162e-01, 1.0287981725535035e-01, 1.0287981733647278e-01 ]
# dakotaMeanErrs = [abs(mean - realMean) for mean in means]
# dakotaStDvErrs = [abs(stdv - realStDv) for stdv in stdvs]
# # SGpp attenuationU 3D
# sgppGridSizes = [10, 26, 70, 188, 499]
# sgppMeans = [0.6150054493496887, 0.6150054493496887, 0.6150042574178122, 0.6150055682637233, 0.6150055735303124 ] 
# sgppStDvs = [0.1028665693980405, 0.1028665693980405, 0.1028832358833377, 0.1028798196682789, 0.1028798173353860]
# sgppMeanErrs = [abs(mean - realMean) for mean in sgppMeans]
# sgppStDvErrs = [abs(stdv - realStDv) for stdv in sgppStDvs]

# attenuationN 3D mu=1, sigma=0.1 
realMean = 0.3684930848010390
realStDv = np.sqrt(0.0004533790571205)
# Dakota attenuationN 3D
gridSizes = [7, 37, 117, 327 ]
means = [3.6849274390687553e-01, 3.6849308473789277e-01, 3.6849308480103860e-01, 3.6849308480103921e-01  ]
stdvs = [2.1251331144322486e-02, 2.1292665690671824e-02, 2.1292699620269046e-02, 2.1292699620303512e-02 ]
dakotaMeanErrs = [abs(mean - realMean) for mean in means]
dakotaStDvErrs = [abs(stdv - realStDv) for stdv in stdvs]
# SGpp attenuationN 3D
sgppGridSizes = [31, 79, 155, 280 , 505]
sgppMeans = [0.3684930858910553, 0.3684931553222115, 0.3684930848181190, 0.3684930846738973 , 0.3684930848011336  ] 
sgppStDvs = [0.0212929405486851, 0.0212924458774328, 0.0212926992490811, 0.0212926999223109, 0.0212926996204060]
sgppMeanErrs = [abs(mean - realMean) for mean in sgppMeans]
sgppStDvErrs = [abs(stdv - realStDv) for stdv in sgppStDvs]

plt.title('mean')
plt.semilogy(gridSizes, dakotaMeanErrs, '-*b', label='dakota')
plt.semilogy(sgppGridSizes, sgppMeanErrs, '-*r', label='SG++')
plt.legend()

plt.figure()
plt.title('stdv')
plt.semilogy(gridSizes, dakotaStDvErrs, '-*b', label='dakota')
plt.semilogy(sgppGridSizes, sgppStDvErrs, '-*r', label='SG++')
plt.legend()

plt.show()

