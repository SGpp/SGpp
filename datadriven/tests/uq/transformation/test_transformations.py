from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder


builder = ParameterBuilder()
dp = builder.defineDeterministicParameters()
up = builder.defineUncertainParameters()

up.new().isCalled('v').withUniformDistribution(0, 10).withInverseCDFTransformation()
dp.new().isCalled('density').hasValue(.3)
up.new().isCalled('theta').withTNormalDistribution(1, 1, -2, 2).withLinearTransformation()
dp.new().isCalled('radius').hasValue(10)
up.new().isCalled('height').withBetaDistribution(3, 3, 0, 2).withInverseCDFTransformation()
params = builder.andGetResult()

ap = params.activeParams()

dist = ap.getIndependentJointDistribution()
trans = ap.getJointTransformation()

for _ in range(1000):
    prob1 = dist.rvs()[0]
    unit = trans.probabilisticToUnit(prob1)
    prob2 = trans.unitToProbabilistic(unit)

    assert all(["%g" % x == "%g" % y for x, y in zip(prob1, prob2)])
