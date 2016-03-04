# define random variables
builder = ParameterBuilder()
up = builder.defineUncertainParameters()

up.new().isCalled('y1').withUniformDistribution(-1, 1)
up.new().isCalled('y2').withUniformDistribution(-1, 1)
up.new().isCalled('y3').withUniformDistribution(-1, 1)

params = builder.andGetResult()
