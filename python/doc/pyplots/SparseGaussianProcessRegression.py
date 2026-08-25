import openturns as ot
import openturns.viewer as otv
from openturns.experimental import SparseGaussianProcessRegression

f = ot.SymbolicFunction(["x"], ["x + x * sin(x)"])
sampleX = [[1.0], [2.0], [3.0], [4.5], [6.0], [7.0], [8.0]]
sampleY = f(sampleX)
covarianceModel = ot.SquaredExponential([1.0])
covarianceModel.setActiveParameter([])
inducingPoints = ot.Sample([[1.5], [3.0], [4.5], [6.0], [7.5]])
algo = SparseGaussianProcessRegression(sampleX, sampleY, covarianceModel, inducingPoints)
algo.run()
result = algo.getResult()

metaModel = result.getMetaModel()
graph = metaModel.draw(0.5, 8.5)
graph.add(ot.Cloud(sampleX, sampleY))
cloud = ot.Cloud(result.getInducingPoints(), metaModel(result.getInducingPoints()))
cloud.setColor("orange")
cloud.setPointStyle("square")
cloud.setLegend("inducing points")
graph.add(cloud)
graph.setLegends(["metamodel", "sample", "inducing points"])
graph.setLegendPosition("upper left")
otv.View(graph, figure_kw={"figsize": (8, 4)})
