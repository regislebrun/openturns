"""
Christoffel subsampling on the square
=====================================
"""

# %%
# Optimal sampling draws the design points from the Christoffel distribution
# associated to the approximation space instead of the reference measure.
# On the square with a tensor Legendre basis everything is explicit, so this
# example reproduces the textbook picture: the optimal density concentrates
# at the exiting corners, and the subsampled weighted least squares needs far
# fewer evaluations than uniform sampling. See Cohen-Dolbeault (Fig. 1, Sec. 6).

# %%
import openturns as ot
import openturns.experimental as otexp
import openturns.viewer as otv

ot.RandomGenerator.SetSeed(0)

# %%
# Reference measure, tensor Legendre basis and total-degree space, degree 3
measure = ot.JointDistribution([ot.Uniform(-1.0, 1.0)] * 2)
basis = ot.OrthogonalBasis(ot.OrthogonalProductPolynomialFactory([ot.Uniform(-1.0, 1.0)] * 2))
spaceDimension = 10
christoffel = otexp.ChristoffelDistribution(basis, spaceDimension)
target = ot.SymbolicFunction(['x0', 'x1'], ['exp(-(x0^2 + x1^2))'])

# %%
# Optimal density heatmap: k_m / m concentrates at the corners
nX = nY = 80
mesh = ot.Box([nX - 2, nY - 2], ot.Interval([-1.0] * 2, [1.0] * 2)).generate()
density = christoffel.computeChristoffel(mesh)
density /= spaceDimension
abscissae = ot.Sample([[-1.0 + 2.0 * i / (nX - 1)] for i in range(nX)])
contour = ot.Contour(abscissae, abscissae, density)
graph = ot.Graph("Christoffel density", "x1", "x2")
graph.add(contour)
view = otv.View(graph)

# %%
# Subsampled design against a uniform design with the same budget
budget = 231
experiment = otexp.SubsampleChristoffelExperiment(basis, budget)
sample, weights = experiment.generateWithWeights()
uniformSample = measure.getSample(budget)


def fit(design, designWeights):
    observations = target(design)
    weights = ot.Point(designWeights)
    expansion = ot.LeastSquaresExpansion(design, weights, observations, measure, basis, spaceDimension)
    expansion.run()
    return expansion.getResult().getMetaModel()


metamodel = fit(sample, weights)
uniformMetamodel = fit(uniformSample, ot.Point(budget, 1.0))

# %%
# Validation error and Gramian conditioning on both designs
validation = measure.getSample(2000)
reference = target(validation)
norm = reference.computeRawMoment(2)[0] ** 0.5
error = ((metamodel(validation) - reference).computeRawMoment(2)[0]) ** 0.5 / norm
uniformError = ((uniformMetamodel(validation) - reference).computeRawMoment(2)[0]) ** 0.5 / norm
print(f"subsampled error={error:.3e} uniform error={uniformError:.3e}")


def condition(design, designWeights):
    rows = design.getSize()
    functions = [basis.build(j) for j in range(spaceDimension)]
    gramian = ot.SymmetricMatrix(spaceDimension)
    for a in range(spaceDimension):
        for b in range(a + 1):
            gramian[a, b] = sum(designWeights[i] * functions[a](design[i])[0] * functions[b](design[i])[0] for i in range(rows)) / rows
    eigenvalues = sorted(gramian.computeEigenValues())
    return eigenvalues[-1] / eigenvalues[0]


print(f"subsampled conditioning={condition(sample, weights):.3f} uniform conditioning={condition(uniformSample, [1.0] * budget):.3f}")

# %%
# The thinned design clusters where the density is large
cloud = ot.Cloud(sample, "blue", "fsquare", "subsample")
uniformCloud = ot.Cloud(uniformSample, "red", "fcircle", "uniform")
designGraph = ot.Graph("Designs", "x1", "x2")
designGraph.add(cloud)
designGraph.add(uniformCloud)
designGraph.setLegends(["subsample", "uniform"])
view = otv.View(designGraph)

# %%
# Display all figures
otv.View.ShowAll()
