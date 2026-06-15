#! /usr/bin/env python

import math
import openturns as ot
import openturns.testing as ott

ot.TESTPREAMBLE()
ot.PlatformInfo.SetNumericalPrecision(5)


# Instantiate one distribution object
dim = 2
copula = ot.EmpiricalBernsteinCopula(ot.Normal(2).getSample(12), 3)
print("Copula ", repr(copula))
print("Copula ", copula)
print("Mean ", repr(copula.getMean()))
print("Covariance ", repr(copula.getCovariance()))

# Is this copula an elliptical distribution?
print("Elliptical distribution= ", copula.isElliptical())

# Is this copula elliptical ?
print("Elliptical copula= ", copula.hasEllipticalCopula())

# Is this copula independent ?
print("Independent copula= ", copula.hasIndependentCopula())

# Test for realization of distribution
oneRealization = copula.getRealization()
print("oneRealization=", repr(oneRealization))

# Test for sampling
size = 10
oneSample = copula.getSample(size)
print("oneSample=", repr(oneSample))

# Test for sampling
size = 10000
anotherSample = copula.getSample(size)
print("anotherSample mean=", repr(anotherSample.computeMean()))
print("anotherSample covariance=", repr(anotherSample.computeCovariance()))

# Define a point
point = [0.2] * dim

# Show PDF and CDF of point
pointPDF = copula.computePDF(point)
pointCDF = copula.computeCDF(point)
print("Point = ", repr(point), " pdf=%.6f" % pointPDF, " cdf=%.6f" % pointCDF)

# Get 50% quantile
quantile = copula.computeQuantile(0.5)
print("Quantile=", repr(quantile))
print("CDF(quantile)=%.6f" % copula.computeCDF(quantile))

# Get 5% quantile
quantile = copula.computeQuantile(0.95, True)
print("Quantile=", repr(quantile))

# Get 95% survival function
inverseSurvival = ot.Point(copula.computeInverseSurvivalFunction(0.95))
print("InverseSurvival=", repr(inverseSurvival))
print(
    "Survival(inverseSurvival)=%.6f" % copula.computeSurvivalFunction(inverseSurvival)
)
print("entropy=%.6f" % copula.computeEntropy())

# Confidence regions
interval, threshold = copula.computeMinimumVolumeIntervalWithMarginalProbability(0.95)
print("Minimum volume interval=", interval)
print("threshold=", ot.Point(1, threshold))
levelSet, beta = copula.computeMinimumVolumeLevelSetWithThreshold(0.95)
print("Minimum volume level set=", levelSet)
print("beta=", ot.Point(1, beta))
interval, beta = copula.computeBilateralConfidenceIntervalWithMarginalProbability(0.95)
print("Bilateral confidence interval=", interval)
print("beta=", ot.Point(1, beta))
interval, beta = copula.computeUnilateralConfidenceIntervalWithMarginalProbability(
    0.95, False
)
print("Unilateral confidence interval (lower tail)=", interval)
print("beta=", ot.Point(1, beta))
interval, beta = copula.computeUnilateralConfidenceIntervalWithMarginalProbability(
    0.95, True
)
print("Unilateral confidence interval (upper tail)=", interval)
print("beta=", ot.Point(1, beta))
print("parameters=", copula.getParameter())
copula.setParameter(copula.getParameter())

# Extract the marginals
for i in range(dim):
    margin = copula.getMarginal(i)
    print("margin=", repr(margin))
    print("margin PDF=%.6f" % margin.computePDF([0.25]))
    print("margin CDF=%.6f" % margin.computeCDF([0.25]))
    print("margin quantile=", repr(margin.computeQuantile(0.95)))
    print("margin realization=", repr(margin.getRealization()))

# Extract a 2-D marginal
indices = [1, 0]
print("indices=", repr(indices))
margins = copula.getMarginal(indices)
print("margins=", repr(margins))
print("margins PDF=%.6f" % margins.computePDF([0.25] * 2))
print("margins CDF=%.6f" % margins.computeCDF([0.25] * 2))
quantile = ot.Point(margins.computeQuantile(0.95))
print("margins quantile=", repr(quantile))
print("margins CDF(qantile)=%.6f" % margins.computeCDF(quantile))
print("margins realization=", repr(margins.getRealization()))

copula6D = ot.EmpiricalBernsteinCopula(ot.Normal(6).getSample(8), 4)
print("Entropy in higher dimension=%.6f" % copula6D.computeEntropy())

dim = 6
x = 0.6
y = [0.2] * (dim - 1)
print("conditional PDF=%.6f" % copula6D.computeConditionalPDF(x, y))
print(
    "conditional PDF ref=%.6f"
    % (copula6D.computePDF(y + [x]) / copula6D.getMarginal(range(5)).computePDF(y))
)
print("conditional CDF=%.6f" % copula6D.computeConditionalCDF(x, y))
print("conditional quantile=%.6f" % copula6D.computeConditionalQuantile(x, y))
pt = [0.05 * (1 + i) for i in range(dim)]
print("sequential conditional PDF=", copula6D.computeSequentialConditionalPDF(pt))
resCDF = copula6D.computeSequentialConditionalCDF(pt)
print("sequential conditional CDF(", pt, ")=", resCDF)
print(
    "sequential conditional quantile(",
    resCDF,
    ")=",
    copula6D.computeSequentialConditionalQuantile(resCDF),
)

# Test isEmpiricalCopulaSample=True with size not multiple of binNumber
ref_sample = ot.Normal(2).getSample(14)
copula_raw = ot.EmpiricalBernsteinCopula(ref_sample, 4, True)
print("isEmpiricalCopulaSample=True: ", copula_raw)
ott.assert_almost_equal(copula_raw.getMean(), [0.5, 0.5])

# Test setCopulaSample
copula_cpy = ot.EmpiricalBernsteinCopula(ot.Normal(2).getSample(12), 3)
copula_cpy.setCopulaSample(ot.Normal(2).getSample(20), True)
print("After setCopulaSample: ", copula_cpy)

# Test setBinNumber
copula_cpy.setBinNumber(5)
print("After setBinNumber: ", copula_cpy)

# Test computeProbability
prob = copula.computeProbability(ot.Interval([0.1, 0.2], [0.6, 0.8]))
print("Probability interval=%.6f" % prob)

# Test computeProbability with interval extending outside [0,1]^d
prob_ext = copula.computeProbability(ot.Interval([-0.5, 0.2], [1.5, 0.8]))
print("Probability interval (extended)=%.6f" % prob_ext)

# Test getSpearmanCorrelation
spearman = copula.getSpearmanCorrelation()
print("Spearman correlation= ", spearman)

# Test 1D case: hasEllipticalCopula and hasIndependentCopula return true
copula1D = ot.EmpiricalBernsteinCopula(ot.Normal(1).getSample(10), 5)
print("1D elliptical= ", copula1D.hasEllipticalCopula())
print("1D independent= ", copula1D.hasIndependentCopula())

# Test boundary values: PDF returns 0 at [0,0] and [1,1]
print("PDF at 0=%.6f" % copula.computePDF([0.0, 0.0]))
print("PDF at 1=%.6f" % copula.computePDF([1.0, 1.0]))
print("LogPDF at 0=%.6f" % copula.computeLogPDF([0.0, 0.0]))

# Test computeLogPDF vs log(computePDF)
point = [0.3, 0.7]
pdf_val = copula.computePDF(point)
logpdf_val = copula.computeLogPDF(point)
print("logPDF=%.6f log(PDF)=%.6f" % (logpdf_val, math.log(pdf_val)))

# Test binNumber=1 edge case (PDF simplifies to 1.0)
copula_b1 = ot.EmpiricalBernsteinCopula(ot.Normal(2).getSample(12), 1)
print("binNumber=1 PDF=%.6f" % copula_b1.computePDF([0.3, 0.7]))
print("binNumber=1 CDF=%.6f" % copula_b1.computeCDF([0.3, 0.7]))

ot.Log.Show(ot.Log.TRACE)
validation = ott.DistributionValidation(copula)
validation.skipGradient()  # computePDFGradient does not handle integer parameters
validation.run()
