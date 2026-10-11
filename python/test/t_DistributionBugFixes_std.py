#! /usr/bin/env python
# Regression tests for the fix_distribution_continuous_A_L branch (continuous laws A-L).
# Each block covers a bug fixed by direct source reading (see doc/distribution_audit.tex).
# Pure asserts, no reference output file needed.

import math as m
import os


import openturns as ot
import openturns.testing as ott

ot.TESTPREAMBLE()
ot.RandomGenerator.SetSeed(0)

studyFile = "myStudyBugFixes.xml"
# --- Arcsine: setAB invalidates cached mean/covariance, kurtosis is 1.5 ---
arcsine = ot.Arcsine(0.0, 1.0)
assert arcsine.getMean()[0] == 0.5
arcsine.setA(-1.0)
ott.assert_almost_equal(arcsine.getMean(), [0.0], 0.0, 0.0)
ott.assert_almost_equal(arcsine.getKurtosis(), [1.5], 0.0, 0.0)
with ott.assert_raises(Exception):
    ot.ArcsineMuSigma().setValues([0.0, -1.0])

# --- Beta: endpoint density, DDF guard, Arcsine shortcut, roughness guard ---
ott.assert_almost_equal(ot.Beta(3.0, 1.0, 0.0, 2.0).computePDF(2.0), 1.5, 1e-12, 0.0)
ott.assert_almost_equal(
    ot.Beta(3.0, 1.0, 0.0, 2.0).computeLogPDF(2.0), m.log(1.5), 1e-12, 0.0
)
ott.assert_almost_equal(ot.Beta(2.0, 1.0, 0.0, 1.0).computeDDF([1.0]), [0.0], 0.0, 0.0)
assert (
    ot.Beta(0.5, 0.5, -1.0, 1.0).getStandardRepresentative().getImplementation().getClassName()
    == "Arcsine"
)
with ott.assert_raises(Exception):
    ot.Beta(0.3, 0.4).getRoughness()
with ott.assert_raises(Exception):
    ot.BetaMuSigma().setValues([0.5, 10.0, 0.0, 1.0])
with ott.assert_raises(Exception):
    ot.BetaMuSigma()([0.5, 10.0, 0.0, 1.0])

# --- Burr: nonexistent moments throw, tail log-density is finite ---
with ott.assert_raises(Exception):
    ot.Burr(2.0, 0.5).getMean()
ott.assert_almost_equal(ot.Burr(2.0, 2.0).getMean()[0], m.pi / 4.0, 1e-12, 0.0)
assert m.isfinite(ot.Burr(2.0, 2.0).computeLogPDF([1.0e150]))

# --- Chi family: gradients valid for small shapes ---
assert m.isfinite(ot.Chi(0.01).computeCDFGradient([0.5])[0])
assert m.isfinite(ot.ChiSquare(0.01).computeCDFGradient([0.5])[0])
with ott.assert_raises(Exception):
    ot.ChiFactory().build(ot.Sample([[1.0], [-2.0], [3.0]]))
with ott.assert_raises(Exception):
    ot.ChiSquareFactory().build(ot.Sample([[1.0], [-2.0], [3.0]]))

# --- CompoundDistribution keeps its link function ---
link1d = ot.SymbolicFunction(["y"], ["y", "1.0"])
compound = ot.CompoundDistribution(ot.Normal(), ot.Exponential(1.0), link1d)
pdfBefore = compound.computePDF([0.0])
compound.setParameter(compound.getParameter())
ott.assert_almost_equal(compound.computePDF([0.0]), pdfBefore, 1e-12, 0.0)

# --- Exponential: gradients agree with density at the atom ---
ott.assert_almost_equal(ot.Exponential().computePDFGradient([0.0]), [0.0, 0.0], 0.0, 0.0)
ott.assert_almost_equal(ot.Exponential().computeCDFGradient([0.0]), [0.0, 0.0], 0.0, 0.0)

# --- Frechet: wrong-dimension inputs are rejected ---
with ott.assert_raises(Exception):
    ot.Frechet().computePDF([0.5, 0.5])

# --- GumbelFactory rejects constant samples, helper prints true names ---
with ott.assert_raises(Exception):
    ot.GumbelFactory().build(ot.Sample([[1.0], [1.0], [1.0]]))
assert "lambda" in str(ot.GumbelLambdaGamma())
with ott.assert_raises(Exception):
    ot.GumbelMuSigma().setValues([0.0, -1.0])
with ott.assert_raises(Exception):
    ot.GumbelLambdaGamma().setValues([-1.0, 0.0])

# --- InverseNormal: CDF stays in [0, 1] over a wide parameter grid ---
for mu, lambd in [(0.01, 100.0), (1.0, 1.0), (2.0, 0.5), (0.5, 50.0)]:
    distIN = ot.InverseNormal(mu, lambd)
    for x in [0.01, 0.1, 1.0, 5.0, 50.0]:
        c = distIN.computeCDF([x])
        assert 0.0 <= c <= 1.0 and m.isfinite(c)
ott.assert_almost_equal(ot.InverseNormal(0.01, 100.0).computeCDF([1.0]), 1.0, 0.0, 0.0)

# --- JointByConditioningDistribution: empty conditioning domain gives 0 ---
link = ot.SymbolicFunction(["y"], ["y", "1.0"])
jbcd = ot.JointByConditioningDistribution(ot.Normal(), ot.Normal(), link)
assert jbcd.computeCDF([-100.0, -100.0]) == 0.0

# --- JointByConditioningDistribution: trailing conditioned state is kept ---
link2d = ot.SymbolicFunction(["y"], ["y", "1.0", "y + 1.0", "2.0"])
conditioned2d = ot.JointDistribution([ot.Normal(), ot.Normal()])
jbcd2d = ot.JointByConditioningDistribution(conditioned2d, ot.Uniform(0.0, 1.0), link2d)
reference2d = ot.JointDistribution(
    [ot.Normal(0.3, 1.0), ot.Normal(1.3, 2.0)], ot.IndependentCopula(2)
)
ott.assert_almost_equal(
    jbcd2d.computeConditionalPDF(1.0, [0.3, 0.5]),
    reference2d.computeConditionalPDF(1.0, [0.5]),
    1e-12,
    0.0,
)

# --- LogNormal helpers validate their inputs ---
with ott.assert_raises(Exception):
    ot.LogNormalMuSigma().setValues([0.0, -1.0, 0.0])
with ott.assert_raises(Exception):
    ot.LogNormalMuSigmaOverMu().setValues([0.0, 0.5, -1.0])
with ott.assert_raises(Exception):
    ot.LogNormalMuErrorFactor().setValues([0.0, 0.5, 0.0])
with ott.assert_raises(Exception):
    ot.GammaMuSigma().setValues([0.0, -1.0, 0.0])
with ott.assert_raises(Exception):
    ot.WeibullMaxMuSigma().setValues([2.0, 1.0, 0.0])
with ott.assert_raises(Exception):
    ot.WeibullMinMuSigma().setValues([0.0, 1.0, 2.0])

# --- LogNormalFactory: invalid shifts are rejected ---
with ott.assert_raises(Exception):
    ot.LogNormalFactory().buildMethodOfLeastSquares(
        ot.Sample([[1.0], [2.0], [3.0]]), 5.0
    )
with ott.assert_raises(Exception):
    ot.ParetoFactory().buildMethodOfLeastSquares(
        ot.Sample([[1.0], [2.0], [3.0]]), 5.0
    )

# --- LogUniform: characteristic function at zero, DDF at the atom ---
ott.assert_almost_equal(
    [ot.LogUniform(1.0, 2.0).computeCharacteristicFunction(0.0).real],
    [1.0],
    0.0,
    0.0,
)
ott.assert_almost_equal(ot.LogUniform(1.0, 2.0).computeDDF([1.0]), [0.0], 0.0, 0.0)

# --- InverseNormalFactory checks positivity on the moments path ---
ot.ResourceMap.SetAsString("InverseNormalFactory-Method", "Moments")
try:
    with ott.assert_raises(Exception):
        ot.InverseNormalFactory().build(ot.Sample([[-1.0], [2.0], [3.0], [4.0]]))
finally:
    ot.ResourceMap.SetAsString("InverseNormalFactory-Method", "MLE")

# --- ChiSquare/FisherSnedecor: persistence round-trips ---
study = ot.Study(studyFile)
study.add("chisquare", ot.ChiSquare(3.0))
study.add("fisher", ot.FisherSnedecor(4.0, 5.0))
study.save()
study2 = ot.Study(studyFile)
study2.load()
chiR = ot.ChiSquare(1.0)
study2.fillObject("chisquare", chiR)
ott.assert_almost_equal(chiR.getParameter(), [3.0], 0.0, 0.0)
ott.assert_almost_equal(chiR.computePDF([2.0]), ot.ChiSquare(3.0).computePDF([2.0]), 1e-12, 0.0)
fisherR = ot.FisherSnedecor(1.0, 1.0)
study2.fillObject("fisher", fisherR)
ott.assert_almost_equal(fisherR.getParameter(), [4.0, 5.0], 0.0, 0.0)
ott.assert_almost_equal(
    fisherR.computePDF([1.0]), ot.FisherSnedecor(4.0, 5.0).computePDF([1.0]), 1e-12, 0.0
)
os.remove(studyFile)

# --- ArcsineFactory/BetaFactory/GumbelFactory: clean errors on bad samples ---
with ott.assert_raises(Exception):
    ot.ArcsineFactory().build(ot.Sample([[1.0], [1.0], [1.0]]))
with ott.assert_raises(Exception):
    ot.GumbelFactory().build(ot.Sample([[1.0], [1.0], [1.0]]))

# --- KernelSmoothing: bandwidth validation and degenerate cache ---
ks = ot.KernelSmoothing()
with ott.assert_raises(Exception):
    ks.build(ot.Normal().getSample(10), [-1.0])
with ott.assert_raises(Exception):
    ks.buildWeighted(
        ot.Sample([[float("nan")], [1.0], [2.0]]), [0.25, 0.25, 0.5], [0.5]
    )
constSample = ot.Sample([[2.0]] * 5)
ks.buildWeighted(constSample, [0.2] * 5, [0.5])
ott.assert_almost_equal(ks.getBandwidth(), [0.5], 0.0, 0.0)

# --- InverseWishart: label count matches parameter count ---
dist = ot.InverseWishart(ot.CovarianceMatrix(2, [1.0, 0.2, 0.2, 1.0]), 5.0)
assert len(dist.getParameterDescription()) == len(dist.getParameter()) == 4
