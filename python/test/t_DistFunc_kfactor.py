#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott

nMin = 2
nMax = 4
pMin = 0.1
pMax = 0.9
nP = 2
alphaMin = 0.1
alphaMax = 0.9
nAlpha = 2
expectedKFactor = {
    (4, 0, 0): 0.0948248,
    (4, 0, 1): 0.159868,
    (4, 0, 2): 0.336537,
    (4, 1, 0): 0.508054,
    (4, 1, 1): 0.852251,
    (4, 1, 2): 1.76562,
    (4, 2, 0): 1.22991,
    (4, 2, 1): 2.03901,
    (4, 2, 2): 4.15708,
    (6, 0, 0): 0.0981609,
    (6, 0, 1): 0.146103,
    (6, 0, 2): 0.246967,
    (6, 1, 0): 0.526396,
    (6, 1, 1): 0.781955,
    (6, 1, 2): 1.31356,
    (6, 2, 0): 1.27862,
    (6, 2, 1): 1.88846,
    (6, 2, 2): 3.14058,
}
for n in range(nMin, nMax):
    for iP in range(nP + 1):
        p = pMin + (pMax - pMin) * iP / nP
        for iAlpha in range(nAlpha + 1):
            alpha = alphaMin + (alphaMax - alphaMin) * iAlpha / nAlpha
            kF = ot.DistFunc.kFactor(2 * n, p, alpha)
            ott.assert_almost_equal(kF, expectedKFactor[(2 * n, iP, iAlpha)], 1e-5, 0.0)
mMin = 1
mMax = 3
expectedKFactorPooled = {
    (4, 1, 0, 0): 0.0948248,
    (4, 1, 0, 1): 0.159868,
    (4, 1, 0, 2): 0.336537,
    (4, 1, 1, 0): 0.508054,
    (4, 1, 1, 1): 0.852251,
    (4, 1, 1, 2): 1.76562,
    (4, 1, 2, 0): 1.22991,
    (4, 1, 2, 1): 2.03901,
    (4, 1, 2, 2): 4.15708,
    (4, 2, 0, 0): 0.101899,
    (4, 2, 0, 1): 0.149319,
    (4, 2, 0, 2): 0.248996,
    (4, 2, 1, 0): 0.546291,
    (4, 2, 1, 1): 0.797623,
    (4, 2, 1, 2): 1.30566,
    (4, 2, 2, 0): 1.32525,
    (4, 2, 2, 1): 1.91438,
    (4, 2, 2, 2): 3.05366,
    (6, 1, 0, 0): 0.0981609,
    (6, 1, 0, 1): 0.146103,
    (6, 1, 0, 2): 0.246967,
    (6, 1, 1, 0): 0.526396,
    (6, 1, 1, 1): 0.781955,
    (6, 1, 1, 2): 1.31356,
    (6, 1, 2, 0): 1.27862,
    (6, 1, 2, 1): 1.88846,
    (6, 1, 2, 2): 3.14058,
    (6, 2, 0, 0): 0.104942,
    (6, 2, 0, 1): 0.140475,
    (6, 2, 0, 2): 0.202791,
    (6, 2, 1, 0): 0.562926,
    (6, 2, 1, 1): 0.752396,
    (6, 2, 1, 2): 1.07853,
    (6, 2, 2, 0): 1.36887,
    (6, 2, 2, 1): 1.82006,
    (6, 2, 2, 2): 2.57114,
}
for n in range(nMin, nMax):
    for m in range(mMin, mMax):
        for iP in range(nP + 1):
            p = pMin + (pMax - pMin) * iP / nP
            for iAlpha in range(nAlpha + 1):
                alpha = alphaMin + (alphaMax - alphaMin) * iAlpha / nAlpha
                kF = ot.DistFunc.kFactorPooled(2 * n, m, p, alpha)
                ott.assert_almost_equal(kF, expectedKFactorPooled[(2 * n, m, iP, iAlpha)], 1e-5, 0.0)

# Non-positive RMax is rejected before bracketing
ot.ResourceMap.SetAsScalar("KFactor-RMax", 0.0)
try:
    with ott.assert_raises(TypeError):
        ot.DistFunc.kFactor(4, 0.5, 0.5)
finally:
    ot.ResourceMap.SetAsScalar("KFactor-RMax", 8.5)

# Non-finite RMax is rejected the same way
ot.ResourceMap.SetAsScalar("KFactor-RMax", float("inf"))
try:
    with ott.assert_raises(TypeError):
        ot.DistFunc.kFactor(4, 0.5, 0.5)
finally:
    ot.ResourceMap.SetAsScalar("KFactor-RMax", 8.5)
