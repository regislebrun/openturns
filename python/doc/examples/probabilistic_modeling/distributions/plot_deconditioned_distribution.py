"""
Create a deconditioned distribution
===================================
"""

# %%
# In this example, we consider :math:`\inputRV` following a distribution whose parameters
# are random variables:
#
# .. math::
#
#      \inputRV \sim  \cL_{\vect{X}}(\vect{\Theta}) \\
#      \vect{\Theta} = g(\vect{Y})\\
#      \vect{Y} \sim \cL_{\vect{Y}}
#
# The objective is to estimate the deconditioned distribution of :math:`\inputRV`.
#
# We illustrate two cases with different regularity of the model:
#
# - Case 1: the support of the distribution of :math:`\inputRV|\vect{\Theta} = \vect{\theta}` does not
#   depend on :math:`\vect{\theta}`.
#   We consider the case of a normal distribution whose mean and standard deviation are random variables:
#
#   .. math::
#
#      \cL_{\vect{X}}(\vect{\Theta}) = \cN(\vect{\Theta})\\
#      \vect{\Theta} = g(\vect{Y})\\
#      g(y) = (y^2, y)\\
#      \cL_{Y} = \cU(1,2)
#
# - Case 2: the support of the distribution of :math:`\inputRV|\vect{\Theta} = \vect{\theta}` depends on
#   :math:`\vect{\theta}`.
#   We consider the case of a uniform distribution whose upper bound :math:`B` is a random variable
#   and whose lower bound is fixed to 0:
#
#   .. math::
#
#      \cL_{\vect{X}}(\vect{\Theta}) = \cU(\vect{\Theta})\\
#      \vect{\Theta} = (A, B) = g(\vect{Y})\\
#      g(y) = (0, 1+y^2) \\
#      \cL_{Y} = \cU(0,1)
#
#   In that case, we randomize  the upper bound only and not the lower bound. As a result,
#   :math:`A` follows the Dirac distribution centered on 0.
#
import openturns as ot
import openturns.viewer as viewer
ot.RandomGenerator.SetSeed(0)

ot.Log.Show(ot.Log.NONE)

# %%.
# Case 1: Regular model
# ^^^^^^^^^^^^^^^^^^^^^
# In that case, we have:
#
# .. math::
#
#      \cL_{\vect{X}}(\vect{\Theta}) = \cN(\vect{\Theta})\\
#      \vect{\Theta} = g(\vect{Y})\\
#      g(y) = (y^2, y)\\
#      \cL_{Y} = \cU(1,2)
#
# The exact deconditioned distribution of :math:`\inputRV` is defined by:
#
# .. math::
#
#   p_{\vect{Y}}(y) = \dfrac{1}{\sqrt{2\pi}} \int_1^2 \dfrac{1}{y}
#            e^{-\dfrac{1}{2}\dfrac{(x-y^2)^2}{y^2}}\, dy
#
# Create the conditioning distribution :math:`\cL_{\vect{Y}}` of :math:`\vect{Y}`:
Y_dist = ot.Uniform(1.0, 2.0)
Y_dist.setDescription(["Y"])

# %%
# Create the link function :math:`g`:
g = ot.SymbolicFunction(["y"], ["y^2", "y"])

# %%
# Create the conditioned distribution  :math:`\cL_{\vect{\Theta}}(\vect{X})`
# of :math:`\inputRV|\vect{\Theta}`:
X_dist_theta = ot.Normal()
X_dist_theta.setDescription(["X"])

# %%
# In this case, the model is regular: the support of the distribution of
# :math:`\inputRV|\vect{\Theta} = \vect{\theta}` does not depend on :math:`\vect{\theta}`.
# We do not need a high number of intergration  nodes to evaluate the deconditioned distribution.
# We fix it to 16 rather than 256 which is the default value.
ot.ResourceMap.SetAsUnsignedInteger( "DeconditionedDistribution-MarginalIntegrationNodesNumber", 16)
ot.ResourceMap.SetAsString("DeconditionedDistribution-ContinuousDiscretizationMethod", "GaussProduct")

# %%
# Build the  deconditioned distribution of :math:`\inputRV`:
X_dist_decond= ot.DeconditionedDistribution(X_dist_theta, Y_dist, g)

# %%
# We define the exact PDF:
formula = '1/sqrt(2*pi_)*(1/y)*exp(-(x-y^2)^2/(2*y^2))'
case1_jointPDF = ot.SymbolicFunction(['x', 'y'], [formula])

def case1_exact_PDF_py(x):
    integrand = ot.ParametricFunction(case1_jointPDF, [0], x)
    return ot.GaussKronrod().integrate(integrand, Y_dist.getRange())

case1_exact_PDF = ot.PythonFunction(1, 1, case1_exact_PDF_py)

# %%
# In order to compute the Kullback-Leibler divergence, we need to compute the exact log-PDF.
case1_exact_LogPDF = ot.ComposedFunction(ot.SymbolicFunction('x', 'log(x)'), case1_exact_PDF)

# %%
# We draw both the exact PDF and the the deconditioned one.
graph_X_dist = X_dist_decond.drawPDF()
xMin = graph_X_dist.getDrawable(0).getData().getMin()[0]
xMax = graph_X_dist.getDrawable(0).getData().getMax()[0]
dr_exact = case1_exact_PDF.draw(xMin, xMax).getDrawable(0)
dr_exact.setLineStyle('dashed')
dr_exact.setLineWidth(2)
graph_X_dist.add(dr_exact)
graph_X_dist.setLegends(['Gauss-Product (16)', 'exact'])
graph_X_dist.setXTitle('x')
graph_X_dist.setTitle(r'Deconditioned distribution of $\mathbb{X}$')
view = viewer.View(graph_X_dist)

# %%
# We can evaluate the  Kullback-Leibler divergence between
# the deconditioned distribution :math:`\cL_d(\vect{X})` and
# the exact distribution :math:`\cL_e(\vect{X})` of :math:`\inputRV`,
# The Kullback-Leibler divergence is defined as:
#
# .. math::
#
#    KL(\cL_d(\vect{X}), \cL_e(\vect{X})) & = \int
#    \mu_{\cL_d(\vect{X})}(\vect{x})
#    \log \dfrac{\mu_{\cL_d(\vect{X})}(\vect{x})}
#    {\mu_{\cL_e(\vect{X})}(\vect{x})} \, d\vect{x}\\
#    & = \Expect{ \log \dfrac{\mu_{\cL_d(\vect{X})}(\vect{X})}{\mu_{\cL_e(\vect{X})}(\vect{X})} }
#
# where the expectation is computed with respect to the :math:`\cL_d(\vect{X})`
# distribution. We estimate it thanks to a Monte Carlo sampling.
#
# We generate a sample from :math:`\cL_d(\vect{X})` and estimate the  Kullback-Leibler
# divergence. The divergence is very low:  we have reached the machine precision with the
# Gauss-Product (GP) method with 16 points. It
# proves that the deconditioned distribution is very near
# the exact one.
sample = X_dist_decond.getSample(10000)
KL_dist = (X_dist_decond.computeLogPDF(sample) - case1_exact_LogPDF(sample)).computeMean()[0]
print("KL divergence GP (16) =%.2e" %KL_dist)

# %%
# At last, we illustrate the influence of the discretization method used to calculate the integrals.
# The Quasi Monte Carlo
# (QMC) sampling method leads to a result much less accurate: the Kullback-Leibler divergence is much
# greater with the QMC method than with the GP one, even with a high number of points (here 256).
ot.ResourceMap.SetAsString("DeconditionedDistribution-ContinuousDiscretizationMethod", "QMC")
ot.ResourceMap.SetAsUnsignedInteger('DeconditionedDistribution-MaximumIntegrationNodesNumber', 2**8)
X_dist_decond= ot.DeconditionedDistribution(X_dist_theta, Y_dist, g)
sample = X_dist_decond.getSample(10000)
KL_dist = (X_dist_decond.computeLogPDF(sample) - case1_exact_LogPDF(sample)).computeMean()[0]
print("KL divergence QMC (256) =%.2e" %KL_dist)

# %%
# Case 2: Irregular model
# ^^^^^^^^^^^^^^^^^^^^^^^
# In that case, we have:
#
# .. math::
#
#      \cL_{\vect{X}}(\vect{\Theta}) = \cU(\vect{\Theta})\\
#      \vect{\Theta} = (A, B) = g(\vect{Y})\\
#      g(y) = (0, 1+y^2) \\
#      \cL_{Y} = \cU(0,1)
#
# The exact deconditioned distribution of :math:`\inputRV` is defined by:
#
# .. math::
#
#     p_{\vect{X}}(x) = \begin{array}{|ll}
#       0              & \mbox{ if } x \leq 0 \\
#       \pi / 4 &  \mbox{ if }0 \leq x \leq 1\\
#       \pi / 4 - \arctan \sqrt{x-1} &  \mbox{ if } 1 \leq x \leq 2\\
#       0              &  \mbox{ if }  x \geq 2
#       \end{array}
#
# Create the conditioning distribution :math:`\cL_{\vect{Y}}` of :math:`\vect{Y}`:
Y_dist = ot.Uniform(0.0, 1.0)
Y_dist.setDescription(["Y"])

# %%
# Create the link function :math:`g`:
g = ot.SymbolicFunction(["y"], ["0.0", "1+y^2"])

# %%
# Create the conditioned distribution  :math:`\cL_{\vect{\Theta}}(\vect{X})` of :math:`\inputRV|\vect{\Theta}`:
X_dist_theta = ot.Uniform()
X_dist_theta.setDescription(["X"])

# %%
# We define the exact PDF:
formula = 'x<0? 0: x< 1? pi_/4: x<2? pi_/4 - atan(sqrt(x-1)):0'
case2_exact_PDF = ot.SymbolicFunction(['x'], [formula])

# %%
# In order to compute the Kullback-Leibler divergence, we need to compute the exact log-PDF.
case2_exact_LogPDF = ot.ComposedFunction(ot.SymbolicFunction('x', 'log(x)'), case2_exact_PDF)

# %%
# We use the Gauss-Product discretization method whose integration nodes number is equal to 16.
ot.ResourceMap.SetAsUnsignedInteger( "DeconditionedDistribution-MarginalIntegrationNodesNumber", 16)
ot.ResourceMap.SetAsString("DeconditionedDistribution-ContinuousDiscretizationMethod", "GaussProduct")

# %%
# Build the  deconditioned distribution of :math:`\inputRV`:
X_dist_decond_GP16 = ot.DeconditionedDistribution(X_dist_theta, Y_dist, g)

# %%
# We draw both the exact PDF and the deconditioned one. We note that the deconditioned
# distribution is not accurate, because of the Gauss-Product discretization method which is not adapted
# to this case.
graph_X_dist = X_dist_decond_GP16.drawPDF()
xMin = graph_X_dist.getDrawable(0).getData().getMin()[0]
xMax = graph_X_dist.getDrawable(0).getData().getMax()[0]
dr_exact_case2 = case2_exact_PDF.draw(xMin, xMax).getDrawable(0)
dr_exact_case2.setLineStyle('dashed')
dr_exact_case2.setLineWidth(2)
graph_X_dist.add(dr_exact_case2)
graph_X_dist.setLegends(['Gauss-Product (16)', 'exact'])
graph_X_dist.setXTitle('x')
graph_X_dist.setTitle(r'Decondiitoned distribution of $\mathbb{X}$')
view = viewer.View(graph_X_dist)

# %%
# We can evaluate the  Kullback-Leibler divergence between
# the deconditioned distribution :math:`\cL_d(\vect{X})` and
# the exact distribution :math:`\cL_e(\vect{X})`.
# We generate a sample from :math:`\cL_d(\vect{X})` and estimate the  Kullback-Leibler
# divergence. The result confirms the bad precision of the deconditioned distribution.
sample = X_dist_decond_GP16.getSample(10000)
KL_dist = (X_dist_decond.computeLogPDF(sample) - case2_exact_LogPDF(sample)).computeMean()[0]
print("KL divergence Gauss-Product (16) =%.2e" %KL_dist)

# %%
# If we increase the integration nodes number to 256 and to 4096, then the Kullback-Leibler divergence
# remains roughly constant of order :math:`1e-5`, even if the graphs of the pdf are superimposed.
ot.ResourceMap.SetAsUnsignedInteger( "DeconditionedDistribution-MarginalIntegrationNodesNumber", 256)
ot.ResourceMap.SetAsString("DeconditionedDistribution-ContinuousDiscretizationMethod", "GaussProduct")
X_dist_decond_GP256= ot.DeconditionedDistribution(X_dist_theta, Y_dist, g)
sample = X_dist_decond_GP256.getSample(10000)
KL_dist = (X_dist_decond_GP256.computeLogPDF(sample) - case2_exact_LogPDF(sample)).computeMean()[0]
print("KL divergence Gauss-Product (256) =%.2e" %KL_dist)

ot.ResourceMap.SetAsUnsignedInteger( "DeconditionedDistribution-MarginalIntegrationNodesNumber", 2**12)
ot.ResourceMap.SetAsString("DeconditionedDistribution-ContinuousDiscretizationMethod", "GaussProduct")
X_dist_decond_GP4096= ot.DeconditionedDistribution(X_dist_theta, Y_dist, g)
sample = X_dist_decond_GP4096.getSample(10000)
KL_dist = (X_dist_decond_GP4096.computeLogPDF(sample) - case2_exact_LogPDF(sample)).computeMean()[0]
print("KL divergence Gauss-Product (4096) =%.2e" %KL_dist)

# %%
# The same experiment with the QMC method rather than the GP method proves that:
#
# - QMC improves the Kullback-Leibler divergence when the sampling size increases,
# - the Kullback-Leibler divergences are better than those obtained with the Gauss-Product method.
#
ot.ResourceMap.SetAsString("DeconditionedDistribution-ContinuousDiscretizationMethod", "QMC")
ot.ResourceMap.SetAsUnsignedInteger('DeconditionedDistribution-MaximumIntegrationNodesNumber', 2**8)
X_dist_decond_QMC256= ot.DeconditionedDistribution(X_dist_theta, Y_dist, g)
sample = X_dist_decond_QMC256.getSample(10000)
KL_dist = (X_dist_decond_QMC256.computeLogPDF(sample) - case2_exact_LogPDF(sample)).computeMean()[0]
print("KL divergence QMC (256) =%.2e" %KL_dist)

ot.ResourceMap.SetAsString("DeconditionedDistribution-ContinuousDiscretizationMethod", "QMC")
ot.ResourceMap.SetAsUnsignedInteger('DeconditionedDistribution-MaximumIntegrationNodesNumber', 2**12)
X_dist_decond_QMC4096= ot.DeconditionedDistribution(X_dist_theta, Y_dist, g)
sample = X_dist_decond_QMC4096.getSample(10000)
KL_dist = (X_dist_decond_QMC4096.computeLogPDF(sample) - case2_exact_LogPDF(sample)).computeMean()[0]
print("KL divergence QMC (4096) =%.2e" %KL_dist)

# %% 
# We compare the PDF built with the QMC method and the Gauss-Product method, each
# of them with 256
# points. There is no visible deviation, even of the Kullback-Leibler divergences in both cases are
# different.
# sphinx_gallery_thumbnail_number = 3
graph_X_dist_2 = X_dist_decond_QMC256.drawPDF()
graph_X_dist_2.add(X_dist_decond_GP256.drawPDF())
graph_X_dist_2.add(dr_exact_case2)

graph_X_dist_2.setLegends(['QMC (256)', 'Gauss-Product (256)', 'exact'])
graph_X_dist_2.setTitle(r'Deconditioned distribution of $\mathbb{X}$')
graph_X_dist_2.setXTitle('x')
view = viewer.View(graph_X_dist_2)

# %%
# In conlusion, if the model is known to be regular, the Gauss-Product discretization is to be
# preferred, while
# if the model is known to be irregular or of its regularity is unknown, it is safer to use a QMC
# discretization with large enough number of points.


# %%
# Display all figures:
viewer.View.ShowAll()

# %%
# Reset default settings:
ot.ResourceMap.Reload()
