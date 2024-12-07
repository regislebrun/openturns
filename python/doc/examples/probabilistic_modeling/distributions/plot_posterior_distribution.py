"""
Create a bayesian posterior distribution
=========================================
"""

# %%
# In this example, we consider :math:`\inputRV` following a normal bivariate distribution with zero
# mean, unit variance and independent components:
#
# .. math::
#
#      \inputRV \sim \cN(\vect{\theta}_v) \\
#      \vect{\theta}_v = (\mu_0, \sigma_0, \mu_1, \sigma_1, \rho) =  (0,1,0,1,0)
#
# where :math:`\mu_i` are the marginal means, :math:`\sigma_i` the marginal standard
# deviations and :math:`\rho` the correlation coefficient.
#
# We assume that we have a sample of observations of :math:`\inputRV`, denoted by
# :math:`(\vect{x}_1, \dots, \vect{x}_\sampleSize)` whith :math:`\sampleSize=25`.
#
# The obective is to estimate :math:`\vect{\theta}_v` using a Bayesian approach. For that, we
# consider the following probabilistic modelisation:
#
# .. math::
#
#      \inputRV|\vect{\Theta} \sim \cL_{\vect{\Theta}}(\vect{X}) = \cN(\vect{\Theta})
#
# where the parameter vector :math:`\vect{\Theta} = (\Theta_0, \Theta_1, \Theta_2, \Theta_3, \Theta_4)` is
# a random vector such that:
#
# .. math::
#
#      \vect{\Theta} = g(\vect{Y})\\
#      \vect{Y} \sim \cL^0_{\vect{Y}}
#
# We illustrate two cases associated to two link functions and two prior distributions:
#
# - Case 1: the prior distribution of :math:`\vect{Y}` is centered on the real value.
#   We consider the following link function and prior distribution:
#
# .. math::
#
#      g(\vect{y}) = (0, y_0, 0, y_1, 0)\\
#      \cL^0_{\vect{Y}} = \cT(0,1,2) \times \cT(0,1,2)
#
# - Case 2: the prior distribution of :math:`\vect{Y}` is not centered on the real value.
#   We consider the following link function and prior distribution:
#
# .. math::
#
#      \cL^0_{\vect{Y}} = \cT(-1, 0, 1) \times \cT(-1, 0,1)\\
#      g(\vect{y}) =  (0,0.5+y_0^2, 0, 0.5+y_1^2, 0)
#
# In both cases, the link function :math:`g` models that the components are independent and with
# zero mean. As a result, :math:`\Theta_0`, :math:`\Theta_2` and :math:`\Theta_4`
# follow the Dirac distribution centered on the value 0.
#
# Furthermore, in both cases, the conditioned distribution family contains the real distribution of
# :math:`\inputRV`.
#
# The objective is to compute the bayesian posterior distribution :math:`\cL^\sampleSize_{\vect{Y}}` of
# :math:`\vect{Y}|\inputRV = (\vect{x}_1, ..., \vect{x}_\sampleSize)`. From that posterior distribution,
# we get its mode :math:`\vect{y}_\sampleSize^m` and we build the distribution
# :math:`\cL_{\vect{\theta^m}}(\vect{X}) = \cN(\vect{\theta^m})` where
# :math:`\vect{\theta}^m = g(\vect{y}_\sampleSize^m)`.
#
import openturns as ot
import openturns.viewer as viewer
import openturns.experimental as otexp

ot.Log.Show(ot.Log.NONE)

# %%
# Case 1:
# ^^^^^^^
# In that case, we have:
#
# .. math::
#
#      g(\vect{y}) = (0, y_0, 0, y_1, 0)\\
#      \vect{Y} \sim \cT(0,1,2) \times \cT(0,1,2)
#
# In this case, the model is regular and the default value of the integration nodes number is too high.
# We fix it to 16 rather than 256, which is sufficient.
ot.ResourceMap.SetAsUnsignedInteger( "DeconditionedDistribution-MarginalIntegrationNodesNumber", 13)
ot.ResourceMap.SetAsString("DeconditionedDistribution-ContinuousDiscretizationMethod", "GaussProduct")
ot.RandomGenerator.SetSeed(0)
# %%
# Create the observations :math:`(\vect{x}_1, ..., \vect{x}_\sampleSize)`:
N_obs = 25
X_dist_real = ot.Normal(2)
obs = X_dist_real.getSample(N_obs)

# %%
# Create the conditioning distribution :math:`\cL^0_{\vect{Y}}` of :math:`\vect{Y}`:
Y_dist_prior = ot.JointDistribution([ot.Triangular(0.0, 1.0, 2.0)]*2)
Y_dist_prior.setDescription(["Y0", "Y1"])

# %%
# Create the link function :math:`g`:
g = ot.SymbolicFunction(["y0", "y1"], ["0.0", "y0", "0.0", "y1", "0.0"])

# %%
# Create the conditioned distribution  :math:`\cL_{\vect{\Theta}}(\vect{X})` of :math:`\inputRV|\vect{\Theta}`:
X_dist_theta = ot.Normal(2)
X_dist_theta.setDescription(["X0", "X1"])

# %%
# Build the  bayesian posterior distribution :math:`\cL^\sampleSize_{\vect{Y}}
# of  :math:`\vect{Y}|\inputRV = (\vect{x}_1, ..., \vect{x}_\sampleSize)`:
X_dist_decond= ot.DeconditionedDistribution(X_dist_theta, Y_dist_prior, g)
Y_dist_posterior = otexp.PosteriorDistribution(X_dist_decond, obs)

# %%
# From the posterior distribution :math:`\cL^\sampleSize_{\vect{Y}}`, we get the mode denoted by
# :math:`\vect{y}_\sampleSize^m`. We compute :math:`\vect{\theta}^m = g(\vect{y}_\sampleSize^m)`.
#
# To compute the mode of a distribution, we first define the following function:
def computeMode(distribution):
    def obj_py(X):
        return distribution.computeLogPDF(X) * (-1.0)

    obj = ot.PythonFunction(distribution.getDimension(), 1, func_sample=obj_py)
    pb = ot.OptimizationProblem(obj)
    pb.setBounds(distribution.getRange())
    algo = ot.Cobyla(pb)
    algo.setStartingPoint(distribution.getMean())
    algo.run()
    return algo.getResult().getOptimalPoint()

theta_bay = g(computeMode(Y_dist_posterior))
print(f"{theta_bay=}")

# %%
# From the posterior distribution :math:`\cL^\sampleSize_{\vect{Y}}`, we compute the bilateral confidence
# interval of level :math:`\alpha`, such that each marginal interval has the same probability
# :math:`\beta`. Note that if the components are independent, then :math:`\beta = \sqrt{\alpha} = 0.9746`.
# Given the particular expression of the link function, this is exactly the bilateral confidence
# interval of level :math:`\alpha` of :math:`(\Theta_1, \Theta_3)`.
alpha = 0.95
interval_bay, beta = Y_dist_posterior.computeBilateralConfidenceIntervalWithMarginalProbability(alpha)
print("Beta = ", beta)
print("IC Bay = \n", interval_bay)

# %%
# We can estimate the parameter :math:`\vect{\theta}` using the maximum likelihood estimator computed on 
# the observations, denoted by :math:`\vect{\theta}^{ML}`. We fix the parameters :math:`\theta_0`,  
# :math:`\theta_2` and :math:`\theta_4` to 0.
#
# We note that the maximum likelihood estimator of
# :math:`\vect{\theta}` and the bayesian estimator based on the mode of the posterior distribution are both
# very good and the bayesian estimator is better than the maximum likelihood one. It is due to the fact
# that the conditioned distribution family contains the real distribution of :math:`\inputRV` and that
# the prior distribution of :math:`\vect{Y}` is centered on the real value!
mv_factory = ot.NormalFactory()
mv_factory.setKnownParameter([0.0, 0.0, 0.0], [0, 2, 4])
theta_ML = mv_factory.buildEstimator(obs).getDistribution().getParameter()
print(f"{theta_ML = }")

# %%
# We can draw the posterior distribution of :math:`(\Theta_1, \Theta_3)`. Given the particular
# expression of the link function, this is exactly the distribution of :math:`\vect{Y}`.
#
ot.ResourceMap.SetAsString("Contour-DefaultColorMapNorm", "rank")
graph_post = Y_dist_posterior.drawPDF([0.5]*2, [1.5]*2)
theta_ML_cloud = ot.Cloud([theta_ML[[1,3]]])
theta_ML_cloud.setColor("red")
theta_ML_cloud.setPointStyle("bullet")
graph_post.add(theta_ML_cloud)
theta_mod_cloud = ot.Cloud([theta_bay[[1,3]]])
theta_mod_cloud.setColor("red")
theta_mod_cloud.setPointStyle("+")
graph_post.add(theta_mod_cloud)
a = interval_bay.getLowerBound()
b = interval_bay.getUpperBound()
interval_bay_cloud = ot.Curve([a, [a[0], b[1]], b, [b[0], a[1]], a])
interval_bay_cloud.setColor("blue")
graph_post.add(interval_bay_cloud)
graph_post.setLegends(["Post. Bay. dist", r"$\boldsymbol{\theta}^{ML}$", r"$\boldsymbol{\theta}^m$", "IC level " + str(int(100*alpha)) + "%"])
graph_post.setXTitle(r"$\theta_1 = \sigma_0$")
graph_post.setYTitle(r"$\theta_3 = \sigma_1$")
graph_post.setTitle(r"Posterior Bayesian dist of $(\Theta_1, \Theta_3)$")
view = viewer.View(graph_post, square_axes=True)

# %%
# We can evaluate the  Kullback-Leibler divergence between
# the distribution :math:`\cL_{\vect{\theta^m}}(\vect{X}) = \cN(\vect{\theta^m})`
# where :math:`\vect{\theta}^m = g(\vect{y}_\sampleSize^m)`  and
# the real distribution of :math:`\inputRV`,
# The Kullback-Leibler divergence is defined as:
#
# .. math::
#
#    KL(\cL_{\vect{\theta^m}}(\vect{X}), \cL_{\vect{\theta_v}}(\vect{X})) & = \int
#    \mu_{\cL_{\vect{\theta^m}}(\vect{X})}(\vect{x})
#    \log \dfrac{\mu_{\cL_{\vect{\theta^m}}(\vect{X})}(\vect{x})}
#    {\mu_{\cL_{\vect{\theta_v}}(\vect{X})}(\vect{x})} \, d\vect{x}\\
#    & = \Expect{\log \dfrac{\mu_{\cL_{\vect{\theta^m}}(\vect{X})}}
#    {\mu_{\cL_{\vect{\theta_v}}(\vect{X})}(\vect{x})}}
#
# where the expectation is computed with respect to the :math:`\cL_{\vect{\theta^m}}(\vect{X})`
# distribution. We compute it thanks to a Monte Carlo sampling.
# 
# We first create the distribution :math:`\cL_{\vect{\theta^m}}(\vect{X}) = \cN(\vect{\theta^m})`:
X_dist_theta_mod =  ot.Distribution(X_dist_theta)
X_dist_theta_mod.setParameter(theta_bay)

# %%
# We generate a sample from :math:`\cL_{\vect{\theta^m}}(\vect{X})` and estimate the  Kullback-Leibler
# divergence. The divergence is very low, which proves that the estimated distribution is very near
# the real one.
sample = X_dist_theta_mod.getSample(10000)
KL_dist = (X_dist_theta_mod.computeLogPDF(sample) - X_dist_theta.computeLogPDF(sample)).computeMean()[0]
print("KL divergence =%.2e" %KL_dist)


# %%
# In this last graph, we draw the PDF of the real distribution of :math:`\inputRV` and the distribution
# :math:`\cL_{\vect{\theta^m}}(\vect{X}) = \cN(\vect{\theta^m})`. This graph confirms that both
# distributions are very near.
graph_X_dist = X_dist_real.drawPDF()
levels = graph_X_dist.getDrawable(0).getLevels()
dr_bay = X_dist_theta_mod.drawPDF().getDrawable(0).getImplementation()
dr_bay.setLevels(levels)
dr_bay.setColorBarPosition("")
dr_bay.setLineStyle("dashed")
graph_X_dist.add(dr_bay)
obs_cloud = ot.Cloud(obs)
obs_cloud.setColor('red')
graph_X_dist.add(obs_cloud)
graph_X_dist.setLegends(["Real dist","Bay. dist (mode)", "observations"])
graph_X_dist.setXTitle(r"$x_0$")
graph_X_dist.setYTitle(r"$x_1$")
graph_X_dist.setTitle("Distribution of $\mathbf{X}$: real one and bayesian one")
view = viewer.View(graph_X_dist, square_axes=True)


# %%
# Case 2:
# ^^^^^^^
# In that case, we have:
#
# .. math::
#
#      g(\vect{y}) =  (0,0.5+y_0^2, 0, 0.5+y_1^2, 0)\\
#      \cL^0_{\vect{Y}} = \cT(-1, 0, 1) \times \cT(-1, 0,1)
#
# Create the conditioning distribution :math:`\cL^0_{\vect{Y}}` of :math:`\vect{Y}`:
Y_dist_prior = ot.JointDistribution([ot.Triangular(-1.0, 0.0, 1.0)]*2)
Y_dist_prior.setDescription(["Y0", "Y1"])

# %%
# Create the link function :math:`g`:
g = ot.SymbolicFunction(["u0", "u1"], ["0.0", "0.5+u0^2", "0.0", "0.5+u1^2", "0.0"])

# %%
# Create the conditioned distribution  :math:`\cL_{\vect{\Theta}}(\vect{X})` of :math:`\inputRV|\vect{\Theta}`:
X_dist_theta = ot.Normal(2)
X_dist_theta.setDescription(["X0", "X1"])

# %%
# Build the  bayesian posterior distribution :math:`\cL^\sampleSize_{\vect{Y}}`
# of  :math:`\vect{Y}|\inputRV = (\vect{x}_1, ..., \vect{x}_\sampleSize)`:
X_dist_decond= ot.DeconditionedDistribution(X_dist_theta, Y_dist_prior, g)
Y_dist_posterior = otexp.PosteriorDistribution(X_dist_decond, obs)

# %%
# From the posterior distribution :math:`\cL^\sampleSize_{\vect{Y}}`, we get the mode denoted by
# :math:`\vect{y}_\sampleSize^m`. We compute :math:`\vect{\theta}^m = g(\vect{y}_\sampleSize^m)`.
theta_bay = g(computeMode(Y_dist_posterior))
print(f"{theta_bay=}")

# %%
# From the posterior distribution :math:`\cL^\sampleSize_{\vect{Y}}`, we can compute the posterior
# distribution of :math:`(\Theta_1, \Theta_3)`. In this case, the link function
# is not as simple as in the first case. To get the posterior distribution of :math:`(\Theta_1, \Theta_3)`,
# we use a Monte Carlo sampling of size :math:`N` generated by the distribution
# :math:`\cL^\sampleSize_{\vect{Y}}` and we build the distribution with
# the kernel smooting technique.
N = 10000
Theta_post_sample = g(Y_dist_posterior.getSample(N)).getMarginal([1, 3])
Theta_post_dist = ot.KernelSmoothing().build(Theta_post_sample)

# %%
# From the posterior distribution of :math:`(\Theta_1, \Theta_3)`, we compute the bilateral confidence
# interval of level :math:`\alpha`, such that each marginal interval has the same probability
# :math:`\beta`. Note that if the components are independent, then :math:`\beta = \sqrt{\alpha} = 0.9746`.
alpha = 0.95
interval_bay, beta = Theta_post_dist.computeBilateralConfidenceIntervalWithMarginalProbability(alpha)
print("Beta = ", beta)
print("IC Bay = \n", interval_bay)

# %%
# We can estimate the parameter :math:`\vect{\theta}` using the maximum likelihood estimator computed on 
# the observations, denoted by :math:`\vect{\theta}^{ML}`. We fix the parameters :math:`\theta_0`,  
# :math:`\theta_2` and :math:`\theta_4` to 0.
#
# We note that the bayesian estimator based on the mode of the
# posterior distribution is less precise than in the first case and becomes worse than the maximum
# likelihood estimator. It is due to the fact even if  the conditioned distribution family contains the real 
# distribution of :math:`\inputRV`,  the prior distribution of :math:`\vect{Y}` is not centered on
# the real value any  longer!
mv_factory = ot.NormalFactory()
mv_factory.setKnownParameter([0.0, 0.0, 0.0], [0, 2, 4])
theta_ML = mv_factory.buildEstimator(obs).getDistribution().getParameter()
print(f"{theta_ML = }")

# %%
# We can draw the posterior distribution of :math:`(\Theta_1, \Theta_3)`. 
graph_post = Theta_post_dist.drawPDF([0.5]*2, [1.5]*2)
theta_ML_cloud = ot.Cloud([theta_ML[[1,3]]])
theta_ML_cloud.setColor("red")
theta_ML_cloud.setPointStyle("bullet")
graph_post.add(theta_ML_cloud)
theta_mod_cloud = ot.Cloud([theta_bay[[1,3]]])
theta_mod_cloud.setColor("red")
theta_mod_cloud.setPointStyle("+")
graph_post.add(theta_mod_cloud)
a = interval_bay.getLowerBound()
b = interval_bay.getUpperBound()
interval_bay_cloud = ot.Curve([a, [a[0], b[1]], b, [b[0], a[1]], a])
interval_bay_cloud.setColor("blue")
graph_post.add(interval_bay_cloud)
graph_post.setLegends(["Post. Bay. dist", r"$\mathbf{\theta}^{ML}$", r"$\mathbf{\theta}^m$", "IC level " + str(int(100*alpha)) + "%"])
graph_post.setXTitle(r"$\theta_1 = \sigma_0$")
graph_post.setYTitle(r"$\theta_3 = \sigma_1$")
graph_post.setTitle(r"Posterior Bayesian dist of $(\Theta_1, \Theta_3)$")
view = viewer.View(graph_post, square_axes=True)

# %%
# We can evaluate the  Kullback-Leibler divergence between
# the distribution :math:`\cL_{\vect{\theta^m}}(\vect{X}) = \cN(\vect{\theta^m})`
# where :math:`\vect{\theta}^m = g(\vect{y}_\sampleSize^m)` and
# the real distribution of :math:`\inputRV`.
# 
# We first create the distribution :math:`\cL_{\vect{\theta^m}}(\vect{X}) = \cN(\vect{\theta^m})`:
X_dist_theta_mod =  ot.Distribution(X_dist_theta)
X_dist_theta_mod.setParameter(theta_bay)

# %%
# We generate a sample from :math:`\cL_{\vect{\theta^m}}(\vect{X})` and estimate the  Kullback-Leibler
# divergence. We note that the divergence is greater in that case than in the first case.
sample = X_dist_theta_mod.getSample(10000)
KL_dist = (X_dist_theta_mod.computeLogPDF(sample) - X_dist_theta.computeLogPDF(sample)).computeMean()[0]
print("KL divergence =%.2e" %KL_dist)


# %%
# In this graph, we draw the PDF of the real distribution of :math:`\inputRV` and the distribution
# :math:`\cL_{\vect{\theta^m}}(\vect{X}) = \cN(\vect{\theta^m})`. This graph confirms that both
# distributions are very near, even if this second modelisation appears to be less good than the first one.
graph_X_dist = X_dist_real.drawPDF()
levels = graph_X_dist.getDrawable(0).getLevels()
dr_bay = X_dist_theta_mod.drawPDF().getDrawable(0).getImplementation()
dr_bay.setLevels(levels)
dr_bay.setColorBarPosition("")
dr_bay.setLineStyle("dashed")
graph_X_dist.add(dr_bay)
obs_cloud = ot.Cloud(obs)
obs_cloud.setColor('red')
graph_X_dist.add(obs_cloud)
graph_X_dist.setLegends(["Real dist","Bay. dist (mode)", "observations"])
graph_X_dist.setXTitle(r"$x_0$")
graph_X_dist.setYTitle(r"$x_1$")
graph_X_dist.setTitle("Distribution of $\mathbf{X}$: real one and bayesian one")
view = viewer.View(graph_X_dist, square_axes=True)


# %%
# In this graph, we draw the posterior distribution  :math:`\cL^\sampleSize_{\vect{Y}}`
# and the prior distribution :math:`\cL^0_{\vect{Y}}` with a sample generated from the posterior
# distribution. We note that the posterior distribution is quadri-modale, due to the parity of the link
# function.
graph_Y_dist_post = Y_dist_posterior.drawPDF()
graph_Y_dist_post.setXTitle(r"$y_0$")
graph_Y_dist_post.setYTitle(r"$y_1$")
graph_Y_dist_post.setTitle("Bayesian posterior distribution of $\mathbf{Y}$")
graph_Y_dist_post.add(ot.Cloud(Y_dist_posterior.getSample(100)))
view = viewer.View(graph_Y_dist_post, square_axes=True)


# %%
# Here we drw the PDF of the priori and the posterior  distribution  of :math:`\vect{Y}`.
graph_Y_dist_post = Y_dist_posterior.drawPDF()
levels = graph_Y_dist_post.getDrawable(0).getLevels()
dr_prior = Y_dist_prior.drawPDF().getDrawable(0).getImplementation()
dr_prior.setLevels(levels)
dr_prior.setColorBarPosition("")
dr_prior.setLineStyle("dashed")
graph_Y_dist_post.add(dr_prior)
graph_Y_dist_post.setXTitle(r"$y_0$")
graph_Y_dist_post.setYTitle(r"$y_1$")
graph_Y_dist_post.setTitle("Bayesian prior and posterior distribution of $\mathbf{Y}$")
view = viewer.View(graph_Y_dist_post, square_axes=True)

# %%
# Display all figures
viewer.View.ShowAll()

# %%
# Reset default settings
ot.ResourceMap.Reload()




