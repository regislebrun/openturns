"""
Create a joint distribution by conditioning
===========================================
"""

# %%
# In this example, we consider :math:`X` such that :math:`X|\vect{\Theta}` follows the
# distribution :math:`\cL_{\vect{\Theta}}(X) = \cN(\vect{\Theta})`, where
# :math:`\vect{\Theta}` is bivariate and is the couple of the
# random mean and the random standard deviation of the normal distribution. 
#
# We assume that :math:`\vect{\Theta} = g(Y)` where :math:`Y` is scalar and follows 
# :math:`\cN(0,1)`and :math:`g: \Rset \rightarrow \Rset^2` is the  link function such that
# :math:`g(y) = (y, 0.1+y^2)`.
#
# The objective is to build the joint distribution of :math:`(\vect{Y}, X)`.

# %%
import openturns as ot
import openturns.viewer as viewer

ot.Log.Show(ot.Log.NONE)

# %%
# Create the conditioning distribution :math:`Y`:
YDist = ot.Normal(0.0, 1.0)

# %%
# Create the link function :math:`g`:
g = ot.SymbolicFunction(["y"], ["y", "0.1 + y^2"])

# %%
# Create the conditioned distribution :math:`X|\vect{\Theta}`:
XgivenThetaDist = ot.Normal()

# %%
# Create the joint distribution of :math:`(\vect{Y}, X)`:
YXDist = ot.JointByConditioningDistribution(XgivenThetaDist, YDist, g)
YXDist.setDescription(["y", "x"])
print(YXDist)

# %%
# Get a sample from the joint distribution:
sample = YXDist.getSample(100)

# %%
# Draw ist PDF with the sample:
ot.ResourceMap.SetAsString("Contour-DefaultColorMapNorm", "rank")
yxMin = sample.getMin()
yxMax = sample.getMax()
delta = 0.05 * (yxMax - yxMin)
graph = YXDist.drawPDF(yxMin - delta, yxMax + delta, [256]*2)
cloud = ot.Cloud(sample)
cloud.setLegend("sample")
graph.add(cloud)
view = viewer.View(graph)

# %%
# Display all figures
viewer.View.ShowAll()

# %%
ot.ResourceMap.Reload()
