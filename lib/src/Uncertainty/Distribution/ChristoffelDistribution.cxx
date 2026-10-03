//                                               -*- C++ -*-
/**
 *  @brief The Christoffel distribution associated to an orthonormal basis
 *
 *  Copyright 2005-2026 Airbus-EDF-IMACS-ONERA-Phimeca
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <cmath>
#include <algorithm>
#include "openturns/ChristoffelDistribution.hxx"
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/ResourceMap.hxx"
#include "openturns/SpecFunc.hxx"
#include "openturns/Exception.hxx"
#include "openturns/OptimizationAlgorithm.hxx"
#include "openturns/RandomGenerator.hxx"
#include "openturns/TBBImplementation.hxx"
#include "openturns/OrthogonalProductPolynomialFactory.hxx"
#include "openturns/OrthogonalProductFunctionFactory.hxx"
#include "openturns/UniVariateFunctionFamily.hxx"
#include "openturns/Uniform.hxx"
#include "openturns/GaussKronrod.hxx"

BEGIN_NAMESPACE_OPENTURNS

/* Unnormalized tensor slice t -> mu_j(t) * S_j(t) for 1D quadrature */
class ChristoffelTensorSliceEvaluation: public EvaluationImplementation
{
public:
  ChristoffelTensorSliceEvaluation(const ChristoffelDistribution * p_distribution,
                                   const Point & partialProducts,
                                   const UnsignedInteger index)
    : EvaluationImplementation()
    , p_distribution_(p_distribution)
    , partialProducts_(partialProducts)
    , index_(index)
  {
    // Nothing to do
  }

  ChristoffelTensorSliceEvaluation * clone() const override
  {
    return new ChristoffelTensorSliceEvaluation(*this);
  }

  Point operator()(const Point & point) const override
  {
    const Distribution marginal(p_distribution_->getMeasure().getMarginal(index_));
    const Scalar density = marginal.computePDF(Point(1, point[0]));
    if (density == 0.0) return Point(1, 0.0);
    return Point(1, density * p_distribution_->computePartialChristoffel(partialProducts_, index_, point[0]));
  }

  UnsignedInteger getInputDimension() const override
  {
    return 1;
  }

  UnsignedInteger getOutputDimension() const override
  {
    return 1;
  }

private:
  const ChristoffelDistribution * p_distribution_;
  Point partialProducts_;
  UnsignedInteger index_;
}; // class ChristoffelTensorSliceEvaluation

CLASSNAMEINIT(ChristoffelDistribution)

static const Factory<ChristoffelDistribution> Factory_ChristoffelDistribution;

/* Default constructor: 1D Legendre basis of size 1, a valid state for plotting */
ChristoffelDistribution::ChristoffelDistribution()
  : DistributionImplementation()
  , orthogonalBasis_()
  , size_(0)
{
  setName("ChristoffelDistribution");
  const OrthogonalProductPolynomialFactory factory(Collection<Distribution>(1, Uniform(-1.0, 1.0)));
  setOrthogonalBasis(OrthogonalBasis(factory), 1);
}

/* Parameters constructor */
ChristoffelDistribution::ChristoffelDistribution(const OrthogonalBasis & basis,
    const UnsignedInteger size)
  : DistributionImplementation()
  , orthogonalBasis_(basis)
  , size_(size)
{
  setName("ChristoffelDistribution");
  if (size == 0) throw InvalidArgumentException(HERE) << "Error: expected a positive space dimension, here size=" << size;
  update();
}

/* Virtual constructor */
ChristoffelDistribution * ChristoffelDistribution::clone() const
{
  return new ChristoffelDistribution(*this);
}

/* Comparison operator: measure and size plus basis values on range probes,
   as Basis has no value-based equality */
Bool ChristoffelDistribution::operator ==(const ChristoffelDistribution & other) const
{
  if (this == &other) return true;
  if (size_ != other.size_) return false;
  if (getDimension() != other.getDimension()) return false;
  if (!(measure_ == other.measure_)) return false;
  const Interval range(getRange());
  Sample probes(0, getDimension());
  probes.add(range.getLowerBound());
  probes.add(range.getUpperBound());
  probes.add((range.getLowerBound() + range.getUpperBound()) * 0.5);
  for (UnsignedInteger p = 0; p < probes.getSize(); ++p)
  {
    const Point probe(probes[p]);
    for (UnsignedInteger i = 0; i < size_; ++i)
      if (basis_[i](probe)[0] != other.basis_[i](probe)[0]) return false;
  }
  return true;
}

Bool ChristoffelDistribution::equals(const DistributionImplementation & other) const
{
  const ChristoffelDistribution * p_other = dynamic_cast<const ChristoffelDistribution *>(&other);
  return p_other && (*this == *p_other);
}

/* String converter */
String ChristoffelDistribution::__repr__() const
{
  OSS oss;
  oss << "class=" << ChristoffelDistribution::GetClassName()
      << " name=" << getName()
      << " dimension=" << getDimension()
      << " size=" << size_
      << " measure=" << measure_;
  return oss;
}

String ChristoffelDistribution::__str__(const String & offset) const
{
  OSS oss;
  oss << getClassName() << "(size = " << size_ << ", measure = " << measure_.__str__(offset) << ")";
  return oss;
}

/* Orthogonal basis and size accessor */
void ChristoffelDistribution::setOrthogonalBasis(const OrthogonalBasis & basis,
    const UnsignedInteger size)
{
  if (size == 0) throw InvalidArgumentException(HERE) << "Error: expected a positive space dimension, here size=" << size;
  orthogonalBasis_ = basis;
  size_ = size;
  update();
}

OrthogonalBasis ChristoffelDistribution::getOrthogonalBasis() const
{
  return orthogonalBasis_;
}

/* Size accessor */
UnsignedInteger ChristoffelDistribution::getSize() const
{
  return size_;
}

/* Reference measure accessor */
Distribution ChristoffelDistribution::getMeasure() const
{
  return measure_;
}

/* Finite basis accessor */
Basis ChristoffelDistribution::getBasis() const
{
  return basis_;
}

/* Christoffel function evaluation on a point */
Scalar ChristoffelDistribution::computeChristoffel(const Point & point) const
{
  if (point.getDimension() != getDimension()) throw InvalidArgumentException(HERE) << "Error: expected a point of dimension=" << getDimension() << ", got dimension=" << point.getDimension();
  Scalar kn = 0.0;
  for (UnsignedInteger i = 0; i < size_; ++i)
  {
    const Scalar value = basis_[i](point)[0];
    kn += value * value;
  }
  return kn;
}

struct ComputeChristoffelPolicy
{
  const Sample & input_;
  Sample & output_;
  const ChristoffelDistribution & distribution_;

  ComputeChristoffelPolicy(const Sample & input,
                           Sample & output,
                           const ChristoffelDistribution & distribution)
    : input_(input)
    , output_(output)
    , distribution_(distribution)
  {
    // Nothing to do
  }

  inline void operator()(const TBBImplementation::BlockedRange<UnsignedInteger> & r) const
  {
    for (UnsignedInteger i = r.begin(); i != r.end(); ++i) output_(i, 0) = distribution_.computeChristoffel(input_[i]);
  }
}; /* end struct ComputeChristoffelPolicy */

/* Christoffel function evaluation on a sample, parallel over points when large */
Sample ChristoffelDistribution::computeChristoffel(const Sample & sample) const
{
  if (sample.getDimension() != getDimension()) throw InvalidArgumentException(HERE) << "Error: expected a sample of dimension=" << getDimension() << ", got dimension=" << sample.getDimension();
  const UnsignedInteger sampleSize = sample.getSize();
  Sample values(sampleSize, 1);
  const ComputeChristoffelPolicy policy(sample, values, *this);
  TBBImplementation::ParallelForIf(isParallel() && sampleSize > 2048, 0, sampleSize, policy, 1024);
  return values;
}

/* Stability factor estimate */
Scalar ChristoffelDistribution::computeKn() const
{
  if (!isAlreadyComputedKn_)
  {
    const UnsignedInteger samplingSize = ResourceMap::GetAsUnsignedInteger("ChristoffelDistribution-KnSamplingSize");
    const Sample sample(measure_.getSample(samplingSize));
    const Sample values(computeChristoffel(sample));
    Scalar knMax = 0.0;
    for (UnsignedInteger i = 0; i < samplingSize; ++i) knMax = std::max(knMax, values(i, 0));
    const Scalar safety = ResourceMap::GetAsScalar("ChristoffelDistribution-KnSafetyFactor");
    knEstimate_ = knMax * safety;
    isAlreadyComputedKn_ = true;
  }
  return knEstimate_;
}

/* PDF evaluation */
Scalar ChristoffelDistribution::computePDF(const Point & point) const
{
  if (size_ == 0) throw InvalidArgumentException(HERE) << "Error: cannot evaluate the PDF of an uninitialized distribution.";
  if (point.getDimension() != getDimension()) throw InvalidArgumentException(HERE) << "Error: expected a point of dimension=" << getDimension() << ", here dimension=" << point.getDimension();
  const Scalar referencePDF = measure_.computePDF(point);
  if (referencePDF == 0.0) return 0.0;
  return referencePDF * computeChristoffel(point) / size_;
}

/* Log-PDF evaluation */
Scalar ChristoffelDistribution::computeLogPDF(const Point & point) const
{
  if (size_ == 0) throw InvalidArgumentException(HERE) << "Error: cannot evaluate the log-PDF of an uninitialized distribution.";
  if (point.getDimension() != getDimension()) throw InvalidArgumentException(HERE) << "Error: expected a point of dimension=" << getDimension() << ", here dimension=" << point.getDimension();
  const Scalar referenceLogPDF = measure_.computeLogPDF(point);
  const Scalar logKn = computeLogChristoffel(point);
  if (!(logKn > SpecFunc::LowestScalar)) return SpecFunc::LowestScalar;
  return referenceLogPDF + logKn - std::log(1.0 * size_);
}

/* Stable log of the Christoffel function: factor out the largest square */
Scalar ChristoffelDistribution::computeLogChristoffel(const Point & point) const
{
  if (point.getDimension() != getDimension()) throw InvalidArgumentException(HERE) << "Error: expected a point of dimension=" << getDimension() << ", here dimension=" << point.getDimension();
  Point values(size_);
  Scalar maxAbs = 0.0;
  UnsignedInteger maxIndex = 0;
  for (UnsignedInteger j = 0; j < size_; ++j)
  {
    const Scalar value = basis_[j](point)[0];
    values[j] = value;
    const Scalar absValue = std::abs(value);
    if (absValue > maxAbs)
    {
      maxAbs = absValue;
      maxIndex = j;
    }
  }
  if (!(maxAbs > 0.0)) return SpecFunc::LowestScalar;
  Scalar reduced = 0.0;
  for (UnsignedInteger j = 0; j < size_; ++j)
    if (j != maxIndex)
    {
      const Scalar ratio = values[j] / maxAbs;
      reduced += ratio * ratio;
    }
  return 2.0 * std::log(maxAbs) + std::log1p(reduced);
}

/* Conditional PDF of Xj | X0..Xj-1, closed slice form in the tensor case */
Scalar ChristoffelDistribution::computeConditionalPDF(const Scalar x,
    const Point & y) const
{
  const UnsignedInteger j = y.getDimension();
  if (j >= getDimension()) throw InvalidArgumentException(HERE) << "Error: cannot compute a conditional PDF with a conditioning point of dimension greater or equal to the distribution dimension.";
  if (!isTensorProduct_ || !measure_.hasIndependentCopula())
    return DistributionImplementation::computeConditionalPDF(x, y);
  // Partial products from the conditioning values, normalizer by orthonormality
  Point partial(size_, 1.0);
  for (UnsignedInteger basisIndex = 0; basisIndex < size_; ++basisIndex)
    for (UnsignedInteger k = 0; k < j; ++k)
    {
      const Scalar factor = computeTensorFactor(basisIndex, k, y[k]);
      partial[basisIndex] *= factor * factor;
    }
  Scalar normalizer = 0.0;
  for (UnsignedInteger basisIndex = 0; basisIndex < size_; ++basisIndex) normalizer += partial[basisIndex];
  if (!(normalizer > 0.0)) return 0.0;
  const Distribution marginalJ(measure_.getMarginal(j));
  const Scalar numerator = marginalJ.computePDF(x) * computePartialChristoffel(partial, j, x);
  return numerator / normalizer;
}

/* Conditional CDF of Xj | X0..Xj-1, slice quadrature in the tensor case */
Scalar ChristoffelDistribution::computeConditionalCDF(const Scalar x,
    const Point & y) const
{
  const UnsignedInteger j = y.getDimension();
  if (j >= getDimension()) throw InvalidArgumentException(HERE) << "Error: cannot compute a conditional CDF with a conditioning point of dimension greater or equal to the distribution dimension.";
  if (!isTensorProduct_ || !measure_.hasIndependentCopula())
    return DistributionImplementation::computeConditionalCDF(x, y);
  const Distribution marginalJ(measure_.getMarginal(j));
  const Interval rangeJ(marginalJ.getRange());
  const Scalar lower = rangeJ.getLowerBound()[0];
  const Scalar upper = rangeJ.getUpperBound()[0];
  if (x <= lower) return 0.0;
  if (x >= upper) return 1.0;
  Point partial(size_, 1.0);
  for (UnsignedInteger basisIndex = 0; basisIndex < size_; ++basisIndex)
    for (UnsignedInteger k = 0; k < j; ++k)
    {
      const Scalar factor = computeTensorFactor(basisIndex, k, y[k]);
      partial[basisIndex] *= factor * factor;
    }
  Scalar normalizer = 0.0;
  for (UnsignedInteger basisIndex = 0; basisIndex < size_; ++basisIndex) normalizer += partial[basisIndex];
  if (!(normalizer > 0.0)) return 0.0;
  const Function slice(ChristoffelTensorSliceEvaluation(this, partial, j));
  Scalar error = -1.0;
  Point ai;
  Point bi;
  Sample fi;
  Point ei;
  const Scalar integral = GaussKronrod().integrate(slice, lower, x, error, ai, bi, fi, ei)[0];
  return SpecFunc::Clip01(integral / normalizer);
}

/* Sequential conditional PDF fanning out to the scalar version */
Point ChristoffelDistribution::computeSequentialConditionalPDF(const Point & x) const
{
  if (x.getDimension() != getDimension()) throw InvalidArgumentException(HERE) << "Error: expected a point of dimension=" << getDimension() << ", got dimension=" << x.getDimension();
  if (!isTensorProduct_ || !measure_.hasIndependentCopula())
    return DistributionImplementation::computeSequentialConditionalPDF(x);
  Point result(getDimension());
  Point prefix;
  for (UnsignedInteger j = 0; j < getDimension(); ++j)
  {
    result[j] = computeConditionalPDF(x[j], prefix);
    prefix.add(x[j]);
  }
  return result;
}

/* Sequential conditional CDF fanning out to the scalar version */
Point ChristoffelDistribution::computeSequentialConditionalCDF(const Point & x) const
{
  if (x.getDimension() != getDimension()) throw InvalidArgumentException(HERE) << "Error: expected a point of dimension=" << getDimension() << ", got dimension=" << x.getDimension();
  if (!isTensorProduct_ || !measure_.hasIndependentCopula())
    return DistributionImplementation::computeSequentialConditionalCDF(x);
  Point result(getDimension());
  Point prefix;
  for (UnsignedInteger j = 0; j < getDimension(); ++j)
  {
    result[j] = computeConditionalCDF(x[j], prefix);
    prefix.add(x[j]);
  }
  return result;
}

/* Univariate factor of a basis function along dimension k at x */
Scalar ChristoffelDistribution::computeTensorFactor(const UnsignedInteger basisIndex,
    const UnsignedInteger k,
    const Scalar x) const
{
  const Indices multi(enumerateFunction_(basisIndex));
  const UnsignedInteger degree = multi[k];
  if (usePolynomialTensor_) return tensorPolynomials_[k][degree](x);
  return tensorFunctions_[k][degree](x);
}

/* Partial Christoffel sum over t in dimension j given partial products */
Scalar ChristoffelDistribution::computePartialChristoffel(const Point & partialProducts,
    const UnsignedInteger j,
    const Scalar t) const
{
  Scalar total = 0.0;
  for (UnsignedInteger basisIndex = 0; basisIndex < size_; ++basisIndex)
  {
    const Scalar factor = computeTensorFactor(basisIndex, j, t);
    total += partialProducts[basisIndex] * factor * factor;
  }
  return total;
}

/* One draw by sequential 1D conditional sampling (tensor independent case) */
Point ChristoffelDistribution::drawSequentialTensor() const
{
  const UnsignedInteger dimension = getDimension();
  const UnsignedInteger gridSize = ResourceMap::GetAsUnsignedInteger("ChristoffelDistribution-SliceGridSize");
  const Scalar safety = ResourceMap::GetAsScalar("ChristoffelDistribution-KnSafetyFactor");
  Point point(dimension);
  Point partial(size_, 1.0);
  for (UnsignedInteger j = 0; j < dimension; ++j)
  {
    const Distribution marginalJ(measure_.getMarginal(j));
    const Interval rangeJ(marginalJ.getRange());
    const Scalar lower = rangeJ.getLowerBound()[0];
    const Scalar upper = rangeJ.getUpperBound()[0];
    Scalar envelope = 0.0;
    for (UnsignedInteger g = 0; g <= gridSize; ++g)
    {
      const Scalar t = lower + (upper - lower) * g / gridSize;
      const Scalar density = marginalJ.computePDF(t) * computePartialChristoffel(partial, j, t);
      if (density > envelope) envelope = density;
    }
    envelope *= safety;
    if (!(envelope > 0.0)) throw InternalException(HERE) << "Error: null slice envelope along dimension " << j;
    // Adaptive rejection with restart on envelope exceedance
    while (true)
    {
      const Scalar proposal = (marginalJ.getRealization())[0];
      const Scalar density = marginalJ.computePDF(proposal) * computePartialChristoffel(partial, j, proposal);
      if (density > envelope)
      {
        envelope = density * safety;
        continue;
      }
      if (RandomGenerator::Generate() * envelope <= density)
      {
        point[j] = proposal;
        break;
      }
    }
    for (UnsignedInteger basisIndex = 0; basisIndex < size_; ++basisIndex)
    {
      const Scalar factor = computeTensorFactor(basisIndex, j, point[j]);
      partial[basisIndex] *= factor * factor;
    }
  }
  return point;
}

/* One draw by rejection from the reference with the Kn envelope */
Point ChristoffelDistribution::drawByRejection() const
{
  const Scalar safety = ResourceMap::GetAsScalar("ChristoffelDistribution-KnSafetyFactor");
  Scalar envelope = computeKn();
  if (!(envelope > 0.0)) throw InternalException(HERE) << "Error: null Christoffel envelope.";
  while (true)
  {
    const Point proposal(measure_.getRealization());
    const Scalar kn = computeChristoffel(proposal);
    if (kn > envelope)
    {
      envelope = kn * safety;
      knEstimate_ = std::max(knEstimate_, envelope);
      continue;
    }
    if (RandomGenerator::Generate() * envelope <= kn) return proposal;
  }
}

/* Realization: tensor sequential, ratio-of-uniforms, then crude rejection */
Point ChristoffelDistribution::getRealization() const
{
  if (size_ == 0) return DistributionImplementation::getRealization();
  if (isTensorProduct_ && measure_.hasIndependentCopula()) return drawSequentialTensor();
  if (getDimension() <= ResourceMap::GetAsUnsignedInteger("ChristoffelDistribution-RatioUniformMaxDimension") && sampler_.isInitialized())
    return sampler_.getRealization();
  return drawByRejection();
}

/* Sample: batched ratio-of-uniforms when selected, sequential loop otherwise */
Sample ChristoffelDistribution::getSample(const UnsignedInteger size) const
{
  if (size_ != 0 && !isTensorProduct_ && getDimension() <= ResourceMap::GetAsUnsignedInteger("ChristoffelDistribution-RatioUniformMaxDimension") && sampler_.isInitialized())
    return sampler_.getSample(size);
  return DistributionImplementation::getSample(size);
}

/* Continuity flags */
Bool ChristoffelDistribution::isContinuous() const
{
  return true;
}

Bool ChristoffelDistribution::isDiscrete() const
{
  return false;
}

/* Parameters accessors: no parametric representation */
ChristoffelDistribution::PointWithDescriptionCollection ChristoffelDistribution::getParametersCollection() const
{
  throw NotYetImplementedException(HERE) << "In ChristoffelDistribution::getParametersCollection() const";
}

void ChristoffelDistribution::setParametersCollection(const PointCollection &)
{
  throw NotYetImplementedException(HERE) << "In ChristoffelDistribution::setParametersCollection(const PointCollection & parametersCollection)";
}

/* Numerical range: the reference range */
void ChristoffelDistribution::computeRange()
{
  setRange(measure_.getRange());
}

/* Rebuild the derived members */
void ChristoffelDistribution::update()
{
  measure_ = orthogonalBasis_.getMeasure();
  if (!measure_.isContinuous()) throw InvalidArgumentException(HERE) << "Error: the reference measure must be continuous.";
  if (measure_.isDiscrete()) throw InvalidArgumentException(HERE) << "Error: the reference measure must not be discrete.";
  const UnsignedInteger dimension = measure_.getDimension();
  Basis::FunctionCollection coll(size_);
  for (UnsignedInteger i = 0; i < size_; ++i)
  {
    const Function phi(orthogonalBasis_.build(i));
    if (phi.getInputDimension() != dimension) throw InvalidArgumentException(HERE) << "Error: basis function " << i << " has input dimension=" << phi.getInputDimension() << ", expected dimension=" << dimension;
    if (phi.getOutputDimension() != 1) throw InvalidArgumentException(HERE) << "Error: basis function " << i << " must be scalar valued.";
    coll[i] = phi;
  }
  const Basis basis(coll);
  basis_ = basis;
  setDimension(dimension);
  setDescription(measure_.getDescription());
  computeRange();
  isAlreadyComputedMean_ = false;
  isAlreadyComputedCovariance_ = false;
  isAlreadyComputedKn_ = false;
  detectTensorProduct();
  // Internal RatioOfUniforms sampler, see https://en.wikipedia.org/wiki/Ratio_of_uniforms
  // Initialized only when a sampling path can use it: skipped for tensor
  // independent references (sequential path) and above the dimension cap
  // (crude rejection path), where its multi-start optimization is pure overhead
  sampler_ = RatioOfUniforms();
  if (!(isTensorProduct_ && measure_.hasIndependentCopula()) && getDimension() <= ResourceMap::GetAsUnsignedInteger("ChristoffelDistribution-RatioUniformMaxDimension"))
  {
    sampler_.setOptimizationAlgorithm(OptimizationAlgorithm::GetByName(ResourceMap::GetAsString("ChristoffelDistribution-OptimizationAlgorithm")));
    sampler_.setCandidateNumber(ResourceMap::GetAsUnsignedInteger("ChristoffelDistribution-RatioUniformCandidateNumber"));
    sampler_.setLogUnscaledPDFAndRange(getLogPDF(), getRange(), true);
  }
}

/* Detect a tensor-product factory and cache univariate evaluators */
void ChristoffelDistribution::detectTensorProduct()
{
  isTensorProduct_ = false;
  usePolynomialTensor_ = false;
  tensorPolynomials_ = Collection<Collection<OrthogonalUniVariatePolynomial>>();
  tensorFunctions_ = Collection<Collection<UniVariateFunction>>();
  const UnsignedInteger dimension = getDimension();
  const OrthogonalProductPolynomialFactory * p_polynomialFactory = dynamic_cast<const OrthogonalProductPolynomialFactory *>(orthogonalBasis_.getImplementation().get());
  if (p_polynomialFactory != nullptr)
  {
    const OrthogonalProductPolynomialFactory::PolynomialFamilyCollection families(p_polynomialFactory->getPolynomialFamilyCollection());
    if (families.getSize() == dimension)
    {
      enumerateFunction_ = p_polynomialFactory->getEnumerateFunction();
      Collection<Collection<OrthogonalUniVariatePolynomial>> cache(dimension);
      for (UnsignedInteger k = 0; k < dimension; ++k)
      {
        UnsignedInteger maxDegree = 0;
        for (UnsignedInteger basisIndex = 0; basisIndex < size_; ++basisIndex) maxDegree = std::max(maxDegree, enumerateFunction_(basisIndex)[k]);
        Collection<OrthogonalUniVariatePolynomial> perDegree(maxDegree + 1);
        for (UnsignedInteger degree = 0; degree <= maxDegree; ++degree) perDegree[degree] = families[k].build(degree);
        cache[k] = perDegree;
      }
      tensorPolynomials_ = cache;
      usePolynomialTensor_ = true;
      isTensorProduct_ = true;
    }
    return;
  }
  const OrthogonalProductFunctionFactory * p_functionFactory = dynamic_cast<const OrthogonalProductFunctionFactory *>(orthogonalBasis_.getImplementation().get());
  if (p_functionFactory != nullptr)
  {
    const OrthogonalProductFunctionFactory::FunctionFamilyCollection families(p_functionFactory->getFunctionFamilyCollection());
    if (families.getSize() == dimension)
    {
      enumerateFunction_ = p_functionFactory->getEnumerateFunction();
      Collection<Collection<UniVariateFunction>> cache(dimension);
      for (UnsignedInteger k = 0; k < dimension; ++k)
      {
        const UniVariateFunctionFamily family(*families[k].getImplementation());
        UnsignedInteger maxDegree = 0;
        for (UnsignedInteger basisIndex = 0; basisIndex < size_; ++basisIndex) maxDegree = std::max(maxDegree, enumerateFunction_(basisIndex)[k]);
        Collection<UniVariateFunction> perDegree(maxDegree + 1);
        for (UnsignedInteger degree = 0; degree <= maxDegree; ++degree) perDegree[degree] = family.build(degree);
        cache[k] = perDegree;
      }
      tensorFunctions_ = cache;
      usePolynomialTensor_ = false;
      isTensorProduct_ = true;
    }
  }
}

/* Method save() stores the object through the StorageManager */
void ChristoffelDistribution::save(Advocate & adv) const
{
  DistributionImplementation::save(adv);
  adv.saveAttribute("orthogonalBasis_", orthogonalBasis_);
  adv.saveAttribute("size_", size_);
}

/* Method load() reloads the object from the StorageManager */
void ChristoffelDistribution::load(Advocate & adv)
{
  DistributionImplementation::load(adv);
  adv.loadAttribute("orthogonalBasis_", orthogonalBasis_);
  adv.loadAttribute("size_", size_);
  if (size_ == 0) throw InvalidArgumentException(HERE) << "Error: cannot reload a distribution with null space dimension.";
  update();
}

END_NAMESPACE_OPENTURNS
