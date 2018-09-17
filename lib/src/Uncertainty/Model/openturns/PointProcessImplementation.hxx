//                                               -*- C++ -*-
/**
 *  @brief An interface for all implementation class of point processes
 *
 *  Copyright 2005-2018 Airbus-EDF-IMACS-Phimeca
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
#ifndef OPENTURNS_POINTPROCESSIMPLEMENTATION_HXX
#define OPENTURNS_POINTPROCESSIMPLEMENTATION_HXX

#include "openturns/PersistentObject.hxx"
#include "openturns/Point.hxx"
#include "openturns/Sample.hxx"
#include "openturns/Function.hxx"
#include "openturns/Indices.hxx"
#include "openturns/Field.hxx"
#include "openturns/TimeSeries.hxx"
#include "openturns/ProcessSample.hxx"
#include "openturns/CovarianceModel.hxx"
#include "openturns/TrendTransform.hxx"
#include "openturns/Pointer.hxx"

BEGIN_NAMESPACE_OPENTURNS

// Forward declaration
class PointProcess;

/**
 * @class PointProcessImplementation
 *
 * An interface for all implementation class of point processes
 */
class OT_API PointProcessImplementation
  : public PersistentObject
{
  CLASSNAME

public:

  /* Some typedefs for easy reading */
  typedef Collection<Interval> IntervalCollection;

  /** Default constructor */
  PointProcessImplementation();

  /** Virtual constructor */
  virtual PointProcessImplementation * clone() const;

  /** String converter */
  String __repr__() const;
  String __str__(const String & offset = "") const;

  /* Here is the interface that all derived class must implement */

  /** Is the point process a Poisson process ? */
  virtual Bool isPoisson() const;

  /** Covariance model accessor */
  virtual CovarianceModel getCovarianceModel() const;

  /** Trend accessor */
  virtual TrendTransform getTrend() const;

  /** Dimension accessor */
  virtual UnsignedInteger getInputDimension() const;
  virtual UnsignedInteger getOutputDimension() const;

  /** TimeGrid accessor */
  virtual RegularGrid getTimeGrid() const;
  virtual void setTimeGrid(const RegularGrid & timeGrid);

  /** Mesh accessor */
  virtual Mesh getMesh() const;
  virtual void setMesh(const Mesh & mesh);

  /** Discrete realization accessor */
  virtual Field getRealization() const;

  /** Continuous realization accessor */
  virtual Function getContinuousRealization() const;

  /** Process sample accessors */
  virtual ProcessSample getSample(const UnsignedInteger size) const;

  /** Continuation of the last realization on a given number of steps */
  virtual TimeSeries getFuture(const UnsignedInteger stepNumber) const;
  virtual ProcessSample getFuture(const UnsignedInteger stepNumber,
                                  const UnsignedInteger size) const;

  /** Get the marginal process corresponding to the i-th marginal component */
  virtual Process getMarginal(const UnsignedInteger i) const;

  /** Get the marginal process corresponding to indices components */
  virtual Process getMarginal(const Indices & indices) const;

  /** Description accessor */
  virtual void setDescription(const Description & description);
  virtual Description getDescription() const;

  /** Method save() stores the object through the StorageManager */
  virtual void save(Advocate & adv) const;

  /** Method load() reloads the object from the StorageManager */
  virtual void load(Advocate & adv);

protected:

  /** Description of each component */
  Description description_;

  /** Dimension of the process */
  UnsignedInteger dimension_;

  /** The collection of intervals on which the point process is defined */
  IntervalCollection intervalCollection_;

}; /* class ProcessImplementation */


END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_POINTPROCESSIMPLEMENTATION_HXX */
