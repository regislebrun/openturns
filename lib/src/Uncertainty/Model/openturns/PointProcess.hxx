//                                               -*- C++ -*-
/**
 *  @brief The interface class that implements all point processes
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
#ifndef OPENTURNS_POINTPROCESS_HXX
#define OPENTURNS_POINTPROCESS_HXX

#include "openturns/TypedInterfaceObject.hxx"
#include "openturns/Pointer.hxx"
#include "openturns/PointProcessImplementation.hxx"

BEGIN_NAMESPACE_OPENTURNS

/**
 * @class PointProcess
 *
 * The class that implements all point processes
 */
class OT_API PointProcess
  : public TypedInterfaceObject<PointProcessImplementation>
{
  CLASSNAME

public:

  /* Some typedefs for easy reading */
  typedef Pointer<PointProcessImplementation>            Implementation;
  typedef PointProcessImplementation::IntervalCollection IntervalCollection;
  
  /** Default constructor */
  PointProcess();

  /** Copy constructors */
  PointProcess(const PointProcessImplementation & implementation);


  /** Constructor from implementation */
  PointProcess(const Implementation & p_implementation);

#ifndef SWIG
  /** Constructor from implementation pointer */
  PointProcess(PointProcessImplementation * p_implementation);
#endif

  /** String converter */
  String __repr__() const;

  /** String converter */
  String __str__(const String & offset = "") const;

  /** Description Accessor */
  void setDescription(const Description & description);
  Description getDescription() const;

  /** Is the underlying a Poisson point process ? */
  Bool isPoisson() const;

  /** Support accessor */
  IntervalCollection getSupport() const;
  void setSupport (const IntervalCollection & support);

  /** Dimension accessor */
  UnsignedInteger getDimension() const;

  /** Discrete realization accessor */
  Sample getRealization() const;

  /** PointProcess sample accessors */
  ProcessSample getSample(const UnsignedInteger size) const;

  /** Continuation of the last realization on a given number of steps */
  TimeSeries getFuture(const UnsignedInteger stepNumber) const;
  PointProcessSample getFuture(const UnsignedInteger stepNumber,
                          const UnsignedInteger size) const;

  /** Get the pointProcess corresponding to the i-th marginal component */
  PointProcess getMarginal(const UnsignedInteger i) const;

  /** Get the marginal pointProcess corresponding to indices components */
  PointProcess getMarginal(const Indices & indices) const;


}; /* class PointProcess */
END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_POINTPROCESS_HXX */
