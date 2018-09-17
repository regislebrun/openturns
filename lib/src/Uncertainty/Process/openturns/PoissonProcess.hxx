//                                               -*- C++ -*-
/**
 *  @brief A class which implements the Poisson point process
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
#ifndef OPENTURNS_POISSONPROCESS_HXX
#define OPENTURNS_POISSONPROCESS_HXX

#include "openturns/ProcessImplementation.hxx"
#include "openturns/Process.hxx"
#include "openturns/Pointer.hxx"
#include "openturns/Distribution.hxx"
#include "openturns/POISSONPROCESSCoefficients.hxx"
#include "openturns/POISSONPROCESSState.hxx"
#include "openturns/WhiteNoise.hxx"

BEGIN_NAMESPACE_OPENTURNS

/**
 * @class POISSONPROCESS
 *
 * An interface class for POISSONPROCESS
 */
class OT_API POISSONPROCESS
  : public ProcessImplementation
{
  CLASSNAME

public:

  /** Some typedefs to ease reading */


  /** Default constructor */
  POISSONPROCESS();

  /** Standard constructor with coefficients and a White Noise */
  POISSONPROCESS(const POISSONPROCESSCoefficients & ARCoefficients,
       const POISSONPROCESSCoefficients & MACoefficients,
       const WhiteNoise & whiteNoise);

  /** Standard constructor with coefficients, a White Noise and a state */
  POISSONPROCESS(const POISSONPROCESSCoefficients & ARCoefficients,
       const POISSONPROCESSCoefficients & MACoefficients,
       const WhiteNoise & whiteNoise,
       const POISSONPROCESSState & state);

  /** Virtual constructor */
  virtual POISSONPROCESS * clone() const;

  /** String converter */
  String __repr__() const;

  /** String converter  - pretty print */
  String __str__(const String & offset = "") const;

  /** Is the underlying a stationary process ? */
  Bool isStationary() const;

  /** Is the underlying a Gaussian process ? */
  Bool isNormal() const;

  /** Realization accessor */
  Field getRealization() const;

  /** Prediction of the N futur iterations of an POISSONPROCESS process */
  using ProcessImplementation::getFuture;
  TimeSeries getFuture(const UnsignedInteger stepNumber) const;

  /** Coefficients accessor : AR & MA */
  POISSONPROCESSCoefficients getARCoefficients() const;
  POISSONPROCESSCoefficients getMACoefficients() const;

  /** State accessor of the POISSONPROCESS process */
  POISSONPROCESSState getState() const;
  void setState(const POISSONPROCESSState & state) const;

  /** WhiteNoise accessor of the POISSONPROCESS process */
  WhiteNoise getWhiteNoise() const;
  void setWhiteNoise(const WhiteNoise & whiteNoise);

  /** Computation of nThermalization */
  UnsignedInteger computeNThermalization(const Scalar epsilon) const;

  /** Nthermalization accessor - Visibility is done */
  UnsignedInteger getNThermalization() const;

  /** Nthermalization accessor - Setting the value */
  void setNThermalization(const UnsignedInteger n);

  /** Get the random vector corresponding to the i-th marginal component */
  Process getMarginal(const UnsignedInteger i) const;

  /** Get the marginal random vector corresponding to indices components */
  Process getMarginal(const Indices & indices) const;

  /** Method save() stores the object through the StorageManager */
  void save(Advocate & adv) const;

  /** Method load() reloads the object from the StorageManager */
  void load(Advocate & adv);


private:


  /** Compute the steps next values of the process starting from the current state.
      The result is the current state extended stepNumber steps further */
  POISSONPROCESSState computeReccurence(const UnsignedInteger stepNumber) const;

  /** thermalize method */
  void thermalize() const;

  /** AR coefficients of the POISSONPROCESS process */
  POISSONPROCESSCoefficients ARCoefficients_;

  /** MA coefficients of the POISSONPROCESS process */
  POISSONPROCESSCoefficients MACoefficients_;

  /** The distribution underlying the White Noise of the process */
  Distribution noiseDistribution_;

  /** Size of AR part */
  UnsignedInteger p_;

  /** Size of MA part */
  UnsignedInteger q_;

  /** Mutable  current state of the POISSONPROCESS process */
  mutable POISSONPROCESSState state_;

  /** Boolean flag - compute once the number of iterations of the thermalize */
  mutable Bool hasComputedNThermalization_;

  /** Number of iterations for the thermalize method */
  mutable UnsignedInteger nThermalization_;

}; /* class POISSONPROCESS */
END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_POISSONPROCESS_HXX */
