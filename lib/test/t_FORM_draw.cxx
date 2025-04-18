//                                               -*- C++ -*-
/**
 *  @brief The test file of FORM class
 *
 *  Copyright 2005-2025 Airbus-EDF-IMACS-ONERA-Phimeca
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
#include "openturns/OT.hxx"
#include "openturns/OTtestcode.hxx"

using namespace OT;
using namespace OT::Test;


int main(int, char *[])
{
  TESTPREAMBLE;
  OStream fullprint(std::cout);

  try
  {

    /* We create a numerical math function */
    /** Analytical construction */
    Description inputFunc(2);
    inputFunc[0] = "x0";
    inputFunc[1] = "x1";
    Description formulas(1);
    formulas[0] = "-(6+x0^2-x1)";
    fullprint << "formulas=" << formulas << std::endl;
    SymbolicFunction myFunction(inputFunc, formulas);

    UnsignedInteger dim = myFunction.getInputDimension();
    /* We create a normal distribution point of dimension 1 */
    Point mean(dim, 0.0);
    mean[0] = 5.0; // x0
    mean[1] = 2.1; // x1
    Point sigma(dim, 0.0);
    sigma[0] = 3.3; // x0
    sigma[1] = 3.0; // x1
    IdentityMatrix R(dim);
    Normal myDistribution(mean, sigma, R);

    /* we name the components of the distribution */
    Description componentDescription(dim);
    componentDescription[0] = "Marginal 1";
    componentDescription[1] = "Marginal 2";
    myDistribution.setDescription(componentDescription);

    /* We create a 'usual' RandomVector from the Distribution */
    RandomVector vect(myDistribution);

    /* We create a composite random vector */
    CompositeRandomVector output(myFunction, vect);
    Description outputDescription(1);
    outputDescription[0] = "Interest Variable 1";
    output.setDescription(outputDescription);

    /* We create an Event from this RandomVector */
    ThresholdEvent myEvent(output, Greater(), 0.0);

    /* We create a NearestPoint algorithm */
    Cobyla myCobyla;
    myCobyla.setStartingPoint(mean);
    myCobyla.setMaximumCallsNumber(200);
    myCobyla.setMaximumAbsoluteError(1.0e-10);
    myCobyla.setMaximumRelativeError(1.0e-10);
    myCobyla.setMaximumResidualError(1.0e-10);
    myCobyla.setMaximumConstraintError(1.0e-10);
    fullprint << "myCobyla=" << myCobyla << std::endl << std::flush;

    /* We create a FORM algorithm */
    /* The first parameter is an OptimizationAlgorithm */
    /* The second parameter is an event */
    /* The third parameter is a starting point for the design point research */
    FORM myAlgo(myCobyla, myEvent);

    fullprint << "FORM=" << myAlgo << std::endl;

    /* Perform the simulation */
    myAlgo.run();

    /* Stream out the result */
    FORMResult result(myAlgo.getResult());
    fullprint << "importance factors=" << result.getImportanceFactors() << std::endl;
    fullprint << "importance factors (classical)=" << result.getImportanceFactors(AnalyticalResult::CLASSICAL) << std::endl;

    /* Graph 1 : Importance Factors graph */
    result.drawImportanceFactors();

    /* Graph 1bis : Importance Factors graph */
    result.drawImportanceFactors(AnalyticalResult::CLASSICAL);

    /* Graph 2 : Hasofer Reliability Index Sensitivity Graphs graph */
    result.drawHasoferReliabilityIndexSensitivity();


    /* Graph 3 : FORM Event Probability Sensitivity Graphs graph */
    result.drawEventProbabilitySensitivity();

  }
  catch (TestFailed & ex)
  {
    std::cerr << ex << std::endl;
    return ExitCode::Error;
  }


  return ExitCode::Success;
}
