/******************************************************************************
 *
 * Copyright (c) 2013-2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 * Written by William Arrighi wjarrighi@llnl.gov
 * CODE-686965
 * All rights reserved.
 *
 * This file is part of libROM.
 * For details, see https://computation.llnl.gov/librom
 * Please also read README_BSD_NOTICE.
 *
 * Redistribution and use in source and binary forms, with or without
 * modifications, are permitted provided that the following conditions are met:
 *
 *    o Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the disclaimer below.
 *    o Redistribution in binary form must reproduce the above copyright
 *      notice, this list of conditions and the disclaimer (as noted below) in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *    o Neither the name of the LLNS/LLNL nor the names of its contributors may
 *      be used to endorse or promote products derived from this software
 *      without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
 * LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OR SUCH DAMAGE.
 *
 *****************************************************************************/

// Description: An abstract class defining the interface for
//              determining the next time at which samples should be
//              taken for hyperreduced model generation using an
//              abstract SVD algorithm for linear basis generation and
//              an abstract HyperreductionAlgorithm for hyperreduced
//              nonlinear basis generation.

#ifndef included_HyperreducedModelSampler_h
#define included_HyperreducedModelSampler_h

#include "SVD.h"
#include "HyperreductionAlgorithm.h"
#include <boost/shared_ptr.hpp>

namespace CAROM {

/**
 * Class HyperreducedModelSampler defines the interface to determine the
 * time at which the next sample is needed when using a basis
 * generation algorithm for reducing the linear part of a PDE/ODE, and
 * a hyperreduction algorithm (containing a second basis generation
 * algorithm) for reducing the nonlinear part of a PDE/ODE.
 */
class HyperreducedModelSampler
{
   public:
      /**
       * @brief Constructor.
       */
      HyperreducedModelSampler();

      /**
       * @brief Destructor.
       */
      ~HyperreducedModelSampler();

      /**
       * @brief Returns true if it is time for the next sample to
       * train the linear basis.
       *
       * @param[in] time Time of interest.
       *
       * @return True if a linear basis sample should be taken at the
       * supplied time.
       */
      virtual
      bool
      isNextLinearSample(
         double time) = 0;

      /**
       * @brief Returns true if it is time for the next sample to
       * train the nonlinear basis.
       *
       * @param[in] time Time of interest.
       *
       * @return True if a nonlinear basis sample should be taken at
       * the supplied time.
       */
      virtual
      bool
      isNextNonlinearSample(
	 double time) = 0;

      /**
       * @brief Sample a new state u_in for the linear basis.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       * @param[in] time The simulation time for the state.
       *
       * @return True if the sampling was successful.
       */
      bool
      takeLinearSample(
         const double* u_in,
         double time)
      {
 	 CAROM_ASSERT(d_linear_basis != NULL);
         return d_linear_basis->takeSample(u_in, time);
      }

      /**
       * @brief Sample a new nonlinear function evaluation f_in for the
       * nonlinear basis state.
       *
       * @pre f_in != 0
       * @pre time >= 0.0
       *
       * @param[in] f_in The nonlinear function evaluation at the
       * specified time.
       * @param[in] time The simulation time for the state.
       *
       * @return True if the sampling was successful.
       */
      bool
      takeNonlinearSample(
         const double* f_in,
         double time)
      {
	 CAROM_ASSERT(d_nonlinear_basis != NULL);
	 return d_nonlinear_basis->takeSample(f_in, time);
      }

      /**
       * @brief Computes next time a linear basis sample is needed.
       *
       * @param[in] u_in The state at the specified time.
       *
       * @param[in] linear_rhs_in The linear part of the right hand
       * side at the specified time.
       *
       * @param[in] rhs_in The nonlinear part of right hand side at
       * the specified time.
       *
       * @param[in] time The simulation time for the state.
       *
       * @return The next time a linear basis sample is needed.
       */
      virtual
      double
      computeNextLinearSampleTime(
         double* u_in,
	 double* linear_rhs_in,
         double* nonlinear_rhs_in,
         double time) = 0;

      /**
       * @brief Computes next time a nonlinear basis sample is needed.
       *
       * @param[in] u_in The state at the specified time.
       *
       * @param[in] linear_rhs_in The linear part of the right hand
       * side at the specified time.
       *
       * @param[in] rhs_in The nonlinear part of right hand side at
       * the specified time.
       *
       * @param[in] time The simulation time for the state.
       *
       * @return The next time a nonlinear basis sample is needed.
       */
      virtual
      double
      computeNextNonlinearSampleTime(
         double* u_in,
         double* linear_rhs_in,
	 double* nonlinear_rhs_in,
         double time) = 0;

      /**
       * @brief Returns the linear basis vectors for the current time
       * interval as a Matrix.
       *
       * @return The linear basis vectors for the current time
       * interval.
       */
      const Matrix*
      getLinearBasis()
      {
	 CAROM_ASSERT(d_linear_basis != NULL);
	 return d_linear_basis->getBasis();
      }

      /**
       * @brief Returns the nonlinear basis vectors for the current
       * time interval as a Matrix.
       *
       * @return The nonlinear basis vectors for the current time
       * interval.
       */
      const Matrix*
      getNonlinearBasis()
      {
	 CAROM_ASSERT(d_nonlinear_basis != NULL);
	 return d_nonlinear_basis->getBasis();
      }

      /**
       * @brief Returns the hyperreduced nonlinear basis vectors for
       * the current time interval as a Matrix.
       *
       * @return The nonlinear basis vectors for the current time
       * interval.
       */
      const Matrix*
      getHyperreducedNonlinearBasis()
      {
	 CAROM_ASSERT(d_nonlinear_basis != NULL);
	 return d_nonlinear_basis->getHyperreducedBasis();
      }

      /**
       * @brief Returns the rows of the nonlinear basis used in the
       * hyperreduced nonlinear basis; also encodes (with the sampled
       * row owners vector) the selection operator used in
       * hyperreduction.
       *
       * @return Rows of the nonlinear basis used in the hyperreduced
       * nonlinear basis.
       */
      const int*
      getHyperreducedNonlinearBasisRows()
      {
	 CAROM_ASSERT(d_nonlinear_basis != NULL);
	 return d_nonlinear_basis->getHyperreducedBasisRows();
      }

      /**
       * @brief Returns the process numbers that own each row of the
       * nonlinear basis used in the hyperreduced nonlinear basis;
       * also encodes (with the sampled row vector) the selection
       * operator.
       *
       * @return Process owners of the nonlinear basis rows used in
       * the hyperreduced basis.
       */
      const int*
      getHyperreducedNonlinearBasisRowOwners()
      {
	 CAROM_ASSERT(d_nonlinear_basis != NULL);
	 return d_nonlinear_basis->getHyperreducedBasisRowOwners();
      }

      /**
       * @brief Returns the singular values for linear basis in the
       * current time interval as a Matrix.
       *
       * @return The singular values for the linear basis in the
       * current time interval.
       */
      const Matrix*
      getLinearBasisSingularValues()
      {
	 CAROM_ASSERT(d_linear_basis != NULL);
         return d_linear_basis->getSingularValues();
      }

      /**
       * @brief Returns the singular values for the nonlinear basis in the
       * current time interval as a Matrix.
       *
       * @return The singular values for the nonlinear basis in the
       * current time interval.
       */
      const Matrix*
      getNonlinearBasisSingularValues()
      {
	 CAROM_ASSERT(d_nonlinear_basis != NULL);
	 return d_nonlinear_basis->getBasisSingularValues();
      }

      /**
       * @brief Returns the number of time intervals on which different sets of
       * linear basis vectors are defined.
       *
       * @return The number of time intervals on which there are
       * linear basis vectors.
       */
      int
      getNumLinearBasisTimeIntervals() const
      {
	 CAROM_ASSERT(d_linear_basis != NULL);
         return d_linear_basis->getNumBasisTimeIntervals();
      }

      /**
       * @brief Returns the number of time intervals on which different sets of
       * nonlinear basis vectors are defined.
       *
       * @return The number of time intervals on which there are
       * nonlinear basis vectors.
       */
      int
      getNumNonlinearBasisTimeIntervals() const
      {
	 CAROM_ASSERT(d_nonlinear_basis != NULL);
	 return d_nonlinear_basis->getNumBasisTimeIntervals();
      }

      /**
       * @brief Returns the start time for the requested linear basis
       * time interval.
       *
       * @pre 0 <= which_interval
       * @pre which_interval < getNumLinearBasisTimeIntervals()
       *
       * @param[in] which_interval Linear basis time interval whose
       * start time is needed.
       *
       * @return The start time for the requested linear basis time
       * interval.
       */
      double
      getLinearBasisIntervalStartTime(
         int which_interval) const
      {
	 CAROM_ASSERT(d_linear_basis != NULL);
         CAROM_ASSERT(0 <= which_interval);
         CAROM_ASSERT(which_interval < getNumLinearBasisTimeIntervals());
         return d_linear_basis->getBasisIntervalStartTime(which_interval);
      }

      /**
       * @brief Returns the start time for the requested linear basis
       * time interval.
       *
       * @pre 0 <= which_interval
       * @pre which_interval < getNumNonlinearBasisTimeIntervals()
       *
       * @param[in] which_interval Linear basis time interval whose
       * start time is needed.
       *
       * @return The start time for the requested nonlinear basis time
       * interval.
       */
      double
      getNonlinearBasisIntervalStartTime(
	 int which_interval) const
      {
	 CAROM_ASSERT(d_nonlinear_basis != NULL);
	 CAROM_ASSERT(0 <= which_interval);
	 CAROM_ASSERT(which_interval < getNumNonlinearBasisTimeIntervals());
	 return d_nonlinear_basis->getBasisIntervalStartTime(which_interval);
      }

      /**
       * @brief Returns true if the next linear basis sample will
       * result in a new linear basis time interval.
       *
       * @return True if the next linear basis sample results in the
       *         creation of a new linear basis time interval.
       */
      bool
      isNewLinearBasisTimeInterval() const
      {
	 CAROM_ASSERT(d_linear_basis != NULL);
         return d_linear_basis->isNewTimeInterval();
      }

      /**
       * @brief Returns true if the next nonlinear basis sample will
       * result in a new linear basis time interval.
       *
       * @return True if the next nonlinear basis sample results in
       *         the creation of a new linear basis time interval.
       */
      bool
      isNewNonlinearBasisTimeInterval() const
      {
	 CAROM_ASSERT(d_nonlinear_basis != NULL);
	 return d_nonlinear_basis->isNewTimeInterval();
      }

      /**
       * @brief Resets linear basis sample time step.
       *
       * @param[in] new_dt New value of linear basis sample time step.
       */
      virtual
      void
      resetLinearBasisDt(
         double new_dt) = 0;

      /**
       * @brief Resets nonlinear basis sample time step.
       *
       * @param[in] new_dt New value of nonlinear basis sample time step.
       */
      virtual
      void
      resetNonlinearBasisDt(
	 double new_dt) = 0;

      /**
       * @brief Returns the dimension of the system on this processor.
       *
       * @return The dimension of the system on this processor.
       */
      int
      getLinearBasisDim()
      {
	 CAROM_ASSERT(d_linear_basis != NULL);
         return d_linear_basis->getDim();
      }

      int
      getNonlinearBasisDim()
      {
	 CAROM_ASSERT(d_nonlinear_basis != NULL);
	 return d_nonlinear_basis->getBasisDim();
      }

      int
      getHyperreducedNonlinearBasisDim()
      {
	 CAROM_ASSERT(d_nonlinear_basis != NULL);
	 return d_nonlinear_basis->getHyperreducedBasisDim();
      }

   protected:
      /**
       * @brief Pointer to the abstract SVD algorithm object.
       */
      boost::shared_ptr<SVD> d_linear_basis;

      /**
       * @brief Pointer to the abstract hyperreduction algorithm
       * object.
       */
      boost::shared_ptr<HyperreductionAlgorithm> d_nonlinear_basis;

   private:
      /**
       * @brief Unimplemented copy constructor.
       */
      HyperreducedModelSampler(
         const HyperreducedModelSampler& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      HyperreducedModelSampler&
      operator = (
         const HyperreducedModelSampler& rhs);
};

}

#endif
