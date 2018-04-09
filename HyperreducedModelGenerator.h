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

// Description: The abstract wrapper class for an abstract SVD algorithm and
//              sampler.  This class provides interfaces to each so that an
//              application only needs to instantiate one concrete
//              implementation of this class to control all aspects of basis
//              vector generation.

#ifndef included_HyperreducedModelGenerator_h
#define included_HyperreducedModelGenerator_h

#include "BasisWriter.h"
#include "HyperreducedModelSampler.h"
#include "SVDSampler.h"

#include <string.h>

namespace CAROM {

class BasisWriter;
class Matrix;

/**
 * Class HyperreducedModelGenerator is an abstract base class defining
 * the interface for the generation of hyperreduced models.  This
 * class wraps the abstract hyperreduction algorithm and sampler and
 * provides interfaces to each so that an application only needs to
 * instantiate one concrete implementation of this class to control
 * all aspects of basis vector generation.
 */
class HyperreducedModelGenerator
{
   public:
      /**
       * @brief Destructor.
       */
      virtual
      ~HyperreducedModelGenerator();

      /**
       * @brief Returns true if it is time for the next sample to
       * train the linear basis.
       *
       * @pre time >= 0.0
       *
       * @param[in] time Time of interest.
       *
       * @return True if it is time for the next sample to be taken.
       */
      bool
      isNextLinearSample(
         double time)
      {
	 CAROM_ASSERT(d_hyperreducedmodelsampler != NULL);
         CAROM_ASSERT(time >= 0.0);
         return d_hyperreducedmodelsampler->isNextLinearSample(time);
      }

      /**
       * @brief Returns true if it is time for the next sample to
       * train the nonlinear basis.
       *
       * @pre time >= 0.0
       *
       * @param[in] time Time of interest.
       *
       * @return True if it is time for the next sample to be taken.
       */
      bool
      isNextNonlinearSample(
         double time)
      {
	 CAROM_ASSERT(d_hyperreducedmodelsampler != NULL);
         CAROM_ASSERT(time >= 0.0);
         return d_hyperreducedmodelsampler->isNextNonlinearSample(time);
      }

      /**
       * @brief Sample the new state, u_in, for the linear basis at
       * the given time.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       * @param[in] time The simulation time for the state.
       * @param[in] dt The current simulation linear basis dt.
       *
       * @return True if the sampling was successful.
       */
      bool
      takeLinearSample(
         const double* u_in,
         double time,
         double dt)
      {
	 CAROM_ASSERT(d_hyperreducedmodelsampler != NULL);
         CAROM_ASSERT(u_in != 0);
         CAROM_ASSERT(time >= 0);

         // Check that u_in is not non-zero.
         Vector u_vec(u_in, d_hyperreducedmodelsampler->getLinearBasisDim(),
		      true);
         if (u_vec.norm() == 0.0) {
            return false;
         }

         if (getNumLinearBasisTimeIntervals() > 0 &&
             d_hyperreducedmodelsampler->isNewLinearBasisTimeInterval()) {
            d_hyperreducedmodelsampler->resetLinearBasisDt(dt);
            if (d_basis_writer) {
               d_basis_writer->writeBasis();
            }
         }
         return d_hyperreducedmodelsampler->takeLinearSample(u_in, time);
      }

      /**
       * @brief Sample the new nonlinear function evaluation, f_in,
       * for the nonlinear basis at the given time.
       *
       * @pre f_in != 0
       * @pre time >= 0.0
       *
       * @param[in] f_in The nonlinear function evaluation at the
       * specified time.
       *
       * @param[in] time The simulation time for the nonlinear
       * function evaluation.
       *
       * @param[in] dt The current simulation nonlinear basis dt.
       *
       * @return True if the sampling was successful.
       */
      bool
      takeNonlinearSample(
         const double* f_in,
         double time,
         double dt)
      {
	 CAROM_ASSERT(d_hyperreducedmodelsampler != NULL);
         CAROM_ASSERT(f_in != 0);
         CAROM_ASSERT(time >= 0);

         // Check that f_in is not non-zero.
         Vector f_vec(f_in, d_hyperreducedmodelsampler->getNonlinearBasisDim(),
		      true);
         if (f_vec.norm() == 0.0) {
            return false;
         }

         if (getNumNonlinearBasisTimeIntervals() > 0 &&
             d_hyperreducedmodelsampler->isNewNonlinearBasisTimeInterval()) {
            d_hyperreducedmodelsampler->resetNonlinearBasisDt(dt);
            if (d_basis_writer) {
               d_basis_writer->writeBasis();
            }
         }
         return d_hyperreducedmodelsampler->takeNonlinearSample(f_in, time);
      }

      /**
       * @brief Signal that the final sample has been taken.
       */
      void
      endSamples()
      {
         if (d_basis_writer) {
            d_basis_writer->writeBasis();
         }
      }

      /**
       * @brief Computes next time a linear basis sample is needed.
       *
       * @pre u_in != NULL
       * @pre linear_rhs_in != NULL
       * @pre nonlinear_rhs_in != NULL
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       *
       * @param[in] linear_rhs_in The linear part of the right hand
       * side at the specified time.
       *
       * @param[in] nonlinear_rhs_in The nonlinear part of the right
       * hand side at the specified time.
       *
       * @param[in] time The simulation time for the state.
       */
      double
      computeNextLinearSampleTime(
         double* u_in,
         double* linear_rhs_in,
	 double* nonlinear_rhs_in,
         double time)
      {
         CAROM_ASSERT(u_in != NULL);
         CAROM_ASSERT(linear_rhs_in != NULL);
	 CAROM_ASSERT(nonlinear_rhs_in != NULL);
         CAROM_ASSERT(time >= 0);

         return d_hyperreducedmodelsampler->computeNextLinearSampleTime
	   (u_in, linear_rhs_in, nonlinear_rhs_in, time);
      }

      /**
       * @brief Computes next time a nonlinear basis sample is needed.
       *
       * @pre u_in != NULL
       * @pre linear_rhs_in != NULL
       * @pre nonlinear_rhs_in != NULL
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       *
       * @param[in] linear_rhs_in The linear part of the right hand
       * side at the specified time.
       *
       * @param[in] nonlinear_rhs_in The nonlinear part of the right
       * hand side at the specified time.
       *
       * @param[in] time The simulation time for the state.
       */
      double
      computeNextNonlinearSampleTime(
         double* u_in,
         double* linear_rhs_in,
	 double* nonlinear_rhs_in,
         double time)
      {
         CAROM_ASSERT(u_in != NULL);
         CAROM_ASSERT(linear_rhs_in != NULL);
	 CAROM_ASSERT(nonlinear_rhs_in != NULL);
         CAROM_ASSERT(time >= 0);

         return d_hyperreducedmodelsampler->computeNextNonlinearSampleTime
	   (u_in, linear_rhs_in, nonlinear_rhs_in, time);
      }

      /**
       * @brief Returns the linear basis vectors for the current
       * linear basis time interval as a Matrix.
       *
       * @return The linear basis vectors for the current linear basis
       * time interval.
       */
      const Matrix*
      getLinearBasis()
      {
         return d_hyperreducedmodelsampler->getLinearBasis();
      }

      /**
       * @brief Returns the nonlinear basis vectors for the current
       * nonlinear basis time interval as a Matrix.
       *
       * @return The nonlinear basis vectors for the current nonlinear basis
       * time interval.
       */
      const Matrix*
      getNonlinearBasis()
      {
         return d_hyperreducedmodelsampler->getNonlinearBasis();
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
	 CAROM_ASSERT(d_hyperreducedmodelsampler != NULL);
	 return d_hyperreducedmodelsampler->getHyperreducedNonlinearBasis();
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
	 CAROM_ASSERT(d_hyperreducedmodelsampler != NULL);
	 return d_hyperreducedmodelsampler->getHyperreducedNonlinearBasisRows();
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
	 CAROM_ASSERT(d_hyperreducedmodelsampler != NULL);
	 return
	   d_hyperreducedmodelsampler->getHyperreducedNonlinearBasisRowOwners();
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
	 CAROM_ASSERT(d_hyperreducedmodelsampler != NULL);
         return d_hyperreducedmodelsampler->getLinearBasisSingularValues();
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
	 CAROM_ASSERT(d_hyperreducedmodelsampler != NULL);
	 return d_hyperreducedmodelsampler->getNonlinearBasisSingularValues();
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
	 CAROM_ASSERT(d_hyperreducedmodelsampler != NULL);
         return d_hyperreducedmodelsampler->getNumLinearBasisTimeIntervals();
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
	 CAROM_ASSERT(d_hyperreducedmodelsampler != NULL);
         return d_hyperreducedmodelsampler->getNumNonlinearBasisTimeIntervals();
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
	 CAROM_ASSERT(d_hyperreducedmodelsampler != NULL);
         CAROM_ASSERT(0 <= which_interval);
         CAROM_ASSERT(which_interval < getNumLinearBasisTimeIntervals());
         return d_hyperreducedmodelsampler->getLinearBasisIntervalStartTime
	   (which_interval);
      }

      /**
       * @brief Returns the start time for the requested nonlinear basis
       * time interval.
       *
       * @pre 0 <= which_interval
       * @pre which_interval < getNumNonlinearBasisTimeIntervals()
       *
       * @param[in] which_interval Nonlinear basis time interval whose
       * start time is needed.
       *
       * @return The start time for the requested nonlinear basis time
       * interval.
       */
      double
      getNonlinearBasisIntervalStartTime(
         int which_interval) const
      {
	 CAROM_ASSERT(d_hyperreducedmodelsampler != NULL);
         CAROM_ASSERT(0 <= which_interval);
         CAROM_ASSERT(which_interval < getNumNonlinearBasisTimeIntervals());
         return d_hyperreducedmodelsampler->getNonlinearBasisIntervalStartTime
	   (which_interval);
      }

   protected:
      /**
       * @brief Constructor.
       *
       * Although all member functions are implemented by delegation to either
       * d_basis_writer or d_hyperreducedmodelsampler, this class is still abstract.  In this
       * context it is not yet known which SVDSampler to instantiate.  Hence an
       * instance of this class may not be constructed.
       *
       * @param[in] basis_file_name The base part of the name of the file
       *                            containing the basis vectors.  Each process
       *                            will append its process ID to this base
       *                            name.
       * @param[in] file_format The format of the file containing the basis
       *                        vectors.
       */
      HyperreducedModelGenerator(
         const std::string& basis_file_name = "",
         Database::formats file_format = Database::HDF5);

      /**
       * @brief Writer of basis vectors.
       */
      BasisWriter* d_basis_writer;

      /**
       * @brief Pointer to the underlying sampling control object.
       */
      boost::shared_ptr<HyperreducedModelSampler> d_hyperreducedmodelsampler;

   private:
      /**
       * @brief Unimplemented copy constructor.
       */
      HyperreducedModelGenerator(
         const HyperreducedModelGenerator& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      HyperreducedModelGenerator&
      operator = (
         const HyperreducedModelGenerator& rhs);
};

}

#endif
