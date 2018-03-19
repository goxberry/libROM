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

// Description: Declaration of an abstract class defining the
//              interface to the generic hyperreduction algorithm.

#ifndef included_HyperreductionAlgorithm_h
#define included_HyperreductionAlgorithm_h

#include "Matrix.h"
#include "SVD.h"
#include <boost/shared_ptr.hpp>

namespace CAROM {

/**
 * Class HyperreductionAlgorithm defines the interface to the generic
 * hyperreduction algorithm.  The API is intentionally small.  One may
 * collect the samples, compute the hyperreduced reduced order model,
 * and get the dimension of the system.
 *
 * For now, force the size of the permutation vector representing the
 * selection operator to be equal to the basis dimension from the SVD.
 * This restriction may be relaxed later; for now, it makes the
 * implementation considerably simpler, because most of the
 * basis-related query methods can be forwarded to the SVD pointer.
 *
 * Note: this class doesn't own its SVD pointer member, which must be
 * managed separately. Semantically, it would be better to declare
 * this member as a reference, but the SVD class does not have a copy
 * constructor or assignment operator.
 */
class HyperreductionAlgorithm
{
   public:
      /**
       * @brief Constructor.
       *
       * @pre svd != NULL
       * @pre svd->dim > 0
       * @pre svd->sample_per_time_interval > 0
       * @pre num_basis_vectors_used >= 0
       *
       * @param[in] svd The SVD algorithm used to generate the basis
       * to be hyperreduced.
       *
       * @param[in] num_basis_vectors_used Upper limit on the number
       *            of basis vectors used to construct a hyperreduced
       *            operator
       *
       * @param[in] debug_algorithm If true results of the algorithm
       * will be printed to facilitate debugging.
       */
      HyperreductionAlgorithm(
	 boost::shared_ptr<SVD> svd,
	 int  max_num_basis_vectors_used,
         bool debug_algorithm = false);

      /**
       * Destructor.
       */
      ~HyperreductionAlgorithm();

      /**
       * @brief Update the hyperreduced basis based on the current
       * state of this object and its SVD pointer member.
       *
       * @pre d_svd != NULL
       *
       */
      virtual
      void
      UpdateHyperreducedBasis() = 0;

      /**
       * @brief Collect a new sample, u_in, for the basis at supplied
       * time, update the basis, and update the hyperreduced basis.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The new sample.
       * @param[in] time The simulation time of the new sample.
       *
       * @return True if the sampling was successful.
       */
      virtual
      bool
      takeSample(
         const double* u_in,
         double time)
      {
	CAROM_ASSERT(d_svd != NULL);
	bool status = d_svd->takeSample(u_in, time);
	UpdateHyperreducedBasis();
	return status;
      }

      /**
       * @brief Returns the dimension of the SVD basis on
       *        this processor.
       *
       * @return The dimension of the SVD basis on this
       *         processor.
       */
      int
      getBasisDim() const
      {
	CAROM_ASSERT(d_svd != NULL);
	return d_svd->getDim();
      }

      /**
       * @brief Returns the number of basis vectors used on this
       * process in hyperreduced basis.
       *
       * @return The number of basis vectors used on this process in
       * the hyperreduced basis.
       */
      int
      getHyperreducedBasisDim() const
      {
	return d_num_basis_vectors_used;
      }

      /**
       * @brief Returns the basis vectors for the current time interval.
       *
       * @return The SVD basis vectors for the current time interval.
       */
      virtual
      const Matrix*
      getBasis()
      {
	CAROM_ASSERT(d_svd != NULL);
	return d_svd->getBasis();
      }

      /**
       * @brief Returns the singular values for the current time interval.
       *
       * @return The singular values for the current time interval.
       */
      virtual
      const Matrix*
      getBasisSingularValues()
      {
	 CAROM_ASSERT(d_svd != NULL);
	 return d_svd->getSingularValues();
      }

      /**
       * @brief Returns the hyperreduced basis for the current time interval.
       *
       * @return The hyperreduced basis for the current time interval.
       */
      virtual
      const Matrix*
      getHyperreducedBasis()
      {
	CAROM_ASSERT(d_svd != NULL);
	return d_basis_sampled;
      }

      /**
       * @brief Returns the rows of the basis used in the hyperreduced
       * basis; also encodes (with the sampled row owners vector) the
       * selection operator.
       *
       * @return Rows of the basis used in the hyperreduced basis.
       */
      virtual
      const int*
      getHyperreducedBasisRows()
      {
	return d_sampled_row;
      }

      /**
       * @brief Returns the process numbers that own each row of the
       * basis used in the hyperreduced basis; also encodes (with the
       * sampled row vector) the selection operator.
       *
       *
       * @return Process owners of the basis rows used in the
       * hyperreduced basis.
       */
      virtual
      const int*
      getHyperrreducedBasisRowOwners()
      {
	return d_sampled_row_owner;
      }

      /**
       * @brief Returns the number of time intervals on which different sets
       * of basis vectors are defined.
       *
       * @return The number of time intervals on which there are basis vectors.
       */
      int
      getNumBasisTimeIntervals() const
      {
	 CAROM_ASSERT(d_svd != NULL);
         return d_svd->getNumBasisTimeIntervals();
      }

      /**
       * @brief Returns the start time for the requested time interval.
       *
       * @pre 0 <= which_interval
       * @pre which_interval < getNumBasisTimeIntervals()
       *
       * @param[in] which_interval The time interval of interest.
       *
       * @return The start time for the requested time interval.
       */
      double
      getBasisIntervalStartTime(
         int which_interval) const
      {
 	 CAROM_ASSERT(d_svd != NULL);
         return d_svd->getBasisIntervalStartTime(which_interval);
      }

      /**
       * @brief Returns true if the next sample will result in a new time
       * interval.
       *
       * @return True if the next sample results in the creation of a
       * new time interval.
       */
      bool
      isNewTimeInterval() const
      {
	 CAROM_ASSERT(d_svd != NULL);
	 return d_svd->isNewTimeInterval();
      }

   protected:
      /**
       * @brief The SVD algorithm used to compute the basis to be
       * hyperreduced.
       */
      boost::shared_ptr<SVD> d_svd;

      /**
       * @brief Number of basis vectors (from d_svd) used in the
       * hyperreduced basis.
       */
      int d_num_basis_vectors_used;

      // NOTE(oxberry1@llnl.gov): d_sampled_row and
      // d_sampled_row_owner could be represented as std::vector or
      // std::array objects; these objects are declared as integer
      // arrays to preserve the existing hyperreduction API
      /**
       * @brief Sampled rows of the basis used in the hyperreduced
       * basis, represented as an array of integers.
       */
      int* d_sampled_row;

      /**
       * @brief Process owner of sampled rows of basis used in the
       * hyperreduced basis, represented as an array of integers.
       *
       */
      int* d_sampled_row_owner;

      /**
       * @brief Possibly large matrix storing the hyperreduced basis.
       */
      Matrix* d_basis_sampled;

      /**
       * @brief Flag to indicate if results of algorithm should be
       * printed for debugging purposes.
       */
      bool d_debug_algorithm;

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      HyperreductionAlgorithm();

      /**
       * @brief Unimplemented copy constructor.
       */
      HyperreductionAlgorithm(
         const HyperreductionAlgorithm& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      HyperreductionAlgorithm&
      operator = (
         const HyperreductionAlgorithm& rhs);
};

}

#endif
