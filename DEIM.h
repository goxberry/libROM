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

// Description: Interface to the DEIM algorithm to determine the rows of the
// rhs to be sampled for the interpolation of the rhs.

#ifndef included_DEIM_h
#define included_DEIM_h

#include "HyperreductionAlgorithm.h"

namespace CAROM {

  class SVD;

/**
 * Class DEIM defines the interface to the discrete empirical
 * interpolation algorithm.
 */
  class DEIM : public HyperreductionAlgorithm {

  public:

    /**
     * @brief Constructor.
     *
     * @pre svd != NULL
     * @pre svd->dim > 0
     * @pre svd->sample_per_time_interval > 0
     * @pre num_basis_vectors_used > 0
     *
     * @param[in] svd The SVD algorithm used to generate the basis
     * to be hyperreduced.
     *
     * @param[in] num_basis_vectors_used The number of basis vectors
     *            used to construct a hyperreduced operator
     *
     * @param[in] debug_algorithm If true results of the algorithm
     * will be printed to facilitate debugging.
     */
    DEIM(
	 boost::shared_ptr<SVD> svd,
	 int num_basis_vectors_used,
	 bool debug_algorithm = false);

    /**
     * Destructor
     */
    ~DEIM();

    /**
     * @brief Update the hyperreduced basis based on the current
     * state of this object and its SVD pointer member.
     *
     * @pre d_svd != NULL
     *
     */
    virtual
    void
    UpdateHyperreducedBasis();

  private:
    /**
     * @brief Unimplemented default constructor.
     */
    DEIM();

    /**
     * @brief Unimplemented copy constructor.
     */
    DEIM(
       const DEIM& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    DEIM&
      operator = (
         const DEIM& rhs);
  };

}

#endif
