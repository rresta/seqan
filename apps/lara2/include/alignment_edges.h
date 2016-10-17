// ==========================================================================
//               LaRAgu - Lagrangian Relaxation Aligner GU
// ==========================================================================
// Copyright (c) 2015-2016, Gianvito Urgese
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Gianvito Urgese nor the names of its contributors
//       may be used to endorse or promote products derived from this software
//       without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL GIANVITO URGESE OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Gianvito Urgese <gianvito.urgese@polito.it>
// ==========================================================================
// This file contains the variable definition and structures of laragu
// application.
// ==========================================================================

#ifndef _INCLUDE_ALIGNMENT_EDGES_H_
#define _INCLUDE_ALIGNMENT_EDGES_H_

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function ()
// ----------------------------------------------------------------------------

template <typename TOption, typename TAlign, typename TScore, typename TRnaAlign>
void checkInterEdgesAndUpdateLambda(TOption const & options, TAlign const & align,
                                    TScore const & alignScore, TRnaAlign & rnaAlign)
{
// TODO if combination of line do not exist create it and fill with 0
    unsigned row0Begin = clippedBeginPosition(row(align, 0));
    unsigned row1Begin = clippedBeginPosition(row(align, 1));
    unsigned gap0 = 0;
    unsigned gap1 = 0;
//TODO this function can be simplifyed using the function toViewPosition(row0, i)
    for(unsigned i=0; i< length(row(align, 0));++i)  //maximum size of this string is length(mapline)
    {
        if(row(align, 0)[i]=='-') // In this choice the assumption that no alignment between - char can exist
        {
            ++gap0;
        }
        else if (row(align, 1)[i]=='-')
        {
            ++gap1;
        }
        else
        {

            rnaAlign.lamb[i+row0Begin-gap0].map[i+row1Begin-gap1] += i; // the default initializer is callet the fist time that set the value to 0
            std::cout << i+row0Begin-gap0 << ":" << i+row1Begin-gap1 << "/" << rnaAlign.lamb[i+row0Begin-gap0].map[i+row1Begin-gap1] << "\t";
        }
    }
    std::cout << align << std::endl;
}

#endif //_INCLUDE_ALIGNMENT_EDGES_H_