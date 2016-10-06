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
// This file contains
// ==========================================================================
#ifndef _INCLUDE_STRUCT_ALIGN_H_
#define _INCLUDE_STRUCT_ALIGN_H_

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

//#include "vienna_rna.h"

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function alignVectorBuild()
// ----------------------------------------------------------------------------
// Used to generate the alignments from a single input file
template <typename TRnaAligns, typename TRnaSeqs, typename TOption>
void alignVectorBuild(TRnaAligns & rnaAligns, TRnaSeqs const & rnaSeqs,
                      TOption const & options)
{
    for (unsigned i = 0; i < length(rnaSeqs) - 1; ++i) {
        TRnaAlign rnaAlign;
        for (unsigned j = i + 1; j < length(rnaSeqs); ++j) {
            if (length(rnaSeqs[i].sequence) <
                length(rnaSeqs[j].sequence)) // in this way the alignment map structure will be always created with the maximum size
            {
                rnaAlign.rna1 = rnaSeqs[j];
                rnaAlign.rna2 = rnaSeqs[i];
            } else {
                rnaAlign.rna1 = rnaSeqs[i];
                rnaAlign.rna2 = rnaSeqs[j];
            }
            if(options.verbose > 2)
            {
                std::cout << rnaAlign.rna1.sequence << std::endl;
                std::cout << rnaAlign.rna2.sequence << std::endl;
            }
            rnaAligns.push_back(rnaAlign);
        }
    }
}

// ----------------------------------------------------------------------------
// Function alignVectorBuild()
// ----------------------------------------------------------------------------
// Used to generate the alignments from two different input files
template <typename TRnaAligns, typename TRnaSeqs, typename TOption>
void alignVectorBuild(TRnaAligns & rnaAligns, TRnaSeqs const & rnaSeqs,
                      TRnaSeqs const & rnaSeqsRef, TOption const & options)
{
    for(unsigned i=0;i<length(rnaSeqs); ++i)
    {
        TRnaAlign rnaAlign;
        for(unsigned j=0;j<length(rnaSeqsRef); ++j)
        {
            if(length(rnaSeqs[i].sequence) < length(rnaSeqsRef[j].sequence)) // in this way the alignment map structure will be always created with the maximum size
            {
                rnaAlign.rna1 = rnaSeqsRef[j];
                rnaAlign.rna2 = rnaSeqs[i];
            }else {
                rnaAlign.rna1 = rnaSeqs[i];
                rnaAlign.rna2 = rnaSeqsRef[j];
            }
            if(options.verbose > 2)
            {
                std::cout << rnaAlign.rna1.sequence << std::endl;
                std::cout << rnaAlign.rna2.sequence << std::endl;
            }
        }
        rnaAligns.push_back(rnaAlign);
    }
}

// ----------------------------------------------------------------------------
// Function setScoreMatrix()
// ----------------------------------------------------------------------------

template <typename TOptions>
void setScoreMatrix(TOptions & options)
{
    options.laraScoreMatrix.data_gap_extend=options.laraGapExtend;
    options.laraScoreMatrix.data_gap_open=options.laraGapOpen;
    if(options.laraScoreMatrixName != "")
    {
        loadScoreMatrix(options.laraScoreMatrix, toCString(getAbsolutePath(toCString(options.laraScoreMatrixName))));
        _V(options, "Provided scoring matrix will be used");
        showScoringMatrix(options.laraScoreMatrix);
    }
    else
    {
        _V(options, "Predefined RIBOSUM matrix will be used");
        setDefaultScoreMatrix(options.laraScoreMatrix, TRibosum());
        showScoringMatrix(options.laraScoreMatrix);
    }
}

// ----------------------------------------------------------------------------
// Function firstSimdAlignsGlobalLocal()
// ----------------------------------------------------------------------------

template <typename TResultsSimd, typename TAlignsSimd, typename TOptions>
void firstSimdAlignsGlobalLocal(TResultsSimd & resultsSimd, TAlignsSimd & alignsSimd, TOptions const & options)
{
    if (!options.globalLocal)  //TODO implement the global-unconstrained alignment using the parameters in the options
    {
        if (options.affineLinearDgs == 0)
            resultsSimd = globalAlignment(alignsSimd, options.laraScoreMatrix, AffineGaps());
        else if (options.affineLinearDgs == 1)
            resultsSimd = globalAlignment(alignsSimd, options.laraScoreMatrix, LinearGaps());
        else
            resultsSimd = globalAlignment(alignsSimd, options.laraScoreMatrix, DynamicGaps());

    } else
    {
        if (options.affineLinearDgs == 0)
            resultsSimd = localAlignment(alignsSimd, options.laraScoreMatrix, AffineGaps());
        else if (options.affineLinearDgs == 1)
            resultsSimd = localAlignment(alignsSimd, options.laraScoreMatrix, LinearGaps());
        else
            resultsSimd = localAlignment(alignsSimd, options.laraScoreMatrix, DynamicGaps());
    }
};

// ----------------------------------------------------------------------------
// Function firstSimdAligns()
// ----------------------------------------------------------------------------

template <typename TResultsSimd, typename TAlignsSimd, typename TRnaAlignVect, typename TOptions>
void firstSimdAligns(TResultsSimd & resultsSimd, TAlignsSimd & alignsSimd, TRnaAlignVect & rnaAligns, TOptions const & options)
{
    for(unsigned i = 0; i < length(rnaAligns); ++i)
    {
        TAlign align;
        resize(rows(align), 2);
        assignSource(row(align, 0), rnaAligns[i].rna1.sequence);
        assignSource(row(align, 1), rnaAligns[i].rna2.sequence);
        appendValue(alignsSimd, align);
    }
    firstSimdAlignsGlobalLocal(resultsSimd, alignsSimd, options);
    if (options.verbose > 2)
    {
        for(unsigned i = 0; i < length(rnaAligns); ++i)
        {
            std::cout << alignsSimd[i];
            std::cout << resultsSimd[i] << std::endl;
        }
    }
}


#endif //_INCLUDE_STRUCT_ALIGN_H_