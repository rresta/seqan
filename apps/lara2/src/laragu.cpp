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
// This file contains the seqan_laragu application.
// ==========================================================================

#define SEQAN_LARAGU

// ============================================================================
// Prerequisites
// ============================================================================

// ----------------------------------------------------------------------------
// std headers
// ----------------------------------------------------------------------------

#include <iostream>
#include <map>
#include <string>
#include <limits>
#include <sstream>
#include <omp.h>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/arg_parse.h>
#include <seqan/store.h>
#include <seqan/stream.h>
#include <seqan/align.h>
#include <seqan/align_rna.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/rna_io.h>


// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

// defines all the constants used in the app
#include "data_types.h"
#include "option.h"
#include "store_seqs.h"
#include "lara_core.h"
#include "alignment.h"
#include "lemon_graph.h"

using namespace seqan;

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    if (argc==1)
        std::cout << "type ./lara_gu --help to get the parameters table (-i option is mandatory)" << std::endl;
    ArgumentParser parser;

    Options options;

    setupArgumentParser(parser, options);
    ArgumentParser::ParseResult res;
    res = parse(options, parser, argc, argv); // Fill the options structure with the standard and the acquired arguments

    if (res != ArgumentParser::PARSE_OK)  // Check the arguments
        return res == ArgumentParser::PARSE_ERROR;

    TRnaVect rnaSeqsRef;
    TRnaVect rnaSeqs; // Define the structure that will store all the input RNAs sequences

    _V(options, "\nOpen the input file and fill the TRnaVect data structure");
    readRnaRecords(rnaSeqs, options, options.inFile);
    _V(options, "\nAnalyse the rnaSeqs structure and check which info are required for the selected analysis");

// Add the weight interaction edges vector map in the data structure
    bppInteractionGraphBuild(rnaSeqs, options);
//  Create the alignment data structure that will be used to store all the alignments
    TRnaAlignVect rnaAligns;
//TODO make a function that do this job
    if( options.inFileRef == "" )
    {
        alignVectorBuild(rnaAligns, rnaSeqs, options);
    }
    else
    {
        readRnaRecords(rnaSeqsRef, options, options.inFileRef);
        bppInteractionGraphBuild(rnaSeqsRef, options);
        alignVectorBuild(rnaAligns, rnaSeqs, rnaSeqsRef, options);
    }
    /*
    for(unsigned i=0; i<length(rnaSeqs); ++i)
    {
        std::cout << rnaSeqs[i].bppMatrGraphs[1].inter << std::endl;
    }
    */
    StringSet<TAlign> alignsSimd;
    String<TScoreValue> resultsSimd;
    // simd vector is created
    createSimdAligns(alignsSimd, rnaAligns);

// first non-structural alignment is computed
    firstSimdAlignsGlobalLocal(resultsSimd, alignsSimd, options);

#pragma omp parallel for num_threads(options.threads)
    for(unsigned i = 0; i < length(alignsSimd); ++i)
    {
        resize(rnaAligns[i].lamb, length(rnaAligns[i].rna1.sequence));
        resize(rnaAligns[i].mask, length(rnaAligns[i].rna2.sequence));
        resize(rnaAligns[i].upperBoundVect, length(rnaAligns[i].rna2.sequence));

// Save the best alignments that give the absolute maximum score
        saveBestAlign(rnaAligns[i], alignsSimd[i], resultsSimd[i]);
// Create the mask of the current alignment to be used for the upper, lower bound computation and the lamb update
        maskCreator(rnaAligns[i], alignsSimd[i]);

        if(options.lowerBoundMethod == LBLEMONMWM) // The MWM is computed to fill the LowerBound
        {
//  Define the datastructure that will be passed to the lemon::MWM function to compute the full lowerBound
            TMapVect lowerBound4Lemon;
            lowerBound4Lemon.resize(rnaAligns[i].maskIndex);
            computeBounds(rnaAligns[i], lowerBound4Lemon);
// Compute the MWM with the Lemon library
            myLemon::computeLowerBound(lowerBound4Lemon, rnaAligns[i]);
            rnaAligns[i].lowerBound = rnaAligns[i].lowerLemonBound.mwmPrimal;
            rnaAligns[i].slm = rnaAligns[i].slm - (rnaAligns[i].lowerLemonBound.mwmCardinality * 2);
        } else if (options.lowerBoundMethod == LBAPPROXMWM) // Approximation of MWM is computed to fill the LowerBound
        {
            computeBounds(rnaAligns[i]);
        } else if (options.lowerBoundMethod == LBMWMTEST) // Function used to test the aproximation of MWM is computed to fill the LowerBound
        {
//  In this branch three different methods are available for the computation: 1) the MWM approx, 2) the lemon MWM, 3) the seqan MWM <to be implemented>
//  The approximation is used while the other structures are computed
//  Define the datastructure that will be passed to the lemon::MWM function to compute the full lowerBound
            TMapVect lowerBound4Lemon;
            lowerBound4Lemon.resize(rnaAligns[i].maskIndex);
            computeBoundsTest(rnaAligns[i], lowerBound4Lemon);
// Compute the MWM with the Lemon library
            myLemon::computeLowerBound(lowerBound4Lemon, rnaAligns[i]);
            std::cout << "Lower bound = " << rnaAligns[i].lowerBound << std::endl;
            std::cout << "Upper bound = " << rnaAligns[i].upperBound << std::endl;
            std::cout << "Slm = " << rnaAligns[i].slm << std::endl;
        }


        unsigned index = 0;
// The alignemnt that give the smallest difference between up and low bound should be saved
        saveBestAlignMinBound(rnaAligns[i], alignsSimd[i], resultsSimd[i], index);
        if (rnaAligns[i].upperBound - rnaAligns[i].lowerBound < options.epsilon)
        {
            std::cout << "computation for this alignment should stopped and the bestAlignMinBounds should be returned" << std::endl;
// FIXME computation for this alignment should stopped and the bestAlignMinBounds should be returned
//            return 0;
        }
//  Compute the step size for the Lambda update
        rnaAligns[i].stepSize = options.my * ((rnaAligns[i].upperBound - rnaAligns[i].lowerBound) / rnaAligns[i].slm);

        _VVV(options, "\nThe step size to be used for Lambda for alignment " << i << " in iteration 0 is " << rnaAligns[i].stepSize);

        updateLambda(rnaAligns[i]);
    }

//    String<TScoringSchemeStruct> alignsSimdLamb;
//    seqan::resize(alignsSimdLamb, length(alignsSimd));
// Add struct scoring scheme pointers to each alignment cell of the alignment vector
    for (unsigned i = 0; i < length(alignsSimd); ++i)
    {
        rnaAligns[i].structScore.lamb = & rnaAligns[i].lamb;
        rnaAligns[i].structScore.score_matrix = options.laraScoreMatrix;
    }

    for(unsigned x = 0; x < options.iterations; ++x)
    {
#pragma omp parallel for num_threads(options.threads)
        for (unsigned i = 0; i < length(alignsSimd); ++i) // TODO replace this function with the SIMD implementation for execute in PARALLEL
        {

            if(options.lowerBoundMethod == LBLEMONMWM) // The MWM is computed to fill the LowerBound
            {
//  Define the datastructure that will be passed to the lemon::MWM function to compute the full lowerBound
                TMapVect lowerBound4Lemon;
                lowerBound4Lemon.resize(rnaAligns[i].maskIndex);
                computeBounds(rnaAligns[i], lowerBound4Lemon);
// Compute the MWM with the Lemon library
                myLemon::computeLowerBound(lowerBound4Lemon, rnaAligns[i]);
                rnaAligns[i].lowerBound = rnaAligns[i].lowerLemonBound.mwmPrimal;
                rnaAligns[i].slm = rnaAligns[i].slm - (rnaAligns[i].lowerLemonBound.mwmCardinality * 2);
            } else if (options.lowerBoundMethod == LBAPPROXMWM) // Approximation of MWM is computed to fill the LowerBound
            {
                computeBounds(rnaAligns[i]);
            } else if (options.lowerBoundMethod == LBMWMTEST) // Function used to test the aproximation of MWM is computed to fill the LowerBound
            {
//  In this branch three different methods are available for the computation: 1) the MWM approx, 2) the lemon MWM, 3) the seqan MWM <to be implemented>
//  The approximation is used while the other structures are computed
//  Define the datastructure that will be passed to the lemon::MWM function to compute the full lowerBound
                TMapVect lowerBound4Lemon;
                lowerBound4Lemon.resize(rnaAligns[i].maskIndex);
                computeBoundsTest(rnaAligns[i], lowerBound4Lemon);
// Compute the MWM with the Lemon library
                myLemon::computeLowerBound(lowerBound4Lemon, rnaAligns[i]);
                std::cout << "Lower bound = " << rnaAligns[i].lowerBound << std::endl;
                std::cout << "Upper bound = " << rnaAligns[i].upperBound << std::endl;
                std::cout << "Slm = " << rnaAligns[i].slm << std::endl;
            }

// The alignemnt that give the smallest difference between up and low bound should be saved
            saveBestAlignMinBound(rnaAligns[i], alignsSimd[i], resultsSimd[i], x);
            if (rnaAligns[i].upperBound - rnaAligns[i].lowerBound < options.epsilon)
            {
                std::cout << "computation for this alignment should stopped and the bestAlignMinBounds should be returned" << std::endl;
// FIXME computation for this alignment should stopped and the bestAlignMinBounds should be returned
//            return 0;
            }

            computeBounds(rnaAligns[i]);

//  Compute the step size for the Lambda update
            rnaAligns[i].stepSize = options.my * ((rnaAligns[i].upperBound - rnaAligns[i].lowerBound) / rnaAligns[i].slm);

            _VVV(options, "\nThe step size to be used for Lambda for alignment " << i << " in iteration" << x << " is " << rnaAligns[i].stepSize);

            updateLambda(rnaAligns[i]);
        }
//        std::cout << "computation " << j << std::endl;
        std::cerr << "|" ;
    }
    std::cerr << std::endl;

    for (unsigned i = 0; i < length(rnaAligns); ++i)
    {
        _VV(options, "Alignment number " << i);
        _VV(options, "Best alignment based on the Score is " << rnaAligns[i].bestAlignScore << "\n" <<  rnaAligns[i].bestAlign);
        _VV(options, "Best alignment based on the Min Bounds is " << rnaAligns[i].bestAlignScoreMinBounds << "\n" <<  rnaAligns[i].bestAlignMinBounds);
        _VV(options, "The step size to be used for Lambda at last iteration is " << rnaAligns[i].stepSize);
        _VV(options, "Best Lower bound is " << rnaAligns[i].lowerMinBound);
        _VV(options, "Best Upper bound is " << rnaAligns[i].upperMinBound);
        _VV(options, "Step size found at iteration " << rnaAligns[i].itMinBounds);
        _VV(options, "Minumum step size is " << rnaAligns[i].stepSizeMinBound << "\n\n");

    }

    return 0;
}

