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
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/align.h>
//#include <seqan/graph_align.h>
//#include <seqan/align_profile.h>
#include <seqan/align_rna.h>
//#include <seqan/vcf_io.h>
//#include <seqan/bpseq_io.h>
#include <seqan/graph_types.h>
#include <seqan/rna_io.h>

// ----------------------------------------------------------------------------
// Boost headers
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

#include "top_data_struct.h" // defined all the constant used in the app
//#include "misc_options.h"
#include "option.h" //
#include "store_seqs.h"
#include "interaction_edges.h"
#include "struct_align.h"

using namespace seqan;
// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    if (argc==1)
        std::cout << "type ./lara_gu --help to get the parameters table (-i option is mandatory)" << std::endl;
    seqan::ArgumentParser parser;

    Options options;

    setupArgumentParser(parser, options);
    ArgumentParser::ParseResult res;
    res = parse(options, parser, argc, argv); // Fill the options structure with the standard and the acquired arguments

    if (res != ArgumentParser::PARSE_OK)  // Check the arguments
        return res == ArgumentParser::PARSE_ERROR;

    TRnaVect rnaSeqsRef;
    TRnaVect rnaSeqs; // Define the structure that will store all the input RNAs sequences

    _V(options, "Open the input file and fill the TRnaVect data structure");
    readRnaRecords(options, options.inFile ,rnaSeqs);
    _V(options, "Analyse the rnaSeqs structure and check which info are required for the selected analysis");

// Add the weight interaction edges vector map in the data structure
    bppInteractionGraphBuild(options, rnaSeqs);
//  Create the alignment data structure that will be used to store all the alignments
    TRnaAlignVect rnaAligns;
    if( options.inFileRef == "" )
    {
        std::cout << "fasta file name of reference is not available"
                "create the alignment data structure on a single file" << std::endl;
        alignVectorBuild(rnaAligns, rnaSeqs, options);
    } else
    {
        std::cout << "fasta file name is " << options.inFileRef << std::endl; //FIXME if input file is not provided the program is stack
        readRnaRecords(options, options.inFileRef, rnaSeqsRef);
        bppInteractionGraphBuild(options, rnaSeqsRef);
        alignVectorBuild(rnaAligns, rnaSeqs, rnaSeqsRef, options);
    }
    if(options.verbose > 2)
        std::cout << rnaSeqs[0].bppMatrGraphs[0].inter << std:: endl;
    StringSet<TAlign> alignsSimd;
    String<TScoreValue> resultsSimd;
    firstSimdAligns(resultsSimd, alignsSimd, rnaAligns, options);
    return 0;
}

