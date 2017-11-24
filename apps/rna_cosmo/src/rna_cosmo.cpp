// ==========================================================================
//               RNA CoSMo - RNA Consensus Structure Module
// ==========================================================================
// Copyright (c) 2015-2017, Gianvito Urgese
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
// Author: Rossella Resta <s222385@studenti.polito.it>
// ==========================================================================
// This file contains the rna_cosmo application.
// ==========================================================================

#define SEQAN_RNA_COSMO

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
#include <ctime>

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
//#include "rna_cosmo_core.h"
#include "rna_cosmo_io.h"
#include "vienna_rna.h"

using namespace seqan;

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main (int argc, char const ** argv)
{
    if (argc == 1)
    {
        std::cout << "type " << argv[0] << " --help to get the parameters table (-i option is mandatory)" << std::endl;
        return 1;
    }

    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser, options);
    ArgumentParser::ParseResult res;
    res = parse(options, parser, argc, argv); // Fill the options structure with the standard and the acquired arguments
    if (res != ArgumentParser::PARSE_OK)  // Check the arguments
        return res == ArgumentParser::PARSE_ERROR;

    // timer start
    std::clock_t begin_clock = std::clock();
    std::chrono::steady_clock::time_point beginChrono = std::chrono::steady_clock::now();

    // Read input files
    RnaStructContents reducedContents;
    unsigned consensusCargo=0;
    double probabilityCargo=0;
    _readMultiStructRnaInputFile(reducedContents, options.inFile, options); //STEP 1
    // add the weight interaction edges vector map in the data structure using Vienna package
    bppInteractionGraphBuild(reducedContents.records, options);
    _VV(options, "Found " << length(reducedContents.records) << " sequences in the input file.");
    for (RnaRecord & currentRecord : reducedContents.records) {
        //std::cout << currentRecord.bppMatrGraphs[0].inter << std::endl;
        unsigned totalConsEdges = 0;  //counter for the edges present in all the fixed graphs
        unsigned notPaired = 0;       //counter for the not paired bases of the Consensus Graph
        unsigned initEdges = 0;       //counter for the edges present only in the first fixed graph
        RnaStructureGraph consensusGraph;
        RnaStructureGraph probabilityConsensusGraph;
        consensusGraph = currentRecord.fixedGraphs[0];              //Creation of the Consensus Graph
        Iterator<Graph<Undirected<double> >, EdgeIterator>::Type consEdgeIt(consensusGraph.inter);
        while (!atEnd(consEdgeIt)) {
            assignCargo(findEdge(consensusGraph.inter, sourceVertex(consEdgeIt), targetVertex(consEdgeIt)),
                        options.firstEdgeWeight);
            goNext(consEdgeIt);
        }
        for (unsigned k = 1; k < length(currentRecord.fixedGraphs); ++k)        //STEP 2
        {
            Iterator<Graph<Undirected<double> >, EdgeIterator>::Type fixedEdgeIt(currentRecord.fixedGraphs[k].inter);
            while (!atEnd(fixedEdgeIt)) {

                if (findEdge(consensusGraph.inter, sourceVertex(fixedEdgeIt), targetVertex(fixedEdgeIt)) != 0)
                {
                    assignCargo(findEdge(consensusGraph.inter, sourceVertex(fixedEdgeIt), targetVertex(fixedEdgeIt)),
                                getCargo(findEdge(consensusGraph.inter, sourceVertex(fixedEdgeIt),
                                                  targetVertex(fixedEdgeIt))) + options.edgeStepWeight);
                } else {
                    addEdge(consensusGraph.inter, sourceVertex(fixedEdgeIt), targetVertex(fixedEdgeIt),
                            options.firstEdgeWeight);
                }
                goNext(fixedEdgeIt);
            }
        }
       _VVV(options,"\nConsensus Graph" << consensusGraph.inter);

        Iterator<Graph<Undirected<double> >, EdgeIterator>::Type edgIt(consensusGraph.inter);
        for (edgIt; !atEnd(edgIt); goNext(edgIt)) {
            double const edgeWeight = getCargo(findEdge(consensusGraph.inter, sourceVertex(edgIt),
                                                        targetVertex(edgIt)));
            if (edgeWeight == options.firstEdgeWeight) {
                ++initEdges;
            }

            if (edgeWeight ==
                options.firstEdgeWeight + (length(currentRecord.fixedGraphs) - 1) * options.edgeStepWeight) {
                ++totalConsEdges;
            }
        }
        _VV(options, "Total Edges: " << numEdges(consensusGraph.inter));
        _VV(options, "Edges with initialization cargo: " << initEdges << ", Total Consensus Edges: " << totalConsEdges);
        probabilityConsensusGraph=consensusGraph;
        _VV(options, "Bpp Matrix Edges found in the consensus graph: ");
        Iterator<Graph<Undirected<double> >, EdgeIterator>::Type bppEdgIt(currentRecord.bppMatrGraphs[0].inter);
        for (bppEdgIt; !atEnd(bppEdgIt); goNext(bppEdgIt)) {                //STEP 3
            if (findEdge(consensusGraph.inter, sourceVertex(bppEdgIt), targetVertex(bppEdgIt)) != 0) {
                consensusCargo=getCargo(findEdge(probabilityConsensusGraph.inter,
                                                 sourceVertex(bppEdgIt),
                                                 targetVertex(bppEdgIt)));
                probabilityCargo=getCargo(findEdge(currentRecord.bppMatrGraphs[0].inter, sourceVertex(bppEdgIt),
                                                   targetVertex(bppEdgIt)));
                assignCargo(findEdge(probabilityConsensusGraph.inter, sourceVertex(bppEdgIt), targetVertex(bppEdgIt)),
                            probabilityCargo+consensusCargo);
                _VV(options, "sourceVertex = " << sourceVertex(bppEdgIt) << " targetVertex = " << targetVertex(bppEdgIt)
                                        << " probability = " << probabilityCargo << " getCargo probConsGraph:"
                                               << getCargo(findEdge(probabilityConsensusGraph.inter,
                                                                    sourceVertex(bppEdgIt),
                                                                    targetVertex(bppEdgIt))));

            }
        }
    }
//////////////////////////////////////////////////////////////////////
    //  RMDB files

    //  test/inputs/rdat_files/GLYCFN_SHP_0004.rdat //
    //  test/inputs/rdat_files/ADD140_1M7_0011.rdat
    //  test/inputs/rdat_files/TRNAPH_SHP_0002.rdat //
    //  test/inputs/rdat_files/TRP4P6_DMS_0008.rdat
    //  test/inputs/rdat_files/5SRRNA_1M7_0008.rdat

    if (isSet(parser, "inFileShape")) {
        std::ifstream rdatFile;
        rdatFile.open(options.inFileShape);
        _V(options, "RDAT FILE: " << options.inFileShape);
        Rna5String RDATseq;
        String<RnaStructureGraph> ifrGraph;
        StringSet<String<char> > stringSetReactivity;
        StringSet<String<char> > stringSetReactivityError;
        rdatContents(RDATseq, ifrGraph, stringSetReactivity, stringSetReactivityError, rdatFile, options);  //STEP 4

        bool seqFlag = false;
        for (RnaRecord &currentRecord : reducedContents.records) {
            if (currentRecord.sequence == RDATseq) {                   //STEP 5
                seqFlag = true;
                currentRecord.reactivity = stringSetReactivity;
                currentRecord.reactError = stringSetReactivityError;
                clear(stringSetReactivity);
                clear(stringSetReactivityError);
/*                typedef Iterator<StringSet<String<float> >, Standard>::Type TIterator;    //iteration on float strings
                for (TIterator it = begin(currentRecord.reactivity, Standard()); it != end(currentRecord.reactivity, Standard());
                     ++it)
                _VVV(options, "REACTIVITY " << position(it, currentRecord.reactivity) << ": " << *it << std::endl);
                for (TIterator it1 = begin(currentRecord.reactError, Standard()); it1 != end(currentRecord.reactError, Standard());
                     ++it1)
                _VVV(options, "REACTIVITY_ERROR " << position(it1, currentRecord.reactError) << ": " << *it1
                                                  << std::endl);*/
            }
        }
        if (!seqFlag) _V(options, "No sequence match Input File - Reference File .");
    } else{
        _V(options, "No file RDAT.")
    }
//////////////////////////////////////////////////////////////////////

    // timer stop
    std::chrono::steady_clock::time_point endChrono= std::chrono::steady_clock::now();
    std::clock_t end_clock = std::clock();
    double elapsed_secs = double(end_clock - begin_clock) / CLOCKS_PER_SEC;

    // Print elapsed time
    _VV(options, "\nTime difference chrono = " << std::chrono::duration_cast<std::chrono::seconds>
            (endChrono - beginChrono).count()); //std::chrono::microseconds
    _VV(options, "\nTime difference = " << elapsed_secs);

/////////////////////////////////////////////////////////////////////////////////////

    _V(options,"\nSTEPS DESCRIPTION at different levels of verbose -v."
            " ['STEP <number>' references present in the code as comments.]");
    _VV(options, "\nSTEP 1: Function that for each RNA record appends all the corresponding Fixed Graphs "
            "of the input file.");
    _VVV(options, "The Graphs are obtained by using different tools. The output will be an RNA Structure "
            "(reducedContents).");
    _VV(options, "\nSTEP 2: Creation of the Consensus Graph. Edge initialization -few, edge step weight -esw.(*)");
    _VVV(options, "The Cargo of the edges will be the result of the consensus.");
    _VV(options, "\nSTEP 3: Creation of the Consensus Probability Graph.(*)");
    _VVV(options, "Add probabilities of the bpp matrix to the cargo of the Consensus Graph.");
    _VV(options, "\nSTEP 4: Function that reads an .rdat file -ifr.");
    _VVV(options, "The output of the function will be the RNA sequence, the secondary structure, the String Set for "
            "Reactivity and the one for Reactivity Error.");
    _VV(options, "\nSTEP 5: Append Reactivty and Reactivity Error for the corresponding sequence in the initial "
            "RNA structure.");
    _VV(options, "\nSTEP 6 (TODO): Filter the bpp matrix.");
    _VV(options, "\nSTEP 7 (TODO): Generate EBPSEQ file format (output).");
    _VV(options, "\nSTEP 8 (TODO): Visualization on jVitz.");
    _VVV(options, "\nRESULTING STRUCTURES:\nConsensus Graph,\nProbability Consensus Graph,\nReactivity and Reactivity "
            "Error(**),\nSecondary Structure of the RDAT file(**).")
    _VV(options, "\nNOTES:\n(*)For each Record.");
    _VVV(options, "(**)Depends on the -ifr.");
    _VVV(options, "INPUT FILE FOR TEST: -i test/inputs/5seq_ipknot_RNAfold_RNAstruct.dbn\nRDAT FILE (example): -ifr "
            "test/inputs/rdat_files/ADD140_1M7_0011.rdat")

    return 0;
}

