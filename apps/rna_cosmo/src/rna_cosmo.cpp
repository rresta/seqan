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
    std::clock_t begin = std::clock();
    std::chrono::steady_clock::time_point beginChrono = std::chrono::steady_clock::now();

    // Read input files
    RnaStructContents reducedContents;
    _readMultiStructRnaInputFile(reducedContents, options.inFile, options);

// add the weight interaction edges vector map in the data structure using Vienna package
// bppInteractionGraphBuild(reducedContents.records, options);

    _VV(options, "Found " << length(reducedContents.records) << " sequences in the input file.");
    for (RnaRecord & currentRecord : reducedContents.records)
    {
        unsigned totalConsEdges=0;  //counter for the edges present in all the fixed graphs
        unsigned totalOutEdges=0;   //counter for all the edges of the Consensus Graph
        unsigned notPaired=0;       //counter for the not paired bases of the Consensus Graph
        unsigned initEdges=0;       //counter for the edges present only in the first fixed graph
        RnaStructureGraph consGraph;                                            //Creation of the Consensus Graph
        consGraph = currentRecord.fixedGraphs[0];
        for (unsigned j = 0; j < currentRecord.seqLen; ++j) {
            RnaAdjacencyIterator adj_it(consGraph.inter, j);           //Access to the base pair of the Consensus Graph
            if (degree(consGraph.inter, j) != 0 && (value(adj_it) > j)) {
                assignCargo(findEdge(consGraph.inter, j, value(adj_it)), options.firstEdgeWeight);
            }
        }
        for (unsigned k = 1; k < length(currentRecord.fixedGraphs); ++k)       //for the other graphs
        {
            for (unsigned l = 0; l < currentRecord.seqLen; ++l) {
                RnaAdjacencyIterator adj_it_fixedGraph(currentRecord.fixedGraphs[k].inter, l);
                if (atEnd(adj_it_fixedGraph) || l < value(adj_it_fixedGraph)) //control for the value of the adj list
                    continue;

                if (findEdge(consGraph.inter, l, value(adj_it_fixedGraph)) != 0)  //if the edge already exists
                {
                    assignCargo(findEdge(consGraph.inter, l, value(adj_it_fixedGraph)),
                                getCargo(findEdge(consGraph.inter, l, value(adj_it_fixedGraph))) +
                                options.edgeStepWeight);
                }
                else
                {
                    addEdge(consGraph.inter, l, value(adj_it_fixedGraph), options.firstEdgeWeight);
                }
            }
        }
        _VVV(options, consGraph.inter);
        _VV(options, "Cargo Edge list:")
        for (unsigned m = 0; m < currentRecord.seqLen; ++m)
        {
            for (RnaAdjacencyIterator adj_it_cons(consGraph.inter, m); !atEnd(adj_it_cons); goNext(adj_it_cons))
            {
                if (value(adj_it_cons)>m)
                {
                double const edgeWeight = getCargo(findEdge(consGraph.inter, m, value(adj_it_cons)));
                _VVV(options,"Source: " << m << ", Target: " << value(adj_it_cons) << ", Cargo: " << edgeWeight);
                if (edgeWeight==options.firstEdgeWeight)
                    ++initEdges;
                if (edgeWeight == options.firstEdgeWeight+(length(currentRecord.fixedGraphs)-1)*options.edgeStepWeight)
                    ++totalConsEdges;
                ++totalOutEdges;
                }

            }
        if (degree(consGraph.inter, m) == 0)
            ++notPaired;
        }
        append(currentRecord.fixedGraphs, consGraph);
        _VV(options,  "Total Edges Count: " << totalOutEdges);
        _VV(options, "Sequence Length: " << currentRecord.seqLen << ", Unpaired bases: " << notPaired << ", "
                "Edges with initialization cargo: " << initEdges << ", Total Consensus Edges: " << totalConsEdges);
        _VVV(options, "Number of graphs after appending the consensus graph: " << length(currentRecord.fixedGraphs)
                                                                               << "\n");

    }

// timer stop
    std::chrono::steady_clock::time_point endChrono= std::chrono::steady_clock::now();
    std::clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

// Print elapsed time
    _VV(options, "\nTime difference chrono = " << std::chrono::duration_cast<std::chrono::seconds>
            (endChrono - beginChrono).count()); //std::chrono::microseconds
    _VV(options, "\nTime difference = " << elapsed_secs);

    return 1;

    //SHAPE
    //Function that controls in RMDB if there are SHAPE of the input sequence

    //if there is a shape

    //if there is a secondary structure: if the edge exist assign SHAPEweight, if it doesn't exist edge creation
    // assign SHAPEweight

    //else modifying the existing ones set a weight to existing edges according to reactivity and reactivity error

    //for all the vertex of the output graph

            //if there are no edges delete the ones of the bppMatrixGraph


// CODE SHOULD BE ADDED HERE
    _V(options, "STEP 1: Create functions that take the RNA seq and "
            "generate many fixed structures using several "
            "combinations of tools and parameters");
    _V(options, "These functions will generate command lines for the various tools");
    _V(options, "STEP 2: Fixed Graph of the outputs of the tools");
    _V(options, "STEP 3: Function that controls in RMDB if there are SHAPE of the input sequence");
    _V(options, "STEP 4: Function that generates a consensus structure module");
    _V(options, "STEP 5: Generate EBPSEQ file format (output)");
    _V(options, "STEP 6: Visualization on jVitz");


    return 0;
}

