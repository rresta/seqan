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
// Graph boost headers
// ----------------------------------------------------------------------------

#include <boost/graph/adjacency_list.hpp> // for customizable graphs
#include <boost/graph/directed_graph.hpp> // A subclass to provide reasonable arguments to adjacency_list for a typical directed graph
#include <boost/graph/undirected_graph.hpp>// A subclass to provide reasonable arguments to adjacency_list for a typical undirected graph
#include <boost/graph/graphviz.hpp>

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

void testGraph();
void AdjacencyList();
void UndirectedGraph();
void DirectedGraph();

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


    // Read input files
    RnaStructContents redContents;
    _readMultiStructRnaInputFile(redContents, options.inFile, options);

// add the weight interaction edges vector map in the data structure using Vienna package
// bppInteractionGraphBuild(redContents.records, options);

    //Creation of the output graph by assigning options.firstEdgeWeight to all the edges of the first graph
    //assignCargo (findEdge(redContents.records[0].fixedGraphs[0], 1, 3), options.firstEdgeWeight);
    /////////////////////////////////////////////////////////////


    //for all the records

        //for all the graphs of the record

            //if it is the first assign options.firstEdgeWeight to all the edges of the first graph
            //assignCargo(findEdge(graph.inter, vertex1, vertex2), newcargo));

            //else (for all the other graphs)
                //if the edge exist add options.edgeStepWeight to the existing weight
                //findEdge
                //else edge creation and assign cargo options.firstEdgeWeight
                //addEdge

    _VVV(options, "reduced Contents records length: " << length(redContents.records));
    for (unsigned i = 0; i < length(redContents.records); ++i)
    {
        unsigned countOrigEdges=0;
        unsigned countConsEdges=0;
        unsigned countNewEdges=0;
        RnaStructureGraph consGraph;
        consGraph = redContents.records[i].fixedGraphs[0];
        for (unsigned j = 0; j < redContents.records[i].seqLen; ++j)
        {
            RnaAdjacencyIterator adj_it(consGraph.inter, j);
            if (degree(consGraph.inter, j) != 0 and (value(adj_it) > j))
            {
                ++countOrigEdges;
                assignCargo(findEdge(consGraph.inter, j, value(adj_it)), options.firstEdgeWeight);
                //std::cout << getCargo(findEdge(consGraph.inter, j, value(adj_it)));
            }
        }
        if (length(redContents.records[i].fixedGraphs) > 1) // in this case we have a consensus graph
        {
            for (unsigned k = 1; k < length(redContents.records[i].fixedGraphs); ++k)   //for the other graphs
            {
                for (unsigned l = 0; l < redContents.records[i].seqLen; ++l) {
                    RnaAdjacencyIterator adj_it(consGraph.inter, l);
                    if (degree(redContents.records[i].fixedGraphs[k].inter, l) != 0 and
                        degree(consGraph.inter, l) != 0 and value(adj_it) > l)
                    {
                        ++countConsEdges;
                        assignCargo(findEdge(consGraph.inter, l, value(adj_it)), getCargo(findEdge(consGraph.inter, l, value(adj_it))) + options.edgeStepWeight);
                    }
                    else if (degree(redContents.records[i].fixedGraphs[k].inter, l) != 0 and value(adj_it) > l)
                    {
                        ++countNewEdges;
                        //addEdge(consGraph.inter, l, value(adj_it), options.firstEdgeWeight);
                    }
                }
            }
        }
            append(redContents.records[i].fixedGraphs, consGraph);
            _VVV(options, "Consensus Edges count: " << countConsEdges << " Original Edges count: " << countOrigEdges << " New Edges count: " << countNewEdges);
            _VVV(options, "Edges Weights for consensus graph for record " << i << " :");
            for (unsigned m = 0; m < redContents.records[i].seqLen; ++m)
            {
                if (degree(consGraph.inter, m) != 0)
                {
                    RnaAdjacencyIterator adj_it(consGraph.inter, m);
                    _VVV(options, getCargo(findEdge(consGraph.inter, m, value(adj_it))));
                }
            }
            if (countConsEdges != 0)
            _VVV(options, "Number of graphs for record " << i << " after appending the consensus graph: " << length(redContents.records[i].fixedGraphs) << "\n");
    }
        _VVV(options, "First edge weight: " << options.firstEdgeWeight);
        _VVV(options, "Added edge weight: " << options.edgeStepWeight);

    return 1;





    /*        for (StructureGraphAdjacencyIterator adj_it(consGraph.inter); !atEnd(adj_it); goNext(adj_it))
        {
            if (degree(consGraph.inter, adj_it) != 0)
            {
            std::cout << value(adj_it) << "\t";
            assignCargo(findEdge(consGraph.inter, i , value(adj_it)), options.firstEdgeWeight);
//            writeValue(adj_it, options.firstEdgeWeight);
            std::cout << value(adj_it) << std::endl;
//            int edgeWeight1 = cargo(findEdge(graph1, line.first, value(adj_it1)));

        }*/

    //SHAPE
    //Function that controls in RMDB if there are SHAPE of the input sequence

    //if there is a shape

            //if there is a secondary structure: if the edge exist assign SHAPEweight, if it doesn't exist edge creation assign SHAPEweight

            //else modifying the existing ones set a weight to existing edges according to reactivity and reactivity error

    //for all the vertex of the output graph

            //if there are no edges delete the ones of the bppMatrixGraph



// timer start
    std::clock_t begin = std::clock();
    std::chrono::steady_clock::time_point beginChrono = std::chrono::steady_clock::now();

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




    testGraph();


// timer stop
    std::chrono::steady_clock::time_point endChrono= std::chrono::steady_clock::now();
    std::clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    double mwmtime = 0.0;
    double lemtime = 0.0;
    double boutime = 0.0;


// Print elapsed time
    _VV(options, "\nTime difference chrono = " << std::chrono::duration_cast<std::chrono::seconds>(endChrono - beginChrono).count()); //std::chrono::microseconds
    _VV(options, "\nTime difference = " << elapsed_secs);

    return 0;
}

void testGraph()
{
    AdjacencyList();
    UndirectedGraph();
    DirectedGraph();
}


void AdjacencyList()
{
    /* Method 1: The most generic
    * The generic class for a graph in Boost is adjacency_list.
    * Up to 7 template parameters can be given, for example:
    * typedef boost::adjacency_list<     boost::listS,             // out-edges stored in a std::list
    *                       boost::vecS,             // vertex set stored here
    *                       boost::undirectedS,    // bidirectional graph.
    *                       boost::no_property,              // vertex properties
    *                       EdgeWeightProperty,       // edge properties
    *                       boost::no_property,       // graph properties
    *                       boost::listS              // edge storage
    *                       > graph_t;
    *
    * The 'S' at the end of the choices (vecS, etc.) stands for 'S'elector.
    */

    {
        // Construct a graph with the vertices container as a vector
        typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> Graph;
        Graph g(3); // Create a graph with 3 vertices.

        // The graph behaves as a new user would expect if the vertex container type is vector. That is, vertices can be indexed with an unsigned int.
        boost::add_edge(0, 1, g);
        boost::add_edge(1, 2, g);
    }

    {
        // Construct a graph with the vertices container as a set
        typedef boost::adjacency_list<boost::vecS, boost::setS, boost::bidirectionalS> Graph;

        // Since the vertex container type is not a vector, the vertices can NOT be indexed with an unsigned int. I.e. the following will not work:
        //Graph g(3); // 3 vertices
        //boost::add_edge(0, 1, g);
        //boost::add_edge(1, 2, g);

        // Instead, you must add vertices individually so that you get a handle to them (a way to reference them, Boost calls this a "vertex_descriptor"),
        // and then add the edges by referencing these descriptors. Note this is a very generic method, so it would work just as well with a vecS vertex container.

        Graph g; // Create a graph.
        Graph::vertex_descriptor v0 = boost::add_vertex(g);
        Graph::vertex_descriptor v1 = boost::add_vertex(g);
        Graph::vertex_descriptor v2 = boost::add_vertex(g);

        boost::add_edge(v0, v1, g);
        boost::add_edge(v1, v2, g);

    }
}

void UndirectedGraph()
{
    // undirected_graph is a subclass of adjacency_list which gives you object oriented access to functions like add_vertex and add_edge, which makes the code easier to understand. However, it hard codes many of the template parameters, so it is much less flexible.

    typedef boost::undirected_graph<> Graph;
    Graph g;
    Graph::vertex_descriptor v0 = g.add_vertex();
    Graph::vertex_descriptor v1 = g.add_vertex();

    g.add_edge(v0, v1);
}

void DirectedGraph()
{
    // directed_graph is a subclass of adjacency_list which gives you object oriented access to functions like add_vertex and add_edge, which makes the code easier to understand. However, it hard codes many of the template parameters, so it is much less flexible.

    typedef boost::directed_graph<> Graph;
    Graph g;
    Graph::vertex_descriptor v0 = g.add_vertex();
    Graph::vertex_descriptor v1 = g.add_vertex();

    g.add_edge(v0, v1);
}
