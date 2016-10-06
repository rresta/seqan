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
// This file contains the variable definition and structures of seqan_laragu
// application.
// ==========================================================================

#ifndef _INCLUDE_INTERACTION_EDGES_H_
#define _INCLUDE_INTERACTION_EDGES_H_

// ----------------------------------------------------------------------------
// std headers
// ----------------------------------------------------------------------------

#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <omp.h>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

#include "vienna_rna.h"

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function bppInteractionGraphBuild()
// ----------------------------------------------------------------------------

template <typename TOption, typename TVect>
void bppInteractionGraphBuild(TOption const & options, TVect  & rnaSeqs)
{
//variable used to plot or not the energy matrix (for now is disabled but the function is tested)
    bool write_dot_plot = 0; //TODO decide if provide this option in the options
// Create a c-style string object for str:
    String<char, CStyle> cStr;
    for(unsigned i=0;i<length(rnaSeqs);++i)  // TODO Execute this part in PARALLEL
    {
        cStr = rnaSeqs[i].sequence;
//		CharString curSeq = rnaSeqs[0].seq;
        std::cout << toCString(cStr) << std::endl;

//TODO once is given support to introduce several BPP from the extended bpseq
// file or from the dot plot file in this position must be placed a condition
// capable to discriminate which bpp matrix to use
        if(length(rnaSeqs[i].bppMatrGraphs) == 0) // if(dotplot or extended bpseq data are not present)
        {
            computeBppMatrix(options, rnaSeqs[i]);
        }else{
// TODO read the BPP matrix from the rnaSeqs[i].bpp_matr_graphs field that contains all the structures acquired by the files
// TODO add the filtering step that involve the biological input and the majority voter of predicted structures
        }
/*

        std::cout << "Input Degree nodo 0 = " << inDegree(rnaSeqs[i].bpProb.interGraph, rnaSeqs[i].bpProb.uVertexVect[0]) << std::endl;
        std::cout << "Output Degree nodo 0 = " << outDegree(rnaSeqs[i].bpProb.interGraph, rnaSeqs[i].bpProb.uVertexVect[0]) << std::endl;
//		std::cout << rnaSeqs[i].bpProb.interGraph << std::endl;

        // Output distances of shortest paths
//	    Iterator<TUgraph, VertexIterator>::Type it(rnaSeqs[i].bpProb.interGraph);
        typedef Iterator<TUgraph, OutEdgeIterator>::Type TOutEdgeIterator;
//TODO place the printing part in an external function printWeightedGraph()
        for(unsigned j=0;j<length(rnaSeqs[i].bpProb.uVertexVect);++j)
        {
            TOutEdgeIterator it(rnaSeqs[i].bpProb.interGraph, rnaSeqs[i].bpProb.uVertexVect[j]);
            for(;!atEnd(it);goNext(it)) {
                //		    std::cout << "Edge id = " << it;
// build the directed graph that will be updated over the iterations
//				addEdge(rnaSeqs[i].bpProb.interGraphUpdated, source(*it), target(*it), cargo(*it));
//				addEdge(rnaSeqs[i].bpProb.interGraphUpdated, target(*it), source(*it), cargo(*it));
                std::cout << "source = " << source(*it) << "\ttarget = " << target(*it) << "\tcargo = " << cargo(*it) << std::endl;
            }
        }

//	    while (!atEnd(it))
//	    {
//	        std::cout << "Distance from 0 to " << getValue(it) << ": ";
//	        std::cout << getProperty(rnaSeqs[i].bpProb.uVertexVect, getValue(it)) << std::endl;
//	        goNext(it);
//	    }

//		dEdge = addEdge(g,v0,v1);
*/
    }
}

#endif //_INCLUDE_INTERACTION_EDGES_H_