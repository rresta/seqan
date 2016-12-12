// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Joerg Winkler <j.winkler@fu-berlin.de>
// ==========================================================================
// Implementation for maximum weighted matching for general graphs.
// The implemented algorithm is the heuristic from Drake/Hougardy (2005),
// which has a performance ratio of 2/3 - epsilon.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_MAXIMUM_WEIGHTED_MATCHING_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_MAXIMUM_WEIGHTED_MATCHING_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function maximumWeightedMatchingGreedy()
// ----------------------------------------------------------------------------

template <typename TCargo>
TCargo maximumWeightedMatchingGreedy(Graph<Undirected<TCargo> > const & graph)
{
    typedef typename Iterator<Graph<Undirected<TCargo> >, EdgeIterator>::Type TEdgeIter;
    typedef typename Iterator<Graph<Undirected<TCargo> >, AdjacencyIterator>::Type TAdjacIterator;
    typedef typename VertexDescriptor<Graph<Undirected<TCargo> > >::Type TVertex;

    Graph<Undirected<TCargo> > tmpGraph (graph); // copy for not changing the original graph
    TCargo maxWeight {};
    TVertex src;
    TVertex trg;

    while (numEdges(tmpGraph) > 0)
    {
        TCargo mw {};
        TEdgeIter edgeIt(tmpGraph);
        for (; !atEnd(edgeIt); goNext(edgeIt))
        {
            if (getCargo(*edgeIt) > mw)
            {
                mw = getCargo(*edgeIt);
                src = getSource(*edgeIt);
                trg = getTarget(*edgeIt);
            }
        }

        TAdjacIterator adj1(tmpGraph, src);
        for (; !atEnd(adj1); goNext(adj1))
            removeEdge(tmpGraph, src, *adj1);
        TAdjacIterator adj2(tmpGraph, trg);
        for (; !atEnd(adj2); goNext(adj2))
            removeEdge(tmpGraph, trg, *adj2);
        removeVertex(tmpGraph, src);
        removeVertex(tmpGraph, trg);
        maxWeight += mw;
    }
    return maxWeight;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_MAXIMUM_WEIGHTED_MATCHING_H_
