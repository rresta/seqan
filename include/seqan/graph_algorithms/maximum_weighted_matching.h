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

// #include <experimental/random>

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

// Compute greedy MWM (performance ratio 1/2)
template <typename TCargo>
TCargo maximumWeightedMatchingGreedy(Graph<Undirected<TCargo> > const & graph)
{
    typedef Graph<Undirected<TCargo> > TUGraph;
    typedef typename EdgeDescriptor<TUGraph>::Type TEdge;
    typedef typename Iterator<TUGraph, EdgeIterator>::Type TEdgeIter;
    typedef typename Iterator<TUGraph, AdjacencyIterator>::Type TAdjacIterator;
    typedef typename VertexDescriptor<TUGraph>::Type TVertex;

    // set up edge vector and bit vector for conflicting edges
    std::vector<TEdgeIter> edges;
    std::vector<bool> conflictFree;
    reserve(edges, numEdges(graph));
    resize(conflictFree, numEdges(graph), true);
    for (TEdgeIter edgeIt(graph); !atEnd(edgeIt); goNext(edgeIt))
    {
        edges.push_back(edgeIt);
    }

    // sort edges with respect to their weight, start with the highest
    std::sort(edges.begin(), edges.end(), [] (auto a, auto b) { return getCargo(*a) >= getCargo(*b); });

    TCargo maxWeight {};
    for (std::size_t idx = 0u; idx < length(edges); ++idx)
    {
        auto const & edge = *edges[idx];
        if (!conflictFree[edge->data_id])  // skip edge if conflict with a previous edge
            continue;

        maxWeight += getCargo(edge);  // edge is contained in the matching

        // mark all adjacent edges
        TVertex const & src = getSource(edge);
        for (TAdjacIterator ai(graph, src); !atEnd(ai); goNext(ai))
        {
            TEdge rmEdge = findEdge(graph, src, *ai);
            conflictFree[rmEdge->data_id] = false;
        }

        TVertex const & trg = getTarget(edge);
        for (TAdjacIterator ai(graph, trg); !atEnd(ai); goNext(ai))
        {
            TEdge rmEdge = findEdge(graph, trg, *ai);
            conflictFree[rmEdge->data_id] = false;
        }

        SEQAN_ASSERT(!conflictFree[edge->data_id]);
    }
    return maxWeight;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_MAXIMUM_WEIGHTED_MATCHING_H_
