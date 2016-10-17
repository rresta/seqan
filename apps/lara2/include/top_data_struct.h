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

#ifndef _INCLUDE_TOP_DATA_STRUCT_H_
#define _INCLUDE_TOP_DATA_STRUCT_H_

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/rna_io.h>
//#include <seqan/align_rna.h>

// ============================================================================
// Prerequisites
// ============================================================================

#define FASTA 0
#define FASTQ 1
#define BPSEQ 2
#define EBPSEQ 3
#define DBN 4
#define EDBN 5
#define UNKNOWN 6

#define LOGARITHMIC 0
#define SCALE       1
#define ORIGINAL    2
#define RIBOSUM     3

// ============================================================================
// Macro utility
// ============================================================================

#define _V(_opt, _str) {if(_opt.verbose>0) std::cerr << _str << std::endl;}
#define _VV(_opt, _str) {if(_opt.verbose>1) std::cerr << _str << std::endl;}

// ============================================================================
// Types used in the program
// ============================================================================

typedef seqan::Rna5String TSequence;
typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;      // align type
typedef unsigned TPosition;
typedef double TScoreValue;
typedef seqan::CharString TString;
typedef Score<double, ScoreMatrix<Rna5, Default> > TScoreMatrix;
typedef seqan::Ribosum65N TRibosum;
typedef float TBioval;
typedef std::map<TPosition, TScoreValue> TMap;
typedef seqan::String<TMap > TMapLine;

/*
typedef float TCargo;
typedef Graph<Directed<TCargo> > TDgraph;
typedef VertexDescriptor<TDgraph>::Type TDVertexDescriptor;
typedef EdgeDescriptor<TDgraph>::Type TDEdgeDescriptor;
typedef Graph<Undirected<TCargo> > TUgraph;
typedef VertexDescriptor<TUgraph>::Type TUVertexDescriptor;
typedef EdgeDescriptor<TUgraph>::Type TUEdgeDescriptor;
*/

//typedef Map<TPosition, TDVertexDescriptor> TMapDGraph;
//typedef seqan::String<TMapDGraph > TMapDGraphStr;

//typedef Map<TPosition, TUVertexDescriptor> TMapUGraph;
//typedef seqan::String<TMapUGraph > TMapUGraphStr;
/*
template <typename TString, typename TPosition>
struct fixedStructElement {
    TString method; // place the method and parameters used to compute the structure
//	seqan::String<unsigned> structure;
    seqan::String<TPosition> seqPos;
    seqan::String<TPosition> interPos;
};
template <typename TString, typename TBioval>
struct bioValStructElement {
    TString method; // place the method and parameters used to compute the structure
//	seqan::String<TBioval> val;
    seqan::String<TBioval> val;
};

struct vectGraphElement {
    std::vector<TUVertexDescriptor> uVertexVect;
    TUgraph interGraph; // this graph represents all the computed interaction edges
    std::vector<TDVertexDescriptor> dVertexVect;
    TDgraph interGraphUpdated;
};
*/

template <typename TSequence, typename TString, typename TPosition, typename TBioval, typename TMapLine>
struct RnaStructSeq {
    TSequence seq;  // from fasta input
    TString qual;  // from fasta input
    TString id;  // from fasta input
    TString info; // raw info from bpseq or ebpseq input
//    seqan::String<fixedStructElement<TString, TPosition> > structPairMate; //TODO use this field to collect all the structural information of this sequence
//    seqan::String<bioValStructElement<TString, TBioval> > structBioVal;
//    vectGraphElement bpProb;
};

typedef RnaStructSeq<TSequence, TString, TPosition, TBioval, TMapLine> TRnaStruct;
//typedef std::vector<TRnaStruct > TRnaVect;
typedef std::vector<seqan::RnaRecord > TRnaVect;

struct lambStruct
{
    // String with size seq1 storing all the aligned lines
    TMap map; //mapLine;
    unsigned indexBest;
};

struct RnaStructAlign
{
//public:
    seqan::RnaRecord rna1; // If we have problems with the memory the index of the TRnaVect can instead saved
    seqan::RnaRecord rna2;
// The best computed alignment is saved in these fields
    TAlign bestAlign;
    TScoreValue bestAlignScore;
// String with size seq1 storing all the aligned lines
    seqan::String<lambStruct > lamb;
};// rnaStructAlign;

typedef RnaStructAlign TRnaAlign;
typedef std::vector<TRnaAlign> TRnaAlignVect;


#endif //_INCLUDE_TOP_DATA_STRUCT_H_
