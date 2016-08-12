// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Lily Shellhammer <>
// ==========================================================================
// This file contains routines to read and write to Rdat format files (.rdat)
// ==========================================================================

#ifndef SEQAN_RDAT_FORMAT_READ_H
#define SEQAN_RDAT_FORMAT_READ_H

#include <seqan/stream.h>


/* IMPLEMENTATION NOTES

Rdat FORMAT example:

=> HEADER START : number of bases in the sequence
=> HEADER END: title of the structure
=> Each line has information about a base pair in the sequence
	Each line is a base, with this order of information:
		- Base number: index n
		- Base (A, C, G, T, U, X)
		- Index n-1
		- Index n+1
		- Number of the base to which n is paired. No pairing is indicated by 0 (zero).
		- Natural numbering. Rdatstructure ignores the actual value given in natural numbering, 
			so it is easiest to repeat n here.

CT Files can hold multiple structures of a single sequence.
This is done by repeating the format for each structure without any blank lines between structures. 

HEADER
 N  SEQUENCE   N-1  	 N+1	J POSITION  N  
 1 	G       	0    	2   	72    		1
 2 	C       	1    	3   	71    		2
 3 	G       	2    	4   	70    		3


*/
namespace seqan{

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// ============================================================================
// Forwards
// ============================================================================
// --------------------------------------------------------------------------
// Tag Rdat
// --------------------------------------------------------------------------

struct Rdat_;
typedef Tag<Rdat_> Rdat;

template <typename T>
struct MagicHeader<Rdat, T> :
    public MagicHeader<Nothing, T> {};

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction FileExtensions
// --------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Rdat, T>
{
    static char const * VALUE[1];
};
template <typename T>
char const * FileExtensions<Rdat, T>::VALUE[1] =
{
    ".rdat"      // default output extension
};


// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function readHeader(); RdatHeader
// ----------------------------------------------------------------------------

template <typename TForwardIter>
inline void
readHeader(RdatHeader & header, RNAIOContext & context, TForwardIter & iter, Rdat const & /*tag*/)
{
    typedef OrFunctor<IsTab, IsNewline> TNextEntry;

    //RDAT_VERSION    0.32
    //NAME        MedLoop

    clear(header);
    clear(context.buffer); 
    skipUntil(iter, NotFunctor<IsWhitespace>()); 
    skipUntil(iter, IsWhitespace());

    readUntil(header.version,  iter, IsNewline());   

    skipUntil(iter, NotFunctor<IsWhitespace>()); 
    skipUntil(iter, IsWhitespace());

    readUntil(header.name,  iter, IsNewline());   
    
}


template <typename TForwardIter>
inline void 
readRecord(RdatRecord & record, RNAIOContext & context, TForwardIter & iter, Rdat const & /*tag*/)
{
//RDAT_VERSION 0.32
//NAME        MedLoop

    clear(context);
    
}


// ----------------------------------------------------------------------------
// Function writeRecord(Rdat);
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target, RdatRecord const & record, Rdat const & /*tag*/)     
{
    writeValue(target, ' ');    //All records start with a space
    appendNumber(target, i+1);
    write(target, record.base[i]);

}

// ----------------------------------------------------------------------------
// Function writeHeader(); RdatHeader
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeHeader(TTarget & target, RdatHeader const & header, Rdat const & /*tag*/)
{
    appendNumber(target, header.amount);
    writeValue(target, ' ');
    write(target, "ENERGY = ");
    writeValue(target, '\t');
    appendNumber(target, header.energy);
    writeValue(target, '\t');
    write(target, header.name);
    writeValue(target, '\n');
}


} //namespace seqan

#endif // SEQAN_RDAT_FORMAT_READ_H