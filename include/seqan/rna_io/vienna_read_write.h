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
// Authors: Joerg Winkler <j.winkler@fu-berlin.de>
// ==========================================================================
// This file contains routines to read/write Vienna format files
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_VIENNA_READ_WRITE_H_
#define SEQAN_INCLUDE_SEQAN_RNA_IO_VIENNA_READ_WRITE_H_

#include <seqan/stream.h>
#include <seqan/rna_io/dot_bracket_read_write.h>  // for bracket-graph transformation
#include <stack>
#include <array>

namespace seqan{

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// --------------------------------------------------------------------------
// Tag Vienna
// --------------------------------------------------------------------------

struct Vienna_;
typedef Tag<Vienna_> Vienna;

// --------------------------------------------------------------------------
// Class Magicheader
// --------------------------------------------------------------------------

template <typename T>
struct MagicHeader<Vienna, T> :
    public MagicHeader<Nothing, T> {};

// ==========================================================================
// Metafunctions
// ==========================================================================

// --------------------------------------------------------------------------
// Metafunction FileExtensions
// --------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Vienna, T>
{
    static char const * VALUE[1];
};

template <typename T>
char const * FileExtensions<Vienna, T>::VALUE[1] =
{
    ".dbv"      // default output extension
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord(); RnaRecord, Vienna
// ----------------------------------------------------------------------------
template <typename TForwardIter>
inline void
readRecord(RnaRecord & record, SEQAN_UNUSED RnaIOContext &, TForwardIter & iter, Vienna const & /*tag*/)
{
    std::string buffer;
    clear(record);

    // read name (and offset)
    skipOne(iter);                                                      // ">" symbol
    readUntil(buffer, iter, IsNewline());
    std::string::size_type pos = buffer.find_last_of('/');
    if (pos == std::string::npos)
    {
        record.name = buffer;
    }
    else
    {
        record.name = buffer.substr(0, pos);
        std::string posStr{buffer.substr(pos + 1)};
        pos = posStr.find('-');
        if (pos == std::string::npos || !lexicalCast(record.offset, posStr.substr(0, pos))) {
            record.name = buffer;
            record.offset = 1;
        }
    }
    clear(buffer);

    // read sequence
    skipOne(iter);
    readUntil(record.sequence, iter, IsNewline());
    record.seqLen = length(record.sequence);

    // read bracket string and build graph
    skipOne(iter);
    readUntil(buffer, iter, IsWhitespace());
    if (length(buffer) != record.seqLen)
        throw std::runtime_error("ERROR: Bracket string must be as long as sequence.");

    TRnaRecordGraph graph;
    for (unsigned idx = 0; idx < length(buffer); ++idx)
        addVertex(graph);

    std::stack<unsigned> stack;
    for(unsigned idx = 0; idx < length(buffer); ++idx)
    {
        if (buffer[idx] == '(')
        {
            stack.push(idx);
        }
        else if (buffer[idx] == ')')
        {
            if (!stack.empty())
            {
                addEdge(graph, idx, stack.top(), 1.);
                stack.pop();
            }
            else
            {
                throw ParseError("Invalid bracket notation: unpaired closing bracket");
            }
        }
    }
    if(!stack.empty())
        throw ParseError("Invalid bracket notation: unpaired opening bracket");
    append(record.fixedGraphs, RnaInterGraph(graph));
    clear(buffer);

    // read energy if present
    if(!atEnd(iter))
    {
        skipUntil(iter, NotFunctor<IsWhitespace>());
        if(*iter == '(')
        {
            skipOne(iter);
            readUntil(buffer, iter, EqualsChar<')'>());
            if (!lexicalCast(record.energy, buffer))
                throw BadLexicalCast(record.energy, buffer);
            clear(buffer);
        }
    }
}


// ----------------------------------------------------------------------------
// Function writeRecord(); RnaRecord, Vienna
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target, RnaRecord const & record, SEQAN_UNUSED RnaIOContext &, Vienna const & /*tag*/)
{
    if (empty(record.sequence) && length(rows(record.align)) != 1)
        throw std::runtime_error("ERROR: Vienna formatted file cannot contain an alignment.");
    if (length(record.fixedGraphs) != 1)
        throw std::runtime_error("ERROR: Vienna formatted file cannot contain multiple structure graphs.");

    Rna5String const sequence = empty(record.sequence) ? source(row(record.align, 0)) : record.sequence;

    // write opening character for new record entry
    writeValue(target, '>');
    // write name
    write(target, record.name);
    // write index beg-end
    writeValue(target, '/');
    write(target, record.offset);
    writeValue(target, '-');
    write(target, record.offset + record.seqLen - 1);
    writeValue(target, '\n');

    // write sequence
    write(target, sequence);
    writeValue(target, '\n');

    // write bracket string
    TRnaRecordGraph const & graph = record.fixedGraphs[0].inter;
    std::string bracketStr;
    resize(bracketStr, numVertices(graph), ' ');
    std::stack<unsigned> stack;

    for (unsigned idx = 0; idx < length(bracketStr); ++idx) // write pairs in bracket notation
    {
        if (degree(graph, idx) == 0)                    // unpaired
        {
            bracketStr[idx] = '.';
            continue;
        }

        TRnaAdjacencyIterator adj_it(graph, idx);
        if (idx < value(adj_it))                        // open bracket
        {
            bracketStr[idx] = '(';
            stack.push(value(adj_it));
        }
        else                                            // close bracket
        {
            bracketStr[idx] = ')';
            if (stack.empty())
                SEQAN_FAIL("Cannot reach here.");
            if (stack.top() == idx)
                stack.pop();
            else
                throw std::runtime_error("ERROR: Vienna format does not allow pseudoknots.");
        }
    }
    write(target, bracketStr);

    // write energy
    if (record.energy != 0.0f)
    {
        write(target, " (");
        write(target, record.energy);
        write(target, ")\n");
    }
}

}
#endif // SEQAN_INCLUDE_SEQAN_RNA_IO_VIENNA_READ_WRITE_H_
