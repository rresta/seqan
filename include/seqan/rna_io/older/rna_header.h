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
// This file contains routines to read RNA sequence containing file headers
// ==========================================================================

#ifndef SEQAN_RNA_FORMAT_HEADER_H
#define SEQAN_RNA_FORMAT_HEADER_H


namespace seqan {
// ----------------------------------------------------------------------------
// Class RNAHeaderRecord
// ----------------------------------------------------------------------------

/*!
 * @class RNAHeaderRecord
 * @headerfile <seqan/RNA_io.h>
 * @brief Store amount of bases/name for RNA header records.
 *
 * @signature class RNAHeaderRecord;
 *
 * @var CharString RNAHeaderRecord::amount;
 * @brief Amount of bases given in the header record.
 *
 * @var CharString RNAHeaderRecord::name;
 * @brief Name of RNA sequence given in the header record.
 */

/*!
 * @fn RNAHeaderRecord::RNAHeaderRecord
 * @brief Constructor
 *
 * @signature RNAHeaderRecord::RNAHeaderRecord();
 * @signature RNAHeaderRecord::RNAHeaderRecord(amount, name);
 *
 * @param[in] amount   Amount of sequences of the header record, @link int32_t @endlink.
 * @param[in] name     Name of sequence of the header record, @link CharString @endlink.
 */

/*!
 * @fn RNAHeaderRecord#clear
 *
 * @brief Clear a RNAHeaderRecord.
 * @signature void clear(record);
 *
 * @param[in,out] record The RNAHeaderRecord to clear.
 */

class RNAHeader
{
public:
    // Amount of records.
    int32_t amount;    
    /////////Do I want a beginning poition??
    float energy;    
    // Record's name.
    CharString name;

    // Default constructor.
    RNAHeader()
    {}

    // Construct directly with amount/name.
    RNAHeader(int32_t const & amount, int32_t const & energy, CharString const & name) :
            amount(amount), energy(energy), name(name)
    {}
};

// ============================================================================
// Functions
// ============================================================================

inline void clear(RNAHeader & record)
{
    clear(record.name);
}


} //namespace seqan
#endif //SEQAN_RNA_FORMAT_HEADER_H