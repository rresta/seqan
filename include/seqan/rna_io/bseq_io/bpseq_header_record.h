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
// Author: Gianvito Urgese <gianvito.urgese@polito.it>
// ==========================================================================

// TODO(holtgrew): Parse more than just the key/value pair.

#ifndef SEQAN_INCLUDE_SEQAN_BPSEQ_IO_BPSEQ_HEADER_RECORD_H_
#define SEQAN_INCLUDE_SEQAN_BPSEQ_IO_BPSEQ_HEADER_RECORD_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class BpseqHeaderRecord
// ----------------------------------------------------------------------------

/*!
 * @class BpseqHeaderRecord
 * @headerfile <seqan/bpseq_io.h>
 * @brief Store key/value pair for Bpseq header records.
 *
 * @signature class BpseqHeaderRecord;
 *
 * @var CharString BpseqHeaderRecord::key;
 * @brief Key of the header record.
 *
 * @var CharString BpseqHeaderRecord::value;
 * @brief Value of the header record.
 */

/*!
 * @fn BpseqHeaderRecord::BpseqHeaderRecord
 * @brief Constructor
 *
 * @signature BpseqHeaderRecord::BpseqHeaderRecord();
 * @signature BpseqHeaderRecord::BpseqHeaderRecord(key, value);
 *
 * @param[in] key   Key of the header record, @link CharString @endlink.
 * @param[in] value Key of the header record, @link CharString @endlink.
 */

/*!
 * @fn BpseqHeaderRecord#clear
 *
 * @brief Clear a BpseqHeaderRecord.
 * @signature void clear(record);
 *
 * @param[in,out] record The BpseqHeaderRecord to clear.
 */

class BpseqHeaderRecord
{
public:
    // Record's key.
    CharString key;
    // Record's value.
    CharString value;

    // Default constructor.
    BpseqHeaderRecord()
    {}

    // Construct directly with key/value.
    BpseqHeaderRecord(CharString const & key, CharString const & value) :
            key(key), value(value)
    {}
};

// ============================================================================
// Functions
// ============================================================================

inline void clear(BpseqHeaderRecord & record)
{
    clear(record.key);
    clear(record.value);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BPSEQ_IO_BPSEQ_HEADER_RECORD_H_
