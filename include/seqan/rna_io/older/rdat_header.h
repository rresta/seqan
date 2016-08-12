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
// This file contains routines to read RDAT header
// ==========================================================================

#ifndef SEQAN_RDAT_FORMAT_HEADER_H
#define SEQAN_RDAT_FORMAT_HEADER_H


namespace seqan {
// ----------------------------------------------------------------------------
// Class RdatHeaderRecord
// ----------------------------------------------------------------------------

/*!
 * @class RdatHeaderRecord
 * @headerfile <seqan/Rdat_io.h>
 * @brief Store amount of bases/name for Rdat header records.
 *
 * @signature class RdatHeaderRecord;
 *
 * @var CharString RdatHeaderRecord::amount;
 * @brief Amount of bases given in the header record.
 *
 * @var CharString RdatHeaderRecord::name;
 * @brief Name of Rdat sequence given in the header record.
 */

/*!
 * @fn RdatHeaderRecord::RdatHeaderRecord
 * @brief Constructor
 *
 * @signature RdatHeaderRecord::RdatHeaderRecord();
 * @signature RdatHeaderRecord::RdatHeaderRecord(amount, name);
 *
 * @param[in] amount   Amount of sequences of the header record, @link int32_t @endlink.
 * @param[in] name     Name of sequence of the header record, @link CharString @endlink.
 */

/*!
 * @fn RdatHeaderRecord#clear
 *
 * @brief Clear a RdatHeaderRecord.
 * @signature void clear(record);
 *
 * @param[in,out] record The RdatHeaderRecord to clear.
 */

class RdatHeader
{
public:
    // RDAT version
    CharString version;
    // Sequence name
    CharString name;

    // Default constructor.
    RdatHeader()
    {}

    // Construct directly with amount/name.
    RdatHeader(CharString const & version, CharString const & name) :
            version(version), name(name)
    {}
};

// ============================================================================
// Functions
// ============================================================================

inline void clear(RdatHeader & record)
{
    clear(record.version);
    clear(record.name);
}


} //namespace seqan
#endif //SEQAN_RDAT_FORMAT_HEADER_H