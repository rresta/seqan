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
// Author: Lily Shellhammer
// ==========================================================================
// Class for reading/writing files in Fasta or Fastq format.
// ==========================================================================

#ifndef SEQAN_RNA_IO_FILE_H_
#define SEQAN_RNA_IO_FILE_H_

namespace seqan {

// ============================================================================
// Classes, Tags
// ============================================================================

// --------------------------------------------------------------------------
// Tag RNA
// --------------------------------------------------------------------------

struct TagRNA_;
typedef Tag<TagRNA_> RNA;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>					//GOT THIS FROM VCF FILE THING, NOT SURE IF NECESSARY BUT SEEMS TO RELATE TO CONTEXT THING WHICH I KNOW IM MISSING
struct FormattedFileContext<FormattedFile<RNA, TDirection, TSpec>, TStorageSpec>
{
    typedef StringSet<CharString>                                  		TNameStore;
    typedef NameStoreCache<TNameStore>                              	TNameStoreCache;
    typedef RNAIOContext<TNameStore, TNameStoreCache, TStorageSpec>     Type;
};

// --------------------------------------------------------------------------
// Metafunction FileExtensions
// --------------------------------------------------------------------------

template <typename T>
struct FileExtensions<RNA, T>
{
    static char const * VALUE[1];
};
template <typename T>
char const * FileExtensions<Raw, T>::VALUE[1] =
{
    ".ct"      // default output extension
};


// ============================================================================
// Typedefs
// ============================================================================

// ----------------------------------------------------------------------------
// Typedef RNAFileIn
// ----------------------------------------------------------------------------

/*!
 * @class RNAFileIn
 * @signature typedef FormattedFile<RNA, Input> ConenctFileIn;
 * @extends FormattedFileIn
 * @headerfile <seqan/RNA_format_io.h>
 * @brief Class for reading RNA files.
 *
 * @see RNAHeader
 * @see RNARecord
 */
typedef FormattedFile<RNA, Input>   RNAFileIn;

// ----------------------------------------------------------------------------
// Typedef VcfFileOut
// ----------------------------------------------------------------------------

/*!
 * @class RNAFileOut
 * @signature typedef FormattedFile<RNA, Output> ConenctFileOut;
 * @extends FormattedFileOut
 * @headerfile <seqan/RNA_format_io.h>
 * @brief Class for writing RNA files.
 *
 * @see RNAHeader
 * @see RNARecord
 */

 typedef FormattedFile<RNA, Output>  RNAFileOut;

// ----------------------------------------------------------------------------
// Function readHeader(); RNAHeader
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readHeader(RNAHeader & header, FormattedFile<RNA, Input, TSpec> & file)
{
    readHeader(header, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function readRecord(); RNAHeader
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readRecord(RNARecord & record, FormattedFile<RNA, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function writeRecord(); RNARecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeRecord(FormattedFile<RNA, Output, TSpec> & file, RNARecord & record)
{
    writeRecord(file.iter, record, context(file), file.format);
}

// ----------------------------------------------------------------------------
// Function writeHeader(); RNAHeader
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeHeader(FormattedFile<RNA, Output, TSpec> & file, RNAHeader & header)
{
    writeHeader(file.iter, header, context(file), file.format);
}



} //seqan namespace

#endif	
//SEQAN_RNA_IO_FILE_H