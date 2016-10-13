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
// Class for Vienna formatted files
// ==========================================================================

#ifndef SEQAN_RNA_IO_VIENNA_FILE_H_
#define SEQAN_RNA_IO_VIENNA_FILE_H_


namespace seqan {

// ============================================================================
// Classes, Tags
// ============================================================================
// ============================================================================
// Typedefs
// ============================================================================

// ----------------------------------------------------------------------------
// Typedef ViennaFileIn
// ----------------------------------------------------------------------------

/*!
 * @class ViennaFileIn
 * @signature typedef FormattedFile<Vienna, Input> ViennaFileIn;
 * @extends FormattedFileIn
 * @recordfile <seqan/rna_io.h>
 * @brief Class for reading Vienna files.
 *
 * @see RnaRecord
 */
typedef FormattedFile<Vienna, Input>   ViennaFileIn;

// ----------------------------------------------------------------------------
// Typedef ViennaFileOut
// ----------------------------------------------------------------------------

/*!
 * @class ViennaFileOut
 * @signature typedef FormattedFile<Vienna, Output> ViennaFileOut;
 * @extends FormattedFileOut
 * @recordfile <seqan/rna_io.h>
 * @brief Class for writing Vienna files.
 *
 * @see RnaRecord
 */
typedef FormattedFile<Vienna, Output>  ViennaFileOut;

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<Vienna, TDirection, TSpec>, TStorageSpec>
{
    typedef RnaIOContext Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<FormattedFile<Vienna, TDirection, TSpec> >
{
    typedef Vienna Type;
};

// ----------------------------------------------------------------------------
// Function writeRecord(); RnaRecord, Vienna File
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeRecord(FormattedFile<Vienna, Output, TSpec> & file, RnaRecord & record)
{
    writeRecord(file.iter, record, file.format);
}

// ----------------------------------------------------------------------------
// Function readRecord(); RnaRecord, Vienna File
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readRecord(RnaRecord & record, FormattedFile<Vienna, Input, TSpec> & file)
{
    readRecord(record, file.iter, file.format);
}



} //seqan namespace

#endif	//SEQAN_RNA_IO_VIENNA_FILE_H_
