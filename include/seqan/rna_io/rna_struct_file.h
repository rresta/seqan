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
// Class for reading/writing RNA structure files.
// ==========================================================================

#ifndef SEQAN_RNA_IO_RNA_STRUCT_FILE_H_
#define SEQAN_RNA_IO_RNA_STRUCT_FILE_H_

namespace seqan {

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// --------------------------------------------------------------------------
// Tag Connect
// --------------------------------------------------------------------------

struct RnaStruct_;
typedef Tag<RnaStruct_> RnaStruct;

// --------------------------------------------------------------------------
// Class RnaStructContents
// --------------------------------------------------------------------------

struct RnaStructContents_
{
    std::vector<RnaRecord> records;
    RnaHeader header;
};
typedef struct RnaStructContents_ RnaStructContents;

// ============================================================================
// Typedefs
// ============================================================================

/*!
 * @class SeqFileIn
 * @signature typedef FormattedFile<Fastq, Input> SeqFileIn;
 * @extends FormattedFileIn
 * @headerfile <seqan/seq_io.h>
 * @brief Class for reading RAW, FASTA, FASTQ, EMBL and GENBANK files containing unaligned sequences.
 */

typedef FormattedFile<RnaStruct, Input> RnaStructFileIn;

/*!
 * @class SeqFileOut
 * @signature typedef FormattedFile<Fastq, Output> SeqFileOut;
 * @extends FormattedFileOut
 * @headerfile <seqan/seq_io.h>
 * @brief Class for writing RAW, FASTA, FASTQ, EMBL and GENBANK files containing unaligned sequences.
 */

typedef FormattedFile<RnaStruct, Output> RnaStructFileOut;

// --------------------------------------------------------------------------
// Tag AutoSeqFormat
// --------------------------------------------------------------------------
// if TagSelector is set to -1, the file format is auto-detected

/*!
 * @class AutoSeqFormat
 * @extends TagSelector
 * @headerfile <seqan/file.h>
 * @brief Auto-detects and stores a file format.
 *
 * @signature typedef TagList<Fastq, TagList<Fasta, TagList<Raw> > > SeqFormats;
 * @signature typedef TagSelector<SeqFormat> AutoSeqFormat;
 */

typedef
TagList<Connect,
TagList<Stockholm,
        TagList<DotBracket,
        TagList<Vienna,
        TagList<Ebpseq,
        TagList<Bpseq
> > > > > >
RnaStructFormats;

typedef TagSelector<RnaStructFormats>   RnaStructFormat;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<RnaStruct, TDirection, TSpec>, TStorageSpec>
{
    typedef RnaIOContext Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TSpec>
struct FileFormat<FormattedFile<RnaStruct, Input, TSpec> >
{
    typedef TagSelector<RnaStructFormats> Type;
};

template <typename TSpec>
struct FileFormat<FormattedFile<RnaStruct, Output, TSpec> >
{
    typedef TagSelector<RnaStructFormats> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord(TagSelector)
// ----------------------------------------------------------------------------

template <typename TFwdIterator>
inline void
readRecord(RnaRecord &, RnaIOContext &, TFwdIterator &, TagSelector<> const &)
{
    SEQAN_FAIL("RnaStructFileIn: File format not specified.");
}

template <typename TFwdIterator, typename TTagList>
inline void
readRecord(RnaRecord & record, RnaIOContext & context, TFwdIterator & iter, TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        readRecord(record, context, iter, TFormat());
    else
        readRecord(record, context, iter, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

template <typename TSpec>
inline void
readRecord(RnaRecord & record, FormattedFile<RnaStruct, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function writeRecord(TagSelector)
// ----------------------------------------------------------------------------

template <typename TFwdIterator>
inline void
writeRecord(TFwdIterator &, RnaRecord const &, RnaIOContext &, TagSelector<> const &)
{
    SEQAN_FAIL("RnaStructFileOut: File format not specified.");
}

template <typename TFwdIterator, typename TTagList>
inline void
writeRecord(TFwdIterator & iter, RnaRecord const & record, RnaIOContext & context, TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        writeRecord(iter, record, context, TFormat());
    else
        writeRecord(iter, record, context, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

template <typename TSpec>
inline void
writeRecord(FormattedFile<RnaStruct, Output, TSpec> & file, RnaRecord const & record)
{
    writeRecord(file.iter, record, context(file), file.format);
}

// ----------------------------------------------------------------------------
// Function readHeader(TagSelector)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readHeader(RnaHeader & header, FormattedFile<RnaStruct, Input, TSpec> & file)
{
    if (isEqual(file.format, Ebpseq()))
        readHeader(header, context(file), file.iter, Ebpseq());
    // other files contain no header
}

// ----------------------------------------------------------------------------
// Function writeHeader(TagSelector)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeHeader(FormattedFile<RnaStruct, Output, TSpec> & file, RnaHeader const & header)
{
    if (isEqual(file.format, Ebpseq()))
        writeHeader(file.iter, header, context(file), Ebpseq());
    // other files contain no header
}

// ----------------------------------------------------------------------------
// Function readRecords
// ----------------------------------------------------------------------------

template <typename TSpec, typename TSize>
inline void readRecords(RnaStructContents & contents, FormattedFile<RnaStruct, Input, TSpec> & file, TSize maxRecords)
{
    readHeader(contents.header, file);
    RnaRecord record;
    TSize numRecords = 0;
    while (!atEnd(file.iter) && ++numRecords <= maxRecords)
    {
        readRecord(record, file);
        append(contents.records, record);
    }
    // for files without header: create pseudo header
    if (empty(contents.header.seqLabels))
        createPseudoHeader(contents.header, contents.records); // in ebpseq_read_write.h
}

// ----------------------------------------------------------------------------
// Function writeRecords
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void writeRecords(FormattedFile<RnaStruct, Output, TSpec> & file, RnaStructContents const & contents)
{
    writeHeader(file, contents.header);
    for (RnaRecord const & record : contents.records)
        writeRecord(file, record);
}

}  // namespace seqan

#endif // SEQAN_RNA_IO_RNA_STRUCT_FILE_H_
