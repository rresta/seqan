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

#ifndef SEQAN_INCLUDE_SEQAN_BPSEQ_READ_BPSEQ_H_
#define SEQAN_INCLUDE_SEQAN_BPSEQ_READ_BPSEQ_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Bpseq
// ----------------------------------------------------------------------------

/*!
 * @tag FileFormats#Bpseq
 * @headerfile <seqan/bpseq_io.h>
 * @brief Variant callinf format file.
 *
 * @signature typedef Tag<Bpseq_> Bpseq;
 */
struct Bpseq_;
typedef Tag<Bpseq_> Bpseq;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()                                            [BpseqHeader]
// ----------------------------------------------------------------------------

//inline void
//_parseBpseqContig(CharString & chromName, CharString const & headerValue)
//{
//    if (length(headerValue) < 3u)
//        return;
//
//    CharString tmp = infix(headerValue, 1, length(headerValue) - 2);
//    StringSet<CharString> tmp2;
//    strSplit(tmp2, tmp, EqualsChar<','>());
//    for (unsigned i = 0; i < length(tmp2); ++i)
//    {
//        if (!startsWith(tmp2[i], "ID="))
//            continue;
//        chromName = suffix(tmp2[i], 3);
//    }
//}

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readHeader(BpseqHeader & header,
           BpseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Bpseq const & /*tag*/)
{
    clear(header);
    CharString buffer;
    BpseqHeaderRecord record;

    while (!atEnd(iter) && value(iter) == '#') // All the information stored in the # lines are saved in a single line
    {
        skipOne(iter);
        clear(buffer);
        // Write header key.
        record.key = "Info";
        // Read header value.
        readLine(record.value, iter);
/*
        if (value(iter) == '#')
        {
            // Is header line.
            skipOne(iter);
            clear(record);

            // Read header key.
            readUntil(record.key, iter, OrFunctor<EqualsChar<'='>, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bpseq> >());

            // Skip '='.
            skipOne(iter);

            // Read header value.
            readLine(record.value, iter);
            appendValue(header, record);

            // Parse out name if headerRecord is a contig field.
            if (record.key == "contig")
            {
                _parseBpseqContig(buffer, record.value);
                appendName(contigNamesCache(context), buffer);
            }
        }
        else
        {
            // Is line "#CHROM\t...".
            readLine(buffer, iter);
            if (!startsWith(buffer, "CHROM"))
                ParseError("Invalid line with samples.");

            // Split line, get sample names.
            StringSet<CharString> fields;
            strSplit(fields, buffer, IsTab());
            if (length(fields) < 9u)
                ParseError("Not enough fields.");

            // Get sample names.
            for (unsigned i = 9; i < length(fields); ++i)
                appendName(sampleNamesCache(context), fields[i]);
        }
*/
    }
    appendValue(header, record);
}

// ----------------------------------------------------------------------------
// Function readRecord()                                            [BpseqRecord]
// ----------------------------------------------------------------------------
// Read record, updating list of known sequences if new one occurs.

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readRecord(BpseqRecord & record,
           BpseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Bpseq const & /*tag*/)
{
    typedef OrFunctor<IsSpace, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bpseq> > NextEntry;
    clear(record);
    CharString &buffer = context.buffer;
    CharString tmpStr="";
    while (!atEnd(iter) && value(iter) != '#')
    {
        // SEQPOS
        clear(buffer);
        readUntil(buffer, iter, NextEntry());
        if (empty(buffer))
            SEQAN_THROW(EmptyFieldError("SEQPOS"));
//        record.seqPos = nameToId(contigNamesCache(context), buffer);
//        std::cout << lexicalCast<__int32>(buffer) << std::endl;//<< "buffer size = " << length(buffer);
//        ciao = lexicalCast<__int32>(buffer);
        appendValue(record.seqPos, lexicalCast<__int32>(buffer));
//        appendValue(record.seqPos, value(iter));
        skipOne(iter);
        // SEQ
        clear(buffer);
        readUntil(buffer, iter, NextEntry());
        appendValue(record.seq,buffer[0]);
//        appendValue(record.seq, value(iter));
        if (empty(record.seq))
            SEQAN_THROW(EmptyFieldError("SEQ"));
//        std::cout << "value(iter)" << value(iter) << std::endl;
        skipOne(iter);
//        std::cout << "value(iter)" << value(iter) << std::endl;
        // INTERPOS
        clear(buffer);
        readUntil(buffer, iter, OrFunctor<IsSpace, IsNewline>());
        if (empty(buffer))
            SEQAN_THROW(EmptyFieldError("INTERPOS"));
//        record.interPos = nameToId(contigNamesCache(context), buffer);
        appendValue(record.interPos, lexicalCast<__int32>(buffer));
//        skipOne(iter);
//        // the following columns are optional
//        if (atEnd(iter) || IsNewline()(value(iter)))
//        {
//            skipLine(iter);
//        }
        // skip line break and optional additional columns
        skipLine(iter);
    }
//    std::cout << "value(iter)" << value(iter) << std::endl;
/*
    // CHROM
    clear(buffer);
    readUntil(buffer, iter, NextEntry());
    if (empty(buffer))
        SEQAN_THROW(EmptyFieldError("CHROM"));
    record.rID = nameToId(contigNamesCache(context), buffer);
    skipOne(iter);

    // POS
    clear(buffer);
    readUntil(buffer, iter, NextEntry());
    if (empty(buffer))
        SEQAN_THROW(EmptyFieldError("POS"));
    record.beginPos = lexicalCast<__int32>(buffer) - 1; // Translate from 1-based to 0-based.
    skipOne(iter);

    // ID
    readUntil(record.id, iter, NextEntry());
    if (empty(record.id))
        SEQAN_THROW(EmptyFieldError("ID"));
    skipOne(iter);

    // REF
    readUntil(record.ref, iter, NextEntry());
    if (empty(record.id))
        SEQAN_THROW(EmptyFieldError("REF"));
    skipOne(iter);

    // ALT
    readUntil(record.alt, iter, NextEntry());
    if (empty(record.id))
        SEQAN_THROW(EmptyFieldError("ALT"));
    skipOne(iter);

    // QUAL
    clear(buffer);
    readUntil(buffer, iter, NextEntry());
    if (empty(buffer))
        SEQAN_THROW(EmptyFieldError("QUAL"));

    if (buffer == ".")
        record.qual = BpseqRecord::MISSING_QUAL();
    else
        lexicalCastWithException(record.qual, buffer);

    skipOne(iter);

    // FILTER
    readUntil(record.filter, iter, NextEntry());
    if (empty(record.filter))
        SEQAN_THROW(EmptyFieldError("FILTER"));
    skipOne(iter);

    // INFO
    readUntil(record.info, iter, OrFunctor<IsTab, IsNewline>());
    if (empty(record.info))
        SEQAN_THROW(EmptyFieldError("INFO"));

    // the following columns are optional
    if (atEnd(iter) || IsNewline()(value(iter)))
    {
        skipLine(iter);
        return;
    }
    skipOne(iter);

    // FORMAT
    readUntil(record.format, iter, NextEntry());
    if (empty(record.format))
        SEQAN_THROW(EmptyFieldError("FORMAT"));
    skipOne(iter);

    // The samples.
    unsigned numSamples = length(sampleNames(context));
    for (unsigned i = 0; i < numSamples; ++i)
    {
        clear(buffer);
        if (i + 1 != numSamples)
        {
            readUntil(buffer, iter, NextEntry());
            skipOne(iter);
        }
        else
        {
            readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());
        }

        if (empty(buffer))
        {
            char buffer[30];    // == 9 (GENOTYPE_) + 20 (#digits in MIN_INT64) + 1 (trailing zero)
            sprintf(buffer, "GENOTYPE_%u", i + 1);
            SEQAN_THROW(EmptyFieldError(buffer));
        }
        appendValue(record.genotypeInfos, buffer);
    }

    // skip line break and optional additional columns
    skipLine(iter);
*/
    return;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BPSEQ_READ_BPSEQ_H_
