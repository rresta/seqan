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

#ifndef SEQAN_INCLUDE_SEQAN_BPSEQ_WRITE_BPSEQ_H_
#define SEQAN_INCLUDE_SEQAN_BPSEQ_WRITE_BPSEQ_H_

namespace seqan {

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function writeHeader()                                           [BpseqHeader]
// ----------------------------------------------------------------------------

template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
writeHeader(TTarget & target,
            BpseqHeader const & header,
            BpseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
            Bpseq const & /*tag*/)
{
    for (unsigned i = 0; i < length(header); ++i)
    {
        write(target, "#");
        write(target, header[i].key);
        writeValue(target, '=');
        write(target, header[i].value);
        writeValue(target, '\n');
    }
//
//    write(target, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
//    for (unsigned i = 0; i < length(sampleNames(context)); ++i)
//    {
//        writeValue(target, '\t');
//        write(target, sampleNames(context)[i]);
//    }
//    writeValue(target, '\n');
}

// ----------------------------------------------------------------------------
// Function writeRecord()                                           [BpseqRecord]
// ----------------------------------------------------------------------------

template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
writeRecord(TTarget & target,
            BpseqRecord const & record,
            BpseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
            Bpseq const & /*tag*/)
{
//	// SEQPOS
//    for (unsigned i = 0; i < length(record.seqPos); ++i)
//    {
//        writeValue(target, '\t');
////        if (empty(record.genotypeInfos[i]))
////            writeValue(target, '.');
////        else
//            write(target, record.seqPos[i]);
//    }
////    write(target, record.seqPos);
//    writeValue(target, '\n');
//
//    // SEQ
//    write(target, record.seq);
//    writeValue(target, '\n');
//
//    // INTERPOS
//    for (unsigned i = 0; i < length(record.interPos); ++i)
//    {
//        writeValue(target, '\t');
////        if (empty(record.genotypeInfos[i]))
////            writeValue(target, '.');
////        else
//            write(target, record.interPos[i]);
//    }
////    write(target, record.interPos);
//    writeValue(target, '\n');

	// SEQPOS
	for (unsigned i = 0; i < length(record.seqPos); ++i)
	{
		write(target, record.seqPos[i]);
		writeValue(target, ' ');
		write(target, record.seq[i]);
		writeValue(target, ' ');
		write(target, record.interPos[i]);
		writeValue(target, '\n');
	}
/*
    // CHROM
    write(target, contigNames(context)[record.rID]);
    writeValue(target, '\t');

    // POS
    appendNumber(target, record.beginPos + 1);
    writeValue(target, '\t');

    // ID
    if (empty(record.id))
        writeValue(target, '.');
    else
        write(target, record.id);
    writeValue(target, '\t');

    // REF
    if (empty(record.ref))
        writeValue(target, '.');
    else
        write(target, record.ref);
    writeValue(target, '\t');

    // ALT
    if (empty(record.alt))
        writeValue(target, '.');
    else
        write(target, record.alt);
    writeValue(target, '\t');

    // QUAL
    if (record.qual != record.qual)  // only way to test for nan
        writeValue(target, '.');
    else
        appendNumber(target, record.qual);

    // FILTER
    writeValue(target, '\t');
    if (empty(record.filter))
        writeValue(target, '.');
    else
        write(target, record.filter);
    writeValue(target, '\t');

    // INFO
    if (empty(record.info))
        writeValue(target, '.');
    else
        write(target, record.info);

    // FORMAT
    writeValue(target, '\t');
    if (empty(record.format))
        writeValue(target, '.');
    else
        write(target, record.format);

    // The samples.
    for (unsigned i = 0; i < length(record.genotypeInfos); ++i)
    {
        writeValue(target, '\t');
        if (empty(record.genotypeInfos[i]))
            writeValue(target, '.');
        else
            write(target, record.genotypeInfos[i]);
    }
    writeValue(target, '\n');
*/
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BPSEQ_WRITE_BPSEQ_H_
