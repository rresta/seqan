// ==========================================================================
//                                   rna_io
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
// Authors: Lily Shellhammer <lily.shellhammer@gmail.com>
//          Joerg Winkler <j.winkler@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_RNA_IO_TEST_RNA_IO_H_
#define TESTS_RNA_IO_TEST_RNA_IO_H_

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/sequence.h>
#include <seqan/rna_io.h>
#include <seqan/graph_types.h>
#include <seqan/align.h>

// ----------------------------------------------------------------------------
// Connect File I/O
// ----------------------------------------------------------------------------

// A test for connect file reading
SEQAN_DEFINE_TEST(test_rna_io_read_connect)
{
    // Path to example.ct
    seqan::CharString rnaPath = SEQAN_PATH_TO_ROOT();
    append(rnaPath, "/tests/rna_io/example.ct");

    seqan::String<char, seqan::MMap<> > mmapString;
    SEQAN_ASSERT(open(mmapString, toCString(rnaPath)));
    seqan::Iterator<seqan::String<char, seqan::MMap<> >, seqan::Rooted>::Type iter = begin(mmapString);

    seqan::RnaRecord rnaRecord;
    readRecord(rnaRecord, iter, seqan::Connect());

    SEQAN_ASSERT_EQ(rnaRecord.recordID, 0u);
    SEQAN_ASSERT_EQ(rnaRecord.seqLen, 73u);
    SEQAN_ASSERT_EQ(rnaRecord.offset, 1u);
    SEQAN_ASSERT_EQ(rnaRecord.energy, -17.50f);
    SEQAN_ASSERT_EQ(rnaRecord.name, "S.cerevisiae_tRNA-PHE");
    seqan::Rna5String base = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA";
    SEQAN_ASSERT_EQ(rnaRecord.sequence, base);
    seqan::TRnaAdjacencyIterator adj_it(rnaRecord.fixedGraphs[0].inter, 0);
    SEQAN_ASSERT_EQ(value(adj_it), 71u);
    SEQAN_ASSERT_EQ(rnaRecord.quality, "");
    // UNUSED: seqID, align, bppMatrGraphs, typeID, reactivity, reactError
}

SEQAN_DEFINE_TEST(test_rna_io_write_connect)
{
    seqan::RnaRecord record{};
    //set values
    record.energy = -17.5f;
    record.name = "S.cerevisiae_tRNA-PHE";
    record.sequence = "GCGGAUUU";
    record.seqLen = static_cast<unsigned>(length(record.sequence));
    seqan::TRnaRecordGraph graph;

    for (unsigned idx = 0; idx < record.seqLen; ++idx)
        addVertex(graph);
    for (unsigned idx = 0; idx < 4; ++idx)
        addEdge(graph, idx, 7u - idx, 1.);
    append(record.fixedGraphs, seqan::RnaInterGraph(graph));

    // Write records to string stream.String<char> out;
    seqan::String<char> outstr;
    writeRecord(outstr, record, seqan::Connect());

    // Compare string stream to expected value.
    seqan::String<char> expected = "8\tENERGY = -17.5\tS.cerevisiae_tRNA-PHE\n"
            " 1\tG\t0\t2\t8\t1\n"
            " 2\tC\t1\t3\t7\t2\n"
            " 3\tG\t2\t4\t6\t3\n"
            " 4\tG\t3\t5\t5\t4\n"
            " 5\tA\t4\t6\t4\t5\n"
            " 6\tU\t5\t7\t3\t6\n"
            " 7\tU\t6\t8\t2\t7\n"
            " 8\tU\t7\t9\t1\t8\n";
    SEQAN_ASSERT_EQ(outstr, expected);
}

// ----------------------------------------------------------------------------
// DotBracket File I/O
// ----------------------------------------------------------------------------

// A test for dot bracket rna file reading.
SEQAN_DEFINE_TEST(test_rna_io_read_dot_bracket)
{
    //Path to example.dt
    seqan::CharString rnaPath = SEQAN_PATH_TO_ROOT();
    append(rnaPath, "/tests/rna_io/example.dt");

    seqan::String<char, seqan::MMap<> > mmapString;
    SEQAN_ASSERT(open(mmapString, toCString(rnaPath)));
    seqan::Iterator<seqan::String<char, seqan::MMap<> >, seqan::Rooted>::Type iter = begin(mmapString);

    seqan::RnaRecord rnaRecord;
    readRecord(rnaRecord, iter, seqan::DotBracket());

    SEQAN_ASSERT_EQ(rnaRecord.recordID, 0u);
    SEQAN_ASSERT_EQ(rnaRecord.seqLen, 73u);
    SEQAN_ASSERT_EQ(rnaRecord.offset, 1u);
    SEQAN_ASSERT_EQ(rnaRecord.energy, -17.50f);
    SEQAN_ASSERT_EQ(rnaRecord.name, "S.cerevisiae_tRNA-PHE M10740");
    seqan::Rna5String base = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA";
    SEQAN_ASSERT_EQ(rnaRecord.sequence, base);
    seqan::TRnaAdjacencyIterator adj_it(rnaRecord.fixedGraphs[0].inter, 0);
    SEQAN_ASSERT_EQ(value(adj_it), 71u);
    SEQAN_ASSERT_EQ(rnaRecord.quality, "");
    // UNUSED: seqID, align, bppMatrGraphs, typeID, reactivity, reactError
}

SEQAN_DEFINE_TEST(test_rna_io_write_dot_bracket)
{
    seqan::RnaRecord record{};
    //set values
    record.energy = -17.5f;
    record.name = "S.cerevisiae_tRNA-PHE";
    record.sequence = "GCGGAUUU";
    record.seqLen = static_cast<unsigned>(length(record.sequence));
    seqan::TRnaRecordGraph graph;

    for (unsigned idx = 0; idx < record.seqLen; ++idx)
        addVertex(graph);
    for (unsigned idx = 0; idx < 3; ++idx)
        addEdge(graph, idx, 7u - idx, 1.);
    append(record.fixedGraphs, seqan::RnaInterGraph(graph));

    // Write records to string stream.String<char> out;
    seqan::CharString outstr;
    writeRecord(outstr, record, seqan::DotBracket());

    // Compare string stream to expected value.
    seqan::String<char> expected = ">S.cerevisiae_tRNA-PHE/1-8\n"
            "GCGGAUUU\n"
            "(((..))) (-17.5)\n";
    SEQAN_ASSERT_EQ(outstr, expected);
}

// ----------------------------------------------------------------------------
// Stockholm File I/O
// ----------------------------------------------------------------------------

// A test for Stockholm rna file reading.
SEQAN_DEFINE_TEST(test_rna_io_read_stockholm)
{
    //Path to example.sth
    seqan::CharString rnaPath = SEQAN_PATH_TO_ROOT();
    append(rnaPath, "/tests/rna_io/example.sth");

    seqan::String<char, seqan::MMap<> > mmapString;
    SEQAN_ASSERT(open(mmapString, toCString(rnaPath)));
    seqan::Iterator<seqan::String<char, seqan::MMap<> >, seqan::Rooted>::Type iter = begin(mmapString);

    seqan::RnaRecord rnaRecord;
    readRecord(rnaRecord, iter, seqan::Stockholm());

    SEQAN_ASSERT_EQ(rnaRecord.recordID, 0u);
    SEQAN_ASSERT_EQ(rnaRecord.seqLen, 74u);
    SEQAN_ASSERT_EQ(rnaRecord.offset, 1u);
    SEQAN_ASSERT_EQ(rnaRecord.energy, 0.0f);
    SEQAN_ASSERT_EQ(rnaRecord.name, "trna");
    seqan::Rna5String base = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA";
    SEQAN_ASSERT_EQ(stringSet(rnaRecord.align)[0], base);
    base = "UCCGUGAUAGUUUAAUGGUCAGAAUGGGCGCUUGUCGCGUGCCAGAUCGGGGUUCAAUUCCCCGUCGCGGAG";
    SEQAN_ASSERT_EQ(stringSet(rnaRecord.align)[2], base);
    SEQAN_ASSERT_EQ(rnaRecord.seqID[0], "DF6280");
    SEQAN_ASSERT_EQ(rnaRecord.seqID[2], "DD6280");

    seqan::TRnaAdjacencyIterator adj_it(rnaRecord.fixedGraphs[0].inter, 10);
    SEQAN_ASSERT_EQ(value(adj_it), 24u);
    // UNUSED: sequence, bppMatrGraphs, quality, typeID, reactivity, reactError
}

SEQAN_DEFINE_TEST(test_rna_io_write_stockholm)
{
    seqan::RnaRecord record{};
    // set values
    record.seqLen = 9u;
    record.name = "trna";
    record.comment = "alignment of 1415 tRNAs";
    seqan::Rna5String seq1 = "GCGGAUUU";
    seqan::Rna5String seq2 = "UCCGAUAUA";
    // create alignment
    seqan::resize(seqan::rows(record.align), 2);
    seqan::assignSource(seqan::row(record.align, 0), seq1);
    seqan::assignSource(seqan::row(record.align, 1), seq2);
    seqan::insertGap(seqan::row(record.align, 0), 4);
    // set sequence identifiers
    seqan::appendValue(record.seqID, seqan::CharString{"seq0"});
    seqan::appendValue(record.seqID, seqan::CharString{"seq1"});

    seqan::TRnaRecordGraph graph;
    for (unsigned idx = 0; idx < record.seqLen; ++idx)
        addVertex(graph);
    for (unsigned idx = 0; idx < 4; ++idx)
        addEdge(graph, idx, 7u - idx, 1.);
    append(record.fixedGraphs, seqan::RnaInterGraph(graph));

    // Write records to string stream.String<char> out;
    seqan::String<char> outstr;
    writeRecord(outstr, record, seqan::Stockholm());

    // Compare string stream to expected value.
    seqan::String<char> expected = "# STOCKHOLM 1.0\n"
            "#=GF ID      trna\n"
            "#=GF DE      alignment of 1415 tRNAs\n\n"
            "seq0        \tGCGG-AUUU\n"
            "seq1        \tUCCGAUAUA\n"
            "#=GC SS_cons\t(((()))).\n"
            "//\n";
    SEQAN_ASSERT_EQ(outstr, expected);
}

// ----------------------------------------------------------------------------
// Bpseq File I/O
// ----------------------------------------------------------------------------

// A test for connect file reading
SEQAN_DEFINE_TEST(test_rna_io_read_bpseq)
{
    // Path to example.bpseq
    seqan::CharString rnaPath = SEQAN_PATH_TO_ROOT();
    append(rnaPath, "/tests/rna_io/example.bpseq");

    seqan::String<char, seqan::MMap<> > mmapString;
    SEQAN_ASSERT(open(mmapString, toCString(rnaPath)));
    seqan::Iterator<seqan::String<char, seqan::MMap<> >, seqan::Rooted>::Type iter = begin(mmapString);

    seqan::RnaRecord rnaRecord;
    readRecord(rnaRecord, iter, seqan::Bpseq());

    SEQAN_ASSERT_EQ(rnaRecord.recordID, 0u);
    SEQAN_ASSERT_EQ(rnaRecord.seqLen, 50u);
    SEQAN_ASSERT_EQ(rnaRecord.offset, 1u);
    SEQAN_ASSERT_EQ(rnaRecord.energy, 0.0f);
    seqan::Rna5String base = "GGGCCGGGCGCGGUGGCGCGCGCCUGUAGUCCCAGCUACUCGGGAGGCUC";
    SEQAN_ASSERT_EQ(rnaRecord.sequence, base);
    seqan::TRnaAdjacencyIterator adj_it(rnaRecord.fixedGraphs[0].inter, 1);
    SEQAN_ASSERT_EQ(value(adj_it), 48u);
    SEQAN_ASSERT_EQ(rnaRecord.quality, "");
    seqan::CharString comment = " A header line beginning with # is for comments not for actual structure information. "
            "PDB ID 1E8O Signal Recognition Particle (SRP) RNA ";
    SEQAN_ASSERT_EQ(rnaRecord.comment, comment);
    // UNUSED: name, seqID, align, bppMatrGraphs, typeID, reactivity, reactError
}

SEQAN_DEFINE_TEST(test_rna_io_write_bpseq)
{
    seqan::RnaRecord record{};
    //set values
    record.name = "S.cerevisiae_tRNA-PHE";
    record.sequence = "GCGGAUUU";
    record.seqLen = static_cast<unsigned>(length(record.sequence));
    seqan::TRnaRecordGraph graph;

    for (unsigned idx = 0; idx < record.seqLen; ++idx)
        addVertex(graph);
    for (unsigned idx = 0; idx < 4; ++idx)
        addEdge(graph, idx, 7u - idx, 1.);
    append(record.fixedGraphs, seqan::RnaInterGraph(graph));

    // Write records to string stream.String<char> out;
    seqan::String<char> outstr;
    seqan::writeRecord(outstr, record, seqan::Bpseq());

    // Compare string stream to expected value.
    seqan::String<char> expected = "# S.cerevisiae_tRNA-PHE\n"
            "1\tG\t8\n2\tC\t7\n3\tG\t6\n4\tG\t5\n5\tA\t4\n6\tU\t3\n7\tU\t2\n8\tU\t1\n";
    SEQAN_ASSERT_EQ(outstr, expected);
}

#endif  // TESTS_RNA_IO_TEST_RNA_IO_H_
