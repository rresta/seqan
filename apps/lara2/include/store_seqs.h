// ==========================================================================
//               LaRAgu - Lagrangian Relaxation Aligner GU
// ==========================================================================
// Copyright (c) 2015-2016, Gianvito Urgese
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
//     * Neither the name of Gianvito Urgese nor the names of its contributors
//       may be used to endorse or promote products derived from this software
//       without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL GIANVITO URGESE OR THE FU BERLIN BE LIABLE
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
// This file contains the SeqStore class.
// ==========================================================================

#ifndef _INCLUDE_STORE_SEQS_H_
#define _INCLUDE_STORE_SEQS_H_

// ----------------------------------------------------------------------------
// std headers
// ----------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/seq_io.h>
#include <seqan/rna_io.h>

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function CheckInstanceFormat();
// ----------------------------------------------------------------------------

template <typename TOption, typename TPath>
unsigned CheckInstanceFormat(TOption const & options, TPath const & inFilePath)
// TODO to be deprecated after all the requested file formats will be available in SeqAn
{
    _V(options, "Check the Instance Format and save the input file code")
    unsigned input_type = UNKNOWN;
    std::ifstream inFile(toCString(inFilePath));
    _V(options, "input file ifstream " << inFilePath);
    std::string cur_line;
    if (inFile.is_open())
    {
        std::getline(inFile, cur_line);
/*
        do
        {
            std::getline(inFile, cur_line);
        } while(cur_line[0] == '#' && !inFile.eof());
        if(cur_line.size() == 0)
            input_type = UNKNOWN;
*/
        if(cur_line[0] == '>' )
        {
            input_type = FASTA;
            do
            {
                std::getline(inFile, cur_line);
                if (cur_line.find("(") != std::string::npos && cur_line[0] != '>')
                {
                    input_type = RNASTRUCT;
                }
            } while(!inFile.eof());
// if we did not find structure, then return FASTA, otherwis we already returned EXTENDED_FASTA
        }
        else if(cur_line[0] == '@' )
        {
            input_type = FASTQ;
        }
        else if(cur_line[0] == '#')
        {
            input_type = RNASTRUCT;
        }
    }
    inFile.close();
    return input_type;
}

// ----------------------------------------------------------------------------
// Function readFastaRecords()
// ----------------------------------------------------------------------------

// This function is able to read Fasta, MultiFasta, FASTQ, EMBL or GenBank formats.
template <typename TOption, typename TSeqVect>
unsigned readFastaRecords(TSeqVect & rnaSeqs, TOption const & options, CharString const & inFilePath)
{
    SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(inFilePath), OPEN_RDONLY))
    {
        std::cerr << "ERROR: Could not open the file " << inFilePath << std::endl;
        return 1;
    }

    StringSet<CharString> ids;
    StringSet<Rna5String> seqs;
    StringSet<CharString> quals;
    try
    {
        readRecords(ids, seqs, quals, seqFileIn);
    }
    catch (Exception const &e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    seqan::resize(rnaSeqs, length(ids));
    if (length(quals) == length(ids))
    {
        for (unsigned i = 0; i < length(ids); ++i)
        {
            rnaSeqs[i].name = ids[i];
            rnaSeqs[i].sequence = seqs[i];
            rnaSeqs[i].quality = quals[i];
            _VV(options, rnaSeqs[i].name << "\t" << rnaSeqs[i].sequence << "\t" << rnaSeqs[i].quality);
        }
    }
    else
    {
        for (unsigned i = 0; i < length(ids); ++i)
        {
            rnaSeqs[i].name = ids[i];
            rnaSeqs[i].sequence = seqs[i];
            _VV(options, rnaSeqs[i].name << "\t" << rnaSeqs[i].sequence);
        }
    }
    return 0;
};

// ----------------------------------------------------------------------------
// Function readRnaStructRecords()
// ----------------------------------------------------------------------------
// This function is able to read RnaStruct formats.
template <typename TOption, typename TSeqVect>
unsigned readRnaStructRecords(TSeqVect & rnaSeqs, TOption const & options, CharString const & inFilePath)
{
    RnaStructFileIn rnaStructFile;
    if (!open(rnaStructFile, toCString(inFilePath), OPEN_RDONLY))
    {
        std::cerr << "ERROR: Could not open the rna struct file " << inFilePath << std::endl;
        return 1;
    }
    RnaStructContents contents;
    readRecords(contents, rnaStructFile, 10000);
    rnaSeqs = contents.records;
    _VV(options, "Read " << length(rnaSeqs) << " records from file.");
    return 0;
};


// ----------------------------------------------------------------------------
// Function readRnaRecords()
// ----------------------------------------------------------------------------
// This function is able to read Fasta, MultiFasta, FASTQ, EMBL, GenBank, dbn, edbn, bpseq and ebpseq formats.
template <typename TOption, typename TFile, typename TSeqVect>
unsigned readRnaRecords(TSeqVect & rnaSeqs, TOption const & options, TFile const & file)
{
    unsigned seqanSupportedFiles = UNKNOWN;
    const char * fileName = toCString(file);
    CharString inFilePath = getAbsolutePath(fileName);
    SeqFileIn seqFileIn;
    seqanSupportedFiles = CheckInstanceFormat(options, inFilePath);
    _V(options,"Reading sequences from file type " << seqanSupportedFiles << " named " << inFilePath );
// the file format should be recognized and a readRecord function with the file type flag should be
// called to acquire the inputs
    switch(seqanSupportedFiles)
    {
        case FASTA:
            _VV(options, "Input file is in Fasta format");
            readFastaRecords(rnaSeqs, options, inFilePath);
            break;
        case FASTQ:
            _VV(options, "Input file is in Fastq format");
            readFastaRecords(rnaSeqs, options, inFilePath);
            break;
        case RNASTRUCT:
            _VV(options, "Input file is RnaStruct");
            readRnaStructRecords(rnaSeqs, options, inFilePath);
            break;
        case UNKNOWN:
            _V(options, "Unable to identify the type of your input instance. Please check the allowed file formats");
            _V(options, "Please check, if you supplied a either valid file.");
            exit(1);
    }
    return 0;
}

#endif // #ifndef _INCLUDE_STORE_SEQS_H_
