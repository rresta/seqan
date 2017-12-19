// ==========================================================================
//               RNA CoSMo - RNA Consensus Structure Module
// ==========================================================================
// Copyright (c) 2015-2017, Gianvito Urgese
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
// Author: Rossella Resta <s222385@studenti.polito.it>
// ==========================================================================
// This file contains functions to read the input and to write the output
// of rna cosmo application.
// ==========================================================================

#ifndef _INCLUDE_LARA_IO_H_
#define _INCLUDE_LARA_IO_H_

// ============================================================================
// Prerequisites
// ============================================================================

// ----------------------------------------------------------------------------
// std headers
// ----------------------------------------------------------------------------

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _readRnaInputFile()
// ----------------------------------------------------------------------------

template <typename TOption>
void _readRnaInputFile(RnaStructContents & filecontents, CharString filename, TOption const & options)
{
    if (empty(filename))
        return;

    RnaStructFileIn rnaStructFile;
    if (open(rnaStructFile, toCString(filename), OPEN_RDONLY))
    {
        _VVV(options, "Input file is RnaStruct.");
        readRecords(filecontents, rnaStructFile, std::numeric_limits<unsigned>::max());
        close(rnaStructFile);
    }
    else
    {
        _VVV(options, "Input file is Fasta/Fastq.");
        SeqFileIn seqFileIn(toCString(filename));
        StringSet<CharString> ids;
        StringSet<IupacString> seqs;
        StringSet<CharString> quals;
        readRecords(ids, seqs, quals, seqFileIn);
        close(seqFileIn);
        resize(filecontents.records, length(ids));
        SEQAN_ASSERT_EQ(length(ids), length(seqs));
        for (typename Size<StringSet<CharString> >::Type idx = 0u; idx < length(ids); ++idx)
        {
            filecontents.records[idx].name = ids[idx];
            filecontents.records[idx].sequence = convert<Rna5String>(seqs[idx]);
        }
        if (length(quals) == length(ids))
        {
            for (typename Size<StringSet<CharString> >::Type idx = 0u; idx < length(ids); ++idx)
                filecontents.records[idx].quality = quals[idx];
        }
    }
}

// ----------------------------------------------------------------------------
// Function _readMultiStructRnaInputFile()
// ----------------------------------------------------------------------------

template <typename TOption>
void _readMultiStructRnaInputFile(RnaStructContents & contentsOut, CharString filename, TOption const & options)
{
    if (empty(filename))
        return;

    RnaStructFileIn rnaStructFile;
    RnaStructContents contentsIn;
    if (open(rnaStructFile, toCString(filename), OPEN_RDONLY))
    {
        readRecords(contentsIn, rnaStructFile, std::numeric_limits<unsigned>::max());
        close(rnaStructFile);
    }
    else
    {
        std::cerr << "Can't open the input file.";
    }

    contentsOut.records.push_back (contentsIn.records[0]); // add current record
    for(unsigned i = 1; i < length(contentsIn.records); ++i)
    {
        bool flag = 0;
        for(unsigned j = 0; j < length(contentsOut.records); ++j)
        {
            if (contentsIn.records[i].sequence == contentsOut.records[j].sequence) {
                append(contentsOut.records[j].fixedGraphs, contentsIn.records[i].fixedGraphs);
                append(contentsOut.records[j].bppMatrGraphs, contentsIn.records[i].bppMatrGraphs);
                flag = 1;
            }
        }
        if (flag == 0)
            contentsOut.records.push_back(contentsIn.records[i]);
    }

    _V(options, "Read " << length(contentsIn.records) << " records from input file.");
    _V(options, "Read " << length(contentsOut.records) << " records from output structure.");
    _VVV(options, contentsOut.header.description);
    for(unsigned i = 0; i < length(contentsOut.records); ++i)
    {
        _VV(options, contentsOut.records[i].name);
        _VV(options, contentsOut.records[i].sequence);
        _VV(options, "Number of Fixed Graphs for the record " << i << " : " << length(contentsOut.records[i].fixedGraphs));
        _VV(options, "Number of bpp Matrices for the record " << i << " : " << length(contentsOut.records[i].bppMatrGraphs));
        for(unsigned j = 0; j < length(contentsOut.records[i].fixedGraphs); ++j)
        {
            _VVV(options, "Parameters of the " << j+1 << " Fixed Graph (specs, energy): ")
            _VVV(options, contentsOut.records[i].fixedGraphs[j].specs);
            _VVV(options, contentsOut.records[i].fixedGraphs[j].energy);
        }
    }
}
//std::ifstream data("/home/osboxes/data_test/sequences/RMDBGLYCFN_SHP_0004.rdat");
// --------------------------------------------------------------------------
// Function rdatContents()
// --------------------------------------------------------------------------
template <typename TOption>
void rdatContents(Rna5String & RDATseq, String<RnaStructureGraph> & ifrGraph,
                  StringSet<String<float> > & stringSetReactivity, StringSet<String<float> > & stringSetReactivityError,
                  std::ifstream & rdatFile, TOption const & options)
{
    unsigned reactivityC = 0;         //number of reactivity values
    unsigned reactivityErrC = 0;      //number of reactivity error values
    bool isFirstString = true;
    bool isFirstStringAgain = true;
    unsigned numReactivities = 0;      //number of reactivity lines in rdat file
    unsigned numReactivitiesErr = 0;      //number of reactivity error in rdat file
    std::string currentLine;
    std::string sequenceStr = "SEQUENCE";
    std::string structureStr = "STRUCTURE";
    std::string reactivityStr = "REACTIVITY:";
    std::string reactivityErrStr = "REACTIVITY_ERROR:";
//    std::string reactivityStr = "DATA:1";

    while (std::getline(rdatFile, currentLine)) {
        std::istringstream iss(currentLine);
        if (std::strstr(currentLine.c_str(), sequenceStr.c_str())) {
            while (iss) {
                std::string lineContent;
                iss >> lineContent;
                if (!std::strstr(lineContent.c_str(), sequenceStr.c_str()) && !lineContent.empty()) {
                    RDATseq = lineContent;
                    _VV(options, "> Sequence\n" << RDATseq);
                }
            }
        }
        if (std::strstr(currentLine.c_str(), structureStr.c_str())) {
            while (iss) {
                std::string lineContent;
                iss >> lineContent;
                if (!std::strstr(lineContent.c_str(), structureStr.c_str()) && !lineContent.empty()) {
                    _VV(options, "> Secondary Structure\n" << lineContent);    //Dot-Bracket String
                    CharString struc = lineContent;
                    //from string to undirected graph
                    bracket2graph(ifrGraph, struc);
                    _VVV(options, "Adjacency list:\n" << "> Graph");
                    Iterator<Graph<Undirected<double> >, EdgeIterator>::Type bppEdgIt(ifrGraph[0].inter);
                    for (bppEdgIt; !atEnd(bppEdgIt); goNext(bppEdgIt)) {
                        _VVV(options,"sourceVertex = " << sourceVertex(bppEdgIt) << " targetVertex = "
                                                       << targetVertex(bppEdgIt));
                    }
                }
            }
        }
        if (std::strstr(currentLine.c_str(), reactivityStr.c_str())) {  //reactivities line
            String<float> reactivityForStruct;
            reactivityC = 0;
            while (iss) {
                std::string lineContent;
                iss >> lineContent;
                if (!std::strstr(lineContent.c_str(), reactivityStr.c_str()) && !lineContent.empty()) {
                    float reactToAppend;
                    reactToAppend=std::strtof((lineContent).c_str(),0);     //from string to float
                    appendValue(reactivityForStruct, reactToAppend);
                    ++reactivityC;
                }
            }
            ///////////////////////////////////////////////
            if (isFirstString) {
                appendValue(stringSetReactivity, reactivityForStruct);
                isFirstString = false;
            } else {
                resize(stringSetReactivity, numReactivities, Exact());
                appendValue(stringSetReactivity, reactivityForStruct);
            }
            ++numReactivities;
        }
        if (std::strstr(currentLine.c_str(), reactivityErrStr.c_str())) {   //reactivity errors line
            String<float> reactivityErrForStruct;
            reactivityErrC = 0;
            while (iss) {
                std::string lineContent;
                iss >> lineContent;
                if (!std::strstr(lineContent.c_str(), reactivityErrStr.c_str()) && !lineContent.empty()) {
                    float reactErrToAppend;
                    reactErrToAppend=std::strtof((lineContent).c_str(),0);     //from string to float
                    appendValue(reactivityErrForStruct, reactErrToAppend);
                    ++reactivityErrC;
                }
            }
            if (isFirstStringAgain) {
                appendValue(stringSetReactivityError,
                            reactivityErrForStruct); // Append string to the end of the string set.
                isFirstStringAgain = false;
            } else {
                resize(stringSetReactivityError, numReactivitiesErr, Exact());
                appendValue(stringSetReactivityError, reactivityErrForStruct);
            }
            ++numReactivitiesErr;
        }
    }

    _VV(options, "> Number of Reactivity : " << length(stringSetReactivity));
    _VV(options, "> Number of Reactivity Error : " << length(stringSetReactivityError));
    if (length(RDATseq) != reactivityC)           //check if the number of react is equal to seq length
    {
        _V(options, "RDAT file:\nNumber of reactivities: " << reactivityC << "\nSequence length: " << length(RDATseq));
        _V(options, "\nRMDB file not valid.");
    }
}


// ----------------------------------------------------------------------------
// Function plotOutput()
// ----------------------------------------------------------------------------

template <typename TOption>
void plotOutput(TOption const & options)
{
    _V(options, "print output using jviz or VARNA");
}
#endif //_INCLUDE_LARA_IO_H_
