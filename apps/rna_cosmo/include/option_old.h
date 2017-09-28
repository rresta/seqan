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
// This file contains the option structures and functions of rna cosmo
// application.
// ==========================================================================

#ifndef INCLUDE_OPTION_H_
#define INCLUDE_OPTION_H_

// ----------------------------------------------------------------------------
// std headers
// ----------------------------------------------------------------------------

#include <iostream>
#include <string.h>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/arg_parse.h>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

// defines all the constants used in the app
#include "data_types.h"

//using namespace std;

// ============================================================================
// Functors
// ============================================================================

typedef EqualsChar<'.'>        IsDot;
typedef EqualsChar<'/'>        IsSlash;
typedef EqualsChar<'\\'>       IsBackSlash;

#ifdef PLATFORM_WINDOWS
typedef IsBackSlash        IsPathDelimited;
#else
typedef IsSlash            IsPathDelimited;
#endif

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------
//TODO this can be re-written as a set of three maps with keys and values but
//if in the future should be embedded in the parsing function is better to
//conserve this style
struct Options
{
// use base pairs or structure
    bool useBasePairs;
// Name of input file
    seqan::CharString inFile;
// Name of input fileRef
    seqan::CharString inFileRef;
// Name of output file (default: stdout)
    seqan::CharString outFile;
// temporary directory where to save intermediate files. Default: use the input file directory.
    seqan::CharString tmpDir;
// Parameter used during the RNAfold execution to select the minimum energy to be considered
    double thrBppm;
// scoring mode, either LOGARITHMIC, SCALE, ORIGINAL, RIBOSUM
    unsigned structureScoring;
// define the weight of _half_ an interaction match for fixed structures
    double fixedStructWeight;
// if structureScoring=SCALING then we have to give a scaling factor
    double scalingFactor;
// time used for an hard timeout
    unsigned timeLimit;
// verbose(0) no outputs,
// verbose(1) Displays global statistics,
// verbose(2) Displays extensive statistics for each batch of reads,
// verbose(3) Debug output.
    unsigned verbose;
// number of threads forced
    unsigned threads;
// number of threads detected
    unsigned threadsCount;
    Options() :
            useBasePairs(true),
            thrBppm(1e-15), // 0.1
            structureScoring(RIBOSUM),
            fixedStructWeight(8.0),
            scalingFactor(1.0),
            timeLimit(-1),
            verbose(0)
    {
#ifdef _OPENMP
        threadsCount = std::thread::hardware_concurrency(); // omp_get_num_threads()
        threads = threadsCount;
#else
        threadsCount = 1;
#endif
    }
};

// ============================================================================
// Functions
// ============================================================================

using namespace seqan;
// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------
template <typename TOption>
void setupArgumentParser(ArgumentParser & parser, TOption const & /* options */)
{
    setAppName(parser, "RNA CoSMo");
    setShortDescription(parser, "RNA Consensus Structure Module");
    setCategory(parser, "RNA Structure Generation");
    setVersion(parser, "1.0");
    setDate(parser, "2017");
    //setDateAndVersion(parser);
    //setDescription(parser);
    addUsageLine(parser, "./rna_cosmo <\\fI-i inFile\\fP> \
            [\\fI-w outFile\\fP] [\\fI -parameters\\fP]");
    addOption(parser, ArgParseOption("v", "verbose", "verbose(0) no outputs, verbose(1) Displays global statistics, "
            "verbose(2) Displays extensive statistics for each batch of reads, verbose(3) Debug output.",
            ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("s", "useBasePairs", "Use structure prediction or fixed structure from extended "
            "input file."));
    addOption(parser, ArgParseOption("tb", "thrBppm", "(Parameter used during the RNAfold execution to select the "
            "minimum energy to be considered (default: 1e-15) ", ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("stsc", "structureScoring", "scoring mode, either LOGARITHMIC, SCALE, ORIGINAL, RIBOSUM",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("fsw", "fixedStructWeight", "define the weight of _half_ an interaction match for "
            "fixed structures", ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("scal", "scalingFactor", "if structurescoring=SCALING then we have to give a "
            "scaling factor", ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("tl", "timeLimit", "some additional option.", ArgParseArgument::INTEGER, "INT"));


    addSection(parser, "Input Options");
    addOption(parser, ArgParseOption("i", "inFile", "Path to the input file", ArgParseArgument::INPUT_FILE, "IN"));
    addOption(parser, ArgParseOption("ir", "inFileRef", "Path to the reference input file",
                                     ArgParseArgument::INPUT_FILE, "IN"));

    addSection(parser, "Output Options");
    addOption(parser, seqan::ArgParseOption( "w", "outFile", "Path to the output file (default: stdout)",
                                             ArgParseArgument::OUTPUT_FILE, "OUT"));
    addOption(parser, ArgParseOption("td", "tmpDir", "Specify a temporary directory where to save intermediate files. \
            Default: use the input file directory.", ArgParseOption::STRING));

    // Setup performance options.
    addSection(parser, "Performance Options");
    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use.", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
#ifdef _OPENMP
    setMaxValue(parser, "threads", std::to_string(std::thread::hardware_concurrency() + 1));
#else
    setMaxValue(parser, "threads", "1");
#endif
//    setDefaultValue(parser, "threads", options.threadsCount);
}

// ----------------------------------------------------------------------------
// Function setEnv()
// ----------------------------------------------------------------------------

template <typename TString, typename TValue>
bool setEnv(TString const & key, TValue const & value)
{
#ifdef PLATFORM_WINDOWS
    return !_putenv_s(toCString(key), toCString(value));
#else
    return !setenv(toCString(key), toCString(value), true);
#endif
}

// ----------------------------------------------------------------------------
// Function getCwd()
// ----------------------------------------------------------------------------

template <typename TString>
void getCwd(TString & string)
{
    char cwd[1000];

#ifdef PLATFORM_WINDOWS
    _getcwd(cwd, 1000);
#else
    ignoreUnusedVariableWarning(getcwd(cwd, 1000));
#endif

    assign(string, cwd);
}

// ----------------------------------------------------------------------------
// Function firstOf()
// ----------------------------------------------------------------------------

template <typename TString, typename TFunctor>
inline typename Iterator<TString const, Standard>::Type
firstOf(TString const & string, TFunctor const & f)
{
    typedef typename Iterator<TString const, Standard>::Type TIter;

    TIter it = begin(string, Standard());
    skipUntil(it, f);

    return it;
}

// ----------------------------------------------------------------------------
// Function lastOf()
// ----------------------------------------------------------------------------

template <typename TString, typename TFunctor>
inline typename Iterator<TString const, Standard>::Type
lastOf(TString const & string, TFunctor const & f)
{
    typedef ModifiedString<TString const, ModReverse>        TStringRev;
    typedef typename Iterator<TStringRev, Standard>::Type    TIterRev;

    TStringRev revString(string);
    TIterRev revIt = firstOf(revString, f);

    return end(string) - position(revIt, revString);
}

// ----------------------------------------------------------------------------
// Function getPath()
// ----------------------------------------------------------------------------

template <typename TString>
inline typename Prefix<TString const>::Type
getPath(TString const & string)
{
    typedef typename Iterator<TString const, Standard>::Type TIter;

    TIter it = lastOf(string, IsPathDelimited());

    if (it != begin(string, Standard())) --it;

    return prefix(string, it);
}

// ----------------------------------------------------------------------------
// Function parseCmd()
// ----------------------------------------------------------------------------
template <typename TOption>
ArgumentParser::ParseResult parse(TOption & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;
    getOptionValue(options.verbose, parser, "verbose");
    getOptionValue(options.useBasePairs, parser, "useBasePairs");
    getOptionValue(options.thrBppm, parser, "thrBppm");
    getOptionValue(options.structureScoring, parser, "structureScoring");
    getOptionValue(options.fixedStructWeight, parser, "fixedStructWeight");
    getOptionValue(options.scalingFactor, parser, "scalingFactor");
    getOptionValue(options.timeLimit, parser, "timeLimit");
    getOptionValue(options.threads, parser, "threads");

    getOptionValue(options.inFile, parser, "inFile");
    if (empty(options.inFile))
        return ArgumentParser::PARSE_ERROR;

    getOptionValue(options.inFileRef, parser, "inFileRef");
    getOptionValue(options.outFile, parser, "outFile");
    if (isSet(parser, "outFile"))
    {
        _V(options, "The specified output file is " << options.outFile);
    }
    else
    {
        CharString tmpDir;
        getOptionValue(tmpDir, parser, "tmpDir");
        if (!isSet(parser, "tmpDir"))
        {
            tmpDir = SEQAN_TEMP_FILENAME();
            // remove "/test_file" suffix
            erase(tmpDir, length(tmpDir) - 10u, length(tmpDir));
        }
        setEnv("TMPDIR", tmpDir);
        options.tmpDir = tmpDir;
        _V(options, "The absolute path where to create the tmpDir is " << tmpDir);
    }
    _V(options, "Initialized the Options structure");
    return ArgumentParser::PARSE_OK;
}

#endif /* INCLUDE_OPTION_SRTUCT_FUNCT_H_ */
