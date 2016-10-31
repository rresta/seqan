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
// This file contains the option structures and functions of seqan_laragu
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

#include "top_data_struct.h" //defined all the constant used in the app
#include "misc_options.h"

//using namespace std;

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
// Use the gap scheme to be used in the N-W alignment (default: affine(0))
    unsigned affineLinearDgs;
// Use the global local or global-Unconstrained algorithm (default: global(0) - local(1) )
    bool globalLocal;
// type used for the global-unconstrained alignment AlignConfig <TTop, Tleft, TRight, TDown>
    bool unTop;
    bool unLeft;
    bool unRight;
    bool unDown;
// Threshold of ratio sizes for automatic choice global and global-Unconstrained algorithm (default: 2/3)
    double thrGlobalUnconstr;
// Threshold of ratio sizes for automatic choice global and Local algorithm (default: 1/2))
    double thrGlobalLocal;
// Parameter used during the RNAfold execution to select the minimum energy to be considered
    double thrBppm;
// number of iterations
    unsigned iterations;
// number of non-decreasing iterations
    unsigned nonDecreasingIterations;
// method to be used for the computation of the Lower bound (MWM or approximation can be chosen)
    unsigned lowerBoundMethod;
// value to be considered for the equality of upper and lower bounds difference
    double epsilon;
// my, necessary for computing appropriate step sizes
    double my; //FIXME to be changed the name
// scoring matrix name that should be used for scoring alignment edges in the actual problem
    seqan::CharString laraScoreMatrixName;
//    Score<double, ScoreMatrix<Rna5, Default> > laraScoreMatrix;
    TScoreMatrix laraScoreMatrix;
//    TScoringSchemeRib laraScoreMatrixRib;
// Gap open and extend costs for generating the alignment edges
    double generatorGapOpen;
    double generatorGapExtend;
    double generatorSuboptimality; // FIXME what means this parameter?
// Gap open and extend costs for generating the alignment edges
    double laraGapOpen;
    double laraGapExtend;
// scaling factor for the scores of the alignment edges
    double sequenceScale;
// gap penalty for RSA
    double rsaGapPenalty;
// scoring mode, either LOGARITHMIC,SCALING or RIBOSUM
    seqan::CharString structureScoring;
// define the weight of _half_ an interaction match for fixed structures
    double fixedStructWeight;
// if structureScoring=SCALING then we have to give a scaling factor
    double scalingFactor;
// specify the location of T-COFFEE
    seqan::CharString tcoffeeLocation;
// Define the id of the sequence that must be splitted
    unsigned splitSequence;
// window size specifies the length of the sliding window when the local alignment algorithm is used (long sequnce vs short)
    unsigned windowSize;
// time used for an hard timeout
    unsigned timeLimit;
// verbose(0) no outputs, verbose(1) Displays global statistics, verbose(2) Displays extensive statistics for each batch of reads.
    unsigned verbose;
// number of threads forced
    unsigned threads;
// number of threads detected
    unsigned threadsCount;
    Options() :
            useBasePairs(true),
            affineLinearDgs(0),
            globalLocal(false),
            unTop(false),
            unLeft(false),
            unRight(false),
            unDown(false),
            thrGlobalUnconstr(0.666667),
            thrGlobalLocal(0.5),
            thrBppm(1e-15), // 0.1 is the value used in the old Lara
            iterations(500),
            nonDecreasingIterations(50),
            lowerBoundMethod(LBMWMTEST),
            epsilon(EPSILON),
            my(1.0),
            laraScoreMatrixName(""), //laraScoreMatrixName("RIBOSUM65"),
            generatorGapOpen(6.0),
            generatorGapExtend(2.0),
            generatorSuboptimality(40),
            laraGapOpen(-12.0),
            laraGapExtend(-5.0),
            sequenceScale(1.0),
            rsaGapPenalty(3.0),
            structureScoring("LOGARITHMIC"),
            fixedStructWeight(8.0),
            scalingFactor(1.0),
            tcoffeeLocation("t_coffee/t_coffee_5.05"),
            splitSequence(1),
            windowSize(100),
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
    setAppName(parser, "LaRA-GU");
    setShortDescription(parser, "Lagrangian Relaxation Structural Alignment Algorithm");
    setCategory(parser, "Structural Alignment Algorithm");
    setVersion(parser, "1.0");
    setDate(parser, "2016");
    //setDateAndVersion(parser);
    //setDescription(parser);
    addUsageLine(parser, "./lara <\\fI-i inFile\\fP> \
            [\\fI-w outFile\\fP] [\\fI -parameters\\fP]");
    addOption(parser, ArgParseOption("v", "verbose", "verbose(0) no outputs, verbose(1) Displays global statistics, verbose(2) Displays extensive statistics for each batch of reads.",
                                     ArgParseArgument::INTEGER, "INT"));

    addSection(parser, "LaRA Alignment Options");
    addOption(parser, ArgParseOption("s", "useBasePairs", "Use structure prediction or fixed structure from extended input file."));
    addOption(parser, ArgParseOption("g", "affineLinearDgs", "Chose the gap scheme affine(0) linear(1) or dynamic(2) to be used in the alignment (default: affine(0))."));
    addOption(parser, ArgParseOption("a", "globalLocal", "Use the local or global algorithm (default: global(false)."));
    addOption(parser, ArgParseOption("ut", "unTop", "type used for the global-unconstrained alignment AlignConfig TTop (default: top(false). "));
    addOption(parser, ArgParseOption("ul", "unLeft", "type used for the global-unconstrained alignment AlignConfig TLeft (default: left(false). "));
    addOption(parser, ArgParseOption("ur", "unRight", "type used for the global-unconstrained alignment AlignConfig TRight (default: right(false). "));
    addOption(parser, ArgParseOption("ud", "unDown", "type used for the global-unconstrained alignment AlignConfig TDown (default: dow(false). "));
    addOption(parser, ArgParseOption("tgu", "thrGlobalUnconstr", "Threshold of ratio sizes for automatic choice global and global-Unconstrained algorithm (default: 2/3)"));
    addOption(parser, ArgParseOption("tgl", "thrGlobalLocal", "Threshold of ratio sizes for automatic choice global and Local algorithm (default: 1/2))"));
    addOption(parser, ArgParseOption("tb", "thrBppm", "(Parameter used during the RNAfold execution to select the minimum energy to be considered (default: 1e-15)"));
    addOption(parser, ArgParseOption("iter", "iterations", "number of iterations.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("nditer", "nonDecreasingIterations", "number of non-decreasing iterations.",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("lbm", "lowerBoundMethod", "method to be used for the computation of the Lower bound (MWM or approximation can be chosen)",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("ep", "epsilon", "value to be considered for the equality of upper and lower bounds difference",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("my", "my", "necessary for computing appropriate step sizes.",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("lsm","laraScoreMatrixName",
                                     "scoring matrix name that should be used for scoring alignment edges in the actual problem",
                                     ArgParseOption::STRING));
    addOption(parser, ArgParseOption("ggo", "generatorGapOpen",
                                     "Gap open costs for generating the alignment edges.",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("gge", "generatorGapExtend",
                                     "Gap extend costs for generating the alignment edges.",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("gso", "generatorSuboptimality",
                                     "suboptimality costs for generating the alignment edges.",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("lgo", "laraGapOpen",
                                     "Gap open costs for generating the alignment edges",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("lge", "laraGapExtend",
                                     "Gap extend costs for generating the alignment edges",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("ssc", "sequenceScale",
                                     "scaling factor for the scores of the alignment edges",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("rsag", "rsaGapPenalty", "gap penalty for RSA",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("stsc", "structureScoring", "scoring mode, either LOGARITHMIC,SCALING or RIBOSUM",
                                     ArgParseOption::STRING));
    addOption(parser, ArgParseOption("fsw", "fixedStructWeight", "define the weight of _half_ an interaction match for fixed structures",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("scal", "scalingFactor", "if structurescoring=SCALING then we have to give a scaling factor", ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("spseq", "splitSequence", "Define the id of the sequence that must be splitted.",
                                     ArgParseArgument::INTEGER, "INT")); // TODO fix the meaning of this parameter
    addOption(parser, ArgParseOption("ws", "windowSize", "window size specifies the length of the sliding window.",
                                     ArgParseArgument::INTEGER, "INT")); // TODO fix the meaning of this parameter
    addOption(parser, ArgParseOption("tl", "timeLimit", "some additional option.", ArgParseArgument::INTEGER, "INT"));


    addSection(parser, "Input Options");
    addOption(parser, ArgParseOption("i", "inFile", "Path to the input file", ArgParseArgument::INPUT_FILE, "IN"));
    addOption(parser, ArgParseOption("ir", "inFileRef", "Path to the reference input file", ArgParseArgument::INPUT_FILE, "IN"));
    addOption(parser, ArgParseOption("tcl", "tcoffeeLocation", "location of T-COFFEE.", ArgParseOption::STRING));

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
    setMaxValue(parser, "threads", "4");
#else
    setMaxValue(parser, "threads", "1");
#endif
//    setDefaultValue(parser, "threads", options.threadsCount);
}

// ----------------------------------------------------------------------------
// Function parseCmd()
// ----------------------------------------------------------------------------
template <typename TOption>
ArgumentParser::ParseResult parse(TOption & options, ArgumentParser & parser, int const argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;
    getOptionValue(options.verbose, parser, "verbose");
    getOptionValue(options.useBasePairs, parser, "useBasePairs");
    getOptionValue(options.affineLinearDgs, parser, "affineLinearDgs");
    getOptionValue(options.globalLocal, parser, "globalLocal");
    getOptionValue(options.unTop, parser, "unTop");
    getOptionValue(options.unLeft, parser, "unLeft");
    getOptionValue(options.unRight, parser, "unRight");
    getOptionValue(options.unDown, parser, "unDown");
    getOptionValue(options.thrGlobalUnconstr, parser, "thrGlobalUnconstr");
    getOptionValue(options.thrGlobalLocal, parser, "thrGlobalLocal");
    getOptionValue(options.thrBppm, parser, "thrBppm");
    getOptionValue(options.iterations, parser, "iterations");
    getOptionValue(options.nonDecreasingIterations, parser, "nonDecreasingIterations");
    getOptionValue(options.lowerBoundMethod, parser, "lowerBoundMethod");
    getOptionValue(options.epsilon, parser, "epsilon");
    getOptionValue(options.my, parser, "my");
    getOptionValue(options.laraScoreMatrixName, parser, "laraScoreMatrixName");
    getOptionValue(options.generatorGapOpen, parser, "generatorGapOpen");
    getOptionValue(options.generatorGapExtend, parser, "generatorGapExtend");
    getOptionValue(options.generatorSuboptimality, parser, "generatorSuboptimality");
    getOptionValue(options.laraGapOpen, parser, "laraGapOpen");
    getOptionValue(options.laraGapExtend, parser, "laraGapExtend");
    getOptionValue(options.sequenceScale, parser, "sequenceScale");
    getOptionValue(options.rsaGapPenalty, parser, "rsaGapPenalty");
    getOptionValue(options.structureScoring, parser, "structureScoring");
    getOptionValue(options.fixedStructWeight, parser, "fixedStructWeight");
    getOptionValue(options.scalingFactor, parser, "scalingFactor");
    getOptionValue(options.tcoffeeLocation, parser, "tcoffeeLocation");
    getOptionValue(options.splitSequence, parser, "splitSequence");
    getOptionValue(options.windowSize, parser, "windowSize");
    getOptionValue(options.timeLimit, parser, "timeLimit");
    getOptionValue(options.threads, parser, "threads");

    getOptionValue(options.inFile, parser, "inFile");
    if( options.inFile == "" )
    {
        std::cout << "fasta file name is " << options.inFile << std::endl; //FIXME if input file is not provided the program is stack
        return ArgumentParser::PARSE_ERROR;
    }
    getOptionValue(options.inFileRef, parser, "inFileRef");
    getOptionValue(options.outFile, parser, "outFile");
    CharString tmpDir, paramFile;
    getOptionValue(tmpDir, parser, "tmpDir");
    if (!isSet(parser, "tmpDir"))
    {
        tmpDir = getPath(options.inFile);
        if (empty(tmpDir))
            getCwd(tmpDir);
    }
    setEnv("TMPDIR", tmpDir);
    _V(options, "The absolute path where to create the tmpDir is " << tmpDir);
    _V(options, "Initialized the Options structure");
    setScoreMatrix(options);
//    showScoringMatrix(options.laraScoreMatrix);
    return ArgumentParser::PARSE_OK;
}

#endif /* INCLUDE_OPTION_SRTUCT_FUNCT_H_ */
