// ----------------------------------------------------------------------------
// std headers
// ----------------------------------------------------------------------------

#include <iostream>
#include <map>
#include <string>
#include <sstream>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/arg_parse.h>
#include <seqan/store.h>
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/align.h>
#include <seqan/graph_align.h>
//#include <seqan/align_profile.h>
#include <seqan/align_rna.h>
#include <seqan/vcf_io.h>
#include <seqan/bpseq_io.h>

#define START 0
#define CONTINUE 1

typedef unsigned TPosition;
typedef float TScoreValue;
typedef seqan::CharString TString;

template <typename TString, typename TPosition>
struct fixedStructElement {
    TString method; // place the method and parameters used to compute the structure
//	seqan::String<unsigned> structure;
//    seqan::String<TPosition> seqPos;
    seqan::String<TPosition> interPos;
};

struct RNARecord
{
    //index (in case we have bases in the middle of a strand, or only a segment of the strand)
    String<int>  index;

    //string of base at each position in RNA strand
    Rna5String seq;

    // Position of n base's pair.
    String<fixedStructElement<TString, TPosition> > pair;  // A list of index can be saved for eac method used to copute the interactions

    //Qual, and information specific to other file formats. Will be set to default value in constructor when I figure out what to set them to.
    CharString qual;

    int32_t energy;
    //RDAT has information such as:
    //Reactivity at each position, reactivity error at each position, XSEl, XSEL_REFINE
    //Comments, annotations
};

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------
//TODO this can be re-written as a set of three maps with keys and values but
//if in the future should be embedded in the parsing function is better to
//conserve this style
struct Options
{
    // Name of input file (default: stdout)
    seqan::CharString inFile;
    Options() :
            inFile("/demo/dotknot.bpseq&/demo/ipknot3seq.bpseq")
    {}
};

template <typename TSeq1, typename TSeq2>
int seqCompare(TSeq1 const & seq1, TSeq2 const & seq2)
{
// Create a c-style string object for str:
    String<char, CStyle> cStr1 = seq1;
    String<char, CStyle> cStr2 = seq2;
// Now use cStyle as char array:
    return strcmp(cStr1, cStr2);
}

template <typename TElement>
void elementClear(TElement & element)
{
    clear(element.method);
//    clear(element.seqPos);
    clear(element.interPos);
}

template <typename TOption, typename TFileList>
void splitBpseqFilenames(TOption const & options, TFileList & fileList)
{
    char fileNames[length(options.inFile)+1]; //= toCString(options.cmd.bpseqFile).str();
    std::cout << options.inFile << std::endl;
    for(unsigned i=0;i<length(options.inFile); ++i ) //the added vector cell is used to save the string terminator
        fileNames[i]=options.inFile[i];
    fileNames[length(options.inFile)]='\0';
    char seps[]   = "&\0";
    char *token;
//  Establish string and get the first token:
    token = strtok( fileNames, seps );
    while( token != NULL )
    {
//  While there are tokens in "string"
        std::cout << token << std::endl;
        fileList.push_back(token);
//  Get next token:
        token = strtok( NULL, seps );
    }
}


template <typename TVect, typename TFileList>
void readBpseqUpdateRnas(TVect & rnaSeqs, TFileList const & fileList)
{
    CharString fileBp;
    CharString tmpStr;
    bool flagSeqUsed=START;
    fixedStructElement<TString, TPosition> tmpStructPairMate;
    RNARecord tmpStructSeq;
    for(unsigned i=0;i<fileList.size();++i)
    {
        std::cout << fileList[i] << std::endl;
// Get path to example file.
        fileBp = getAbsolutePath(toCString(fileList[i]));
// Open Bpseq input file.
        BpseqFileIn bpseqIn(toCString(fileBp));
// Get file extension
        std::vector<std::string> formats = getFileExtensions(bpseqIn);
        std::cout << "file format in bpseq = " << formats[0] << std::endl;
// Attach Bpseq to standard output.
        BpseqFileOut bpseqOut(bpseqIn);
        open(bpseqOut, std::cout, Bpseq());
//Open an output file
        BpseqFileOut bpseqOutB(bpseqIn);
        tmpStr = "/demo/"+std::to_string(i)+"-seq.bpseq";
        open(bpseqOutB, toCString(getAbsolutePath(toCString(tmpStr)))); //TODO for now files are called with numbers but they can be called with the original file name
// Check the output file format
        formats = getFileExtensions(bpseqOut);
        std::cout << "file format out bpseq = " << formats[0] << std::endl;
        std::cout << std::endl;
// Write all the output on the std::out and on the files
        BpseqHeader headerBp;
        BpseqRecord recordBp;
        while (!atEnd(bpseqIn))
        {
// Copy over headerBp.
            readHeader(headerBp, bpseqIn);
//          writeHeader(bpseqOut, headerBp);  // Print on video the file contents
            writeHeader(bpseqOutB, headerBp);  // Print on an output file the file contents
// Copy the Bpseq file record by record.
            readRecord(recordBp, bpseqIn);
//          writeRecord(bpseqOut, recordBp); // Print on video the file contents
            writeRecord(bpseqOutB, recordBp); // Print on an output file the file contents
            flagSeqUsed=START;
            for(unsigned i=0; i<length(rnaSeqs);++i)
            {
                if(!seqCompare(rnaSeqs[i].seq, recordBp.seq))
                {
//                  std::cout << "seq recordBp.seq" << " = " << recordBp.seq << std::endl;
//                  std::cout << "seq rnaSeqs" << i << " = " << rnaSeqs[i].seq << std::endl;
//                  std::cout << "header " << headerBp[0].key << tmpStr << std::endl;
//                  tmpStructPairMate.method = fileList[i];
                    tmpStructPairMate.method = headerBp[0].value; // TODO in this position must be added the tool name and the parameters and not the sequence name
//                    tmpStructPairMate.seqPos = recordBp.seqPos;
                    tmpStructPairMate.interPos = recordBp.interPos;
                    appendValue(rnaSeqs[i].pair, tmpStructPairMate);
                    flagSeqUsed=CONTINUE;
                    elementClear(tmpStructPairMate);
                }
            }
            if(flagSeqUsed==START)
            {
// create the new item in the sequence structure
                tmpStructPairMate.method = headerBp[0].value;  // TODO in this position must be added the tool name and the parameters and not the sequence name
//              tmpStructPairMate.method = fileList[i];
                tmpStructSeq.index = recordBp.seqPos;
                tmpStructPairMate.interPos = recordBp.interPos;
// Add the new reported sequence in the rnaSeqs vector
                tmpStructSeq.seq = recordBp.seq;
                appendValue(tmpStructSeq.pair, tmpStructPairMate);
                appendValue(rnaSeqs,tmpStructSeq);
//              std::cout << "seq recordBp.seq has not been used" << " = " << recordBp.seq << std::endl;
                flagSeqUsed=CONTINUE;
                elementClear(tmpStructPairMate);
            }
        }
// Print the last cached sequences from .bpseq file
//      for (unsigned i = 0; i < length(recordBp.seqPos); ++i)
//      {
//          std::cout << recordBp.seqPos[i] << "\t";
//      }
//      std::cout << std::endl;
//      for (unsigned i = 0; i < length(recordBp.seq); ++i)
//      {
//          std::cout << recordBp.seq[i] << "\t";
//      }
//      std::cout << std::endl;
//      for (unsigned i = 0; i < length(recordBp.interPos); ++i)
//      {
//          std::cout << recordBp.interPos[i] << "\t";
//      }
//      std::cout << std::endl;
//      std::cout << recordBp.seq << std::endl;

// Print all the content of rnaSeqs vector

// Fill the RNA fields
//      appendValue(header, record);
// Check the rnaSeqs[i].seq vector and in case you find the right sequence append the structure info
    }
}
//template <typename TVect>
//void createRnaBppGraph(TVect & rnaSeqs_i, unsigned const & seqSize)
//{
//  TUVertexDescriptor uVertex[seqSize];
//  TUEdgeDescriptor uEdge[seqSize];
//  rnaSeqs_i.seq;
//}



int main(int argc, char const ** argv)
{
    if (argc==1)
        std::cout << "type ./bpseq_io --help to get the parameters table" << std::endl;

    seqan::ArgumentParser parser;

    Options options;

//  Read file list produced by Ipknot, dotknot and other tools
    std::vector<std::string> fileList;
    splitBpseqFilenames(options, fileList);

    String<RNARecord> rnaSeqs;

//  Read .bpseq file list and update the data structure
    readBpseqUpdateRnas(rnaSeqs, fileList);

//  TODO add the .bpseq writer
    return 0;
}
