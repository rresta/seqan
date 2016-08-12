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
// Author: Lily Shellhammer <>
// ==========================================================================
// This file contains routines to read and write to RNA format files (.ct)
// ==========================================================================

#ifndef SEQAN_RNA_FORMAT_READ_H
#define SEQAN_CONNEXT_FORMAT_READ_H

/* IMPLEMENTATION NOTES

RNA FORMAT example:

=> HEADER START : number of bases in the sequence
=> HEADER END: title of the structure
=> Each line has information about a base pair in the sequence
	Each line is a base, with this order of information:
		- Base number: index n
		- Base (A, C, G, T, U, X)
		- Index n-1
		- Index n+1
		- Number of the base to which n is paired. No pairing is indicated by 0 (zero).
		- Natural numbering. RNAstructure ignores the actual value given in natural numbering, 
			so it is easiest to repeat n here.

CT Files can hold multiple structures of a single sequence.
This is done by repeating the format for each structure without any blank lines between structures. 

HEADER
 N  SEQUENCE   N-1  	 N+1	J POSITION  N  
 1 	G       	0    	2   	72    		1
 2 	C       	1    	3   	71    		2
 3 	G       	2    	4   	70    		3


*/
 namespace seqan
{

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function _readHeader(); RNAHeader
// ----------------------------------------------------------------------------

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
_readHeader(RNAHeader & header,
           RNAIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           RNA const & /*tag*/)
{
    clear(header);
    CharString buffer;
    RNAHeaderRecord record;

    clear(buffer);  ////////I never used buffer, do I need it at all??
    clear(record);

                        ///////// Need to read energy amount after = and figure out how to store that value and use the Assert Functor??
    readUntil(record.amount, IsWhitespace());
                            /////Read until Energy??? how to do??
    readUntil(record.amount, iter, OrFunctor<EqualsChar<'='>, AssertFunctor<NotFunctor<IsNewline>, ParseError, RNA> >());      

    // Skip '='.
    skipOne(iter);
    skipUntil(iter, NotFunctor<IsWhitespace>());

    readUntil(record.energy, iter, IsWhitespace()); //read energy
    skipUntil(iter, NotFunctor<IsWhitespace>());    //skip tab

    readUntil(name, IsNewLine());   //read name until end of header
    appendValue(header, record);    //append values captured in record

}


// ----------------------------------------------------------------------------
// Function _readRecord(); RNAHeader
// ----------------------------------------------------------------------------
template <typename TNIndex, typename TSeqChar, typename TBasePairPos, typename TFwdIterator>
inline void _readRecord(TNIndex & index, TSeqChar & base, TBasePairPos & pair, TFwdIterator & iter, RNA const & RNA)
{
    typedef OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, RNA> > NextEntry;               //This was in Vcf io stuff and I don't know if I need it here

    /*
    typedef typename FastaIgnoreOrAssertFunctor_<TSeqAlphabet>::Type        TSeqIgnoreOrAssert;         //I see how ^ is related to assert fn but still lost about how to implement and use
    typedef typename FastaIgnoreFunctor_<TQualAlphabet>::Type               TQualIgnore;
    typedef typename FastaIgnoreOrAssertFunctor_<TQualAlphabet>::Type       TQualIgnoreOrAssert;
	*/

    clear(index);
    clear(base);
    clear(pair);

    CharString buffer;  ///////how to use (reiteration of quesiton above^^)
    clear(buffer);

    skipOne(iter);  //All records start with a space char
    readUntil(index, iter, IsWhitespace());             //read index position
    readUntil(seq, iter, IsWhitespace());               //read base 
    skipUntil(iter, NotFunctor<IsWhitespace>());
    skipUntil(iter, IsWhitespace());                    //skip both redundant indices
    skipUntil(iter, NotFunctor<IsWhitespace>());
    readUntil(pair, iter, IsWhitespace());              //read pair position for base
    skipUntil(iter, IsNewLine());

    ///////Do I use IsNewLine or do I use NextEntry??


    /////////What is this for reading and how does it work 

    /*if (IsBlank()(value(iter)))
    {
        skipUntil(iter, NotFunctor<IsBlank>());
        readLine(val, iter);
    }

    while (IsBlank()(value(iter)))
    {
        appendValue(val, '\n');
        readLine(val, iter);
    }
    */
}



            ////// This function has lots of parameters. In most IO files, you only pass the filetypeOut and the record. Is there a function elesewhere
            //that directs it here? And if so, does it work for the templated type I have or do i have to write code elsewehre to direct it here
            //and pass the correct parameters?

// ----------------------------------------------------------------------------
// Function _writeRecord(RNA);
// ----------------------------------------------------------------------------

template <typename TNIndex, typename TSeqChar, typename TBasePairPos>
inline void
_writeRecord(TTarget & target,
            TIdString const & index,
            TSeqString const & base,
            TBasePairPos const & pair,
            RNA const & /*tag*/,
            SequenceOutputOptions const & options = SequenceOutputOptions())            ///////RNA const& is this jus to direct writeRecord to this fn with tag of RNA file
{
    writeValue(target, ' ');    //All records start with a space
    writeValue(target, '\t');
    write(target, index);
    //Need to implement way to increse/decrease index by one and also print that to ct file
    writeValue(target, '\t');
    write(target, base);
    writeValue(target, '\t\t');
    //IMPLEMENT: n-1 tab n+1
    write(target, pair);

                                            /////////////////////what is this?? it was in fasta_fastq file at the end of writing characters and values
   writeWrappedString(target, seq, (options.lineLength < 0)? 70 : options.lineLength); // 70bp wrapping, by default
}

// ----------------------------------------------------------------------------
// Function _writeHeader(); RNAHeader
// ----------------------------------------------------------------------------

template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
_writeHeader(TTarget & target,
            RNAHeader const & header,
            RNAIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
            RNA const & /*tag*/)
{
                                ///////Can I do this or how do I access members of header?? With . operator?
    write(target, amount);
    writeValue(target, '\t')
    writeValue(target, "ENERGY = ");
    writeValue(target, '\t')
    write(target, energy);
    writeValue(target, '\t')
    write(target, name);
    writeValue(target, '\n');
}



} //namespace seqan

#endif 
//SEQAN_RNA_FORMAT_READ_H