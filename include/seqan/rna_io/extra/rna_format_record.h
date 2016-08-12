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
// Author: Lily Shellhammer
// ==========================================================================

// TODO(holtgrew): Get out more than just strings...

#ifndef SEQAN_RNA_RECORD_H_
#define SEQAN_RNA_RECORD_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class VcfRecord
// ----------------------------------------------------------------------------
//
/*!
 * @class RNARecord
 * @implements FormattedFileRecordConcept
 * @headerfile <seqan/RNA_io.h>
 * @brief Information for one RNA record.
 *
 * @signature class RNARecord;
 *
 * @section Remarks
 *
 * _______________?
 *
 * //COPIED FROM VCF NOT MY WRITING//Invalid qualities are stored as a float <tt>NaN</tt> (not a number).  To test a float quality <tt>q</tt> for being
 * ////<tt>NaN</tt>, test for <tt>q != q</tt>.  Only <tt>NaN</tt> has the property that <tt>NaN != NaN</tt>.
 */

/*!
 * @fn VcfRecord::MISSING_QUAL                    //Again, don't know how to implement this yet
 * @brief Return IEEE <tt>NaN</tt> float value.
 *
 * @signature float VcfRecord::MISSING_QUAL();
 *
 * @section Remarks
 *
 * This cannot be implemented portably as a constant in a header-only library.  In C++11, there is a function
 * <tt>std::nanf()</tt> that could be used instead.  For C++98, this is the best we can do.
 */

/*!
 * @fn RNARecord::RNARecord
 * @brief Default constructor.
 *
 * @signature RNARecord::RNARecord();
 */

/*!
 * @var VariableType RNARecord::index
 * @brief Position of base.
 *
 * @var VariableType RNARecord::base
 * @brief Base at n index 
 *
 * @var VariableType RNARecord::INVALID_REFID
 * @brief Static member as marker for invalid reference (<tt>int32_t</tt>)
 */

class RNARecord
{
public:
    
    // Constant for invalid index
    static const int32_t INVALID_INDEX = -1;        //HOW TO IMPLEMENT SOMETHING LIKE THIS/IS IT NECESSARY FOR .CT?
    // Constant for invalid pair to index's base
    static const int32_t INVALID_PAIR = -1;
    

    // Index of n base in sequence.
    int32_t index;                                          //IS THIS THE CORRECT DATA TYPE
    // Base at n index in one-letter notation.
    CharString base;
    // Position of n base's pair.
    int32_t pair;
    
    
    // Default constructor.
    RNARecord() : index(INVALID_INDEX), pair(INVALID_PAIR)      ///I understand what is happening here conceptually
    {}                                                                                      //but I don't know how to recreate it for a simpler file

    // Actually this is IEEE NaN.
                                    //CALLED IN CONSTRUCTOR LIKE THIS qual(MISSING_QUAL() ????
    /*
    static float MISSING_QUAL()
    {
        union {
            uint32_t u;
            float f;
        } tmp;
        tmp.u = 0x7F800001;
        return tmp.f;
    }
    */
   
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn RNARecord#clear
 * @brief Clear a RNARecord.
 *
 * @signature void clear(record);
 *
 * @param[in,out] record The RNARecord to clear.
 */

inline void clear(RNARecord & record)
{
    clear(record.index);
    clear(record.base);
    clear(record.pair);

}

}  // namespace seqan

#endif  
//SEQAN_RNA_RECORD_H
