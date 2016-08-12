///////////////////////////
///////////////////////////
//ALL (INCOMPLETE) PROTOTYPES FOR THE FUNCTIONS IN THE DIFFERENT FILES I'M USING. FOR REFERENCE AS TO WHAT FUNCTIONS I HAVE AND WHERE I DEFINE THEM
///////////////////////////
///////////////////////////

/*
RNA_FORMAT_FILE_H
*////////////////////////////
struct TagRNA_;
struct MagicHeader<Bed, T> {};

struct FormattedFileContext<FormattedFile<RNA, TDirection, TSpec>, TStorageSpec> {};

struct FileExtensions<RNA, T> {};

typedef FormattedFile<RNA, Input>   RNAFileIn;
typedef FormattedFile<RNA, Output>  RNAFileOut;

readHeader(RNAHeader & header, FormattedFile<RNA, Input, TSpec> & file);
readRecord(RNARecord & record, FormattedFile<RNA, Input, TSpec> & file);
writeRecord(FormattedFile<RNA, Output, TSpec> & file, RNARecord & record);
writeHeader(FormattedFile<RNA, Output, TSpec> & file, RNAHeader & header);

/*
RNA_FORMAT_HEADER_H
*////////////////////////////
class RNAHeaderRecord {};
inline void clear(RNAHeaderRecord & record);
typedef String<RNAHeaderRecord> RNAHeader;

/*
RNA_FORMAT_RECORD_H
*////////////////////////////
class RNARecord {};
inline void clear(RNARecord & record)
// Do I need a similar function to the third shown above in header.h file

/*
RNA_FORMAT_IO_H
*////////////////////////////
struct SequenceOutputOptions {};

readHeader(RNAHeader & header,
           RNAIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           RNA const & /*tag*/);

readRecord(TNIndex & index, TSeqChar & base, TBasePairPos & pair, TFwdIterator & iter, RNA const & RNA);

writeRecord(TTarget & target,
            TIdString const & index,
            TSeqString const & base,
            TBasePairPos const & pair,
            RNA const & /*tag*/,
            SequenceOutputOptions const & options = SequenceOutputOptions());

writeHeader(TTarget & target,
            RNAHeader const & header,
            RNAIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
            RNA const & /*tag*/);

