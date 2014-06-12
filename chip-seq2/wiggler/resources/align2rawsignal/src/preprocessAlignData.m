function [alignStruct] = preprocessAlignData(alignStruct,chr,uMap,logFile)
% Reads a tagalign/BAM file , filters all non-unique tags and shifts reads in the 5' to 3' direction by half the fragment length
% function [alignStruct] = preprocessAlignData(alignStruct,chr,uMap,logFile)
% --------------------------------------------------------------------------------------------------
% INPUT ARGUMENTS:
% -------------------
% alignStruct<struct>: Structure represents an Alignment file
%   .fName<string>: name of tagAlign/BAM file
%   .format<string>: format of alignment file (tagalign/bam)
%   .chrNames{string}: cell array containing names of chrosomosomes
%   .readStart{[double]}: cell array of vectors of read start positions for each chromosome
%        reads are shifted fragLen/2 in the 5' to 3' direction in a strand dependent manner               
%   .maxReadLen<double>: maximum read length
%   .minReadLen<double>: minimum read length
%   .modeReadLen<double>: mode of read length
%   .fragLen<double>: tag extension length
%   .nReads[double]: array containing number of reads per chromosome
%   .mgSize<double>: mappable genome size
%
% chr<struct>
%     .dir<string>: directory containing chromosome files
%     .fNames{string}: names with paths of chromosome fasta files
%     .mmNames{string}: names with paths of chromosome memory mapped files
%     .names{string}: names of chromosomes
%     .len[double]: length of chromosomes
%     .nChr<double>: number of chromosomes
%
% uMap<struct>
%   .minReadLen<double>: min read length covered in uniqueness map
%   .maxReadLen<double>: max read length covered in uniqueness map
%   .mmNames{string}: cell array of names with paths of uniqueness map files
%
% logFile: path and name of log file
% -------------------
% OUTPUT ARGUMENTS:
% -------------------
% alignStruct: same as input argument
% -------------------
% NOTES:
% -------------------
% AVOID CONFUSION ABOUT READ STARTS
% [*] The readTagAlign/readBAM function will read in a tagAlign/BAM file.
% The read start positions returned are the true 5'-end position of the reads
% w.r.t. the + strand. So for a read on the + strand, the readStart value
% returned is the <chromStart> value in the tagAlign file. For a read on the
% - strand, the readStart value returned is -<chromEnd>. Note the -ve sign.
% [*] preprocessAlignData will remove multimapping reads and also shift all
% reads in the 5' to 3' direction by fragLen/2. Also the - sign is now removed.
% So all readStart positions are now positive EXCEPT for - strand reads that are 
% at the beginning of the chromosome. These will get left shited beyond 0 if they 
% are in the range [1:fragLen/2]. 
% So note that the meaning of the - sign is totally different from before.
% All - strand reads DO NOT have negative read start values.
% --------------------------------------------------------------------------------------------------

% -----------------------------
% Get format of alignment file 
% alignStruct.format
% -----------------------------
if regexpi( alignStruct.fName , '\.bam(\.gz)?$' )
    alignStruct.format = 'bam';
elseif regexpi( alignStruct.fName , '\.tagalign(\.gz)?$' )
    alignStruct.format = 'tagalign';
else
    error( '%s Unsupported file format' , alignStruct.fName );
end

% --------------------------------------------------------------------------------------------------
% Read alignment file
% alignFileData[<struct>]: array of structures, each structure represents one chromosome
%   .chrName<string>: name of the chromosome
%   .readStart[double]: start position of tag (if +ve then + strand, else - strand)
%   .readSeqLen[uint8]: lengths of reads
% --------------------------------------------------------------------------------------------------
alignStruct.chrNames = chr.names;
switch alignStruct.format
    case 'bam'
        alignFileData = readBAM( alignStruct.fName , alignStruct.chrNames , logFile );
    case 'tagalign'        
        alignFileData = readTagAlign( alignStruct.fName , alignStruct.chrNames , logFile );        
end

% --------------------------------------------------------------------------------------------------
% Get predominant read length and check if uMap supports read length
% --------------------------------------------------------------------------------------------------
alignStruct.maxReadLen = double(max(cell2mat( ...
    transpose(arrayfun(@(x) max(x.readSeqLen), alignFileData,'UniformOutput',0)) ...
    )));
alignStruct.minReadLen = double(min(cell2mat( ...
    transpose(arrayfun(@(x) min(x.readSeqLen), alignFileData,'UniformOutput',0)) ...
    )));
alignStruct.modeReadLen = double(mode(cell2mat( ...
    transpose(arrayfun(@(x) mode(double(x.readSeqLen)), alignFileData,'UniformOutput',0)) ...
    )));
assert( ((alignStruct.minReadLen >= uMap.minReadLen) && (alignStruct.maxReadLen <= uMap.maxReadLen)), ...
    '%s - Read lengths not supported by mappability data\n' , alignStruct.fName );

% --------------------------------------------------------------------------------------------------
% For each chromsome filter reads using mapability tracks
% --------------------------------------------------------------------------------------------------
for iChr = 1 : chr.nChr
    
    writeLogFile( logFile , sprintf( '-Filtering chromosome %s\n' , chr.names{iChr} ) );    
    
    % ------------------------------------
    % Remove illegal read starts > chrLen
    % ------------------------------------
    % logical indicating if abs(readstart) > chrlen
    rmIdx = ( abs( alignFileData(iChr).readStart ) > chr.len(iChr) ); 
    if nnz(rmIdx)        
        alignFileData(iChr).readStart(rmIdx) = [];
        alignFileData(iChr).readSeqLen(rmIdx) = [];
        disp('Removing reads here');
    end
    clear rmIdx;
    
    % ------------------
    % Update alignStruct
    % ------------------            
    alignStruct.readStart{iChr} = alignFileData(iChr).readStart; % Set read start positions    
    readSeqLen = alignFileData(iChr).readSeqLen; % Get read lengths    
    nReads = length(readSeqLen); % number of reads
    % Save memory by erasing data for current chromosome from alignFileData
    alignFileData(iChr).readStart = [];    
    alignFileData(iChr).readSeqLen = [];
    
    % ---------------------------
    % Get appropriate uMap file
    % ---------------------------
    tmp_uMap = fopen(uMap.mmNames{iChr},'r');
    uMapData = fread(tmp_uMap,'*uint8');
    fclose(tmp_uMap);
    
    % --------------------------------
    % Filter out non-unique positions
    % IMPLEMENTATION NOTE: uMapData has mapability values only for +
    % strand. For mapability at a position 'x' on - strand, you need to
    % look at uMapData(x-readSeqLen+1)
    % --------------------------------
    writeLogFile( logFile , sprintf( '--#Reads before uniqueness filter = %d\n' , nReads ) );
    
    % mapVals will hold mapability values corresponding to each readStart    
    mapVals = zeros(size(readSeqLen),'uint8');
    % + strand logical indices into readStart
    pIdx = (alignStruct.readStart{iChr} >= 0);
    % get mapVals for + strand positions
    mapVals(pIdx) = uMapData(alignStruct.readStart{iChr}(pIdx));
    % get mapVals for - strand positions (refer to imp. note above)
    mapVals(~pIdx) = uMapData( -alignStruct.readStart{iChr}(~pIdx) - double(readSeqLen(~pIdx)) + 1 ); 
    % valid indices
    keepIdx = ( (mapVals > uint8(0)) & (mapVals <= readSeqLen) ); 
    clear mapVals readSeqLen;
    % prune readStart
    alignStruct.readStart{iChr} = alignStruct.readStart{iChr}(keepIdx); 
    % prune + strand logical indices
    pIdx = pIdx(keepIdx); 
    clear keepIdx;
    % update number of reads
    nReads = length(alignStruct.readStart{iChr});         
    alignStruct.nReads(iChr) = nReads;
    
    writeLogFile(logFile,sprintf('--#Reads after uniqueness filter = %d\n',nReads));
        
    % ---------------------------------------------------------------------
    % Shift + strand positions to the right by floor(fragLen/2) bp
    % Shift - strand positions to the left by floor(fragLen/2) bp
    % ---------------------------------------------------------------------
    alignStruct.readStart{iChr}(pIdx) = alignStruct.readStart{iChr}(pIdx) + floor(alignStruct.fragLen/2);
    alignStruct.readStart{iChr}(~pIdx) = -alignStruct.readStart{iChr}(~pIdx) - floor(alignStruct.fragLen/2);
end
end