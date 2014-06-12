function [signalVals] = computeSignalVals(readStart , winSize , kernel , tukeyRatio , chunk)
% Extracts, aggregates and tag extends reads that belong to a chunk
% function [signalVals] = computeSignalVals(readStart , winSize , kernel , tukeyRatio , chunk)
% -------------------------------------------------------------------------
% INPUT ARGUMENTS:
% --------------------
% readStart[double]: start indices of the reads (can be -ve values because - strand
%            locations are shifted to the left by winSize-1)
% winSize: smoothing window size
% kernel: smoothing kernel
% tukeyRatio: kernel tukey Ratio (set to NaN if not used)
% chunk: details about the chunk
%           .chrName: name of chromosome
%           .start: start coordinate
%           .stop: stop coordinate
%           .len: length of chunk
%           .isFirst: if 1 then it is the first chunk
%           .isLast: if 1 then it is the last chunk
% --------------------
% OUTPUT ARGUMENTS:
% --------------------
% signalVals[double]: signal values (column vector)
% --------------------
% DEPENDENCIES:
% --------------------
% cumsumwin()
% --------------------
% NOTES:
% --------------------
% AVOID CONFUSION ABOUT READ STARTS
% [*] readTagAlign()/reaBAM() will read in a tagAlign/BAM file. The read
% start positions returned are the 5'-end position of the reads w.r.t. the
% + strand. So for a read on the + strand, the readStart value returned is
% the chromStart value in the tagAlign file. For a read on the - strand,
% the readStart value returned is -chromEnd. Note the -ve sign.
% [*] preprocessAlignData() will remove multimapping reads and also move
% all - strand reads to the left by the tag extension value. Also the -
% sign is now removed. So all readStart positions are now positive EXCEPT
% for - strand reads that are at the beginning of the chromosome. These
% will get left shited beyond 0 if they are in the range [1:winSize]. So
% note that the meaning of the - sign is totally different from before. All
% - strand reads DO NOT have negative read start values.
% -------------------------------------------------------------------------

% ---------------------------------------------------------------------
% smoothing extension on either side of each position
% (i-smoothLen:i:i+smoothLen)
% ---------------------------------------------------------------------
smoothLen = floor(winSize/2);

% -------------------------
% Initialize signalVals
% -------------------------
chunkValLen = chunk.len + 2*smoothLen; % extended length of signalVals
signalVals = zeros( chunkValLen , 1 , 'double' ); % initialize signalVals

% -------------------------
% Chunk boundaries
% -------------------------
tmp_start = chunk.start - smoothLen;
tmp_stop = chunk.stop + smoothLen;
validIdx = (readStart >= tmp_start) & (readStart <= tmp_stop); % reads that fall within chunk

if ~isempty(validIdx)
    % -------------------------
    % Offset the readStarts
    % -------------------------
    tmp_readStart = readStart(validIdx) - tmp_start + 1;
    
    % -------------------------
    % Aggregate reads
    % -------------------------
    signalVals = accumarray(tmp_readStart,1,size(signalVals)); % accumulate reads
    
    % -------------------------
    % Smooth window sum of read counts
    % -------------------------    
    signalVals = kdesmooth( signalVals , smoothLen , kernel , 0 , tukeyRatio);
end

% ---------------------------------------------------------------------
% Remove the first (smoothLen) positions and the last (smoothLen)
% positions
% ---------------------------------------------------------------------
signalVals = signalVals( (smoothLen+1) : (end-smoothLen) );

end