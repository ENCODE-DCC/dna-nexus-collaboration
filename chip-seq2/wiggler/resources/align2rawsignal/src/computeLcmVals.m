function [lcmVals] = computeLcmVals( uMapFileData , readLen , fragLen , winSize , kernel , tukeyRatio , chunk )
% Calculates local cumulative mappability values for a chunk of a mappability file
% function [lcmVals] = computeLcmVals( uMapFileData , readLen , fragLen , winSize , kernel , tukeyRatio, chunk )
% -------------------------------------------------------------------------
% INPUT ARGUMENTS:
% ------------------
% uMapFileData: uint8 mappability values for a chromosome
% readLen: mode length of reads
% fragLen: fragment length
% winSize: smoothing window size
% kernel: smoothing kernel
% tukeyRatio: taper to fixed ratio for tukey kernel (set to NaN if not used)
% chunk: details about the chunk
%           .chrName: name of chromosome
%           .start: start coordinate
%           .stop: stop coordinate
%           .len: length of chunk
%           .isFirst: if 1 then it is the first chunk
%           .isLast: if 1 then it is the last chunk
% ------------------
% OUTPUT ARGUMENTS:
% ------------------
% lcmVals[double]: max tags output for that chunk
% ------------------
% DEPENDENCIES:
% ------------------
% cumsumwin()
% -------------------------------------------------------------------------

% uMapFileData = memmapfile(uMapFile,'format','uint8');

% ---------------------------------------------------------------------
% smoothing window on either side of each position
% (i-smoothLen:i:i+smoothLen)
% ---------------------------------------------------------------------
smoothLen = floor(winSize/2);
tagShift = floor(fragLen/2);

% ---------------------------------------------------------------------
% extended length of lcmVals
% ---------------------------------------------------------------------
extension = max(smoothLen,tagShift);
chunkValLen = chunk.len + 2*extension;

% ---------------------------------------------------------------------
% Get uMap values from + strand
% tmp_start:tmp_stop index into the uniqueness map
% ---------------------------------------------------------------------

% initialize + strand uniqueness map
tmp_vals = zeros( chunkValLen , 1 , 'uint8' );

chunkStatus = ...
    1*(chunk.isFirst && chunk.isLast) + ...
    2*(chunk.isFirst && ~chunk.isLast) + ...
    3*(~chunk.isFirst && chunk.isLast);

switch chunkStatus
    case 1 % first and last chunk (done)
        tmp_start = chunk.start;
        tmp_stop = chunk.stop - (readLen-1);
        tmp_vals( ( 1 + extension + tagShift ) :  ( chunkValLen - extension + tagShift - (readLen-1) ) ) = uMapFileData(tmp_start:tmp_stop);
    case 2 % first chunk (done)
        tmp_start = chunk.start;
        tmp_stop = chunk.stop + extension - tagShift;        
        tmp_vals( ( 1 + extension + tagShift ) : chunkValLen) = uMapFileData(tmp_start:tmp_stop);
    case 3 % last chunk (done)
        tmp_start = chunk.start - extension - tagShift;
        tmp_stop = chunk.stop - (readLen-1);
        tmp_vals(1 : ( chunkValLen - extension + tagShift - (readLen-1) ) ) = uMapFileData(tmp_start:tmp_stop);
    case 0 % middle chunk (done)
        tmp_start = chunk.start - extension - tagShift;
        tmp_stop = chunk.stop + extension - tagShift;
        tmp_vals = uMapFileData(tmp_start:tmp_stop);
    otherwise
        error('Invalid chunk type');
end

% set coordinates that have valid uMap values
lcmVals = double( (tmp_vals > uint8(0)) & (tmp_vals <= uint8(readLen)) );

% ---------------------------------------------------------------------
% Get uMap values from - strand
% tmp_start:tmp_stop index into the uniqueness map
% ---------------------------------------------------------------------

% initialize - strand uniqueness map
tmp_vals = zeros(chunkValLen,1,'uint8');

switch chunkStatus
    case 1 % first and last chunk (done)
        tmp_start = chunk.start;
        tmp_stop = chunk.stop - (readLen - 1);
        tmp_vals( ( 1 + extension - tagShift + (readLen - 1) ) : ( chunkValLen - extension - tagShift ) ) = uMapFileData(tmp_start:tmp_stop);
    case 2 % first chunk (done)
        tmp_start = chunk.start;
        tmp_stop = chunk.stop + extension + tagShift - (readLen - 1);
        tmp_vals( ( 1 + extension - tagShift + (readLen - 1) ) : chunkValLen) = uMapFileData(tmp_start:tmp_stop);
    case 3 % last chunk (done)
        tmp_start = chunk.start - extension + tagShift - (readLen - 1);
        tmp_stop = chunk.stop - (readLen - 1);
        tmp_vals( 1 : ( chunkValLen - extension - tagShift ) ) = uMapFileData(tmp_start:tmp_stop);
    case 0 % middle chunk (done)
        tmp_start = chunk.start - extension + tagShift - (readLen - 1);
        tmp_stop = chunk.stop + extension + tagShift - (readLen - 1);
        tmp_vals = uMapFileData(tmp_start:tmp_stop);
    otherwise
        error('Invalid chunk type');
end

% add to coordinates that have valid uMap values
lcmVals = lcmVals + double((tmp_vals > uint8(0)) & (tmp_vals <= uint8(readLen)));
clear tmp_vals;

% ----------------------------------
% Smooth window sum of local mappability values
% ----------------------------------
lcmVals = kdesmooth( lcmVals , smoothLen , kernel , 0 , tukeyRatio);

% ---------------------------------------------------------------------
% Remove the first (smoothLen) positions and the last (smoothLen)
% positions
% ---------------------------------------------------------------------
lcmVals = lcmVals( (extension+1) : (end-extension) );
end