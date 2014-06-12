function [] = writeWiggle(fileParams,chunk,signal)
% Converts a signal vector into a wiggle file
% function [] = writeWiggle(fileParams,signal,chunk)
% -------------------------------------------------------------------------
% INPUT ARGUMENTS:
% ------------------
% fileParams: structure containing details about the output wiggle file
%           .oFname: full path and name of output wiggle file
%                   if '' then print to STDOUT
%           .oAppend: if set to 1 the file will be appended else a new file
%           will be created
% chunk: details about the current chunk of signal
%           .chrName: name of chromosome
%           .start: start coordinate
%           .stop: stop coordinate
%           .len: length of chunk
%           .isFirst: if 1 then it is the first chunk
%           .isLast: if 1 then it is the last chunk
% signal[double]: column vector of signal values (invalid positions have
%                 value NaN)
% ------------------
% OUTPUT ARGUMENTS:
% ------------------
% None
% -------------------------------------------------------------------------

signal = round( signal / 1e-1 ) * 1e-1;
% --------------------------------------------------------------
% Convert chunk structure fields into scalars to speed up loops
% --------------------------------------------------------------
chunk_start = chunk.start;
chunk_chrName = chunk.chrName;

% -------------------------------------------------------------------------
% Find valid parts (finite signal) inside the chunk
% -------------------------------------------------------------------------
validIdx = diff( [ 0 ; isfinite(signal) ; 0 ] );
part_start = find( validIdx == 1 ); % start idx
part_stop = find( validIdx == -1 ) - 1; % stop idx
nParts = numel(part_start); % number of valid contiguous chunks
clear validIdx;

% -------------------------------------------------------------------------
% Open wiggle file
% -------------------------------------------------------------------------
% file permission
fPerm = char( fileParams.oAppend * 'a' + ~fileParams.oAppend * 'w' );
switch lower(fileParams.oFname)
    case {'stdout',''}
        fp = 1; % STDOUT
    case 'stderr'
        fp = 2; % STDERR
    otherwise
        fp = fopen(fileParams.oFname,fPerm);    
end
assert( (fp~=-1) , 'Unable to open wiggle file for writing\n' );

% -------------------------------------------------------------------------
% Write each valid part into wiggle file
% WIGGLE FILE FORMAT
% fixedStep  chrom=chrN  start=position(1-based)  step=stepInterval  [span=windowSize]
% dataValue1
% dataValue2
% OR
% variableStep  chrom=chrN  [span=windowSize]
% chromStartA(1-based)  dataValueA
% chromStartB  dataValueB
% -------------------------------------------------------------------------

for iPart = 1:nParts
    fprintf( fp , 'fixedStep  chrom=%s  start=%d  step=1\n', ...
        chunk_chrName , (chunk_start - 1) + part_start(iPart) );
    fprintf( fp , '%.1f\n' , signal( part_start(iPart) : part_stop(iPart) ) );
end

% -------------------------------------------------------------------------
% Close file
% -------------------------------------------------------------------------
if ~any( fp == [1,2] )
    fclose(fp);
end
end
