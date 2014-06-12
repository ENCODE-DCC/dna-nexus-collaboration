function [] = writeMatFile(fileParams,chunk,varargin)
% Saves signal vectors as chunks into a .mat file
% function [] = writeMatFile(fileParams,chunk,varargin)
% -------------------------------------------------------------------------
% Converts a number of signal vectors in varargin (invalid positions have
% value NaN) into variables in a .mat file. The signal values are broken
% into chunks of size fileParams.oChunkLen and stored as individual column
% vectors of size (chunksize X 1). 
% Chunks are named <prefix>_<chrName>_<chunkId>
% -------------------------------
% INPUT ARGUMENTS:
% -------------------------------
% fileParams: structure containing details about the output wiggle file
%           .oFname: full path and name of output .mat file
%           .oAppend: if set to 1 the file will be appended else a new file
%                     will be created
%           .oChunkLen: chunking size to use for .mat file
%           .oType: cell array of prefixes to be used for each signal
%                   vector 
% chunk: details about the current chunk of signal
%           .chrName: name of chromosome
%           .start: start coordinate
%           .stop: stop coordinate
%           .len: length of chunk
%           .isFirst: if 1 then it is the first chunk
%           .isLast: if 1 then it is the last chunk
% varargin{}[double]: each varargin is a column vector of signal values
%                     (invalid positions have value NaN)
% -------------------------------
% OUTPUT ARGUMENTS:
% -------------------------------
% None
% -------------------------------------------------------------------------

% ----------------------------------
% Find number of signal vectors
% ----------------------------------
nSig = nargin - 2;
if (nSig == 0)
    error('No signal vector to save in the mat file');
end

% --------------------------------------------
% Find number of parts and write each part
% --------------------------------------------
% length of each part
part_len = fileParams.oChunkLen;
% number of parts
nParts = ceil(chunk.len / part_len);
% global part offset in the chromosome for the chunk
part_offset = floor(chunk.start / part_len); 

for iPart = 1:nParts
    part_start = (iPart-1)*part_len + 1; % start position in the chunk
    part_stop = min(part_len*iPart , chunk.len); % stop position in the chunk
    % Create dynamic variable <prefix>_<chrName>_<absoluteChunkId>
    for iSig = 1:nSig
        tmp_var.(sprintf('%s_%s_%d', fileParams.oType{iSig}, chunk.chrName, part_offset+iPart)) = ...
            single(round(varargin{iSig}(part_start:part_stop)/1e-1)*1e-1);         
    end
    % save variable
    save(fileParams.oFname,'-struct','tmp_var','-append');
    clear tmp_var;
end
end

