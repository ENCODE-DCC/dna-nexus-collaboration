function [mgSize] = getMappableGenomeSize(readLen,uMap,readsPerChr)
% Gets the total mappable size of the genome (excluding chromosomes that
% have no reads)
% function [mgSize] = getMappableGenomeSize(readLen,uMap,readsPerChr)
% --------------------------------------------------------------------------------------------------
% INPUT ARGUMENTS
% -----------------
% readLen[<double>]: array of read lengths
% uMap<struct>
%   .minReadLen<double>: min read length covered in uniqueness map
%   .maxReadLen<double>: max read length covered in uniqueness map
%   .mmNames{string}: cell array of names with paths of uniqueness map files
% readsPerChr[double]: number of reads per chromosome
% ------------------
% OUTPUT ARGUMENTS
% ------------------
% mgSize[<double>]: array of mappable size of the genome corresponding to each read length
% --------------------------------------------------------------------------------------------------

assert( ( numel(uMap.mmNames) == numel(readsPerChr) ) , ...    
    'Number of chromosomes in uMap does not match numel(readsPerChr)');

% chr with > 0 read counts
validChr = find(readsPerChr > 0);
% get unique read length values
[uniqueReadLen,~,sameReadLen] = unique(readLen);
% initialize total genome size
mgSize = zeros(size(readLen));

for iChr = 1:length(validChr)
    
    curr_uMapFile = uMap.mmNames{validChr(iChr)};
    
    % -------------
    % Read uMap
    % -------------
    tmp_uMap = fopen(curr_uMapFile,'r');
    uMapdata = fread(tmp_uMap,'*uint8');
    fclose(tmp_uMap);
    
    % --------------------------------------------------
    % For each unique ReadLen get mappable genome size
    % --------------------------------------------------
    
    for irl = 1:numel(uniqueReadLen)
        tmp_readLen = uniqueReadLen(irl); % current value of readLen
        tmp_idx = find(sameReadLen==irl); % which idx in mgSize have same readLen
        chrMgSize = nnz((uMapdata > uint8(0)) & (uMapdata <= uint8(tmp_readLen)));
        mgSize(tmp_idx) = mgSize(tmp_idx) + chrMgSize;
    end
    
end

mgSize = 2*mgSize; % both strands

end