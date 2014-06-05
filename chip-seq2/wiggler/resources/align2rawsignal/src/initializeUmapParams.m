function [uMap] = initializeUmapParams(iParams,chr)
% Create datastructure uMap containing uniqueness map related information
% function [uMap] = initializeUmapParams(uMapDir,chrNames,logFile)
% -------------------------------------------------------------------------
% INPUT ARGUMENTS:
% ---------------------
% iParams<struct>
%     .logFile<string> : log file
%     .alignFname{<string>}: alignment file names
%     .seqDir<string>: sequence directory
%     .uMapDir<string>: mappability tracks directory
%     .outFile.name<string>: name of output file
%     .outFile.format<string>: format of output file
%     .outFile.lcmFile<string>: local cumulative mappability file
%     ,normFlag<double>: normalization flag
%     .smooth.fragLen[double]: fragment length
%     .smooth.winLen<double>: smoothing window size
%     .smooth.kernel<string>: smoothing kernel
%     .mapFilter<double>: local cumulative mappability filter
%     .maxMem<double>: approximate maximum memory
%     .outChunk<double>: output chunk length
%     .processChunk<double>: processing chunk length
%
% chr<struct>
%     .dir<string>: directory containing chromosome files
%     .fNames{string}: names with paths of chromosome fasta files
%     .mmNames{string}: names with paths of chromosome memory mapped files
%     .names{string}: names of chromosomes
%     .len[double]: length of chromosomes
%     .nChr<double>: number of chromosomes
% ---------------------
% OUTPUT ARGUMENTS:
% ---------------------
% uMap<struct>
%   .minReadLen<double>: min read length covered in uniqueness map
%   .maxReadLen<double>: max read length covered in uniqueness map
%   .mmNames{string}: cell array of names with paths of uniqueness map files
% -------------------------------------------------------------------------

writeLogFile( iParams.logFile , sprintf('Preprocessing uniqueness maps\n') );
minmax = regexp( iParams.uMapDir , ['globalmap_k(\d+)tok(\d+)' , regexptranslate( 'escape' , filesep() ) , '?$' ] , 'tokens' );
try
    uMap.minReadLen = str2double(minmax{1}{1}); % min read length covered in uniqueness map
    uMap.maxReadLen = str2double(minmax{1}{2}); % max read length covered in uniqueness map
catch ME
    error('Mappability track directory MUST end in globalmap_k[kmin]tok[kmax]');
end

assert( (uMap.minReadLen <= uMap.maxReadLen) && all( isfinite( [uMap.minReadLen , uMap.maxReadLen ] ) ) , ...
    'Illegal values of mink and maxk in Mappability track directory name' );

nChr = length(chr.names);

for iChr = 1:nChr
    tmpName = dir( fullfile( iParams.uMapDir , [ chr.names{iChr} , '.*.unique' ] ) );
    tmpName = tmpName.name;
    uMap.mmNames{iChr} = fullfile( iParams.uMapDir , tmpName );
    assert( logical( exist( uMap.mmNames{iChr} , 'file' ) ) , ...
        'Mappability track %s does not exist\n', ...
        uMap.mmNames{iChr});    
    % -----------------------------------------------
    % Check that uMaps have same length as sequences
    % -----------------------------------------------    
    assert( ( chr.len(iChr) == getUMapLength(uMap.mmNames{iChr}) ), ...        
        'Sequence length and mappability track length do not match for chromosome %s' , chr.names{iChr} );
end

end

% ##################################################################################################
% AUXILIARY FUNCTION
% ##################################################################################################

function [uMapLen] = getUMapLength(uMapFile)
% Gets the length of a mappability track
% function [uMapLen] = getUMapLength(uMapFile)
% --------------------------------------------------------------------------------------------------
% INPUT ARGUMENTS
% -----------------
% uMapFile<string>: Path and name of a memory mappable mappability/uniqueness track (uint8)
% -----------------
% OUTPUT ARGUMENTS
% -----------------
% uMapLen<double>: length of mappability track
% --------------------------------------------------------------------------------------------------

uMap = memmapfile(uMapFile,'format','uint8');
uMapLen = numel(uMap.data);
end