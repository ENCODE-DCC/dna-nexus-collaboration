function [chr] = initializeChromosomeParams(iParams)
% Create datastructure chr containing chromosome related information,
% convert fasta files to memory mapped files
% function [chr] = initializeChromosomeParams(iParams)
% -------------------------------------------------------------------------
% INPUT ARGUMENTS:
% ----------------
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
% ------------------
% OUTPUT ARGUMENTS:
% ------------------
% chr<struct>
%     .dir<string>: directory containing chromosome files
%     .fNames{string}: names with paths of chromosome fasta files
%     .mmNames{string}: names with paths of chromosome memory mapped files
%     .names{string}: names of chromosomes
%     .len[double]: length of chromosomes
%     .nChr<double>: number of chromosomes
% -------------------------------------------------------------------------

writeLogFile( iParams.logFile , sprintf('Preprocessing chromosome files\n') );

% set chr.dir
chr.dir = iParams.seqDir;

% names of chromosomes
tmpNames = dir( fullfile( chr.dir , 'chr*.fa' ) ); % get names of chromosome files
tmpNames = {tmpNames.name};
chr.names = regexprep( tmpNames , '(chr.+)\..*' , '$1' ); 

% number of chromosomes
chr.nChr = length(tmpNames); 
assert( (chr.nChr > 0) , 'Chromosome *.fa files are missing' );

% names of fasta file names (with path)
chr.fNames = cellfun( @(x) fullfile(chr.dir,x) , tmpNames , 'UniformOutput' , 0 ); 

% convert the fasta files to memory mapped files
chr.mmNames = cellfun( @fastaToMemmap , chr.fNames , repmat( {iParams.logFile} , size(chr.fNames) ) , ...
    'UniformOutput' , 0 ); 

% initialize lengths of chromosomes
chr.len = zeros(chr.nChr,1); 
for ichr = 1:chr.nChr
    mm = memmapfile( chr.mmNames{ichr} , 'format' , 'uint8' );
    chr.len(ichr) = length(mm.data);
end

end

% *************************************************************************
% AUXILIARY FUNCTION
% *************************************************************************
function [mmFile] = fastaToMemmap(fastaFile,logFile)
% The function converts a FASTA file into memory mappable uint8 binary file
% function [mmFile] = fastaToMemmap(fastaFile,logFile)
% -------------------------------------------------------------------------
% INPUT ARGUMENTS:
% -----------------
% fastaFile: Input FASTA file name
% logFile: name of logFile
% -----------------
% OUTPUT ARGUMENTS:
% -----------------
% mmFile: Output memory map file name ( it is the same name as the fasta
% file with an .mm extension 
% -------------------------------------------------------------------------

[fullPath, filename, extension] = fileparts(fastaFile);
mmFile = fullfile(fullPath, [filename,extension,'.mm']); % Output file name
if exist(mmFile,'file')
    return;
end
tmp_mmFile = tempname(); % Write the mm file to the temp directory
fidIn = fopen(fastaFile,'r');
header = fgetl(fidIn); %#ok<NASGU>
fidOut = fopen(tmp_mmFile,'w');

% ----------------------------------------------
% Read the FASTA file in blocks of 50MB ...
% Remove new line characters
% Convert to uint8
% Write to the MM file.
% ----------------------------------------------
newLine = sprintf('\n');
blockSize = 50 * 2^20;
writeLogFile(logFile,sprintf('Converting fasta to memory mappable file %s\n',fastaFile));
while ~feof(fidIn)
    % Read in the data
    charData = fread(fidIn,blockSize,'*char')';
    % Remove new lines
    charData = strrep(charData,newLine,'');
    % Convert to integers
    intData = nt2int(charData);
    % Write to the new file
    fwrite(fidOut,intData,'uint8');
    writeLogFile(logFile,sprintf('\tDone converting %d MB of data\n',blockSize/(2^20)));
end

% Close the files.
fclose(fidIn);
fclose(fidOut);

% Move the tmp_mmFile to mmFile is required
if ~exist(mmFile,'file')
    writeLogFile(logFile,sprintf('\tMoving from TMP directory\n'));
    [mvStatus,mvMsg,mvMsgId] = movefile(tmp_mmFile,mmFile,'f');
    assert(mvStatus,mvMsgId,'Unable to move mmFile %s to %s: %s',tmp_mmFile,mmFile,mvMsg);
end
end