function [bamData] = readBAM(bamFname,chrNames,logFile)
% Reads a BAM file
% function [bamData] = readBAM(bamFname,chrNames,logFile)
% --------------------------------------------------------------------------------------------------
% INPUT ARGUMENTS
% ---------------
% bamFname<string>: BAM file name
% chrNames{<string>}: cell array of chromosome names
% logFile<string>: name of log file
% ---------------
% OUTPUT ARGUMENT
% ---------------
% bamData[struct]: array of structures, each element represents a chromosome
%       .chrName<string>: name of the chromosome
%       .readStart[double]: start position of tag (if +ve then + strand else - strand. For - strand
%                           tags this is the negative of the STOP value in the tagalign file)
%       .readSeqLen[uint8]: lengths of reads
% -------
% NOTES:
% -------
% AVOID CONFUSION ABOUT READ STARTS
% readBAM() will read in a BAM file. The read start positions returned are the 5'-end position
% of the reads w.r.t. the + strand. So for a read on the + strand, the readStart value returned is
% the <chromStart> value in the BAM file. For a read on the - strand, the readStart value returned
% is -<chromEnd>. NOTE the -ve sign.
% --------------------------------------------------------------------------------------------------
% BAM FILE FORMAT
% Field2 : name of chromosome
% Field3 : start position (1-based)
% Field10: read sequence
% --------------------------------------------------------------------------------------------------
[tmp_var.path , tmp_var.name , tmp_var.ext ] = fileparts(bamFname);

% Create temporary directory
tmpDir = tempname();
while exist(tmpDir , 'dir')
    tmpDir = tempname();
end
[st,msg,msgId] = mkdir(tmpDir); % Create temporary directory
assert( st , msgId , 'Unable to create temporary directory %s : %s', tmpDir, msg);
pause(2);

tmp_BamFname = fullfile( tmpDir , [tmp_var.name,tmp_var.ext] );

% -----------------------------------------------------------
% Use Samtools to extract relevant information
% -----------------------------------------------------------
% FORMAT: chr[tab]+/-start[tab]readLen
tmp_BamExtractFname = [tmp_BamFname , '.ext']; 
tmp_BamExtractErrorSamtoolsFname = [tmp_BamFname , '.sterr']; % error file for samtools
tmp_BamExtractErrorAwkFname = [tmp_BamFname , '.awkerr']; % error file for awk

[samStatus,samPath] = system('which samtools'); % get path to samtools
samPath = strtrim(samPath); % remove trailing new lines or blank characters
assert( (samStatus==0) , 'Samtools not found');

writeLogFile( logFile , sprintf('Extracting relevant fields from BAM file %s\n',bamFname) );
if regexp( tmp_BamFname , '\.gz$' )
    samCommand = sprintf( ...
        'zcat %s | %s view -F 0x0204 - 2> %s | awk ''{if (and($2,16) > 0) {print $3,-($4-1+length($10)),length($10)} else {print $3,$4,length($10)}}'' 1> %s 2> %s', ...
        bamFname , samPath , tmp_BamExtractErrorSamtoolsFname , tmp_BamExtractFname , tmp_BamExtractErrorAwkFname );    
else
    samCommand = sprintf( ...
        '%s view -F 0x0204 %s 2> %s | awk ''{if (and($2,16) > 0) {print $3,-($4-1+length($10)),length($10)} else {print $3,$4,length($10)}}'' 1> %s 2> %s', ...
        samPath , bamFname , tmp_BamExtractErrorSamtoolsFname , tmp_BamExtractFname , tmp_BamExtractErrorAwkFname);
end
[samCommandStatus, samCommandResult] = system(samCommand);
assert( (samCommandStatus==0) , 'Unable to extract data from BAM file: %s' , samCommandResult);
pause(2);

if exist( tmp_BamExtractErrorSamtoolsFname , 'file' )
    [~,errResult] = system( sprintf( 'cat %s 2> /dev/null' , tmp_BamExtractErrorSamtoolsFname ) );    
    assert( isempty( strtrim(errResult) ) , ...
        'Samtools error while extracting reads from BAM file %s : %s', bamFname , errResult);
end

if exist(tmp_BamExtractErrorAwkFname,'file')
    [~,errResult] = system( sprintf( 'cat %s 2> /dev/null' , tmp_BamExtractErrorAwkFname ) );    
    assert( isempty( strtrim(errResult) ) , ...
        'Awk error while extracting reads from BAM file %s : %s' , bamFname , errResult);
end
writeLogFile( logFile , sprintf('Done extracting relevant fields  ...\n') );

% -----------------------------------------------------------
% Initialize Output Data
% -----------------------------------------------------------
nChr = length(chrNames); % number of chromosomes
bamData = struct( ...
    'chrName' , chrNames , ...
    'readStart', [] , ...
    'readSeqLen' , [] ); % output data structure

% -----------------------------------------------------------
% Read BAM extract File in chunks
% -----------------------------------------------------------
chunkSize = 1e6;
totReads = 0;
writeLogFile( logFile , sprintf('Reading BAM file %s\n',bamFname) );
fp = fopen( tmp_BamExtractFname , 'r' );
assert( (fp ~= -1) , 'Unable to read extracted data from BAM file %s' , bamFname );
while ~feof(fp)
    bamRawData = textscan( fp , '%s%n%u8' , chunkSize ); % chrName, (+/-)readStart, readLen
    
    % ----------------------------------
    % Filter out illegal read Starts i.e. (-readLen < readStart <= 0)
    % ----------------------------------
    rmIdx = (bamRawData{2} <= 0) & (bamRawData{3} > -bamRawData{2}) ;    
    if nnz(rmIdx)
        bamRawData{1}(rmIdx) = [];
        bamRawData{2}(rmIdx) = [];
        bamRawData{3}(rmIdx) = [];
    end
    clear rmIdx;
    nReads = numel(bamRawData{1});
    if (nReads==0)
        continue;
    else
        totReads = totReads + nReads; % Update total number of reads
        writeLogFile(logFile,sprintf('%g..',totReads));
        % -------------------------------------
        % Process chromosome id and group reads
        % -------------------------------------
        [~,grpId] = ismember(bamRawData{1},chrNames);
        for c = 1:nChr
            % read indices that belong to chromosome 'c'
            tmp_idx = (grpId==c); 
            % Update bamData(c)
            bamData(c).readStart = [ bamData(c).readStart ; bamRawData{2}(tmp_idx) ];
            bamData(c).readSeqLen = [ bamData(c).readSeqLen ; bamRawData{3}(tmp_idx) ];
        end
    end    
end
writeLogFile( logFile , sprintf('\nDone reading BAM file\n') );
fclose(fp);

rmdir(tmpDir,'s'); % Remove temporary directory

assert( (totReads ~= 0) ,'BAM file %s seems to have no valid tags' , bamFname );

end
