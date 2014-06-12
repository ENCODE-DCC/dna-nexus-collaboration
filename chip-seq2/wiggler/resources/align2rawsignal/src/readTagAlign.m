function [taData] = readTagAlign(taFname,chrNames,logFile)
% Reads a tagAlign file
% function [taData] = readTagAlign(taFname,chrNames,logFile)
% --------------------------------------------------------------------------------------------------
% INPUT ARGUMENTS
% ----------------
% taFname: tagAlign File name
% chrNames{string}: names of all chromosomes (should be the same naming convention as the tagAlign
%                   file) 
% logFile: path and name of log file
% ----------------
% OUTPUT ARGUMENT
% ----------------
% taData[struct]: array of structures, each element represents a chromosome
%       .chrName: name of the chromosome
%       .readStart[double]: start position of tag (if +ve then + strand else - strand. For - strand
%                           tags this is the negative of the STOP value in the tagalign file)
%       .readSeqLen[uint8]: lengths of reads
% ----------------
% NOTES:
% ----------------
% AVOID CONFUSION ABOUT READ STARTS
% readTagAlign() will read in a tagAlign file. The read start positions returned are the 5'-end
% position of the reads w.r.t. the + strand. So for a read on the + strand, the readStart value
% returned is the <chromStart> value in the tagAlign file. For a read on the - strand, the readStart
% value returned is -<chromEnd>. Note the -ve sign.
% --------------------------------------------------------------------------------------------------
% TAGALIGN FILE FORMAT
% (1) chrom 	string 	Name of the chromosome
% (2) chromStart 	int 	Starting position on the chromosome (0-based)
% (3) chromEnd 	int 	Ending position on the chromosome (1-based)
% (4) sequence 	string 	Sequence of this read
% (5) score 	int 	Indicates uniqueness or quality (preferably 1000/alignmentCount).
% (6) strand 	char 	Orientation of this read (+ or -)
% --------------------------------------------------------------------------------------------------

[tmp_var.path,tmp_var.name,tmp_var.ext] = fileparts(taFname);

if regexp( taFname , '\.gz$' )
    isGz = 1;
    tmpDir = tempname();
    while exist(tmpDir , 'dir')
        tmpDir = tempname();
    end
    [st,msg,msgId] = mkdir(tmpDir); % Create temporary directory
    assert( st , msgId , 'Unable to create temporary directory %s : %s', tmpDir, msg);
    pause(2);    
    tmp_taFname = fullfile( tmpDir , [tmp_var.name,tmp_var.ext] );    
    tmp_taExtractFname = [tmp_taFname , '.ext']; % unzipped tagAlign file
    tmp_taExtractErrorFname = [tmp_taFname , '.exterr']; % unzip error file
    writeLogFile( logFile , sprintf('Decompressing TagAlign file %s\n',taFname) );
    unzipCommand = sprintf( 'zcat %s 1> %s 2> %s' , taFname, tmp_taExtractFname , tmp_taExtractErrorFname );
    [gzStatus,gzResult] = system(unzipCommand);
    pause(2);
    assert( (gzStatus==0) , 'Error unzipping tagAlign File %s: %s' , taFname , gzResult );
    if exist( tmp_taExtractErrorFname , 'file' )
        [~,errResult] = system( sprintf( 'cat %s 2> /dev/null' , tmp_taExtractErrorFname ) );        
        assert( isempty( strtrim(errResult) ) , ...
            'Error while decompressing tagAlign file %s : %s', taFname , errResult);        
    end
else
    isGz = 0;
    tmp_taExtractFname = taFname;
end

% ----------------------------------
% Check that file is readable
% ----------------------------------
[fp,errmsg] = fopen(tmp_taExtractFname,'r');
assert( (fp ~= -1) , 'Error reading tagAlign file %s : %s' , taFname, errmsg );

% ----------------------------------
% Initialize Output Data
% ----------------------------------
nChr = length(chrNames); % number of chromosomes
taData = struct( ...
    'chrName' , chrNames , ...
    'readStart' ,[] , ...
    'readSeqLen' , [] ); % output data structure

% ----------------------------------
% Read in data in chunks and filter
% ----------------------------------
chunkSize = 1e6;
totReads = 0;
writeLogFile( logFile , sprintf( 'Reading TagAlign file %s\n' , taFname) );
while ~feof(fp)
    
    taRawData = textscan(fp,'%s%n%n%*s%*f%c',chunkSize);   
    nReads = length(taRawData{1});
    
    % ----------------------------------
    % Process Reads
    % ----------------------------------
    if (nReads==0)        
        continue; % if no reads pass filter continue        
    else
        % Update total number of reads
        totReads = totReads + nReads; 
        writeLogFile(logFile,sprintf('%g..',totReads));
        
        [readChr, readStart, readStop] = taRawData{1:3}; % First 3 fields chr,star,stop
        readStrand = taRawData{end}; % Last field is strand                        
        clear taRawData;        
        % ----------------------------------
        % Filter out illegal readStart < 0
        % ----------------------------------
        rmIdx = (readStart<0);
        if nnz(rmIdx)
            readChr(rmIdx) = [];
            readStart(rmIdx) = [];
            readStop(rmIdx) = [];            
            readStrand(rmIdx) = [];
        end
        % ----------------------------------
        % Process readLength
        % ----------------------------------
        readSeqLen = uint8(readStop - readStart); % length of reads
        % ----------------------------------
        % Process start,stop and strand
        % ----------------------------------        
        readStrand = (readStrand == '+'); 
        % convert startpos from 0-based to 1-based index        
        readStart(readStrand) = readStart(readStrand)+1;
        % for reads on - strand replace start with -stop
        readStart(~readStrand) = -readStop(~readStrand);
        clear readStrand readStop;               
        % -------------------------------------
        % Process chromosome id and group reads
        % ----------------------------------
        [~,grpId] = ismember(readChr,chrNames);
        for c = 1:nChr
            % read indices that belong to chromosome 'c'
            tmp_idx = (grpId==c); 
            % Update taData(c)
            taData(c).readStart = [taData(c).readStart ; readStart(tmp_idx)];
            taData(c).readSeqLen = [taData(c).readSeqLen ; readSeqLen(tmp_idx)];
        end
    end
end
writeLogFile( logFile , sprintf('\nDone reading TagAlign file\n') );
fclose(fp);

if isGz
    rmdir(tmpDir,'s');
end
assert( (totReads ~= 0) , 'TagAlign file %s seems to have no valid tags' , taFname );

end