function [] = align2rawsignal(varargin)
%  Converts tagAlign/BAM files into normalized signal
% --------------------------------------------------------------------------------------------------
% USAGE:
% -------
% align2rawsignal -i=<alignFname> -s=<seqDir> <OPTIONAL ARGUMENTS>
%
% --help/-h print usage information and quit
% ----------------
% INPUT OPTIONS:
% ----------------
% -i=<alignFname> (MANDATORY, MULTIPLE ALLOWED)
% One or more tagAlign/BAM files (replicates) as input. 
% The tagAlign files can have extensions (gz,tagAlign). BAM files MUST have extension (bam,bam.gz)
%
% -s=<seqDir> (MANDATORY)
% Full path to directory containing chromosome fasta files (eg. chr1.fa ..)
% The file names MUST match the chromosome names used in the tagAlign files
% e.g: /seq/hg19
%
% -u=<uMapDir> (MANDATORY)
% Full path to directory containing binary mappability tracks. 
% The directory name must be of the form [PATH]/globalmap_k<min>tok<max> 
% e.g: /umap/hg19/globalmap_k20tok54
%
% ----------------
% OUTPUT OPTIONS:
% ----------------
% -o=<oFname> (OPTIONAL)
% Full path and name of output signal file.
% Set to stdout if you want to print to stdout which is also the default behavior 
% Default: stdout
%
% -of=<outputFormat> (OPTIONAL)
% Output signal file format
% wiggle (wig) or bedGraph (bg) or matfile (mat)
% Default: mat
%
% -m=<localCumMapFile> (OPTIONAL)
% Calculate, for each position 'i' in the genome, the maximum number of uniquely mappable
% surrounding positions that contribute to the signal value at position 'i'.
% This is a function of the mappability of the surrounding positions, 
% smoothing length and number of replicates.
% The local cumulative mappability is output in <maxTagsFile>.
% If -of=mat then, <localCumMapFile> can have the same name as <oFname>.
% In this case, the local cummap output is stored as a separate set of variables with prefix 'lcm'
% in the .mat file.
% Default: local cumMap is not output to a file
%
% -v=<logFile> (OPTIONAL)
% verbose mode.
% <logFile> Full path and name of file for logging.
% Set to stdout or stderr if you want to output logging info to stdout/stderr 
% Default: off
%
% -n=<normFLag> (OPTIONAL)
% a flag indicating whether the signal output should be normalized
% <normalization_flag> = 0,1,2,3,4,5
% Default: 5 (fold change wrt. expected value from a uniform distribution of reads)
% 0: no normalization
% 1: normSignal(i) = signal(i) * (1e9 / #reads)
% 2: normSignal(i) = (signal(i)/winsize) * (1e9 / #reads)
% 3: normSignal(i) = (signal(i)/localCumMap(i)) * (1e9 / #reads)
% 4: normSignal(i) = (signal(i)/winsize) * (#total_mappable_bases / #reads)
% 5: normSignal(i) = (signal(i)/localCumMap(i)) * (#total_mappable_bases / #reads)
%
% ----------------
% PARAMETERS:
% ----------------
% -l=<fragLen> (OPTIONAL, MULTIPLE ALLOWED)
% Fragment-length / 2*Tag-shift
% Default: 1 (no extension)
% Tags are shifted by floor(fragLen/2) in a 3' direction, relative to the strand of the tag
% NOTE: If a single fragLen is specified then it is applied to all tagAlign/BAM files.
% NOTE: Multiple arguments of this type are allowed. In such a case, 
%       number of fragLen arguments MUST BE == no. of Align files.
%       The ORDER of these arguments is important. 
%       e.g. The first <fragLen> is matched with the first tagAlign/BAM file and so on.
%
% -w=<smoothingWindow> (OPTIONAL)
% Smoothing window size for signal
% Default: mean(1.5*fragLen)
%
% -k=<smoothingKernel> (OPTIONAL)
% Smoothing kernel to use 
% Valid kernels (rectangular,triangular,epanechnikov,biweight,triweight,cosine,gaussian,tukey)
% Default: tukey (with taper ratio of max( 0.25 , min (0.5,max(w-mean(l),0)/(2*w)) ) )
%
% -f=<localCumMapFilter> (OPTIONAL)
% Will nullify positions that have localCumMap <= <mappability_threshold>
% <mappability_threshold> <= 1 implies threshold is in terms of percentage localCumMap
%           e.g. 0.1 means atleast 10% of the positions in the smoothing window MUST be mappable 
%                > 1 implies threshold is on actual localCumMap values (per replicate)
%           e.g. 30 means atleast 30 positions (per replicate) in the extension/smoothing window 
%                MUST be mappable 
% Default: 0.25
%
% -mm=<memory> (OPTIONAL)
% Total memory to use in GB
% Default: 2
% --------------------------------------------------------------------------------------------------

allArgs = varargin;

try
    iParams = processInputArgumentsDeployed( allArgs );    
catch ME
    if strcmp(ME.identifier,'normalExit:normalExit')
        return;
    else
        rethrow(ME);
    end    
end

%% Print input parameters
printIParams( iParams );

%% Run align2rawsignal
align2rawsignalMain( iParams );

end

% #############################################################################################
% AUXILIARY FUNCTIONS
% #############################################################################################

function iParams = processInputArgumentsDeployed( allArgs )
% Process all input arguments

% --------------
% Generate help output
% --------------
helpLine = getUsageHelp();
minArgin = 3; % minimum number of mandatory arguments
if ( numel(allArgs) < minArgin ) || any( ismember( {'--help','-h'} , allArgs ) )
    fprintf( 2 , helpLine );
    error('normalExit:normalExit','This is a normal exit');
end

% --------------
% ARGUMENTS: -v -i -s -u -o -of -m -n -l -w -k -f -mm
% --------------

% --------------
% Verbose flag -v=<logFile>
% default: ''
% --------------
iParams.logFile = processStringArg( allArgs , {'-v='} , '' , false , {} , '' , false );
% Check if logFile directory exists
if ~any( strcmp( iParams.logFile , {'stdout','stderr',''} ) )
    [logDir,~,~] = fileparts(iParams.logFile);
    assert( logical(exist(logDir , 'dir')) || isempty(logDir) , 'Directory for log file output does not exist');
end
% If logfile already exists replace it
if exist( iParams.logFile , 'file' )
    movefile( iParams.logFile , [ iParams.logFile , '.old' ] );    
end

% --------------
% TagAlign/BAM files -i=
% --------------
iParams.alignFname = cellstr( processStringArg( allArgs , {'-i='} , '' , true , {} , 'file' , true ) );
% If any of the align files are BAM then samtools must be in the path
if any( cellfun( @(x)( ~isempty(x) ) , regexpi( iParams.alignFname , '\.bam(\.gz)?$' ) ) )
    [stStatus,~] = system('which samtools');
    assert( ~stStatus, 'Samtools not found. Alignment files are in BAM format. Samtools must be installed and in the path $PATH');
end

% --------------
% SeqDir -s=
% --------------
iParams.seqDir = processStringArg( allArgs , {'-s='} , '' , true , {} , 'dir' , false );

% --------------
% umapDir -u=
% --------------
iParams.uMapDir = processStringArg( allArgs , {'-u='} , '' , true , {} , 'dir' , false );

% --------------
% output File name -o=
% default: stdout
% --------------
iParams.outFile.name = processStringArg( allArgs , {'-o='} , 'stdout' , false , {} , '' , false );
if ~any( strcmp( iParams.outFile.name , {'stdout','stderr',''} ) )
    [outDir,~,~] = fileparts(iParams.outFile.name);
    assert( logical(exist(outDir , 'dir')) || isempty(outDir) , 'Output directory does not exist' );
end
% If output file already exists rename it
if exist( iParams.outFile.name , 'file' )
    movefile( iParams.outFile.name , [iParams.outFile.name,'.old'] );
end

% --------------
% output File format -of=
% default: mat
% --------------
iParams.outFile.format = processStringArg( allArgs , {'-of='} , 'mat' , false , {'mat','bg','wig'} , '' , false );
% if .mat output format then output file cannot be '' or 'stdout' or 'stderr'
assert( ~( ...
    strcmp( iParams.outFile.format , 'mat' ) && ...
    any( strcmp( iParams.outFile.name , {'stdout','stderr',''} ) ) ) , ...
    'If output format is .mat then output file cannot be stdout or stderr');

% --------------
% local cumumlative mappability file -m=
% default: ''
% --------------
iParams.outFile.lcmFile = processStringArg( allArgs , {'-m='} , '' , false , {} , '' , false );
if ~any( strcmp( iParams.outFile.lcmFile , {'stdout','stderr',''} ) )
    [outDir,~,~] = fileparts(iParams.outFile.lcmFile);
    assert( logical(exist(outDir , 'dir')) || isempty(outDir) , 'Output directory for local cumulative mappability file does not exist' );
end
% If lcmFile already exists rename it
if exist( iParams.outFile.lcmFile , 'file' )
    movefile( iParams.outFile.lcmFile , [ iParams.outFile.lcmFile , '.old' ] );    
end
% If lcmFile == outFile.name BUT oFormat is NOT 'mat' then error
assert( ~( ...
    strcmpi( iParams.outFile.name , iParams.outFile.lcmFile ) && ...
    ~strcmpi( iParams.outFile.format , 'mat' ) ) , ...
    'Local CumMap Filename and Signal file name can only be the same if output format is mat');

% --------------
% Normalization flag -n=
% default: 5
% --------------
iParams.normFlag = str2double( processStringArg( allArgs , {'-n='} , '5' , false , {'0','1','2','3','4','5'} , '' , false ) );

% --------------
% fragment length  -l=
% default: 1
% --------------
iParams.smooth.fragLen = ceil( processNumericArg( allArgs , {'-l='} , 1 , false , [1 Inf] , true ) );
assert( ...
    ( numel(iParams.smooth.fragLen) == 1 ) || ...
    ( numel(iParams.alignFname) == numel(iParams.smooth.fragLen) ) , ...
    'Number of fragment length values do not match number of alignment files');
if numel(iParams.smooth.fragLen) == 1
    iParams.smooth.fragLen = repmat(iParams.smooth.fragLen , numel(iParams.alignFname) , 1);
end

% --------------
% Smoothing window  -w=
% default: NaN (will be computed later)
% --------------
iParams.smooth.winLen = ceil( processNumericArg( allArgs , {'-w='} , NaN , false , [1 Inf] , false ) );
if rem(iParams.smooth.winLen,2) == 0
    iParams.smooth.winLen = iParams.smooth.winLen + 1; % make it odd
end

% --------------
% Smoothing kernel  -k=
% default: 'tukey'
% --------------
validKernels = {...
    'rectangular','triangular','epanechnikov','biweight', ...
    'triweight','cosine','gaussian','tukey'}; 
iParams.smooth.kernel = processStringArg( allArgs , {'-k='} , 'tukey' , false , validKernels , '' , false );

% --------------
% Local Mappability filter  -f=
% default: 0.25
% --------------
iParams.mapFilter = processNumericArg( allArgs , {'-f='} , 0.25 , false , [0 Inf] , false );

% --------------
% Approximate memory usage  -mm=
% default: 2
% --------------
iParams.maxMem = processNumericArg( allArgs , {'-mm='} , 2 , false , [2 Inf] , false );

% ----------------
% Chunking parameters
% ----------------
iParams.outChunk = 1e6;
iParams.processChunk = iParams.outChunk * round( 25 * iParams.maxMem / 4 );

end

% ===============================================================================================

function [value] = processNumericArg( allArgs , attributeFlags , defaultVal , mandatory , valRange , multiArgs )
% Will process numeric arguments
% ------------------------------------------------------------------------------------------
% INPUT ARGUMENTS
% -----------------
% attributeFlags{<string>}/<string>: flag preceeding argument value eg. {'--mem='}
% defaultVal<double>: default value to use
% mandatory[0/1]: set to 1 if the field is mandatory
% valRange[double,double]: [min max] legal values
% multiArgs[0/1]: if set to 1 more than 1 instance of the attribute value pair is allowed
% -----------------
% OUTPUT ARGUMENTS
% -----------------
% value[double]: value of the numeric argument
% ------------------------------------------------------------------------------------------

attributeRegExp = strcat( '^' , cellstr(attributeFlags) , '|' ); % Convert cell array of attribute Flags into a combined regular expression
attributeRegExp = regexprep( [ attributeRegExp{:} ] , '\|$' , '' ); % Remove extra trailing |

matchIdx = find( cellfun( @length, regexp( allArgs , attributeRegExp ) ) ); % find index of attribute in allArgs
matchIdx = sort(matchIdx); % sort so that indices are in left to right order
nMatchIdx = numel(matchIdx); % number of matches

if (nMatchIdx==0)
    
    if mandatory % if mandatory throw error else set to default value
        error( 'ERROR: Missing mandatory argument %s\n' , attributeFlags{1} );
    else
        value = defaultVal;
    end
    
elseif ( ~multiArgs && (nMatchIdx > 1) )
    
    error( 'ERROR: Too many arguments of type %s\n', attributeFlags{1} );
    
else
    
    value = regexprep( allArgs(matchIdx) , attributeRegExp , '' ); % remove attributeFlag part
    value = regexprep( value , '^"|"$' , '' ); % remove leading and lagging quotes
    value = str2double(value); % convert to double
    % throw error if illegal value
    assert( all( ~isnan(value) ) , 'ERROR: Illegal value for attribute %s: Value MUST be numeric\n', attributeFlags{1} );
    % Check for illegal range of value
    if ~isempty(valRange)
        assert( (all( value >= valRange(1) ) && all( value <= valRange(2) ) ), 'ERROR: Illegal value for attribute %s: Legal range is [%f %f]\n', attributeFlags{1} , valRange(1) , valRange(2) );
    end
    
end

end

% ===============================================================================================

function [value] = processStringArg( allArgs, attributeFlags , defaultVal , mandatory , valRange , existCond , multiArgs )
% Will process string arguments
% ------------------------------------------------------------------------------------------
% INPUT ARGUMENTS
% -----------------
% attributeFlags{<string>}/<string>: flag preceeding argument value e.g. {'--mem='}
% defaultVal<string>: default value to use
% mandatory[0/1]: set to 1 if the field is mandatory
% valRange{}: cell arrray of legal string values. Set to {} if not used
% existCond<string>: (OPTIONAL) use the exist() function with the flag existCond. Useful for
%            checking if a file or directory exists. Set to '' if not used
% multiArgs[0/1]: if set to 1 more than 1 instance of the attribute value pair is allowed
% -----------------
% OUTPUT ARGUMENTS
% -----------------
% value{string}/<string>: value of the string argument
% ------------------------------------------------------------------------------------------

attributeRegExp = strcat( '^' , cellstr(attributeFlags) , '|' ); % Convert cell array of attribute Flags into a combined regular expression
attributeRegExp = regexprep( [attributeRegExp{:}] , '\|$' , '' ); % Remove extra trailing |

matchIdx = find( cellfun( @length , regexp( allArgs , attributeRegExp ) ) ); % find index of attribute in allArgs
matchIdx = sort(matchIdx); % sort so that indices are in left to right order
nMatchIdx = numel(matchIdx); % number of matches

if (nMatchIdx==0)
    
    if mandatory % if mandatory throw error else set to default value
        error( 'ERROR: Missing mandatory argument %s\n', attributeFlags{1} );
    else
        value = defaultVal;
    end
    
elseif ( ~multiArgs && (nMatchIdx > 1) )
    
    error( 'ERROR: Too many arguments of type %s\n', attributeFlags{1} );
    
else
    
    value = regexprep( allArgs(matchIdx) , attributeRegExp , '' ); % remove attributeFlag part
    value = regexprep( value , '^"|"$' , '' ); % remove leading and lagging quotes
    % throw error if illegal value
    if ~isempty(valRange)
        assert( all( ismember(value,valRange) ) , 'ERROR: Illegal value for attribute %s: Check help for valid values\n', attributeFlags{1} );
    end
    % check any exist() conditions if applicable
    if ~isempty(existCond)
        assert( all ( logical ( cellfun( @(x) exist(x,existCond), value, 'UniformOutput', true ) ) ) , ...
            'ERROR: Illegal value for attribute %s: %s does not exist\n', attributeFlags{1} , existCond );
    end
    
end

% Convert single cell into a string
if ( iscell(value) && ( numel(value)==1 ) )
    value = [ value{:} ];
end

end

% ===============================================================================================

function [] = printIParams(iParams)
% Prints run parameters
% function [] = printIParams(iParams)
% --------------------------------------------------------------------------------------------------
% INPUT ARGUMENTS
% ------------------
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
% --------------------------------------------------------------------------------------------------

writeLogFile(iParams.logFile,date());
writeLogFile(iParams.logFile, sprintf('\n============\nRun parameters\n============\n'));

for i = 1 : numel(iParams.alignFname)
    writeLogFile(iParams.logFile, sprintf('Alignment filename %d=%s\tTag shift=%d\n', ...
        i , iParams.alignFname{i} , floor(iParams.smooth.fragLen(i)/2) ) );
end

writeLogFile(iParams.logFile, sprintf('Chromosome sequence directory=%s\n',iParams.seqDir));
writeLogFile(iParams.logFile, sprintf('Uniqueness map directory=%s\n',iParams.uMapDir));
writeLogFile(iParams.logFile, sprintf('Signal filename=%s\n',iParams.outFile.name));
writeLogFile(iParams.logFile, sprintf('Signal file format=%s\n',iParams.outFile.format));
writeLogFile(iParams.logFile, sprintf('Logfile=%s\n',iParams.logFile));
writeLogFile(iParams.logFile, sprintf('Max Tags filename=%s\n',iParams.outFile.lcmFile));
writeLogFile(iParams.logFile, sprintf('Normalize Flag = %d\n',iParams.normFlag));
writeLogFile(iParams.logFile, sprintf('Mappability filter = %g\n',iParams.mapFilter));
writeLogFile(iParams.logFile, sprintf('Chunk Length = %g\n',iParams.outChunk));
writeLogFile(iParams.logFile, sprintf('Smoothing window = %g\n',iParams.smooth.winLen));
writeLogFile(iParams.logFile, sprintf('Smoothing kernel = %s\n\n',iParams.smooth.kernel));

end

% ===============================================================================================

function helpLine = getUsageHelp()
%GETUSAGEHELP generates help/usage information
% function helpLine = getUsageHelp()
% hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = 0;
hlnum = hlnum + 1; helpLine{hlnum} = '--------------------------------------------------------------\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Program: align2rawsignal (Converts tagAlign/BAM files into normalized signal)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Version: 2.0\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Contact: Anshul Kundaje (akundaje@stanford.edu)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '--------------------------------------------------------------\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = 'USAGE: align2rawsignal -i=<alignFname> -s=<seqDir> <OPTIONAL ARGUMENTS>\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '--help/-h print usage information and quit\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '----------------\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'INPUT OPTIONS:\n';
hlnum = hlnum + 1; helpLine{hlnum} = '----------------\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-i=<alignFname> (MANDATORY, MULTIPLE ALLOWED)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'One or more tagAlign/BAM files (replicates) as input.\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'The tagAlign files can have extensions (gz,tagAlign). BAM files MUST have extension (bam,bam.gz)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-s=<seqDir> (MANDATORY)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Full path to directory containing chromosome fasta files (eg. chr1.fa ..)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'The file names MUST match the chromosome names used in the tagAlign files\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'e.g: /seq/hg19\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-u=<uMapDir> (MANDATORY)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Full path to directory containing binary mappability tracks.\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'The directory name must be of the form [PATH]/globalmap_k<min>tok<max>\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'e.g: /umap/hg19/globalmap_k20tok54\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '----------------\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'OUTPUT OPTIONS:\n';
hlnum = hlnum + 1; helpLine{hlnum} = '----------------\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-o=<oFname> (OPTIONAL)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Full path and name of output signal file.\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Set to stdout if you want to print to stdout which is also the default behavior\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Default: stdout\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-of=<outputFormat> (OPTIONAL)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Output signal file format\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'wiggle (wig) or bedGraph (bg) or matfile (mat)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Default: mat\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-m=<localCumMapFile> (OPTIONAL)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Calculate, for each position ''i'' in the genome, the maximum number of uniquely mappable\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'surrounding positions that contribute to the signal value at position ''i''.\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'This is a function of the mappability of the surrounding positions, \n';
hlnum = hlnum + 1; helpLine{hlnum} = 'tag extension/smoothing length and number of replicates.\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'The local cumulative mappability is output in <maxTagsFile>.\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'If -of=mat then, <localCumMapFile> can have the same name as <oFname>.\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'In this case, the local cummap output is stored as a separate set of variables with prefix maxTags\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'in the .mat file.\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Default: local cumMap is not output to a file\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-v=<logFile> (OPTIONAL)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'verbose mode.\n';
hlnum = hlnum + 1; helpLine{hlnum} = '<logFile> Full path and name of file for logging.\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Set to stdout/stderr if you want to output logging info to stdout/stderr \n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Default: off\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-n=<normFLag> (OPTIONAL)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'a flag indicating whether the signal output should be normalized\n';
hlnum = hlnum + 1; helpLine{hlnum} = '<normalization_flag> = 0,1,2,3,4,5\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Default: 5 (fold change wrt. expected value from a uniform distribution of reads)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '0: no normalization\n';
hlnum = hlnum + 1; helpLine{hlnum} = '1: normSignal(i) = signal(i) * (1e9 / #reads)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '2: normSignal(i) = (signal(i)/winsize) * (1e9 / #reads)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '3: normSignal(i) = (signal(i)/localCumMap(i)) * (1e9 / #reads)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '4: normSignal(i) = (signal(i)/winsize) * (#total_mappable_bases / #reads)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '5: normSignal(i) = (signal(i)/localCumMap(i)) * (#total_mappable_bases / #reads)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '----------------\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'PARAMETERS:\n';
hlnum = hlnum + 1; helpLine{hlnum} = '----------------\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-l=<fragLen> (OPTIONAL, MULTIPLE ALLOWED)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Fragment-length / 2*Tag-shift\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Default: 1 (no extension)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Tags are shifted by floor(fragLen/2) in a 3'' direction, relative to the strand of the tag\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'NOTE: If a single fragLen is specified then it is applied to all tagAlign/BAM files.\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'NOTE: Multiple arguments of this type are allowed. In such a case, \n';
hlnum = hlnum + 1; helpLine{hlnum} = '      number of fragLen arguments MUST BE == no. of Align files.\n';
hlnum = hlnum + 1; helpLine{hlnum} = '      The ORDER of these arguments is important. \n';
hlnum = hlnum + 1; helpLine{hlnum} = '      e.g. The first <fragLen> is matched with the first tagAlign/BAM file and so on.\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-w=<smoothingWindow> (OPTIONAL)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Smoothing window size for signal\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Default: mean(1.5*fragLen)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-k=<smoothingKernel> (OPTIONAL)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Smoothing kernel to use \n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Valid kernels (rectangular,triangular,epanechnikov,biweight,triweight,cosine,gaussian,tukey)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Default: tukey (with taper ratio of max( 0.25 , min (0.5,max(w-mean(l),0)/(2*w)) ) if w is specified)\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-f=<localCumMapFilter> (OPTIONAL)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Will nullify positions that have localCumMap <= <mappability_threshold>\n';
hlnum = hlnum + 1; helpLine{hlnum} = '<mappability_threshold> <= 1 implies threshold is in terms of percentage localCumMap\n';
hlnum = hlnum + 1; helpLine{hlnum} = '          e.g. 0.1 means atleast 10 percent of the positions in the smoothing window MUST be mappable \n';
hlnum = hlnum + 1; helpLine{hlnum} = '               > 1 implies threshold is on actual maxtags values\n';
hlnum = hlnum + 1; helpLine{hlnum} = '          e.g. 30 means atleast 30 positions (per replicate) in the extension/smoothing window \n';
hlnum = hlnum + 1; helpLine{hlnum} = '               MUST be mappable \n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Default: 0.25\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

hlnum = hlnum + 1; helpLine{hlnum} = '-mm=<memory> (OPTIONAL)\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Total memory to use in GB\n';
hlnum = hlnum + 1; helpLine{hlnum} = 'Default: 2\n';
hlnum = hlnum + 1; helpLine{hlnum} = '--------------------------------------------------------------------------------------------------\n';
hlnum = hlnum + 1; helpLine{hlnum} = '\n';

helpLine = cell2mat(helpLine);
end

% ===============================================================================================