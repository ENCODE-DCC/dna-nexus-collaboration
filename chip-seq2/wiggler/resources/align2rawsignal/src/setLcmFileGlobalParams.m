function [globalParameters] = setLcmFileGlobalParams(iParams, globalParameters)
% Sets global maxtags file parameters
% function [globalParameters] = setLcmFileGlobalParams(iParams, globalParameters)
% ------------------------------------------------------------------------------------------
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
%
% globalParameters<struct>
%   .oFname<string>: output file name
%   .oFormat<string>: output file format (wig,mat,bg)
%   .oType{<string>}: type of output files e.g. signal, lcm, peak etc. 
%                     is used as a prefix to name the variables in the mat files. Hence the order is
%                     very important. It must the same as the order of the variables that are passed
%                     to writeMatFile(...,varargin)  
%   .oAppend<logical>: flag indicating if output file should be appended
%   .oChunkLen<double>: chunking size for .mat file output
%   .iFname{<string>}: input file names used
%   .iFormat{<string>}: input file formats
%   .iType<string>: type of input files e.g. align, signal, peak etc.
%   .program<string>: program used to create the output file
% ------------------
% OUTPUT ARGUMENT
% ------------------
% globalParameters<struct>
%   Same as input
% ------------------------------------------------------------------------------------------

% Change relevant parameters
globalParameters.oFname = iParams.outFile.lcmFile;
% Get file extension of maxTagsFname
[~, ~, tmp_ext] = fileparts(iParams.outFile.lcmFile);
switch lower(tmp_ext)
    case {'.mat','.wig'}
        globalParameters.oFormat = lower(tmp_ext(2:end));
    otherwise
        globalParameters.oFormat = 'bg';
end
globalParameters.oType = {'lcm'};

end
