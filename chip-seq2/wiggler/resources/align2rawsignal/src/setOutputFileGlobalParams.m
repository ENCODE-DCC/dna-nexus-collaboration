function [globalParameters] = setOutputFileGlobalParams(iParams,alignData)
% Sets global output file parameters
% function [globalParameters] = setOutputFileGlobalParams(iParams,alignData)
% --------------------------------------------------------------------------------------------------
% INPUT ARGUMENTS
% ----------------
% iParams<struct>
%     .logFile<string> : log file
%     .alignFname{<string>}: alignment file names
%     .seqDir<string>: sequence directory
%     .uMapDir<string>: mappability tracks directory
%     .outFile.name<string>: name of output file
%     .outFile.format<string>: format of output file
%     .outFile.lcmFile<string>: local cumulative mappability file
%     .normFlag<double>: normalization flag
%     .smooth.fragLen[double]: fragment length
%     .smooth.winLen<double>: smoothing window size
%     .smooth.kernel<string>: smoothing kernel
%     .mapFilter<double>: local cumulative mappability filter
%     .maxMem<double>: approximate maximum memory
%     .outChunk<double>: output chunk length
%     .processChunk<double>: processing chunk length
%
% alignData[<struct>]: array of structures. Each structure represents an Align file
%   .fName<string>: name of tagAlign/BAM file
%   .format<string>: format of alignment file (tagalign/bam)
%   .chrNames{string}: cell array containing names of chrosomosomes
%   .readStart{[double]}: cell array of vectors of read start positions for each chromosome
%        reads are shifted fragLen/2 in the 5' to 3' direction in a strand dependent manner               
%   .maxReadLen<double>: maximum read length
%   .minReadLen<double>: minimum read length
%   .modeReadLen<double>: mode of read length
%   .fragLen<double>: tag extension length
%   .nReads[double]: array containing number of reads per chromosome
%   .mgSize<double>: mappable genome size
% ----------------
% OUTPUT ARGUMENT
% ----------------
% globalParameters<struct>
%   .oFname<string>: output file name
%   .oFormat<string>: output file format (wig,mat,bed)
%   .oType{<string>}: type of output files e.g. signal, lcmap, peak etc. 
%                     is used as a prefix to name the variables in the mat files. Hence the order is
%                     very important. It must the same as the order of the variables that are passed
%                     to writeMatFile(...,varargin)  
%   .oAppend<logical>: flag indicating if output file should be appended
%   .oChunkLen<double>: chunking size for .mat file output
%   .iFname{<string>}: input file names used
%   .iFormat{<string>}: input file formats
%   .iType<string>: type of input files e.g. align, signal, peak etc.
%   .program<string>: program used to create the output file
% --------------------------------------------------------------------------------------------------

globalParameters.oFname = iParams.outFile.name;
globalParameters.oFormat = iParams.outFile.format;
if strcmp( iParams.outFile.name , iParams.outFile.lcmFile )
    globalParameters.oType = {'signal','lcm'};
else
    globalParameters.oType = {'signal'};
end
globalParameters.oAppend = 1;
globalParameters.oChunkLen = iParams.outChunk;
globalParameters.iFname = iParams.alignFname;
globalParameters.iFormat = {alignData.format};
globalParameters.iType = 'align';
globalParameters.program = 'align2rawsignal';

end
