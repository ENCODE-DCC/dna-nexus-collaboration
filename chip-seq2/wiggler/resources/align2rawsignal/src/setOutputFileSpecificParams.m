function [specificParameters] = setOutputFileSpecificParams(iParams,chr,alignData)
% Sets specific output file parameters
% function [specificParameters] = setOutputFileSpecificParams(iParams,chr,alignData)
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
% specificParamets<struct>
%      .nReps<double>: number of replicates
%      .chrNames{<string>}: names of chromosomes
%      .chrLen[<double>]: lengths of chromosomes
%      .nReadsPerChrPerRep[<double>,<double>]: reads per chromosome per replicate
%      .nReadsPerChr[<double>]: reads per chrosomome
%      .nReadsPerRep[<double>]: total reads per replicate
%      .nTotalReads<double>: total number of reads
%      .mapGenomeSize<double>: mappable size of the genome
%      .modeReadLen[<double>]: mode of read lengths in each file
%      .seqDir<string>: chromosome sequence directory
%      .uMapdir<string>: uniqueness map directory
%      .fragLen[<double>]: tag extension lengths
%      .winLen<double>: smoothing window size
%      .kernel<string>: smoothing kernel
%      .normFlag<double>: normalization flag
%      .mapFilter<double>: mappability filter
%      .mergeSignalWithMaxTags<logical>: merge maxtags and signal flag
% --------------------------------------------------------------------------------------------------

specificParameters.nReps = numel( iParams.alignFname );
specificParameters.chrNames = chr.names; 
specificParameters.chrLen = chr.len;
specificParameters.nReadsPerChrPerRep = [alignData.nReads]; 
specificParameters.nReadsPerChr = sum(specificParameters.nReadsPerChrPerRep,2);
specificParameters.nReadsPerRep = sum(specificParameters.nReadsPerChrPerRep,1);
specificParameters.nTotalReads = sum(specificParameters.nReadsPerRep);
specificParameters.mapGenomeSize = [alignData.mgSize];
specificParameters.modeReadLen = [alignData.modeReadLen]; 
specificParameters.seqDir = iParams.seqDir; 
specificParameters.uMapdir = iParams.uMapDir; 
specificParameters.fragLen = iParams.smooth.fragLen; 
specificParameters.winLen = iParams.smooth.winLen;
specificParameters.kernel = iParams.smooth.kernel;
specificParameters.normFlag = iParams.normFlag; 
specificParameters.mapFilter = iParams.mapFilter;
specificParameters.mergeSignalWithMaxTags = strcmpi( iParams.outFile.name , iParams.outFile.lcmFile);

end
