function alignData = initializeAlignData(iParams,chr)
% initializes the alignData datastructure
% alignData = initializeAlignData(iParams,chr)
% --------------------------------------------------------------------------------------------------
% INPUT ARGUMENTS
% -----------------
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
% -----------------
% OUTPUT ARGUMENT
% -----------------
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
% --------------------------------------------------------------------------------------------------

[ alignData( 1 : numel(iParams.alignFname) ).fName ] = iParams.alignFname{:};
[alignData.format] = deal('');
[alignData.chrNames] = deal(chr.names);
[alignData.readStart] = deal(cell(chr.nChr,1));
[alignData.maxReadLen] = deal(0);
[alignData.minReadLen] = deal(0);
[alignData.modeReadLen] = deal(0);
tmp_extLen = num2cell(iParams.smooth.fragLen);
[alignData.fragLen] = tmp_extLen{:};
[alignData.nReads] = deal(zeros(chr.nChr,1));
[alignData.mgSize] = deal(NaN);

end
