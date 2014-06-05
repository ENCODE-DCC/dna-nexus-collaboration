function [] = align2rawsignalMain( iParams )
%  Reads in a set of tagAlign/BAM files and creates a consolidated normalized signal
%  wiggle/bedGraph/mat file 
% --------------------------------------------------------------------------------------------------
% INPUT ARGUMENTS
% ---------------------
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

% --------------------------------------------------------------------------------------------------
%% Set chromosome related parameters
% chr<struct>
%     .dir<string>: directory containing chromosome files
%     .fNames{string}: names with paths of chromosome fasta files
%     .mmNames{string}: names with paths of chromosome memory mapped files
%     .names{string}: names of chromosomes
%     .len[double]: length of chromosomes
%     .nChr<double>: number of chromosomes
% --------------------------------------------------------------------------------------------------
chr = initializeChromosomeParams(iParams);

% --------------------------------------------------------------------------------------------------
%% Set uniqueness map related parameters
% uMap<struct>
%   .minReadLen<double>: min read length covered in uniqueness map
%   .maxReadLen<double>: max read length covered in uniqueness map
%   .mmNames{string}: cell array of names with paths of uniqueness map files
% --------------------------------------------------------------------------------------------------
uMap = initializeUmapParams( iParams , chr );

% --------------------------------------------------------------------------------------------------
%% Initialize alignData datastructure
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
alignData = initializeAlignData( iParams , chr );

% --------------------------------------------------------------------------------------------------
%% Read and preprocess each align file
% --------------------------------------------------------------------------------------------------
for ita = 1 : numel( iParams.alignFname )
    alignData(ita) = preprocessAlignData( alignData(ita) , chr , uMap , iParams.logFile );
    writeLogFile( iParams.logFile , sprintf( 'Mode of read lengths = %d\n' , alignData(ita).modeReadLen ) );
end

% --------------------------------------------------------------------------
%% Get mappable size of the genome (excludes chromosomes that have no reads)
% --------------------------------------------------------------------------
readsPerChr = sum( [alignData.nReads] , 2 ); % number of reads per chromosome
writeLogFile( iParams.logFile , sprintf('\nComputing mappable genome size\n') );
mgSize = getMappableGenomeSize( [alignData.modeReadLen] , uMap , readsPerChr );
mgSize = num2cell(mgSize);
[alignData.mgSize] = mgSize{:};
clear mgSize;
writeLogFile( iParams.logFile , sprintf( 'Mappable genome size = %g bases\n' , [alignData.mgSize] ) );

% --------------------------------------------------------------------------
%% Compute window size if NaN
% --------------------------------------------------------------------------
if isnan( iParams.smooth.winLen )
    iParams.smooth.winLen = round(1.5*mean(iParams.smooth.fragLen));
end
% Make window size odd number
if rem(iParams.smooth.winLen , 2) == 0
    iParams.smooth.winLen = iParams.smooth.winLen + 1;
end
writeLogFile( iParams.logFile , sprintf( 'Smoothing window size = %d\n' , iParams.smooth.winLen ) );

% Compute tukeyRatio max( 0.25 , min( 0.5 , max(w-mean(l),0)/(2*mean(l)) ) )
tukeyRatio = NaN;
if strcmp(iParams.smooth.kernel , 'tukey')
    tukeyRatio = max( 0.25 , ...
                             min( 0.5 , ...
                                        max( 0 , (iParams.smooth.winLen - mean(iParams.smooth.fragLen)) ) ...
                                        / (2 * mean(iParams.smooth.fragLen) ) ...
                                ) ...
                    );
    writeLogFile( iParams.logFile , sprintf( 'Tukey kernel taper ratio = %g\n' , tukeyRatio ) ); 
end

% --------------------------------------------------------------------------
%% Set output file Parameters
% oFileParams<struct>
%  .globalParameters<struct>
%  .specificParamets<struct>
% --------------------------------------------------------------------------------------------------
% Global parameters
% globalParameters<struct>
%       .oFname<string>: output file name
%       .oFormat<string>: output file format (wig,mat,bg)
%       .oType{<string>}: type of output files e.g. signal, maxtags, peak etc.
%                         is used as a prefix to name the variables in the mat files. Hence the
%                         order is very important. It must the same as the order of  the variables 
%                         that are passed to writeMatFile(...,varargin)
%       .oAppend<logical>: flag indicating if output file should be appended
%       .oChunkLen<double>: chunking size for .mat file output
%       .iFname{<string>}: input file names used
%       .iFormat{<string>}: input file formats
%       .iType<string>: type of input files e.g. align, signal, peak etc.
%       .program<string>: program used to create the output file
% --------------------------------------------------------------------------------------------------
% Specific parameters
% specificParameters<struct>
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
oFileParams.globalParameters = setOutputFileGlobalParams(iParams,alignData);
oFileParams.specificParameters = setOutputFileSpecificParams(iParams,chr,alignData);
% Initialize output file if .mat
if strcmpi(oFileParams.globalParameters.oFormat,'mat')
    save(oFileParams.globalParameters.oFname,'-struct','oFileParams');
end

% --------------------------------------------------------------------------------------------------
%% Set lcMap file parameters if we are not merging lcMap and signal
% lcmFileParams<struct>
%  .globalParameters<struct>
%  .specificParamets<struct>
% SAME AS ABOVE
% --------------------------------------------------------------------------------------------------
if ( ~isempty(iParams.outFile.lcmFile) && ~strcmp(iParams.outFile.name , iParams.outFile.lcmFile) )
    
    lcmFileParams.globalParameters = setLcmFileGlobalParams( iParams , oFileParams.globalParameters );
    lcmFileParams.specificParameters = oFileParams.specificParameters;
      
    % Initialize lcmFile if .mat
    if strcmpi(lcmFileParams.globalParameters.oFormat,'mat')
        save( lcmFileParams.globalParameters.oFname , '-struct' , 'lcmFileParams' );
    end    
end

% ==================================================================================================
%% MAIN PROCESSING OF EACH CHROMOSOME (READS TO SIGNAL)
% ==================================================================================================
% *IMPORTANT VARIABLES*
% lcmVal[double]: maxtags values
% signalVal[double]: total signal value after tag extension
% totalReads<double>; total number of reads in all files
% readsPerChr[double]; number of reads per chromosome
% ------------------------------------------------
% Check if all files have same tagext and readlen
% ------------------------------------------------
isUniqueReadLen = ( numel( unique( [alignData.modeReadLen] ) ) == 1 );
isUniqueFraglen = ( numel( unique( [alignData.fragLen] ) ) == 1 );

for iChr = 1:chr.nChr
    
    % number of chunks
    nChunks = ceil( chr.len(iChr) / iParams.processChunk );
    writeLogFile( iParams.logFile , sprintf( '\nProcessing chromosome %s in %d chunks\n' , chr.names{iChr} , nChunks ) );
    
    % ----------------------------------------------------------------------------------------------
    % Concatenate read starts from all align files    
    % catReadStart[double]: concatenated read starts
    % ----------------------------------------------------------------------------------------------
    catReadStart = concatenateReads(iChr);

    % Read current mappability track
    fu = fopen( uMap.mmNames{iChr} , 'r' );
    currUmap = fread( fu , '*uint8' );
    fclose(fu);

    
    % ----------------------------------------------------------------------------------------------
    %% PROCESS EACH CHROMOSOME CHUNK
    % ----------------------------------------------------------------------------------------------
    for iChunk = 1:nChunks
        
        writeLogFile( iParams.logFile , sprintf( '-Chunk %d of %d\n' , iChunk , nChunks ) );
        % ------------------------------------------------------------------------------------------
        % Set chunk boundaries
        % chunk<struct>
        %  .chrName<string> : name of chromosome
        %  .start<double> : start position on chromosome
        %  .stop<double> : stop position on chromosome
        %  .len<double> : length of chunk
        %  .isFirst<logical> : flag indicating if it is the first chunk for the chromosome
        %  .isLast<logical> : flag indicating if this is the last chunk for the chromosome
        % ------------------------------------------------------------------------------------------                
        chunk.chrName = chr.names{iChr};
        chunk.start = ( (iChunk-1) * iParams.processChunk ) + 1;
        chunk.stop = min( iChunk * iParams.processChunk , chr.len(iChr) );
        chunk.len = chunk.stop - chunk.start + 1;
        chunk.isFirst = (iChunk == 1);
        chunk.isLast = (iChunk == nChunks);
        
        % ------------------------------------------------------------------------------------------
        % Calculate max tags and signal values for that chunk
        % lcmVal = NaN if position is unmappable (after extension/smoothing)
        % ------------------------------------------------------------------------------------------
        if (isUniqueReadLen && isUniqueFraglen) 
            % If there is a single value of fragLen and modeReadLen, then maxTags needs to be
            % computed for only one of the align files and signalVal is computed from the
            % concatenated readStarts            
            [signalVal,lcmVal] = catReads2SignalPerChunk( 0.5 );                       
        else
            % If fragLen or modeReadLen are different for the different align files, then maxTags and
            % signalVal need to computed separatedly for each align file            
            [signalVal,lcmVal] = perRepReads2SignalPerChunk( 0.5 );            
        end
               
        % ------------------------------------------------------------------------------------------
        % Post process signal values
        % signalVal = NaN if position is unmappable OR does not pass mappability threshold
        % ------------------------------------------------------------------------------------------
        mapFilterSignal();
        
        % ------------------------------------------------------------------------------------------        
        % lcmVal = 0 if position is mappable after ext/smoothing BUT lies in sequence gaps (Ns) 
        % ------------------------------------------------------------------------------------------
        markSeqGapsInLcmVals();
        
        % ------------------------------------------------------------------------------------------
        % Write signal file
        % ------------------------------------------------------------------------------------------
        writeLogFile( iParams.logFile, sprintf( '--Writing Signal %s file\n' , iParams.outFile.format ) );
        switch iParams.outFile.format
            case 'wig'
                writeWiggle( oFileParams.globalParameters , chunk , signalVal );
            case 'bg'
                writeBedGraph( oFileParams.globalParameters , chunk , signalVal );
            case 'mat'
                % Merge signal output with lcmVals output if required
                if (oFileParams.specificParameters.mergeSignalWithMaxTags)
                    writeMatFile( oFileParams.globalParameters , chunk , signalVal , lcmVal );
                else
                    writeMatFile( oFileParams.globalParameters , chunk , signalVal );
                end
            otherwise
                error('Illegal output format\n');
        end
        
        % ------------------------------------------------------------------------------------------
        % Write local cumulative mappability file if required
        % ------------------------------------------------------------------------------------------
        if ( ~isempty(iParams.outFile.lcmFile) && ~strcmp( iParams.outFile.name , iParams.outFile.lcmFile ) )
            writeLogFile( iParams.logFile , sprintf( '--Writing MaxTags %s file\n' , lcmFileParams.globalParameters.oFormat ) );
            switch lcmFileParams.globalParameters.oFormat
                case 'wig'
                    writeWiggle( lcmFileParams.globalParameters , chunk , round(lcmVal) );
                case 'bg'
                    writeBedGraph( lcmFileParams.globalParameters , chunk , round(lcmVal) );
                case 'mat'
                    writeMatFile( lcmFileParams.globalParameters , chunk , lcmVal );
            end
        end
        
    end
    
end

writeLogFile(iParams.logFile, sprintf('============\n'));

% **************************************************************************************************
%% NESTED FUNCTIONS
% **************************************************************************************************
    
    function catReadStart = concatenateReads(tmp_iChr)
        % Concatenates readStarts from all the align Files for a particular chromosome
        % function catReadStart = concatenateReads(tmp_iChr)
        % ------------------------------------------------------------------------------------------
        % INPUT ARGUMENT
        % -------------------
        % tmp_iChr<double>: chromosome index
        % -------------------
        % OUTPUT ARGUMENT
        % -------------------
        % catReadStart[double]: readStarts from all align Files are concatenated
        % -------------------
        % VARIABLES USED
        % -------------------
        % alignData (modified), iParams
        % ------------------------------------------------------------------------------------------
        
        % cumsum of #reads for chromosome tmp_iChr in each align file
        tmp_totReads = cumsum( arrayfun( @(x) ( x.nReads(tmp_iChr) ) , alignData ) );
        % initialize catReadStart
        catReadStart = zeros( tmp_totReads(end) , 1 , 'double' );
        % get reads from first tagAlign file
        catReadStart( 1 : tmp_totReads(1) ) = alignData(1).readStart{tmp_iChr};
        % save memory
        alignData(1).readStart{tmp_iChr} = [];
        % concatenate readStarts from remaining align files
        for tmp_ita = 2 : numel(alignData)
            catReadStart( ( tmp_totReads(tmp_ita - 1) + 1 ) : tmp_totReads(tmp_ita) ) = ...
                alignData(tmp_ita).readStart{tmp_iChr};
            % save memory
            alignData(tmp_ita).readStart{tmp_iChr} = [];
        end
    end

% **************************************************************************************************

    function [signalVal,lcmVal] = catReads2SignalPerChunk( minMaxLcmVal )
        % Normalizes signalVal for concatenated readStart
        % function [signalVal,lcmVal] = catReads2SignalPerChunk(tmp_iChr, minMaxLcmVal)
        % ------------------------------------------------------------------------------------------
        % INPUT ARGUMENTS
        % ----------------                
        % minMaxLcmVal: lcmVal <= minMaxLcmVal*nReps are capped to minMaxLcmVal*nReps before
        %             normalization (Used only for norm3 and 5)
        %             If < 1 it is assumed to be in terms of percentage
        % ----------------
        % OUTPUT ARGUMENTS
        % ----------------
        % signalVal[double]: normalized signal values
        % lcmVal[double]: maxTags values
        % ----------------
        % VARIABLES USED
        % ----------------
        % chunk, catReadStart, iParams, oFileParams, currUmap
        % ------------------------------------------------------------------------------------------
               
        % --------------------------------------------------------------------------------------
        % Compute signal values
        % --------------------------------------------------------------------------------------
        winSize = iParams.smooth.winLen;
        writeLogFile( iParams.logFile , sprintf('--Aggregating reads and extending tags\n') );
        signalVal = computeSignalVals( catReadStart , winSize , iParams.smooth.kernel , tukeyRatio , chunk ); % Aggregate and extend        
                                             
        % --------------------------------------------------------------------------------------
        % Compute local cumulative mappability values
        % --------------------------------------------------------------------------------------
        tmp_modeReadLen = oFileParams.specificParameters.modeReadLen(1);
        writeLogFile( iParams.logFile , sprintf('--Calculating Local Cummulative mappability\n' ) );
        lcmVal = numel(iParams.alignFname) * computeLcmVals( currUmap , tmp_modeReadLen , iParams.smooth.fragLen(1) , winSize , iParams.smooth.kernel , tukeyRatio , chunk );
                
        % --------------------------------------------------------------------------------------
        % Normalize signal values if required
        % Optimal normalization (norm5)
        % NS = (s1 + s2 + ...) /(av1 + av2 + ....)
        %       PER REPLICATE
        %       av = lcmVals X p                         : poisson rate
        %       p  = #reads / mappable_genome_size       : probability of read at each pos.
        % --------------------------------------------------------------------------------------

        % ----------------------------------------------------
        % Convert fractional minMaxLcmVal into absolute value
        % ----------------------------------------------------
        if (minMaxLcmVal < 0)
            minMaxLcmVal = 0;
        end

        minMaxLcmVal = ...
            (minMaxLcmVal < 1)  * minMaxLcmVal * 2 * sum( generateKernel( iParams.smooth.kernel , floor(winSize/2) , 0 , tukeyRatio) ) * numel(iParams.alignFname) + ...
            (minMaxLcmVal >= 1) * minMaxLcmVal * numel(iParams.alignFname); % both strands
        
        totalReads = oFileParams.specificParameters.nTotalReads; % total number of reads from all replicates                
        mgSize = oFileParams.specificParameters.mapGenomeSize(1); % mappable genome size      
        
        switch iParams.normFlag
            
            case 0 % no normalization
                signalVal = round(signalVal);                
            case 1
                % normSignal = (1e9 / totalReads) * signalVal
                writeLogFile(iParams.logFile, sprintf('--Normalizing signal\n'));
                signalVal = (1e9 / totalReads) * signalVal;
                
            case 2
                % normSignal = (1e9 / totalReads) * (signalVal / winSize)
                writeLogFile(iParams.logFile, sprintf('--Normalizing signal\n'));
                signalVal = (1e9 / (winSize * totalReads)) * signalVal;
                
            case 3 % normSignal = (1e9 / totalReads) * (signalVal / (maxtags/nReps))
                writeLogFile(iParams.logFile, sprintf('--Normalizing signal\n'));
                signalVal = (1e9 * numel(iParams.alignFname) / totalReads) * ( signalVal ./ lcmVal );
                truncateIdx = (lcmVal <= minMaxLcmVal);
                signalVal(truncateIdx) = signalVal(truncateIdx) .* lcmVal(truncateIdx) / minMaxLcmVal;
                signalVal(~isfinite(signalVal)) = 0; % replace Inf/Nan values with 0
                
            case 4 % normSignal = (mgSize / totalReads) * (signalVal / winSize)
                writeLogFile(iParams.logFile, sprintf('--Normalizing signal\n'));
                signalVal = (mgSize / (winSize * totalReads)) * signalVal;
                
            case 5 % normSignal = (mgSize / totalReads) * (signalVal / (maxtags/nReps))
                writeLogFile(iParams.logFile, sprintf('--Normalizing signal\n'));
                signalVal = (mgSize * numel(iParams.alignFname) / totalReads) * (signalVal ./ lcmVal);
                truncateIdx = (lcmVal <= minMaxLcmVal);
                signalVal(truncateIdx) = signalVal(truncateIdx) .* lcmVal(truncateIdx) / minMaxLcmVal;
                signalVal(~isfinite(signalVal)) = 0; % replace Inf/Nan values with 0
                
            otherwise
                error('Illegal normalization flag\n');
                
        end
        
        % Set lcmVals==0 to NaN (unmappable positions)
        lcmVal(lcmVal==0) = NaN;        
        
    end

% **************************************************************************************************

    function [signalVal,lcmVal] = perRepReads2SignalPerChunk( minMaxLcmVal )
        % Normalizes signalVal for concatenated readStart (unique readLen & fragLen across all reps)
        % function [signalVal,lcmVal] = perRepReads2SignalPerChunk(tmp_iChr,minMaxLcmVal)
        % ------------------------------------------------------------------------------------------
        % INPUT ARGUMENTS
        % ----------------                
        % minMaxLcmVal: lcmVal <= minMaxLcmVal are capped to minMaxLcmVal before normalization
        %             (Used only for norm3 and 5)
        %             If < 1 it is assumed to be in terms of percentage
        % ----------------
        % OUTPUT ARGUMENTS
        % ----------------
        % signalVal[double]: normalized signal values
        % lcmVal[double]: maxTags values
        % ----------------
        % VARIABLES USED
        % ----------------
        % chunk, alignData, iParams, oFileParams, currUmap
        % ------------------------------------------------------------------------------------------

        % --------------------------------------------------------------------------------------
        % Compute signal values
        % --------------------------------------------------------------------------------------
        winSize = iParams.smooth.winLen;
        writeLogFile( iParams.logFile , sprintf('--Aggregating reads and extending tags\n') );
        signalVal = computeSignalVals( catReadStart , winSize , iParams.smooth.kernel , tukeyRatio , chunk ); % Aggregate and extend        
        if iParams.normFlag == 0
            signalVal = round(signalVal);
        end
        
        % ----------------------------------------------------
        % Convert fractional minMaxLcmVal into absolute value
        % ----------------------------------------------------
        if (minMaxLcmVal < 0)
            minMaxLcmVal = 0;
        end
        
        minMaxLcmVal = ...
            (minMaxLcmVal < 1)  * minMaxLcmVal * 2 * sum( generateKernel( iParams.smooth.kernel , floor(winSize/2) , 0 , tukeyRatio) ) + ...
            (minMaxLcmVal >= 1) * minMaxLcmVal;
               
        % ------------------------------------------
        % total number of reads from all replicates
        % ------------------------------------------
        totalReads = oFileParams.specificParameters.nTotalReads;
                
        % --------------------------------------------------------------------------------------
        % For each replicate compute lcmVal and expectedRate
        % --------------------------------------------------------------------------------------
        lcmVal = 0;        
        expectedRate = 0;

        for tmp_ita = 1:numel(iParams.alignFname)
        
            % ------------------------------------------
            % Number of reads in current replicate
            % ------------------------------------------
            readsPerRep = oFileParams.specificParameters.nReadsPerRep(tmp_ita);             
            
            % -------------------------------------------
            % mappable genome size for current replicate
            % -------------------------------------------
            mgSize = oFileParams.specificParameters.mapGenomeSize(tmp_ita);
                                   
            % ------------------------
            % Compute local cumulative mappability values
            % ------------------------
            tmp_modeReadLen = oFileParams.specificParameters.modeReadLen(tmp_ita);
            writeLogFile( iParams.logFile , sprintf('--Calculating MaxTags\n') );
            tmp_lcmVal = computeLcmVals( currUmap , tmp_modeReadLen , iParams.smooth.fragLen(tmp_ita) , winSize , iParams.smooth.kernel , tukeyRatio, chunk );            
            lcmVal = lcmVal + tmp_lcmVal;
            
            % -----------------------
            % Compute expected rate
            % -----------------------           
            switch iParams.normFlag
                
                case 0 % no normalization
                    expectedRate = 1;
                    
                case 1
                    % normSignal = (1e9 / totalReads) * signalVal
                    writeLogFile(iParams.logFile, sprintf('--Normalizing signal\n'));
                    expectedRate = totalReads / 1e9;                    
                    
                case 2
                    writeLogFile(iParams.logFile, sprintf('--Normalizing signal\n'));
                    expectedRate = expectedRate + (readsPerRep * winSize) / 1e9;
                    
                case 3
                    writeLogFile(iParams.logFile, sprintf('--Normalizing signal\n'));
                    tmp_lcmVal( (tmp_lcmVal <= minMaxLcmVal) & (tmp_lcmVal > 0) ) = minMaxLcmVal;
                    expectedRate = expectedRate + (readsPerRep / 1e9) * tmp_lcmVal;                    
                    
                case 4
                    writeLogFile(iParams.logFile, sprintf('--Normalizing signal\n'));
                    expectedRate = expectedRate + (readsPerRep * winSize) / mgSize;
                    
                case 5
                    writeLogFile(iParams.logFile, sprintf('--Normalizing signal\n'));
                    tmp_lcmVal( (tmp_lcmVal <= minMaxLcmVal) & (tmp_lcmVal > 0) ) = minMaxLcmVal;
                    expectedRate = expectedRate + (readsPerRep / mgSize) * tmp_lcmVal;                    
                    
                otherwise
                    error('Illegal normalization flag\n');
                    
            end
            
        end
        
        % ------------------------------------------------------------------------------------------        
        % Normalize signal values if required        
        % NS = (s1 + s2 + ...) /(av1 + av2 + ....)        
        % ------------------------------------------------------------------------------------------                
        signalVal = signalVal ./ expectedRate;
        signalVal( ~isfinite(signalVal) ) = 0; % replace Inf/Nan values with 0
        
        % Set lcmVals==0 to NaN (unmappable positions)
        lcmVal(lcmVal==0) = NaN;
        
    end

% **************************************************************************************************

    function [] = mapFilterSignal()
        % Set signalVal=NaN if pos is unmappable OR does not pass mappability threshold
        % signalVal = 0 if position passed mappability threshold BUT no extended/smoothed reads
        %               overlap the position 
        % function [] = mapFilterSignal()
        % ------------------------------------------------------------------------------------------
        % VARIABLES USED
        % ----------------
        % signalVal, lcmVal, iParams
        % ------------------------------------------------------------------------------------------
        
        tmp_mapFilter = iParams.mapFilter;
        
        % Compute maximum value of maxTags (sum of smooth winSize over all replicates, both strands)               
        maxLcmVal = 2 * numel(iParams.alignFname) * sum( generateKernel( iParams.smooth.kernel , floor(iParams.smooth.winLen / 2) , 0 , tukeyRatio) ) ; 
        
        % if mapFilter < 1 convert to absolute scale
        tmp_mapFilter = ...
            (tmp_mapFilter < 1)  * tmp_mapFilter * maxLcmVal + ...
            (tmp_mapFilter >= 1) * tmp_mapFilter * numel(iParams.alignFname);

        % Set signalVal=NaN if position is unmappable OR does not pass mappability threshold
        signalVal( isnan(lcmVal) | ( lcmVal <= tmp_mapFilter ) ) = NaN;
        
    end

% **************************************************************************************************
    
    function [] = markSeqGapsInLcmVals()
        % Set lcmVal=0 if position is mappable after ext/smoothing but lies in sequence gaps(Ns)
        % function [] = markSeqGapsInLcmVals()
        % ------------------------------------------------------------------------------------------
        % VARIABLES USED
        % -----------------
        % lcmVal (modified), chr, chunk
        % ------------------------------------------------------------------------------------------               
        
        % find positions that are 'Ns'
        tmp_chrMap = memmapfile( chr.mmNames{iChr} , 'format' , 'uint8' );
        unMapIdx = ( tmp_chrMap.data( chunk.start : chunk.stop ) == nt2int('N') );
        
        % Set positions that are 'N' AND have finite maxtags values to 0
        unMapIdx = ( unMapIdx & isfinite(lcmVal) );
        lcmVal(unMapIdx) = 0;
        
    end

end