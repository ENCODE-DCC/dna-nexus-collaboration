function [] = writeBedGraph( fileParams , chunk , signal )
% Converts a signal vector into a bedgraph file
% function [] = writeBedGraph(fileParams,signal,chunk)
% -------------------------------------------------------------------------
% INPUT ARGUMENTS:
% ------------------
% fileParams: structure containing details about the output wiggle file
%           .oFname: full path and name of output file
%                    OR 'stdout','stderr',''
%                    if '' then print to STDOUT
%           .oAppend: if set to 1 the file will be appended else a new file
%                     will be created
% chunk: details about the current chunk of signal
%           .chrName: name of chromosome
%           .start: start coordinate
%           .stop: stop coordinate
%           .len: length of chunk
%           .isFirst: if 1 then it is the first chunk
%           .isLast: if 1 then it is the last chunk
% signal[double]: column vector of signal values (invalid positions have
%                 value NaN)
% ------------------
% OUTPUT ARGUMENTS:
% ------------------
% None
% -------------------------------------------------------------------------

signal = round( signal / 1e-1 ) * 1e-1;
% --------------------------------------------------
% Convert chunk fields into scalars to speed up loops
% --------------------------------------------------
chunk_start = chunk.start;
chunk_chrName = chunk.chrName;

% --------------------------------------------------
% Find valid parts (finite signal) inside the chunk
% --------------------------------------------------
validIdx = diff([ 0 ; isfinite(signal) ; 0 ]);
% start idx relative to signal(1)
part_start = find( validIdx == 1 );
% stop idx relative to signal(1)
part_stop = find( validIdx == -1 ) - 1;
% number of valid contiguous valid chunks
nParts = numel(part_start); 
clear validIdx;

% --------------------------------------------------
% Open bedgraph file
% --------------------------------------------------
fPerm = char( fileParams.oAppend * 'a' + ~fileParams.oAppend * 'w' ); % file permission
switch lower(fileParams.oFname)
    case {'stdout',''}
        fp = 1; % STDOUT
    case 'stderr'
        fp = 2; % STDERR
    otherwise
        fp = fopen(fileParams.oFname,fPerm);    
end
assert( (fp ~= -1) ,'Unable to open bedgraph file for writing\n');

% -------------------------------------------------------------------------
% Write valid parts into a bedgraph file
% BEDGRAPH FORMAT
% chromA  chromStartA(0-based)  chromEndA(1-based)  dataValueA
% -------------------------------------------------------------------------

% format string to use for printing
formatString = [ chunk_chrName , '\t%d\t%d\t%.1f\n' ]; 

for iPart = 1:nParts
    
    % --------------------------
    % if length of part is 1 bp
    % --------------------------
    if ( part_start(iPart) == part_stop(iPart) ) 
        fprintf( fp , formatString , ...
            ( chunk_start + part_start(iPart) - 2 ) , ...
            ( chunk_start + part_stop(iPart) - 1 ) , ...
            signal( part_start(iPart) ) );

    % --------------------------
    % if length of part is > 1 bp
    % --------------------------        
    else
       
        % -----------------------------------------------------------------
        % Get subParts of contigous repeated values
        % -----------------------------------------------------------------
        
        % Mark contiguous repeat values with 0 (using diff). Then set those
        % all repeat value positions to 1 and all non-repeat values to 0. A
        % second diff() will mark all start/stop positions of contiguous
        % repeated regions with +1/-1
        validIdx = diff( [ -1 ; signal( part_start(iPart) : part_stop(iPart) ) ; -1 ] );        
        validIdx = diff(validIdx == 0);
        
        % #contiguous repeat chunks
        nVals = nnz(validIdx == 1); 
        subPart = [];        
        if nVals
            % (chromStart,chromEnd,signalval) repeated ...
            subPart = zeros(nVals,3); 
            % starts of contiguous chunks relative to signal(1)
            tmp_start = (part_start(iPart) - 1) + find(validIdx==1); 
            % signal values
            subPart(:,3) = signal(tmp_start);
            % absolute start of contiguous chunk (0-based)
            subPart(:,1) = (chunk_start - 2) + tmp_start; 
            clear nVals tmp_start;
            % absolute stop of contiguous chunk (1-based)
            subPart(:,2) = (chunk_start - 1 + part_start(iPart) - 1) + find(validIdx==-1); 
        end
        
        % -----------------------------------------------------------------
        % Get subParts of contigous NON-repeated values
        % -----------------------------------------------------------------
        
        % length of current part
        part_len = part_stop(iPart) - part_start(iPart) + 1; 
        % start of subPart relative to part_start(iPart)
        tmp_start = [ 1 ; find( validIdx == -1 ) + 1 ]; 
        % stop of subPart relative to part_start(iPart)
        tmp_stop = [ find( validIdx == 1 ) - 1 ; part_len ]; 
        % find valid start <= stop
        keepIdx = (tmp_start <= tmp_stop); 
        if any(keepIdx)
            tmp_start = tmp_start(keepIdx); % remove invalid starts
            tmp_stop = tmp_stop(keepIdx); % remove invalid stops
            clear keepIdx;
            
            % Mark non-repeated positions with 1s
            validIdx = zeros( part_len + 1 , 1 );
            % mark start of non-repeated regions with 1
            validIdx(tmp_start) = 1; 
            clear tmp_start;
            % mark stop+1 of non-repeated regions with -1
            validIdx(tmp_stop + 1) = -1; 
            clear tmp_stop;
            % cumsum will mark all positions in the non-repeated regions with 1
            validIdx = find(cumsum(validIdx)); 
            % #non-repeated positions
            nVals = numel(validIdx); 
            if nVals
                % (chromStart,chromEnd,signalval) repeated ...
                subPartNonRepeat = zeros( nVals , 3 );
                % absolute start (0-based)
                subPartNonRepeat(:,1) = (chunk_start - 2 + part_start(iPart) - 1) + validIdx;
                % absolute stop (1-based)
                subPartNonRepeat(:,2) = (chunk_start - 1 + part_start(iPart) - 1) + validIdx;
                % signal value at that position
                subPartNonRepeat(:,3) = signal( ( part_start(iPart) - 1 ) + validIdx ); 
                clear validIdx;
                % merge subPart and subPartNonRepeat
                subPart = [ subPart ; subPartNonRepeat ];  %#ok<AGROW>
                clear subPartNonRepeat;
                % Sort the subParts by chromStart
                [ subPart(:,1) , sidx ] = sort( subPart(:,1) );
                subPart( : , [2,3] ) = subPart( sidx , [2,3] );
                clear sidx;
            end            
        end
        
        % -----------------------------------------------------------------
        % Print the subParts
        % -----------------------------------------------------------------
        if ~isempty(subPart)
            fprintf( fp , formatString , subPart' );
            clear subPart;
        end
    end
end

% -------------
% Close file
% -------------
if ~any( fp == [1,2] )
    fclose(fp);
end
end