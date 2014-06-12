function [] = writeLogFile(logFile,logStatement)
% Write debugging/verbose statements to a logFile
% function [] = writeLogFile(logFile,logStatement)
% -------------------------------------------------------------------------
% INPUT PARAMETERS:
% -------------------
% logFile: full path and name of logFile. If '' no output
% logStatement: statement to be logged
% -------------------------------------------------------------------------

if isempty(logFile)
    return;
end

switch lower(logFile)
    case {'stdout',''}
        lp = 1; % stdout
        fprintf(lp,'%s',logStatement);        
    case 'stderr'
        lp = 2; % stderr
        fprintf(lp,'%s',logStatement);        
    otherwise
        lp = fopen(logFile,'a');
        fprintf(lp,'%s',logStatement);
        fclose(lp);        
end
end