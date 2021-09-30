memoryForMac()


function memoryForMac()
% This function will return the memory used by MATLAB on the MAC
%
%% First get the version of MATLAB
curVer = version('-release');
%% get the PID for MATLAB
sysStr = ['ps -o ppid,command |grep ',curVer,'.app'];
[status,verStr] = system(sysStr);
if (status ~= 0)
    error('MATLAB was not found: That is odd since you are in MATLAB');
end
%% Get where the string is located
% Format looks like: interested in PPID
%  PPID COMMAND
%  4151 /Applications/MATLAB_R2019b.app/bin/maci64/matlab_helper /dev/ttys000 yes
slash = findstr('/',verStr);
pidStr = verStr(1:slash(1)-1);
%% Now get the memory string
sysStr = ['top -pid ',pidStr,' -l 1 -stats MEM'];
[status,info] = system(sysStr);
if (status ~= 0)
    error('Invalid PID provided')
else
    % now parse out the memory
    memLoc = findstr(info,'MEM');
    MEM = info(memLoc+5:end-1);
    fprintf('Total memory used: %s\n',MEM);
end

end