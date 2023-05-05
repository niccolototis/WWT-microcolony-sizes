function varargout = loadData(varNames,isBigFile,label)
global T_sim folderInput
% n = nargin
% The first n-3 inputs have to be the variables that need to be written
% The n-2th input has to be a cell array of character vectors indicating
% the names of the variables
% The n-1th input is a boolean array that says if the format for big files needs to
% be used
% The nth input is a char variable that provide an additional label put as
% a suffix to the variables

nVars = length(varNames);
if ~iscellstr(varNames) || length(varNames)~= (nVars)
 ME = MException('Provide ALL the names of the nVars variables to be written to file');
 throw(ME);
end
if ~islogical(isBigFile) || nVars~= length(isBigFile)
 ME = MException('Say what variables need to be treated as big files, use boolean([0,1,1,0]) for instance');
 throw(ME);
end
if ~ischar(label)
 label = ''; 
else
 label = strcat([label,'_']);
end
for i = 1:nVars
 if isBigFile(i)
    load(strcat([folderInput,varNames{i},'_',label,num2str(T_sim),'.mat']));
    varargout{i} = aVar; %aVar is the name that has been used in writeData
 else
    varargout{i} = readmatrix(strcat([folderInput,varNames{i},'_',label,num2str(T_sim),'.txt']));
 end
end
return
