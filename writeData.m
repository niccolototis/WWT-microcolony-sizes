function [] = writeData(varargin)
    global T_sim folderData
    
    n = nargin;
    nVars = n - 3;
    varNames = varargin{n - 2};
    isBigFile = varargin{n - 1};
    label = varargin{n};

    if ~iscellstr(varNames) || length(varNames) ~= nVars
        ME = MException('Provide ALL the names of the nVars variables to be written to file');
        throw(ME);
    end
    
    if ~islogical(isBigFile)
        ME = MException('Say what variables need to be treated as big files, use boolean([0,1,1,0]) for instance');
        throw(ME);
    end
    
    if ~ischar(label)
        label = ''; 
    else
        label = strcat([label, '_']);
    end
    
    for i = 1:nVars
        aVarName = varNames{i};
        aVar = varargin{i};
        
        if isBigFile(i)
            save(strcat([folderData, aVarName, '_', label, num2str(T_sim), '.mat']), 'aVar', '-v7.3');
        else
            writematrix(aVar, strcat([folderData, aVarName, '_', label, num2str(T_sim), '.txt']));
        end
    end
    return
end
