function P = changeParam(P,idx,parChangingNames)
% idx is the index of the parameter set I am using in this run. P.(str) can
% be a vector or a matrix
% filename = strcat([P.folderOutput,'debug.txt']);
% fileID = fopen(filename, 'a');
% fprintf(fileID, 'run %d: ', idx);
 for jj = 1:length(parChangingNames) % for each of the parameters that I change
    str = strcat(['p' num2str(jj)]);       
    P.(char(parChangingNames(jj))) = P.(str)(idx); %I assign the ii-th par value to its specific parameter name 
    % fprintf(fileID, strcat([char(parChangingNames(jj)),' = %.3f']), P.(str)(idx));
 end
% fprintf(fileID, '\n');
% fclose(fileID);
end
