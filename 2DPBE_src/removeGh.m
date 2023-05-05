function clean = removeGh(g) %g is either a matrix or a vector
 global nCellsGh
 noGhostInd = ones(nCellsGh,1);
 [noGhostInd(1),noGhostInd(end)] = deal(0);
 noGhostInd = logical(noGhostInd);
 if ~any(size(g) == nCellsGh)
        ME = MException('Wrong dimension of g in order to remove ghost cells it has to have at leas one dimension equal to nCellsGh');
        throw(ME)
 end 
 if (size(g,1) ==1 || size(g,2)==1)
    clean = g(noGhostInd);
 else
    if size(g,1) ==nCellsGh
        clean = g(noGhostInd,:);
    else
        clean = g(:,noGhostInd);
    end
 end
end
