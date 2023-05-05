function x = ghostToZero(x,t)
global ghostIdx 
if isvector(x)
numericalThr = 0; %1e-10
 if (abs(x(1))>numericalThr || abs(x(end))>numericalThr)
%     t
    strano = 3; % I should not have positive stuff in these ghost cells!
 end
 x(ghostIdx) = zeros(size(ghostIdx)); 
else
 if size(x,1) ==ghostIdx(2)
    x(ghostIdx,:) = 0;
 else
    if size(x,2) ==ghostIdx(2)
        x(:,ghostIdx) = 0;
    else
        ME = MException('Wrong matrix dimensions');
        throw(ME)
    end
 end
end
return
