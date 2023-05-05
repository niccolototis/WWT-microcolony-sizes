function [BrRateKer_b, BrRateKer_v] = BrekageKernel(v_all,v_lb,v_ub,nCellsGh)

% BrRateKer_b is a square probability matrix of obtaining praticle with size vD from mother sized vM. vD
% on rows vM on columns. Bb = 2*BrRateKer_b*(BrRate.*N) returns a column vector of size(N)
% BrRateKer_v is a square matrix that returns the particle volume appearing in cell vD, originating from mother sized vM. Vb = 2*BrRateKer_v*(BrRate.*N) returns a column vector of size(N)
% 
% Check analytical derivation of the breakage kernel! int{over_z}(Dirac(z)dz) = 1, z=(2*xD-xM) has to hold. as 2 is the coefficient that multplies the integrating variable then a 2 factor has to be inserted if we integrate over xD int{over_xD}(2*Dirac(2*xD-xM)dxD)=1 
% This means that the sum of the matrix vertically along the columns (integration over the rows) needs to
% equal 1

% for aScale = scales
%     sC = string(aScale);
 BrRateKer_b = zeros(nCellsGh); % adding ghost cells
 for i = 1 : nCellsGh
    if(i ==1)
        lb = 0;
    else
        lb = v_lb(i);
    end
  for k = i : nCellsGh
        vM = v_all(k);
        vD = vM/2;
        if (i == k)
            ub = v_all(i);
        else
            ub = v_ub(i);
        end
        if (vD>= lb && vD<ub) % NB the extreme is included in just one side
            BrRateKer_b(i,k) = 1; % IMPO!!! the breakage kernel must integrate to 1 over xD (== summing vertically to 1)
        end % remains zero otherwise
  end
  
 end
 BrRateKer_b = BrRateKer_b./sum(BrRateKer_b,1); % normalization condition, the breakage kernel must integrate to 1 over xD (== summing vertically to 1)
 BrRateKer_b(isnan(BrRateKer_b)) = 0; % resolving division by zero
 BrRateKer_v = BrRateKer_b.*v_all;
% end
% [0, 0, v_all'; v_lb,v_ub, BrRateKer_b]
return 
