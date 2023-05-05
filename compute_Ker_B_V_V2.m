function [Ker_B, Ker_V] = compute_Ker_B_V_V2(aSt)
% BrRateKer_b is a square probability matrix of obtaining praticle with size vD from mother sized vM. vD
% on rows vM on columns. Bb = 2*BrRate*Ker_B.*N returns a column vector of size(N)
% Ker_V is a square matrix that returns the total size/volume of newly generated particles in cell vD, originating from mother sized vM.
global noGhostInd nCellsGh noGhostInd speciesPBE spNames
for aSp = speciesPBE
 sp = char(spNames{aSp});
 % aSt.(sp).x2SCell is the lowest particle size created through erosion. This
 %is the lower end of the particle sizes generated from the redistribution/kernel function
 % NB: It is important that this DOES NOT coincide with the
 % lower boundary of the first cell of the grid. If this was the
 % case, because of the U-shaped distribution, I would have an
 % accumulation of new particles with size close to the
 % boundary, which cannot be represented correctly by the first
 % grid point, cannot be reassigned to a lower grid point (first grid point is the first!)
 % and this would cause loss of the consistency of the first
 % moment. In this way, instead, the smallest particle size
 % created through erosion is the single cell, whose size
 % coincide with the first grid point.
 if ~strcmp(string(aSt.(sp).x2SCell),string(aSt.(sp).x2_all(noGhostInd(1))))
    ME = MExeption('x2_all(1) should be centered at x2SCell. Check in DefineGrid inside PBEParams that things is done correctly');
    throw(ME)
 end
 % adding ghost cells
 Ker_B.(sp) = zeros(nCellsGh); 
 Ker_V.(sp) = zeros(nCellsGh); 
 for k = noGhostInd %i index of mother Cell with vM. I neglect the first cell as I would have vM=x2SCell and vD=x2SCell, thus I would have the particle breaking just into itself
        vM = aSt.(sp).x2_all(k);
        all_i = noGhostInd(1):k; %k-long vector of all possible daughter cell indexes for this k mother (up to its own size)
        lbs = aSt.(sp).x2_lb(all_i);
        % In the first (noGhost) cell vD can take a minimal value of x2SCell ->
        % thus the lb for the first (noGhost) cell is redefined
        lbs(1) = aSt.(sp).x2SCell;
        ubs = aSt.(sp).x2_ub(all_i);% all cell upper bounds
        ubs(end) = vM;% for the cell of the mother itself the daugter can have size just up to the reference (vM)
        
        % Compute the zeroth moments (B terms) for all_i cells
        % intBrKerWWT(vD,vM,vMin,BorV)
        distFnB_ub = intBrKerWWT_V2(ubs,vM,aSt.(sp).x2SCell,'B');
        distFnB_lb = intBrKerWWT_V2(lbs,vM,aSt.(sp).x2SCell,'B');
        distFnB = distFnB_ub-distFnB_lb;
        Ker_B.(sp)(all_i,k) = distFnB;
        
        % Compute the first moments (V terms) for all_i cells
        distFnx2_ub = intBrKerWWT_V2(ubs,vM,aSt.(sp).x2SCell,'V');
        distFnx2_lb = intBrKerWWT_V2(lbs,vM,aSt.(sp).x2SCell,'V');
        distFnx2 = distFnx2_ub-distFnx2_lb;
        Ker_V.(sp)(all_i,k) = distFnx2;

        % Check that the zeroth and the first moments of newly
        % generated particles in each cell is not negative
        if ~isempty(distFnB(distFnB<0)) || ~isempty(distFnx2(distFnx2<0)) 
            disp(string([lbs ubs]))
            disp(string([distFnx2_lb distFnx2_ub]))
            ME = MExeption('Integral evaluated at the ub should always be bigger than the one evaluated at the lb');
            throw(ME);
        end
 end
 
 % recheck negativities
 if ~isempty(Ker_B.(sp)(Ker_B.(sp)<0)) || ~isempty(Ker_V.(sp)(Ker_V.(sp)<0))
    ME = MExeption('There should not be negative values in the breakage kernels');
    throw(ME);
 end
 
% cumulative density functions for each vM have to reach 1 (i.e. I have to have 100% probability of finding an allocation for a daughter given a specific mother)
 cums = sum(Ker_B.(sp),1);
 if ~isequal(string(cums(noGhostInd(2:end)))',string(ones(size(noGhostInd(2:end))))')
    ME = MExeption('cumulative density functions for each vM have to reach 1 ');
    throw(ME);
 end
end
end


