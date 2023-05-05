 function [Ker_B, Ker_V] = defineBrekageKernel(aSt)
% Ker_B is a square probability matrix of obtaining praticle with size vD
% from mother sized vM. vD on rows vM on columns. Bb = 2*Ker_B*(BrRate.*N)
% returns a column vector of size(N) Ker_V is a square matrix that returns
% the particle volume appearing in cell vD, originating from mother sized
% vM. Vb = 2*Ker_V*(BrRate.*N) returns a column vector of size(N)

global spNames FirstRefNoBreak speciesPBE nCells x2MinIsSCell
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
 if ~strcmp(string(aSt.(sp).x2SCell),string(aSt.(sp).x2_all(1)))
    ME = MExeption('x2_all(1) should be centered at x2SCell. Check in DefineGrid inside PBEParams that things is done correctly');
    throw(ME)
 end
 Ker_B.(sp) = zeros(nCells); % adding ghost cells
 Ker_V.(sp) = zeros(nCells); 
 if FirstRefNoBreak || x2MinIsSCell
    motherIdx = 2:nCells;
 else
    motherIdx = 1:nCells;
 end
 vDMin = aSt.(sp).x2Min;
 if any((aSt.(sp).x2_all(motherIdx)-aSt.(sp).x2_lb(motherIdx))<vDMin)
    % in this case I know some vDMax will fall below the lb of the cell
    % that contains vM -> I will have that in the cell that contains vM
    % I put lb ==ub to have
 end
 for k = motherIdx %i index of mother Cell with vM. I neglect the first cell as I would have vM=x2SCell and vD=x2SCell, thus I would have the particle breaking just into itself
        vM = aSt.(sp).x2_all(k);
        all_i = (1):k; %k-long vector of all possible daughter cell indexes for this k mother (up to its own size)
        lbs = aSt.(sp).x2_lb(all_i);
        % lbs(1) = aSt.(sp).x2SCell;
        % NB potentially I can have vM<lbs(k)
        ubs = aSt.(sp).x2_ub(all_i);% all cell upper bounds
        % -
        vDMax = vM-vDMin;
        % - I want to have the breakage kernel to be zero for any
        % v>vDmax. I do so by putting ubs = vDMax for the cell that
        % contains vDmax, and ubs = lbs for all cells above that one up
        % till vM
        ubs(ubs>vDMax) = max(vDMax,lbs(ubs>vDMax));
        [m0_vD,m1_vD] = deal(zeros(size(lbs)));
        isE = ubs==lbs;
        %  m0_vD(isE) = 0; m1_vD(isE)=0; % REDUNDANT
        if any(lbs>ubs)
            ME = MException('they should not even be equal. It might be that because of ubs(end)=vM-vDMin, the ubs(end) has fallen below the lbs(end) ');
            throw(ME);
        end
        isIN = ~isE;
        % NB if any(isE) then I have that the size of lbs_x m0_x is
        % smaller than lbs and m0_vD
        lbs_x = varSubst(lbs(isIN),vM,vDMin,'vD_to_x');
        ubs_x = varSubst(ubs(isIN),vM,vDMin,'vD_to_x');
        
        % - Compute Ker_B
        m0_x = inKerWWT_x(ubs_x,'B')-inKerWWT_x(lbs_x,'B');
        % -the zeroth moment gives me the expected value E[r^0], so does not need to be transformed
        m0_vD(isIN) = m0_x; 
        Ker_B.(sp)(all_i,k) = m0_vD; 
        
        % - Compute Ker_V
        m1_x = inKerWWT_x(ubs_x,'V')-inKerWWT_x(lbs_x,'V');
        % - See paper for details
        m1_vD(isIN) = (vM-2*vDMin).*m1_x+vDMin.*m0_x;
        Ker_V.(sp)(all_i,k) = m1_vD;
 end
 if ~isempty(Ker_B.(sp)(Ker_B.(sp)<0)) || ~isempty(Ker_V.(sp)(Ker_V.(sp)<0))
    ME = MExeption('There should not be negative values in the breakage kernels');
    throw(ME);
 end
 
% cumulative density functions for each vM have to reach 1 (i.e. I have to
% have 100% probability of finding an allocation for a daughter given a
% specific mother) Sum of the differences in cumulative density function is
% already the integral itself
 cums = sum(Ker_B.(sp),1);
 indNZ = find(cums~=0);
 if ~isequal(string(cums(indNZ))',string(ones(size(cums(indNZ))))')
    ME = MException('cumulative density functions for each vM have to reach 1 ');
    throw(ME);
 end 
 
 cums = sum(Ker_V.(sp),1);
 expectedV = aSt.(sp).x2_all/2;
 if ~isequal(string(cums(indNZ))',string(expectedV(indNZ')))
    ME = MException('cumulative density functions for each vM have to reach 1 ');
    throw(ME);
 end 

% - check that the average volumes in each cell is inbetween the boundaries
vAV = Ker_V.(sp)./Ker_B.(sp);
for ind_vD = 1:nCells
    if any(vAV(ind_vD,:)<aSt.(sp).x2_lb(ind_vD)) || any(vAV(ind_vD,:)>aSt.(sp).x2_ub(ind_vD))
         ind_l = find(vAV(ind_vD,:)<aSt.(sp).x2_lb(ind_vD));
         ind_u = find(vAV(ind_vD,:)>aSt.(sp).x2_ub(ind_vD));
         ind_vM = sort([ind_l;ind_u]);       
         disp([ind_vD ind_vM])
    end
end
end
return

function varOUT = varSubst(varIN,vM,vDMin,sense)
% - Variable substitution, check the paper for details
global sC
switch sense
 case 'vD_to_x'
    vD = varIN;
    x = (vD-vDMin)./(vM-2*vDMin);
    
    if any(x>1) || any(x<0)
        %     ME = MException('Look at the paper. The pdf(x) and cdf(x) are defined only for 0<=x<=1');
        %     throw(ME);
        % - Check that this does not mess things
        ind_u = find(x>1);
        if ~isempty(ind_u)
            for j = ind_u
                [~,A,~] = reallyEqual(x(j),1,'~','Make sure that ubs(end)==vDMax');
                x(j) = A;
            end
        end
        ind_l = find(x<0);
        ind_l(~isempty(ind_l))
        if ~isempty(ind_l)
            for j = ind_l
                [~,A,~] = reallyEqual(x(j),0,'~','Make sure that x2_lb_all(1)==x2SCell');
                x(j) = A;
            end
        end
    end    
    
    varOUT = x;
end
return

