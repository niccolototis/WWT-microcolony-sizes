function B_CA_aSp = computeRearrangeCAT_new(B,V,aStp,P,t) 
global out_vAV nCells
% [1] J. Kumar, M. Peglow, G. Warnecke, and S. Heinrich, “An efficient
% numerical technique for solving population balance equation involving
% aggregation, breakage, growth and nucleation,” Powder Technology, vol.
% 182, no. 1, pp. 81–104, 2008.

x2_all = aStp.x2_all;

M_NisTilde = false;
zeroth_out_vAVs = 0;
integralOfTheMarg = computeMarg(B,M_NisTilde,[],aStp.x2_widths,[],'x1'); % with this I obtain the marginal for x1 Marg_x1 (no tilde)

[B_CA_aSp,f,fp,fm] = deal(zeros(size(x2_all)));
% checkFirstMomentConserved(B,V,x2_all,t)
vAV = V./B;
vAV(isnan(vAV)) = 0;
if ~isempty(vAV(vAV<0))
 ME = MException('It should never happen that the average volume inside a cell is negative');
 throw(ME);
end
ind_nn0 = find(vAV~=0);

% - MEMO: if x2MinIsSCell ==true Cells of the first grid point do not break
% and when I have negative growth the vAV would become
% vAV((1))<x2_lb((1)) and vAV((1))<x2_all((1))
% = > Here I make an adjustment that keeps the consistency of the first
% moment by removing a number of particles from the first grid point adding
% up to a volume that accounts for the difference vAV((1))<x2_all((1))
% Only after this is done, I can check that all vAVs are inbetween lbs and
% ubs

% - NB these adjustments are made to keep the consistency of the FIRST
% moment, whereas the zeroth moment is affected
if ismember(1,ind_nn0) && vAV((1))<x2_all((1))
 zeroth_out_vAVs = 1;
 % - How many cells of volume x2_all(1ref) can I add at x2_all(1ref) to come up
 % with the same volume B(1ref)*vAV(1ref)?
 % B(1ref)*vAV(1ref) = Badj(1ref)*x2_all(1ref) -> Badj(1ref)=(B(1ref)*vAV(1ref))/x2_all(1ref)
 B((1)) = (B((1))*vAV((1)))/x2_all((1));
 % - reset the new vAV
 vAV((1)) = x2_all((1));
 % - IMPO: still I go trough the CAT because the first reference can
 % still receive a fraction from the second cell
end
if ismember(P.nCells,ind_nn0) && vAV(P.nCells)>x2_all(P.nCells)
 zeroth_out_vAVs = 1;
 % - How many cells of volume x2_all(lastref) can I add at x2_all(lastref) to come up
 % with the same volume B(lastref)*vAV(lastref)?
 % B(lastref)*vAV(lastref) = Badj(lastref)*x2_all(lastref) -> Badj(lastref)=(B(lastref)*vAV(lastref))/x2_all(lastref)
 B(P.nCells) = (B(P.nCells)*vAV(P.nCells))/x2_all(P.nCells);
 % - reset the new vAV
 vAV(P.nCells) = x2_all(P.nCells);
 % - IMPO: still I go trough the CAT because the lastref reference can
 % still receive a fraction from the bottom cell
end

if any(vAV(ind_nn0)<aStp.x2_lb(ind_nn0)) || any(vAV(ind_nn0)>aStp.x2_ub(ind_nn0))
 if ~P.identify_lambda && ~P.sweep_lambda
 ind = find(vAV(ind_nn0)<aStp.x2_lb(ind_nn0));
 out_vAV = max(out_vAV,vAV(ind_nn0(ind))-aStp.x2_lb(ind_nn0(ind)));
 [aStp.x2_lb(ind_nn0(ind)) vAV(ind_nn0(ind)) aStp.x2_ub(ind_nn0(ind))]
 end
 disp(string([ vAV V B ]));
 scsME = MException('The average particle volume in any cell needs to stay between the lb and ub of the cell itself. If not, either B or V are not created accurately');
 throw(ME);
end

% - In cells where average volume vAV falls exactly on the reference, Births
% are not reassigned
ind_E = ind_nn0(vAV(ind_nn0)==x2_all(ind_nn0));

% - Here If I want to consider that smaller differences do not matter I can
% use the string (not recommended though)
% ind_E = ind_nn0(strcmp(string(vAV(ind_nn0)),string(x2_all(ind_nn0))));
% if isempty(ind_E) && isempty(ind_e) && ~isequal(ind_E,ind_e)
% ME = MException('This can create problem ind_E, ind_U, ind_L not being disjucted');
% throw(ME)
% end

B_CA_aSp(ind_E) = B_CA_aSp(ind_E)+B(ind_E);

% - Find cells where average volume vAV falls above (U = up) the reference
ind_U = ind_nn0(vAV(ind_nn0)>x2_all(ind_nn0));
if ismember(ind_U,P.nCells)
 ME = MException('This should not happen because of the adjustment made before');
 throw(ME)
end

i = ind_U;
% - For these i-cells I will have a fraction of B(i) staying at x2_all(i)
% (called f(i), is the same of a1/B(i) in [1]) and a fraction of B(i) to be
% reassigned to x2_all(i+1) (called fp(i) , is the same of a2/B(i) in [1])
f(i) = (vAV(i)-x2_all(i+1))./(x2_all(i)-x2_all(i+1));
fp(i) = 1-f(i); %(vAV(i)-x2_all(i))./(x2_all(i+1)-x2_all(i));
in01(f(i))
in01(fp(i))
if ~isequal((f(i)+fp(i)),ones(size(f(i))))
 ME = MException('Fractions are calculated wrongly');
 throw(ME);
end

% - Here doing either one of the following is the same
% B_CA_aSp(ind_U) = B(ind_U).*f(ind_U);
B_CA_aSp(ind_U) = B_CA_aSp(ind_U)+B(ind_U).*f(ind_U);

% - Here I need to sum up
B_CA_aSp(ind_U+1) = B_CA_aSp(ind_U+1)+B(ind_U).*fp(ind_U);


% - Find cells where average volume vAV falls above (L = low) the reference
ind_L = ind_nn0(vAV(ind_nn0)<x2_all(ind_nn0));
if ismember(ind_L,(1))
 ME = MException('This should not happen because of the adjustment made before');
 throw(ME)
end
i = ind_L;
% - For all these i-cells I will have a fraction of B(i) staying at x2_all(i) (called f(i)) 
% and a fraction of B(i) to be reassigned to x2_all(i-1) (called fm(i)) 
f(i) = (vAV(i)-x2_all(i-1))./(x2_all(i)-x2_all(i-1));
fm(i) = 1-f(i); %(vAV(i)-x2_all(i))./(x2_all(i-1)-x2_all(i));
in01(f(i))
in01(fm(i))
if ~isequal((f(i)+fm(i)),ones(size(f(i))))
 ME = MException('Fractions are calculated wrongly');
 throw(ME);
end
% If I have the vAV in the first cell below the first reference, I need to
% remove a number of particles from the first reference which is consistent
% with the volume of this fraction of particles


% - Here I need to sum up
B_CA_aSp(ind_L) = B_CA_aSp(ind_L)+B(ind_L).*f(ind_L);

% - Here I need to sum up
B_CA_aSp(ind_L-1) = B_CA_aSp(ind_L-1)+B(ind_L).*fm(ind_L);

ind_union = sort(cat(1,ind_E,ind_U,ind_L));
if ~isequal(ind_union,unique(ind_union)) || ~isequal(sort(cat(1,ind_E,ind_U,ind_L)),ind_nn0)
 ME = MException('These should be disjucted sets');
 throw(ME)
end

integralOfTheMarg_CA = computeMarg(B_CA_aSp,M_NisTilde,[],aStp.x2_widths,[],'x1'); % with this I obtain the marginal for x1 Marg_x1 (no tilde)
% if ~zeroth_out_vAVs && ~isequal(integralOfTheMarg,integralOfTheMarg_CA)
% ME = MException('If I dont make the out_vAVs to keep consistency of the first moment then the zeroth moment should not out_vAV');
% throw(ME);
% end
return

function [] = in01(v)
 if ~(isempty(v(v<0)) && isempty(v(v>1)))
    ME = MException('fractions are only defined between 0 and 1');
    throw(ME);
 end
return
% 
% function [] = checkFirstMomentConserved(B,V,x2_all,t)
% global noGhostInd lostFirstMoment_lb lostFirstMoment_ub tFlaglb tFlagub approx_lb_max approx_ub_max
% %     If these conditions are not met I am loosing conservation of first moment after CAT
% %     reassignement (too many particles appearing below the first reference or above the last)
% vAV_lb = (V(1)+V((1)))/(B(1)+B((1)));% average volume of particles born in the (ghost1 + first) 
% vAV_ub = (V(nCells)+V(end))/(B(nCells)+B(end));% average volume of particles born in the (end + ghostEnd) 
% approx_lb = x2_all((1))-vAV_lb;
% approx_ub = vAV_ub-x2_all(nCells);
% approx_lb_max = max(approx_lb,approx_lb_max);
% approx_ub_max = max(approx_ub,approx_ub_max);
% flaglb = (lostFirstMoment_lb==false)&&(approx_lb>0);
% flagub = (lostFirstMoment_ub==false)&&(approx_ub>0);
% if flaglb
%     tFlaglb = t;
%     lostFirstMoment_lb = true;   
% end
% if flagub
%     tFlagub = t;
%     lostFirstMoment_ub = true;   
% end
% return
