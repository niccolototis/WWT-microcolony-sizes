function B_CA_aSp = computeRearrangeCAT(B,V,aStp,t,ghost1,ghostEnd) 
global noGhostInd 
x2_all = aStp.x2_all;
[B_CA_aSp, am_vm,am_v,ap_v,ap_vp] = deal(zeros(size(x2_all)));
checkFirstMomentConserved(B,V,x2_all,t)
vAV = V./B;
am_vm(noGhostInd) = (vAV(noGhostInd-1)-x2_all(noGhostInd-1))./(x2_all(noGhostInd)-x2_all(noGhostInd-1)); % alpha_minus(vAV_i_minus) == portion of vAV_i-1 which gets assigned to v(i) from below (vAV_i-1 goes to v(i) from minus)
am_v(noGhostInd) = (vAV(noGhostInd)-x2_all(noGhostInd-1))./(x2_all(noGhostInd)-x2_all(noGhostInd-1)); % alpha_minus(vAV_i) == portion of vAV_i which gets assigned to v(i) from below (vAV_i goes to v(i) from minus)
ap_v(noGhostInd) = (vAV(noGhostInd)-x2_all(noGhostInd+1))./(x2_all(noGhostInd)-x2_all(noGhostInd+1)); % alpha_plus(vAV_i) == portion of vAV_i which gets assigned to v(i) from above (vAV_i goes to v(i) from plus)  
ap_vp(noGhostInd) = (vAV(noGhostInd+1)-x2_all(noGhostInd+1))./(x2_all(noGhostInd)-x2_all(noGhostInd+1)); % alpha_plus(vAV_i_plus) == portion of vAV_i+1 which gets assigned to v(i) from above (vAV_i+1 goes to v(i) from plus)   
%   I further need to calculate: 
%   ap_v and ap_vp for ghost1, >>> reassigned to B_CA_aSp(noGhostInd(1))
%   am_vm e am_v for ghostEnd, >>> reassigned to B_CA_aSp(noGhostInd(end))
%   this assures the conservation of zeroth moment
ap_vp(1) = (vAV(1+1)-x2_all(1+1))./(x2_all(1)-x2_all(1+1)); % proportion of births of the first cell which would fall below the first lower bound
ap_v(1) = (vAV(1)-x2_all(1+1))./(x2_all(1)-x2_all(1+1)); % I am sure that particles cannot appear below x2_all(1)        
am_v(end) = (vAV(end)-x2_all(end-1))./(x2_all(end)-x2_all(end-1)); % proportion of births of the last cell which would fall above the max upper bound
am_vm(end) = (vAV(end-1)-x2_all(end-1))./(x2_all(end)-x2_all(end-1)); % I am sure that particles cannot appear above x2_all(end)         
%   disp([am_vm am_v ap_v ap_vp])

%   IMPO need to set also all the negative values that might be created to zero! these are due to the ghost
%   cells added. I want to keep fractions referring to ghost cells only when they are positive
am_vm(isnan(am_vm)) = 0; 
am_v(isnan(am_v)) = 0;
ap_v(isnan(ap_v)) = 0;
ap_vp(isnan(ap_vp)) = 0;    
vAV(isnan(vAV)) = 0; % do this only AFTER the computation of a_minus and a_plus, because having zeros before would cause a_minus and a_plus to have abnormal nonzero values

%   check and check2 should be all 1 for positions noGhostInd
check = ap_v(noGhostInd(1):noGhostInd(end)) + am_vm(noGhostInd(2):noGhostInd(end)+1);
check2 = ap_vp((noGhostInd(1)-1):noGhostInd(end-1)) + am_v(noGhostInd(1):noGhostInd(end));

for j = noGhostInd 
 c_extra = 0;
 c1 = B(j-1).*am_vm(j)*heaviside(vAV(j-1)-x2_all(j-1));
 c2 = B(j)*am_v(j)*heaviside(x2_all(j)-vAV(j));
 c3 = B(j)*ap_v(j)*heaviside(vAV(j)-x2_all(j));
 c4 = B(j+1)*ap_vp(j)*heaviside(x2_all(j+1)-vAV(j+1));

 switch j
    case noGhostInd(1)
        c1 = 0;
     case noGhostInd(end) 
        c4 = 0;
 end

%         switch j
%             case noGhostInd(1)
%                 JJ = j-1; % calculate the component of births in the first cell that CAT would reassign to ghost1
%                 c4ghost1 = B(JJ+1)*ap_vp(JJ)*heaviside(x2_all(JJ+1)-vAV(JJ+1));
%                 c_extra = B(JJ)+c4ghost1; % All births by breakage B(JJ) appearing in the ghost1 cell have to be reassigned to the first normal cell
%             case noGhostInd(end) 
%                 JJ = j+1; % calculate c1 and c2 components for ghostEnd
%                 c1extra = B(JJ-1).*am_vm(JJ)*heaviside(vAV(JJ-1)-x2_all(JJ-1));
% %                 c2extra = B(JJ)*am_v(JJ)*heaviside(x2_all(JJ)-vAV(JJ)); % this will actually be used ONLY IF THERE IS AGGREGATION, because growth does not shift particles across cell boundaries!
%                 c_extra = c1extra+B(JJ); % The B(JJ) is >0 only if I have birth by aggregation of particles greater than the maximum size. All of these births reassigned to the last normal cell
%         end
 if (B(noGhostInd(1)-1)>0&&(heaviside(vAV(noGhostInd(1)-1)-x2_all(noGhostInd(1)-1))~= 1))||(B(noGhostInd(end)+1)>0&&((heaviside(x2_all(noGhostInd(end)+1)-vAV(noGhostInd(end)+1))~=1))) % If I have births in ghost cells I wanna be sure they are reassigned
    ME = MException('Ghost cells should be set so that these always equal 1');
    throw(ME);
 end
 if (c1<0 || c2<0 || c3<0 || c4<0 || c_extra<0)
    ME = MException('negative values of alphas coefficients should be canceled by heaviside functions');
    throw(ME)
 end
 B_CA_aSp(j) = c1+c2+c3+c4+c_extra;       
end 
return

function [] = checkFirstMomentConserved(B,V,x2_all,t)
 global noGhostInd lostFirstMoment_lb lostFirstMoment_ub tFlaglb tFlagub approx_lb_max approx_ub_max
%     If these conditions are not met I am loosing conservation of first moment after CAT
%     reassignement (too many particles appearing below the first reference or above the last)
 vAV_lb = (V(1)+V(noGhostInd(1)))/(B(1)+B(noGhostInd(1)));% average volume of particles born in the (ghost1 + first) 
 vAV_ub = (V(noGhostInd(end))+V(end))/(B(noGhostInd(end))+B(end));% average volume of particles born in the (end + ghostEnd) 
 approx_lb = x2_all(noGhostInd(1))-vAV_lb;
 approx_ub = vAV_ub-x2_all(noGhostInd(end));
 approx_lb_max = max(approx_lb,approx_lb_max);
 approx_ub_max = max(approx_ub,approx_ub_max);
 flaglb = (lostFirstMoment_lb==false)&&(approx_lb>0);
 flagub = (lostFirstMoment_ub==false)&&(approx_ub>0);
 if flaglb
    tFlaglb = t;
    lostFirstMoment_lb = true;   
 end
 if flagub
    tFlagub = t;
    lostFirstMoment_ub = true;   
 end
return
