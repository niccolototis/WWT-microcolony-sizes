function [dist_ss_from_IC] = compare_and_sum(X_c_end,S_c_end,X_c_0,S_c_0)
 diffs_X = abs((X_c_end(3:end)-X_c_0(3:end))./X_c_0(3:end)); % Excluding the Het and AUt biomass
 diffs_X(isinf(diffs_X)) = 0;
 dist_X = sum(diffs_X);
 diffs_S = abs((S_c_end-S_c_0)./S_c_0); % Excluding the Het and AUt biomass
 diffs_S(isinf(diffs_S)) = 0;
 dist_S = sum(diffs_S);
 dist_ss_from_IC = dist_X+dist_S;
return
