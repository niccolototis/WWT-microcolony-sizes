function P = initializeSmatrix(P)
%Initialize the stoichiometric matrix
S = zeros(8,12); % zeros(8,13); leaving ALK out
%process 1: aerobic growth of heterotrophs
S(1,2) = -1/P.YH;
S(1,5) = 1;
S(1,8) = -(1-P.YH)/P.YH;
S(1,10) = -P.iXB;
% S(1,13) = -iXB/14; % leaving ALK out
%process 2: anoxic growth of heterotrophs
S(2,2) = -1/P.YH;
S(2,5) = 1;
S(2,9) = -(1-P.YH)/(2.86*P.YH);
S(2,10) = -P.iXB;
% S(2,13) = ((1-YH)/(14*2.86*YH)) - iXB/14; % leaving ALK out
%process 3: aerobic growth of autotrophs
S(3,6) = 1;
S(3,8) = -(4.57-P.YA)/P.YA;
S(3,9) = 1/P.YA;
S(3,10) = -P.iXB - (1/P.YA);
% S(3,13) = -(iXB/14) - (1/(7*YA)); % leaving ALK out
%process 4: decay of heterotrophs
S(4,4) = 1-P.fp; 
S(4,5) = -1;
S(4,7) = P.fp;
S(4,12) = P.iXB - P.fp*P.iXP; 
%process 5: decay of autotrophs
S(5,4) = 1-P.fp;
S(5,6) = -1;
S(5,7) = P.fp;
S(5,12) = P.iXB - P.fp*P.iXP;
%process 6: ammonification
S(6,10) = 1;
S(6,11) = -1; 
% S(6,13) = 1/14; % leaving ALK out
%process 7: hydrolysis of entrapped organics
S(7,2) = 1;
S(7,4) = -1; 
%process 8: hydrolysis of entrapped nitrogen
S(8,11) = 1;
S(8,12) = -1;
P.Smat = S;
return
