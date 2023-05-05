function P = runMorris(P,comb)
% This function executes the Morris sensitvity analysis.
% INPUTS:
% - lb: lower bound on parameter values
% - ub: upper bound on parameter values
% - M: number of Morris trajectories generated
% - p: number of grid points
% - r: number of chosen Morris trajectories
% OUTPUTS:
% - mu: a vector containing all averages of the elementary effects of a
% parameter
% - sigma: a vector containing all standard deviations of the elementary 
% effects of a parameter

lb = P.EE.lb;
ub = P.EE.ub;

k = P.EE.k;
p = P.EE.p;
M = P.EE.M;
r = P.EE.r;

filename_X = strcat(['outputs/parameter_Morris/EE_X_k',num2str(P.EE.k),'_p',num2str(P.EE.p),'_M',num2str(P.EE.M),'_r',num2str(P.EE.r),'.mat']);
if isfile(filename_X)
 load(filename_X)
 disp('Retrieved previously generated Morris trajectories')
else

%% Generate required constants
% Create the stepsize of the elementary effects
delta = p/(2*(p-1));

% Create base values 
% @NT these are the ones that can be considered for the given delta
base = 0;
v = 1;
while v/(p-1) < 1-delta
 base = [base, v/(p-1)];
 v = v+1;
end

% Create a matrix B with the following properties:
% - size k+1 by k witk k the number of parameters
% - all elements are either 0 or 1
% - for every column index j, j = 1, , k, there are two rows of B that
% differ only in the jth entry

B = zeros(k+1,k);
for j = 2:size(B,1)
 B(j,1:j-1) = ones(1,size(1:j-1,2));
end

%% Step 1: Generate M trajectories

% X = zeros(M*(k+1),k);
X = zeros((k+1),k,M);

%  __lb/ub___
%  __________
% |\         \
% | \M        \
% |  \         \
% |   \ ____k___\
%  \   |        |
%   \  |k+1     |
%    \ |        |
%     \|________|

for i = 1:M
 % Generate a random point xi in the hypercube on the p-level grid
 xi = zeros(1,k);
 for z = 1:k
    ind = randperm(size(base,2));
    xi(z) = base(ind(1));
 end
 
 % Generate a diogonal matrix D with on its diagonal only 1 or -1 with
 % equal probability
 d = zeros(1,k);
 for z = 1:k
    element = -1+2*rand;
    if element < 0
        element = -1;
    elseif element >= 0
        element = 1;
    end
    d(z) = element;
 end
 D = diag(d);
 
 % Generate a random permutation matrix P
 Pmat = zeros(k,k);
 indx = randperm(k);
 for t = 1:length(indx)
    Pmat(t,indx(t)) = 1;
 end
 
 % Calculate a sampling matrix Bi
 J = ones(k+1,k);
 Bi = ones(k+1,1)*xi + ((delta/2).*((2.*B-J)*D+J))*Pmat;
 
 % Store the sampling matrix in X with size M*(k+1) by k @NT all
 % matrices == trajectories Bi are stored one below the other inside X
% X((i-1)*(k+1)+1:i*(k+1),:) = Bi;
 X(:,:,i) = Bi;

 end

%% Step 2: Choose r trajectories with the highest dispersion in x-space

% Calculate all pairwise distances between the trajectories
distances = zeros(M,M); 
for i = 1:M
 for j = i:M
    if i == j
        d = 0;
    else
        %@NT here scrolling down along X the distances are computed
        %only with Bj further down in X wrt Bi, because the others have already
        %been computed in previous rounds
        d = determinedistance(X(:,:,i),X(:,:,j)); 
    end
    % in the end it is a symmetric matrix
    distances(i,j) = d;
    distances(j,i) = d;
 end
end

% ﻿Campolongo et al. (2007) proposed an improvement of the sampling strategy just described that facilitates a better scanning of the input domain without increasing the number of model executions needed. The idea is to select the r trajectories in such a way as to maximize their spread in the
% 
% Choose the ones with the highest "spread Ds" (r trajectories are selected)
% comb = combnk(1:M,r); %@NT I consider all the possible combinations of the set of r trajectories inside the M generated
Ds = zeros(size(comb,1),1); %@NT Ds store the square distances for each combination comb. NB distance() calculates the distance among the M trajectories. Here instead, these are summed up to come up with distances of combos of trajectories
for i = 1:size(comb,1)
 comb1 = nchoosek(comb(i,:),2);
 for j = 1:size(comb1,1)
    Ds1 = distances(comb1(j,1),comb1(j,2)); % this takes the distances of single trajectories inside each combo
    Ds(i) = Ds(i) + Ds1^2;
 end
 Ds(i) = sqrt(Ds(i)); 
end
[~,ind] = sort(Ds,'descend');
indmax = ind(1); %@NT takes the index of the item with the max Ds
trajIndxs = comb(indmax,:); %@NT Run is a row index which contains the indexes of the r trajectories selected among the M
if P.EE.r~= length(trajIndxs)
 ME = MException('r should represent the number of trajectories');
 throw(ME)
end
X = X(:,:,trajIndxs); % @NT I am estracting the r trajectories out of the M generated that correspond to the set with the higher distance 
disp('Generated and chosen Morris trajectories')
save(filename_X,'X');
end
P.EE.nTraj = P.EE.r;

P.EE.X = X;
if size(ub,1)~= 1 || size(lb,1)~=1
 ME = MException('ub an lb need to be given as row vectors!');
 throw(ME);
end
P.EE.X_Params = X.*(ub-lb)+lb; 
P.EE.nModelRuns = (k+1)*P.EE.nTraj;

%  __________
% |\         \
% | \nTraj    \
% |  \         \
% |   \ ____k___\
%  \   |        |
%   \  |k+1     |
%    \ |        |
%     \|        |
if P.EE.confirmMorrisLikeLHS
P.EE.AllParCombRun = [];
for kk = 1:size(P.EE.X_Params,3)
 P.EE.AllParCombRun = [P.EE.AllParCombRun;P.EE.X_Params(:,:,kk)];
end
end

% %% Step 3: Calculate elementary effects
% 
% %@NT each row of EEtot contains the EEs for one parameter
% EEtot = zeros(k,length(runs));
% for t = 1:length(runs)
% Xi = X((runs(t)-1)*(k+1)+1:runs(t)*(k+1),:);
% EE = determineelementaryeffects(model,Xi,lb,ub);
% EEtot(:,t) = EE;
% end
% 
% %% Step 4: Calculate mu and sigma
% 
% % Note that the absolute value of the elementary effects is taken while
% % calculating mu!
% mu = zeros(size(EEtot,1),1);
% sigma = zeros(size(EEtot,1),1);
% for i = 1:size(EEtot,1)
% mu(i) = mean(abs(EEtot(i,:)));
% sigma(i) = std(EEtot(i,:));
% end

end
