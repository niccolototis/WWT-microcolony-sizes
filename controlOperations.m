function P = controlOperations(P)
% . @@@@ NB all the control is done in minutes!! @@@@. Outside this function, all
% the rest is in hours
P.onAer = false; P.onFeed=false; P.onSettling=false; P.onEfflux=false; P.onPurge=false; %baseline
P.pastPhase = P.phaseNames{P.phaseCtrl};
P.phaseCtrl = P.phaseSeq_all(P.indPhase+1);
% phases 2,3,4 have to be repeated 3 times per 10h cycle
switch P.phaseCtrl
%     NB I am using the condition with greater than because otherwise I am
%     imprecise and I tend to miss switches
 case 1
        P.onAer = true;  % phase 1: aeration  
        P.thisPhase = 'onAer';
 case 2 
        % do nothing 
        P.thisPhase = '';
 case 3
        P.onFeed = true;% phase 3: feed
        P.thisPhase = 'onFeed';
 case 4
        P.onAer = true;% phase 4: aeration
        P.thisPhase = 'onAer';
 case 5
        P.onSettling = true;% phase 5: settling does nothing actually
        P.thisPhase = 'onSettling';
 case 6
        P.onPurge = true;
        P.thisPhase = 'onPurge';
 case 7
        P.onEfflux = true;% phase 6: remove effluent
        P.thisPhase = 'onEfflux';
end
end
