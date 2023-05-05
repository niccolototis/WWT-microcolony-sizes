function P = setOperationParams(P)
global  T_sim iT  cycleModStr parametersEffectOnMSD gsToKgs

if iT>0 && parametersEffectOnMSD
    P=changeImprovements(P);
else
    P.QIN=210; % Influent Flow rate [m^3/d]  
    P.SRT=20; % IZICO Sludge Retention Time [d]
    P.O2_c_ref=1.5; % Reference 02 concentration maintainted in the reactor [mg/L=g/m^3]
    %
    cycleModStr='';
    % NB Purge window is considered as the 10 final minutes of the settling phase
    P.phaseSeq_1C=[1 2 3 4 2 3 4 2 3 4 5 6 7]; % refers to the index of phaseNames
    P.nPhases_1C=length(P.phaseSeq_1C);
    P.phaseLenght_1C=[50,15,25,117,15,25,117,15,25,118,130+3,10,55]; % [min]
end

% -- constants
P.hoursOneDay=24;
P.minsOneHour=60; 
P.kgsToGrams=1000;
gsToKgs=1/1000;

% -- Operation parameters
P.V_0=1000; % reactor volume at IZICO=1000 m^3  
P.HRT=P.V_0/P.QIN; % Hidraulic Retention Time [d]
P.sludgeDensity=1010*1000; % [mg/L] , used to compute the volume removed during the purge phase

% -- Phase description, sequence and timings
P.phaseNames={'onAer','','onFeed','onAer','onSettling','onPurge','onEfflux'}; 
if sum(P.phaseLenght_1C)~=12*60
    error('The cycle should be 12 hours long')
end

%
wherethisPhase=find(strcmp(P.phaseNames,'onFeed'));
lengthsThisPhase=find(P.phaseSeq_1C==wherethisPhase);
influxMins_1C=sum(P.phaseLenght_1C(lengthsThisPhase)); % Feeding phase lasts 25 min, which is repeated 3 times per cycle [min/cycle]
%  % [d/cycle]
wherethisPhase=find(strcmp(P.phaseNames,'onEfflux'));
lengthsThisPhase=find(P.phaseSeq_1C==wherethisPhase);
effluxMins_1C=sum(P.phaseLenght_1C(lengthsThisPhase)); % lenght of Efflux phase [min/cycle]
% 
wherethisPhase=find(strcmp(P.phaseNames,'onPurge'));
lengthsThisPhase=find(P.phaseSeq_1C==wherethisPhase);
purgeMins_1C=sum(P.phaseLenght_1C(lengthsThisPhase));% lenght of Purge phase [min/cycle]
% 
P.influx_tDays=influxMins_1C/P.minsOneHour/P.hoursOneDay; % [d/cycle]
P.efflux_tDays=effluxMins_1C/P.minsOneHour/P.hoursOneDay; % [d/cycle]
P.purge_tDays=purgeMins_1C/P.minsOneHour/P.hoursOneDay; % [d/cycle]
%
P.phaseLenght_1C=P.phaseLenght_1C/P.minsOneHour/P.hoursOneDay; % [d]
P.switchTimes_1C=cumsum(P.phaseLenght_1C); % [d]
P.CycleLength=P.switchTimes_1C(end); % [d/cycle]
P.cyclesOneDay=1/P.CycleLength; % [cycle/d]
%
P.W1=P.influx_tDays*P.cyclesOneDay; % [(d/cycle)*(cycles/day)=dimensionless] Fraction of 1 d actually used for influx
P.W2=P.efflux_tDays*P.cyclesOneDay; % [(d/cycle)*(cycles/day)=dimensionless] Fraction of 1 d actually used for efflux
P.W3=P.purge_tDays*P.cyclesOneDay; % [(d/cycle)*(cycles/day)=dimensionless] Fraction of 1 d actually used for purge
P.QInfl=P.QIN/P.W1; 
% This is the continuous flow rate sustained in the time window W1 to match
% QIN, the flow rate defined as if the flow was continuous thorgh the whole
% day
P.V_in_1Cy=P.QIN*P.CycleLength; % Volume influx per cycle [m^3/cycle]
switch P.modelWaterSludgeSegreg 
    case 'segregated'
        % do nothing, QEffl and QPurge flow rates depend each time from the amount of
        % sludge actually present in the system
    case 'stirred' 
        % Purge happens together with the removal of the effluent, in phase 7, not in phase 6
        P.QEffl=P.QIN/P.W2; % QOUT=QIN
end


% -- Time span of the simulation
tMins=0.1; % time step in minutes [min]
P.nPointsMin=1/tMins; % This is later used to thin out the data before plotting
P.dtWholebag=tMins/P.hoursOneDay/P.minsOneHour;
P.tspan_all_ASM1seq=0:P.dtWholebag:T_sim;

% -- Tspan is cut into integration time intervals, one for each phase
P.switchTimes_all=0; % Need the zero time point!
P.phaseSeq_all=0;
intSegmEnd=0;
while intSegmEnd<T_sim && ~reallyEqual(intSegmEnd,T_sim) 
    P.switchTimes_all=[P.switchTimes_all intSegmEnd+P.switchTimes_1C];
    P.phaseSeq_all=[P.phaseSeq_all P.phaseSeq_1C];
    intSegmEnd=P.switchTimes_all(end); %take last element to stop when T_sim is reached
end

% If by adding cycles I go past the T_sim, I resize switchTimes_all
final=find(strcmp(num2str(T_sim),(string(num2cell(P.switchTimes_all)))));
if length(final)>1
    final=min(final);
end
P.switchTimes_all=P.switchTimes_all(1:final);
P.switchTimes_all(1,end)=T_sim;
P.phaseSeq_all=P.phaseSeq_all(1:final);

% - define time intervals/fragments and store them in structure tspans_ASM1seq
P.nPhases_all=final-1;
P.tFragm_lengths=zeros(1,P.nPhases_all);
P.tspallin=[];
for nf = 1 : P.nPhases_all     
    P.tspans_ASM1seq.(strcat(['f',num2str(nf)]))=P.switchTimes_all(nf):P.dtWholebag:P.switchTimes_all(nf+1);
    P.tspallin=[P.tspallin P.tspans_ASM1seq.(strcat(['f',num2str(nf)]))];
    P.tFragm_lengths(nf)=length(P.tspans_ASM1seq.(strcat(['f',num2str(nf)])));
end
P.tFragm_lengths=[0, cumsum(P.tFragm_lengths)]; % Adding the zero is important for the main
P.ntpts_ASM1seq=P.tFragm_lengths(end);
return

function P=changeImprovements(P)
global cycleModified cycleModStr iT
    % overwrite the values just defined with the ones defined in
    % paramsToModify.m
    fN=fieldnames(P.impro_chosen);
    for q=1:length(fN)
        P.(fN{q})=P.impro_chosen.(fN{q})(iT);
    end

    if cycleModified
        cycleModStr=['_cycle_v' num2str(P.cycleVariant)];
    else
        cycleModStr='';
    end
    switch P.cycleVariant
        case 0
            cycleModStr='';
            % NB Purge window is considered as the 10 final minutes of the settling phase
            P.phaseSeq_1C=[1 2 3 4 2 3 4 2 3 4 5 6 7]; % refers to the index of phaseNames
            P.nPhases_1C=length(P.phaseSeq_1C);
            P.phaseLenght_1C=[50,15,25,117,15,25,117,15,25,118,130+3,10,55]; % [min]
        case 1 
            % only two feeding phases instead of 3        
            P.phaseSeq_1C=[1 2 3 4 2 3 4 2 4 5 6 7]; 
            P.nPhases_1C=length(P.phaseSeq_1C);
            P.phaseLenght_1C=[50,15,37,117,15,38,117,15,118,130+3,10,55]; 
        case 2
            % larger feeding phases, shorter aeration phases
            P.phaseSeq_1C=[1 2 3 4 2 3 4 2 3 4 5 6 7]; 
            P.nPhases_1C=length(P.phaseSeq_1C);
            P.phaseLenght_1C=[50,15,25+15,117-15,15,25+15,117-15,15,25+15,118-15,130+3,10,55]; 
        case 3
            % larger feeding phases, much shorter aeration phases, more nothing
            P.phaseSeq_1C=[1 2 3 4 2 3 4 2 3 4 5 6 7]; 
            P.nPhases_1C=length(P.phaseSeq_1C);
            P.phaseLenght_1C=[50,15+40,25+15,117-15-40,15+40,25+15,117-15-40,15+40,25+15,118-15-40,130+3,10,55]; 
        case 4
            % larger feeding phases, Minimal aeration phases, more nothing
            P.phaseSeq_1C=[1 2 3 4 2 3 4 2 3 4 5 6 7]; 
            P.nPhases_1C=length(P.phaseSeq_1C);
            P.phaseLenght_1C=[50,110,25+15,7,110,25+15,7,110,25+15,8,130+3,10,55]; 
        case 5
            % larger feeding phases, NO aeration phases, more nothing
            P.phaseSeq_1C=[1 2 3 2 3 2 3 5 6 7]; 
            P.nPhases_1C=length(P.phaseSeq_1C);
            P.phaseLenght_1C=[50,117,25+15,117,25+15,118,25+15,130+3,10,55]; 
    end      
return


% NB in odeFunASM1_whole I use dilutions, defined in this way:
% dilInfl=Q1/V_0; % for the influx --> HRT/W1 * this is what I am doing combining the previous
% dilEffl=Q2/V_0; % for the efflux --> HRT/W2
% dilPurge=Q3/V_0; % for the purge --> SRT/W3
