function P = settingsASM1seq(P)
global myPlotOff runs gcs ltys LatinHS nnsos

% - what model running
P.isRunning = 'ASM1seq';
P.parallelMode = 'parfor';
P.runningParallel = false;

% - Save to file;
P.writeASM1seqtoFile = true;

% -- Plotting
myPlotOff = true;
P.figHET = figure('visible','off');
P.figAUT = figure('visible','off');

% -- Initializing reactor control counters and flags
P.phaseCtrl = 1;P.timeLastCycleStart=0;P.repeat234=1; P.allChecks=1; P.step=0; P.indSwitchTimes=1;

% -- Initializing reactor control counters and flags
runs = true; % Refers to the variable runASM1seq in runASM1seqModel_ASM1.m 
nnsos = [1,2];
gcs = {'g','c'}; %{'grams','concs'};
ltys = {'-','o'};

% -- Flags to spot negativities 
P.thr = 0;      % When variables go below this thr, they are reset to thr within odeFunASM1_whole
P.minNeg = 0;

% -- Spot numerical inaccuracies
P.epsyFeed = 0;
P.epsyPurge = 0;
return
