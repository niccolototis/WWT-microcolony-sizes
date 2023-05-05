function [] = printInconsistencies(P)
global epsyBreak epsyGrow nInconsGrowth nInconsBreak lostFirstMoment_lb lostFirstMoment_ub approx_lb_max approx_ub_max out_vAV
global tFlaglb tFlagub minNeg epsylon
% Print inconsistencies found  
if ~P.identify_lambda && ~P.sweep_lambda

 disp(['epsyBreak:', num2str(epsyBreak),' epsyGrow:', num2str(epsyGrow),' nInconsGrowth:', num2str(nInconsGrowth),'  nInconsBreak:', num2str(nInconsBreak)])
 
 if lostFirstMoment_lb || lostFirstMoment_ub
    if lostFirstMoment_lb 
        disp(['@@@@@@@@@@ At t = ',num2str(tFlaglb),'  lost conservation first moment at lb @@@@@@@@@@@@@'])
        disp(['Max number of particles lost at lb per integration step: ', num2str(approx_lb_max)])
        disp('NB: lost of first moment at lower bound: the average volume of particles born at the lower bound is below the first grid point: these particles cannot be represented')
    end
    if lostFirstMoment_ub
        disp(['At t = ',num2str(tFlagub),'  lost conservation first moment at ub'])
        disp(['Max number of particles lost at ub per integration step: ', num2str(approx_ub_max)])
        disp('NB: lost of first moment at upper bound: the average volume of particles born at the upper bound is above the first grid point: these particles cannot be represented')
    end
    disp('Need a grid point closer to the bound. Suggestion: change how the grid is defined (try logarithmic) or change growth/breakage rates so less particles get so close to the bound')
 end
 disp(['### Max numerical inaccuracy introduced per cycle',num2str(epsylon)]); %In order to assess this you have only breakage no flowin no washout no growth
 disp(['### Minimal negative value re-adjusted to zero: ',num2str(minNeg)]);
 
 disp('')
 disp('')
 disp(strcat(['Average volumes vAVs were found outside the cell boundaries for a max value of ',num2str(out_vAV)]))
end
end
