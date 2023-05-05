function [] = checkEqualStr_ForLoop(A,B,origin,P)
global epsyBreak epsyGrow nInconsGrowth nInconsBreak
if ~P.identify_lambda && ~P.sweep_lambda
 if ~isequal(num2str(A),num2str(B)) 
    switch origin
        case 'growth'
            if abs(A-B)>epsyGrow
                epsyGrow = abs(A-B);
            end
            nInconsGrowth = nInconsGrowth+1;
        case 'breakage'
            if abs(A-B)>epsyBreak
                epsyBreak = abs(A-B);
            end
            nInconsBreak = nInconsBreak+1;
    end
%     disp( [num2str(A),num2str(B)]);
%     ME = MException('`B and D by growth should coincide see paper: Kumar, Jitendra, et al. "An efficient numerical technique for solving population balance equation involving aggregation, breakage, growth and nucleation." Powder Technology (2008)');
%     throw(ME)
 end
end
return
