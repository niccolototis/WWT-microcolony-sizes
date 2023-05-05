function [varargout] = reallyEqual(A,B,varargin)
%type1 output: [true] = reallyEqual(A,B) if A and B are approxequal
%          [true,A,B] = reallyEqual(A,B,'~','error') if A and B are approxequal
%          [false] = reallyEqual(A,B) if A and B are different
%          error = reallyEqual(A,B,'error') if A and B are different
%          error = reallyEqual(A,B,'~','error') if A and B are different
%type2 output: [true] = reallyEqual(A,B,'=') if A and B are identical
%          [true,false] = reallyEqual(A,B,'=') if A and B are approxequal not identical
%          false = reallyEqual(A,B,'=') if A and B are different
%          error = reallyEqual(A,B,'=','error') if A and B are different
%type1 output: [true,true] = reallyEqual(A,B,'=','~','error') if A and B are identical
%           [true,false] = reallyEqual(A,B,'=','error') if A and B are approxequal not identical
%           [true,false,A,B] = reallyEqual(A,B,'=','~','error') if A and B are approxequal not identical
%           error = reallyEqual(A,B,'=','error') if A and B are different
% when I ask for forceMatch and (A,B) are approxequal, then A is adjusted to the value of B: A<-B 
checkEqual = false; forceMatch=false; incorrectI=false; MEout=false; MEstr=[]; additNargin=nargin-2;
thr = 1e-6;
if additNargin>3
 incorrectI = true;
else
switch additNargin
 case 1
    switch varargin{1}
        case ' = '
            checkEqual = true;
        case '~'
            forceMatch = true;
        otherwise
            MEout = true;
            MEstr = varargin{1};
    end
 case 2
    switch varargin{1}
        case ' = '
            checkEqual = true;
        case '~'
            forceMatch = true;
        otherwise
            incorrectI = true;
    end
    switch varargin{2}
        case ' = '
            incorrectI = true; % I want '=' always at the first position
        case '~'
            forceMatch = true;
        otherwise
            MEout = true;
            MEstr = varargin{2};
    end
 case 3
    switch varargin{1}
        case ' = '
            checkEqual = true;
        case '~'
            forceMatch = true;
        otherwise
            incorrectI = true;
    end
    switch varargin{2}
        case ' = '
            checkEqual = true;
        case '~'
            forceMatch = true;
        otherwise
            incorrectI = true;
    end
    MEout = true;
    MEstr = varargin{3};
end
if MEout && ~ischar(MEstr)
 incorrectI = true;
end
if incorrectI
 error('Incorrect input');
end

% ---

if checkEqual
 if isequal(A,B) 
    varargout{1} = true;
 else
    varargout{2} = false;        
 end
end

if strcmp(string(A),string(B)) 
 if abs(A-B)<thr
    varargout{1} = true;
    A = B;
 else
    error(strcat(['Strings the same but difference is larger than threshold thr = ',num2str(thr)]));
 end
else
 varargout{1} = false;
 if forceMatch
    error('Values are different');
 end
 if MEout
    error(MEstr);
    throw(ME);  
 end % else do nothing
end

if forceMatch
 if checkEqual
    varargout{3} = A;
    varargout{4} = B;
 else
    varargout{2} = A;
    varargout{3} = B;
 end
end
end
