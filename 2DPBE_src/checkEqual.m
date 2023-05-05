function [] = checkEqual(A,B,str)      
if ~strcmp(num2str(A),num2str(B))
 ME = MException(strcat([str,' should be equal']));
 throw(ME) 
end 
end
