function A = setNNG(A)
 if ~isempty(A(A<0)>0) % they should sum up to one aprt from where I have no births
        ME = MException('Should never be negative !');
        throw(ME)
 end
 A(A<= realmin)=0;
return
