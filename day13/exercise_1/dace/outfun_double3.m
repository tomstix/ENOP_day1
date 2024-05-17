function stop = outfun_double3(x, optimValues, state, varargin)
   stop = false;
t=optimValues.fval
k=optimValues.constrviolation
     if optimValues.fval+optimValues.constrviolation<= 0.3162
%    if optimValues.fval<= 0.3162
     stop = true;
   end 
