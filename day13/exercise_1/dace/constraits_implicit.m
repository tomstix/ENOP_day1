function [c,ceq]=constraits(x)
ceq=[];
c(1)=-(x(1)-1.5);
c(2)=-(3-x(1));
c(3)=-(x(2)-0.5);
c(4)=-(1.5-x(2));