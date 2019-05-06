function [x_long] = long(x_short)
global load1

x_long=zeros(length(load1.dofFree) +  ...
             length(load1.dofCC),1);
x_long(load1.dofFree) = x_short;
x_long(load1.dofCC) = load1.fixedvalues;
 