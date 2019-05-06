function [x_short] = short(x_long)
global load1

x_short=x_long;
x_short(load1.dofCC) = [];
 