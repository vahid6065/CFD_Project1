function [Re] = Reynolds_Number(Reynolds)

if Reynolds == 1
    Re = 1;
elseif Reynolds == 2
    Re = 10;
elseif Reynolds == 3
    Re = 100;
elseif Reynolds == 4
    Re = 500;
end
end
