clc
clear all
close all

Delta_1 = 0.2;

for i = [03 05 08 10]
    sec{i} = fprintf(['section' num2str(i, '%02.f') '_ref.dat']);
end
