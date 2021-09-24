% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Copyright 2021, Liheng Zheng
%
% This file is part of UBER.
%
%    UBER is free software: you can redistribute it and/or modify it under the
%    terms of the MIT License as published by Massachusetts Institute of
%    Technology. UBER is distributed in the hope that it will be useful, but
%    WITHOUT ANY WARRANTY, without even the implied warranty of MERCHANTABILITY or
%    FITNESS FOR A PARTICULAR PURPOSE. See the MIT License for more details.
%
%    You should have received a copy of the MIT License along with UBER. If not,
%    see <https://opensource.org/licenses/MIT>.
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

% this script plots the UBER solutions for Example2
clear all;

datfile = './example2_out.dat';
ostr = uber_read_output_file(datfile);

indices = [1, 2, 3, 5, 8, 13, 21];

nt = 7;
tt = ostr.tt(indices);
nx = ostr.nx1;
xx = ostr.xx1;
sol = ostr.sol(:, indices);

canvas = [18, 15];
fig = figure(1);
clf;
clcarr = lines(nt);
for i=1:nt
   plot(xx, sol(:,i), 'linestyle', '-', 'color', clcarr(i,:) ,'linewidth', 1.5);
   hold on;
end
xlabel('x');
ylabel('f(x)');

% plot colorbar
colormap(clcarr);
hcb = colorbar('eastoutside','direction','reverse');
cticklabels = {};
for i=1:nt
   cticklabels{i} = num2str(tt(i));
end
set(hcb,'ticks',linspace(0,1,nt),'ticklabels',cticklabels);
hcb.Label.String = 'Time';
setupfig(fig, canvas);
printfig(fig, 'example2.png', canvas);




