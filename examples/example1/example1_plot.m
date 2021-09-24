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

% this script plots the UBER solutions for Example 1 against the reference solutions
% obtained using staggered grid finite difference.
clear all;

% load reference solution file
matfile = '../data/ex1_ref_solution.mat';
mymat = load(matfile);
refs = mymat.solstr;

% load UBER solution file
datfile = './example1_out.dat';
uber = uber_read_output_file(datfile);

% make plot
nt = uber.nt;
tt = uber.tt;

fig = figure(1);
canvas = [18, 15];
clf;
bottom = 1.d-10;
top = 1.d1;
clcarr = lines(nt);

h1 = semilogy(refs.xx, refs.uu(1,:), 'color', clcarr(1,:), 'linewidth', 2);
hold on;

for i=2:nt
   semilogy(refs.xx, refs.uu(i,:), 'color', clcarr(i,:), 'linestyle', '-.',...
            'linewidth', 2);
   logErrorPlot(uber.xx1, uber.sol(:,i)', uber.err(:,i)', bottom,...
      {'color', clcarr(i,:), 'marker', 'o', 'markersize', 8, 'linestyle', 'none',...
      'linewidth', 1.5});
end

dummy = 0.d0*refs.xx + 1.d3;
h2 = semilogy(refs.xx, dummy, 'color', clcarr(1,:), 'linestyle', '-.',...
              'linewidth', 2);
h3 = semilogy(refs.xx, dummy, 'color', clcarr(1,:), 'linestyle', 'none',...
              'marker', 'o', 'markersize', 8, 'linewidth', 1.5);
lid = legend([h1, h2, h3], 'Initial condition', 'Finite difference', 'UBER');
set(lid, 'box', 'off', 'location', 'southwest');
hold off;

xlim([0.d0, 1.d0]);
xticks([0.0:0.2:1.0]);
ylim([bottom, top]);
xlabel('r');
ylabel('f(r)');

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
printfig(fig, 'example1.png', canvas);


