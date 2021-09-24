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

% this script plots the UBER solutions for Example3
clear all;

Po180 = pi/180.d0;
E = 0.304d0; % MeV

datfile = 'example3_out.dat';
ostr = uber_read_output_file(datfile);

nt = ostr.nt - 1;
tt = ostr.tt(2:end)/3600;
nx1 = ostr.nx1;
nx2 = ostr.nx2/2;
aa0 = ostr.xx1/Po180;
pphi = ostr.xx2(nx2+1:2*nx2)/Po180;
psd = ostr.sol(:, nx2+1:2*nx2, 2:end);
psd(psd<0.d0) = NaN;
lg_psd = log10(psd);

canvas = [18, 15];
for i=1:nt
   fig = figure(i);
   clf;
   colormap(jet(256));
   pcolor(pphi, aa0, lg_psd(:, :, i))
   shading flat;
   caxis([-10, -7.5]);
   cid = colorbar('eastoutside');
   cid.Label.String = 'log_{10}[psd (c^{3}cm^{-3}MeV^{-3})]';
   cid.Label.FontSize = 12;
   tstr = ['T = ',num2str(tt(i)),' hr'];
   title(tstr);
   ylim([30, 90]);
   xlabel('Geomagnetic longitude (deg)');
   ylabel('Equatorial pitch angle (deg)');
   setupfig(fig, canvas);
   printfig(fig, ['example3_T',num2str(tt(i)),'hr.png'], canvas);
end





