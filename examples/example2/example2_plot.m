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




