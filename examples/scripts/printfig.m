function [] = printfig(varargin)
% Syntax:
%    printfig(varargin)

switch nargin
   case 0
      fid = gcf;
      figname = 'printfig_tmp.eps';
   case 1
      fid = gcf;
      figname = varargin{1};
   case 2
      fid = gcf;
      figname = varargin{1};
      dimarr = varargin{2};
   otherwise
      fid = varargin{1};
      figname = varargin{2};
      dimarr = varargin{3};
end

set(fid,'resize','off')
set(fid,'paperpositionmode','auto');

if( nargin>1 )
   set(fid,'paperunits','centimeters');
   set(fid,'paperposition',[1 1 dimarr]);
   set(fid,'papersize',dimarr);
end

if( strcmp(figname(end-3:end),'.eps') )
   print(fid,'-depsc','-painters',figname);
elseif( strcmp(figname(end-3:end),'.pdf') )
   figname2 = [figname(1:end-4),'.eps'];
   print(fid,'-depsc','-painters',figname2);
elseif( strcmp(figname(end-3:end),'.png') )
   set(findall(fid,'type','axes'),'xtickmode','manual','ytickmode','manual');
   saveas(fid,figname,'png');
else
end

return
end % function printfig
