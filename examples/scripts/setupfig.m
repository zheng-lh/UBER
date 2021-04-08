function [] = setupfig(varargin)
% Syntax:
%    setupfig(varargin)

switch nargin
   case 0
      fid = gcf;
   case 1
      fid = gcf;
      dimarr = varargin{1};
   otherwise
      fid = varargin{1};
      dimarr = varargin{2};
end

if( nargin>0 )
   set(fid,'units','centimeters');
   set(fid,'position',[1 1 dimarr]);
end

set(findall(fid,'Type','text'), 'fontsize', 20);
set(findall(fid,'Type','axes'), 'fontsize', 20, 'tickdir', 'out', 'xminortick',...
'on', 'yminortick', 'on');

if( nargin>2 )
   for i=3:nargin
      thisvar = varargin{i};
      % e.g., set(findall(fid,'type','line'),'linewidth',2);
      set(findall(fid,'type',thisvar{1}),thisvar{2},thisvar{3});
   end
end

return
end % function setupfig
