function save_figure(h,fname)
% SAVE_FIGURE   Save a figure both as 'fig' file and as 'eps'.
%
% Yoel Shkolnisky, December 2013.

if nargin==1    
    fname=h;
    h=gcf;
end

[pathstr,name,ext]=fileparts(fname);
figname=sprintf('%s.%s',[pathstr filesep name],'fig');
epsname=sprintf('%s.%s',[pathstr filesep name],'eps');

hgsave(h,figname);
print(h,epsname,'-depsc');