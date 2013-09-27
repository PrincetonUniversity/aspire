function numdrawn=ImagicDisplay(ImageStack, scale, startlabel, MenusOn)
% function numdrawn=ImagicDisplay(ImageStack,scale, startlabel, MenusOn)
% Make an EMAN-like display of the image stack in the current
% figure window (or create one if no figures exist). This version does not
% do automatic resizing.  Only the first
% argument is required. startlabel gives the starting number labeling
% an image (default is 1).  Use MenusOn=1 (default) if you
% want to have the ability to print or save the figure.  The returned
% variable tells how many images were drawn when the routine was first
% called; however resizing works as long as the window exists.

if nargin<4
    MenusOn=0;
end;
if nargin<3
    startlabel=1;
end;
offset=startlabel-1;
LabelsOn=(offset>=0);
if nargin<2
    scale=1;
end;
h=gcf;
SetGrayscale;
set(h,'Units','pixels','Toolbar','none');  % no redrawing
if MenusOn
    set(h,'MenuBar','figure');
end;
% set(h,'color','k');
index=0;
Redraw(0,0);
numdrawn=index-1;

    function Redraw(src,evt)
        h0=gcf;
%         set(h0,'color','k');
        clf;
        %     iters=iters+1
        P=get(h0,'Position');
        xsize=P(3);
        ysize=P(4);
        [nx ny ni]=size(ImageStack);
        nx=nx*scale;
        ny=ny*scale;

        mx=floor((xsize-1)/(nx+1));  % Number of panels horizontally.
        my=floor((ysize-1)/(ny+1));  % We leave 1 pixel of border everywher.

        index=1;
        stx=2;  % Starting X value
        sty=ysize-ny;  % Starting y value.
        for iy=1:my
            for ix=1:mx
                if index>ni
                    break;
                end;
                x=stx+(ix-1)*(nx+1);
                y=sty+(1-iy)*(ny+1);
                P1=[x y nx ny];
                h(index)=axes('Units','pixels','Position',P1,'Visible','off');
                imacs(ImageStack(:,:,index));
                axis off;
                if LabelsOn
                    text(x+nx,y,num2str(index+offset),'units','pixels','Position',[nx-2 1],...
                    'horizontalalignment','right','VerticalAlignment','bottom',...
                    'color','w','BackgroundColor',[.2 0 .2],'FontSize',9);
                end;
                index=index+1;
            end;
            if index>ni
                break;
            end;
        end;
    end % Redraw

end % ImagicDisplay
