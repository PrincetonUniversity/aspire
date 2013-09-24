function ImagicDisplay2(ImageStack, scale, ShowLabels, LabelOffset)
% function ImagicDisplay2(ImageStack,scale, ShowLabels, LabelOffset)
% Make an EMAN-like, resizable, scrollable display of the image stack in the current
% figure window (or create one if no figures exist).  Only the first
% argument is required. Scale is the magnification of the displayed image;
% 2 means each stack pixel is shown as 2x2 pixels on the screen.
% If ShowLabels=1 (default) only first one in each 
% row is labeled.  If ShowLabels=0, no labels are drawn (runs fastest);
% ShowLabels=2 causes every panel to be labeled.
% LabelOffset is an offset added to each panel number; by default this is
% zero, so the first panel is panel 1.  Resizing works as
% long as the window exists.
% Yunhui Liu, Aug 2011 based on an earlier version by fs

if nargin<4
    LabelOffset=0;
end;
if nargin<3
    ShowLabels=1;    % Labels only at the start of columns
end;
if nargin<2
    scale=1;
end;
h=gcf;         

SetGrayscale;

set(h,'Units','pixels','ResizeFcn',@redraw_boxes,'Toolbar','none');

% if MenusOn
%     set(h,'MenuBar','figure');
% end;
% set(h,'color','k');
% fsize=9;  % font size: this seems to be the smallest that is readable.
index=0;
%Redraw(0,0);
redraw_boxes(0,0);
numdrawn=index-1;


  function  redraw_boxes(src,evt)
    h0=gcf;
%     set(h0,'color','k');
    clf;
    P=get(h0,'Position');
    [boxsize, ~ ,num]=size(ImageStack);
    
    nboxsize= boxsize+1;            % add  1 pixel as border  
    showboxsize=nboxsize*scale;
%     fsize=max(7,min(10,round(showboxsize/4)));  % font size
fsize=9;  % smallest readable size
    win_width = P(3); win_height=P(4);
    col_max= floor(win_width/showboxsize);  % each line can contain as max as col_max number image
 
    if (mod(num,col_max)==0) 
        cols=col_max ;
    else
        cols=mod(num,col_max);
    end
    rows = ceil(num/col_max);

    bigmatrix =zeros(rows*nboxsize,col_max*nboxsize); % the whole image
    index = 1;

    % reconstruce the bigmatrix( reconstruct the small images as one image) 
    if (rows ==1 )
        for i=1:num
            bigmatrix(1:boxsize,(i-1)*nboxsize+1:i*nboxsize-1 )=rot90(ImageStack(:,:,i));
        end
    else
        for j=1: rows-1
            for i= 1: col_max
                bigmatrix((j-1)*nboxsize+1:j*nboxsize-1,(i-1)*nboxsize+1:i*nboxsize-1) = rot90(ImageStack(:,:,index));
                index =index+1;
            end
         end
 
        for i=1:num-(index-1)
            bigmatrix((rows-1)*nboxsize+1:rows*nboxsize-1,(i-1)*nboxsize+1:i*nboxsize-1) = rot90(ImageStack(:,:,index+i-1));
        end
    end
    
    apos1=[0 win_height-rows*showboxsize showboxsize*col_max showboxsize*rows];
 
    axes('Units','pixels','Position',apos1);
    set(gca, 'color', [0 0 0],'xtick',[],'ytick',[]); 

    hIm = imshow(bigmatrix,'InitialMagnification',100*scale,'DisplayRange',[min(min(bigmatrix)) max2d(bigmatrix)]);
 
    %%% label on
    if ShowLabels>0  % Show some labels
        index=1;
        if (rows==1)
            if (ShowLabels==1)     % only label the first one
                text(showboxsize-3,5,num2str(index+LabelOffset),'units','pixels',...
                      'horizontalalignment','right','VerticalAlignment','bottom',...
                    'color','w','BackgroundColor',[.2 0 .2],'FontSize',fsize);
            end    
            
            for j=1:num 
                if (ShowLabels>1)
                    text(showboxsize-3+(j-1)*showboxsize,5,num2str(index+LabelOffset),'units','pixels',...
                      'horizontalalignment','right','VerticalAlignment','bottom',...
                    'color','w','BackgroundColor',[.2 0 .2],'FontSize',fsize);
                end
                index=index+1;
                % here we need not use showboxsize,because the magnification factor
                 % in this figure takes effect on patch
                patch([(j-1)*nboxsize j*nboxsize j*nboxsize (j-1)*nboxsize],[0 0 nboxsize nboxsize],[0 0 0] ,'EdgeColor',[0 0 0],'FaceColor','none');
            end
        else
            for i=rows-1:-1:1
                if (ShowLabels==1)     % only label the first one
                    text(showboxsize-3 ,i*showboxsize+5 ,num2str(index+LabelOffset),'units','pixels',...
                        'horizontalalignment','right','VerticalAlignment','bottom',...
                        'color','w','BackgroundColor',[.2 0 .2],'FontSize',fsize);
                end

                for j= 1: col_max
                    if (ShowLabels>1)
                        text(showboxsize-3+(j-1)*showboxsize,i*showboxsize+5 ,num2str(index+LabelOffset),'units','pixels',...
                        'horizontalalignment','right','VerticalAlignment','bottom',...
                        'color','w','BackgroundColor',[.2 0 .2],'FontSize',fsize);
                    end
                    index =index+1;
                    patch([(j-1)*nboxsize j*nboxsize j*nboxsize (j-1)*nboxsize],[(i-1)*nboxsize (i-1)*nboxsize i*nboxsize i*nboxsize],[0 0 0] ,'EdgeColor',[0 0 0],'FaceColor','none');
                end
            end
            
            %% lable the last row
            if (ShowLabels==1)     % only label the first one
                text(showboxsize-3 , 5 ,num2str(index+LabelOffset),'units','pixels',...
                'horizontalalignment','right','VerticalAlignment','bottom',...
                'color','w','BackgroundColor',[.2 0 .2],'FontSize',fsize);
            end

            for i=1:cols
                if (ShowLabels>1)  
                    text(showboxsize-3+(i-1)*showboxsize,5,num2str(index+LabelOffset),'units','pixels', ...
                    'horizontalalignment','right','VerticalAlignment','bottom',...
                    'color','w','BackgroundColor',[.2 0 .2],'FontSize',fsize);
                end
                index =index+1;
                patch([(i-1)*nboxsize i*nboxsize i*nboxsize (i-1)*nboxsize],[(rows-1)*nboxsize (rows-1)*nboxsize rows*nboxsize rows*nboxsize],[0 0 0] ,'EdgeColor',[0 0 0],'FaceColor','none');
            end 
        end
    end
 
    if (rows*showboxsize > win_height)
        hSP = imscrollpanel(h0,hIm);
        set(hSP,'Units','normalized','Position',[0 0 1 1])
        api = iptgetapi(hSP);
        api.setMagnification( 1*scale);
        api.setVisibleLocation(0, 0) ;
    end
  end    % redraw_boxes
 

end % ImagicDisplay
