function h=cryo_mask_radius_3d_auxplot(V,H,r,cntr)

h=gcf;
clf;

% Plot x,y,z projections of the volume and the corresponding bounding
% circle.

th = 0:pi/50:2*pi;
p=size(V,1);

% If H==0, the our hypothesis was not rejected, that is, we
% think that the ring contains only noise pixels, and so our
% object is enclosed in it. In such a case the circle is green.
% Otherwise, be think the circle intersects the object, and so
% the circle is red.
if H==0
    clr='g';
else
    clr='r';
end

% X projection
subplot(1,3,1)
imagesc(squeeze(sum(V,1)))
colormap(gray);
axis image;
hold on;
xunit = r * cos(th) + cntr(2);
yunit = r * sin(th) + cntr(3);

% Don't plot parts of the circle that fall outside the image.
xunit(xunit<=0)=1; xunit(xunit>p)=p;
yunit(yunit<=0)=1; yunit(yunit>p)=p;

plot(xunit, yunit,'LineWidth',2,'Color',clr);
hold off;
title('X')

% Y projection
subplot(1,3,2)
imagesc(squeeze(sum(V,2)))
colormap(gray);
axis image;
hold on;
xunit = r * cos(th) + cntr(1);
yunit = r * sin(th) + cntr(3);

% Don't plot parts of the circle that fall outside the image.
xunit(xunit<=0)=1; xunit(xunit>p)=p;
yunit(yunit<=0)=1; yunit(yunit>p)=p;

plot(xunit, yunit,'LineWidth',2,'Color',clr);
hold off;
title('Y')

% Y projection
subplot(1,3,3)
imagesc(squeeze(sum(V,3)))
colormap(gray);
axis image;
hold on;
xunit = r * cos(th) + cntr(1);
yunit = r * sin(th) + cntr(2);

% Don't plot parts of the circle that fall outside the image.
xunit(xunit<=0)=1; xunit(xunit>p)=p;
yunit(yunit<=0)=1; yunit(yunit>p)=p;

plot(xunit, yunit,'LineWidth',2,'Color',clr);
hold off;
title('Z')