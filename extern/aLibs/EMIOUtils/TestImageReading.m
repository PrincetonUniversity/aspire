% Test image reading
% Compare local image buffering, i.e. BufferedReadImage() with reading
% images with fread().  The conclusion is that fread() is about 30% faster,
% so the file system is quite well cached.

name='start.hed';

  [h s]=ReadImagic(name,1,-1);
tic
for i=1:10000
      m1=fread(h,s.nx*s.ny,s.string);
      m1=reshape(m1,s.nx,s.ny);
%       imacs(m1); title(i); drawnow;
  end;
toc;
fclose(h);
%%

name='FlCtrRot.mrc';

  [h s]=ReadMRC(name,1,-1);
tic
for i=1:10000
      m1=fread(h,s.nx*s.ny,s.string);
      m1=reshape(m1,s.nx,s.ny);
%        imacs(m1); title(i); drawnow;
  end;
toc;
fclose(h);

%%
name='FlCtrRot.mrc';
name='start.hed';
tic
for i=1:10000
    m1=BufferedReadImage(name,i);
%        imacs(m1); title(i); drawnow;
end;
toc
