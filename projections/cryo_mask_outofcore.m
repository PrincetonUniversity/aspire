function cryo_mask_outofcore(instackname,outstackname,r,risetime)
% XXX FILL
% Yoel Shkolnisky, April 2017.

instackReader=imagestackReader(instackname);
szprojs=instackReader.dim;
outstackWriter=imagestackWriter(outstackname,szprojs(3));
for k=1:szprojs(3)
    p=instackReader.getImage(k);
    p=cryo_mask(p,1,r,risetime);
    outstackWriter.append(p);
end
outstackWriter.close;
