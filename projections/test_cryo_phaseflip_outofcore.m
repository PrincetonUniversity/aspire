% Test the the output of the functions cryo_phaseflip and
% cryo_phaseflip_outofcore match
%
% Yoel Shkolnisky, October 2017

CTFtestdata=readSTAR('phaseflip_testdata.star');
FPprojs=cryo_phaseflip(CTFtestdata,'',-1,3);
tmpfile=tempmrcname;
cryo_phaseflip_outofcore(CTFtestdata,'',tmpfile,-1,3);
FPprojs2=ReadMRC(tmpfile);
delete(tmpfile);

err=norm(FPprojs(:)-FPprojs2(:)); % Should be zero.
fprintf('diff between two methods = %5.3e\n',err);
if err<1.0e-14
    fprintf('test ok\n');
else
    fprintf('test FAILED\n');
end
