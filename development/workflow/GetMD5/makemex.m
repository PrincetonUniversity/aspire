mexFile = ['GetMd5', '.', mexext];

%Check for existing mex in path
whichMex = which(mexFile);
if ~isempty(whichMex)
   fprintf('  Mex exist in already:  %s\n', whichMex);
else
    try
		GetMD5();
    catch
		fprintf('Error in compling GetMD5 function\n');
    end		
end

% %Check only in this directory
% if ( ~exist(fullfile(pathstr, 'development/workflow/GetMD5', mexFile),'file') )
% 	try
% 		GetMD5();
%     catch
% 		fprintf('Error in compling GetMD5 function\n');
%     end		
% end