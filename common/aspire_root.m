% ASPIRE_ROOT Obtain root directory of ASPIRE package
%
% Usage
%    root = aspire_root();
%
% Output
%    root: A string containing the path to the ASPIRE package.

function root = aspire_root()
    full_path = mfilename('fullpath');

    [root, ~, ~] = fileparts(full_path);

    ind = find(root == filesep, 1, 'last');

    root = root(1:ind-1);
end
