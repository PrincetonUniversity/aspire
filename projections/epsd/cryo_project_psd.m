% CRYO_PROJECT_PSD Project PSD estimate onto a subspace
%
% Usage
%    x_proj_psd = cryo_project_psd(x_psd, V);
%
% Input
%    x_psd: The PSD estimates to be projected in an array of size sig_sz-by-n.
%    V: An array of basis vectors of size sig_sz-by-r.
%
% Output
%    x_proj_psd: The PSD estimates projected onto the linear subspace spanned
%       by the basis vectors in V.

function x_proj_psd = cryo_project_psd(x_psd, V)
    n = size(x_psd, 3);
    r = size(V, 3);

    sig_sz = size(x_psd);
    sig_sz = sig_sz(1:2);

    sig_len = prod(sig_sz);

    x_psd = reshape(x_psd, [sig_len n]);
    V = reshape(V, [sig_len r]);

    [Q, ~] = qr(V, 0);

    x_proj_psd = Q*(Q'*x_psd);

    x_proj_psd = reshape(x_proj_psd, [sig_sz n]);
end
