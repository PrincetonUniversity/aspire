
function [corrs_out]=...
    scoreCandidatesDn(scls_lookup_data,pf,max_shift,shift_step,doFilter)

cl_matrix=[scls_lookup_data.scls_lookup1;scls_lookup_data.scls_lookup2];
%eq_class=[scls_lookup_data.eq_class1;scls_lookup_data.eq_class2]';
%eq_class=repmat(eq_class,scls_lookup_data.ntheta,1);
%eq_class=eq_class(:);
tic
[corrs_out]=...
    cryo_clmatrix_scl_ML_Dn2(pf,cl_matrix,max_shift,shift_step,doFilter,scls_lookup_data);
toc



