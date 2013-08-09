% test_all_systems

list = {
	'aspire_buff2mat test'
	'fatrix2_tests test'
	'Fatrix_test_basic test'
	'fatrix_one test'
	'block_fatrix_test'
	'Gblock_test'
...
	'Gblur test'
	'Gcone test'
	'Gdelaysum1 test'
	'Gdft test'
	'Gdsft test'
	'Gdown test'
%	'Gembed test'
	'Glinear test'
%	'Gnearest test' % obsolete
	'Gnufft test' % do before Gmri; calls 'Gnufft_test0' 'Gnufft_test'
...
	'Gmri test'
	'Gsparse_test'
%	'Gtomo2_dscex test' % within Gwtmex_test
	'Gtomo2_moj_test'
	'Gtomo2_strip test'
	'Gtomo2_table_test'
	'Gtomo2_wtmex_test'
	'Gtomo3_test'
	'Gtomo3_test_adj'
	'Gtomo_nufft_test'
	'Gtomo_nufft_test_adj'
	'Gtranslate test'
	'wtf_read test'
};

im nan-fail
run_mfile_local(list)
