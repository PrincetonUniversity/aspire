function  demean_y_batch( F, EstMean, index, nbatch, nim)

% Subtract estimated "mean image" from stack of images
% y - A \mu : Demeaned y
% New y_i = y_i - A_i \mu
% Do in place
% Tejal, Oct 2015

for nb=1:nbatch
	filename=fullfile('/scratch/tbhamre/cwf_batch/', sprintf('set%d',nb));
	load (filename);
	offset=(nb-1)*(nim/nbatch);
	for i=1:size(curr_batch,3)
		curr_batch(:,:,i) = curr_batch(:,:,i)-(F(:,:,index(i+offset)).*EstMean);
	end
	filename=fullfile('/scratch/tbhamre/cwf_batch/', sprintf('demean_set%d',nb));
	save(filename,'curr_batch','-v7.3');
end
