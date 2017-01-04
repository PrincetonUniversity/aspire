

M = [1,5];

for N = 1000%200
    for nn = 2%[7,5]
        for m = 1%[0,1]
            for j = 0%[0,1]
                fprintf('N=%d, nn%02d, Sw=%d, Jw=%d\n', N, nn, m, j);

                d = get_reconstruction_data(89,'',N,nn,m,j,1);
                
                load(sprintf(...
                    '../../aspire_merge/output/80s_89_complete/N%d_nn%02d_g1_3N_method%d_jscores%d_itr1/data_final_%d.mat',...
                    N,nn,M(m+1),j,N));
                
                fprintf('projs\t%e\n', norm(d.projs.projs(:)-data.xxN1.projs(:)) );
                fprintf('cls\t%e\n', norm(d.cls.clstack(:)-data.N1N1.clstack(:)) );
                fprintf('Rij0\t%e\n', norm(d.rij.Rij0(:)-data.xxN2.Rij0(:)) );
                fprintf('jsync\t%e\n', norm(d.jsync.J_sync(:)-data.N2.J_sync(:)) );
                fprintf('Rij\t%e\n', norm(d.jsync.Rij(:)-data.xxN2.Rij(:)) );
                
                if m; fprintf('Pij\t%e\n', norm(d.w.Pij(:)-data.N2.Pij(:)) ); end
                
                fprintf('Ri\t%e\n', norm(d.ri.rotations(:)-data.xxN1.R(:)) );
                R = data.xxN1.R;
                s = sign(d.ri.rotations(:,1,1).*data.xxN1.R(:,1,1));
                for i = 1:3
                    R(i,:,:) = s(i)*R(i,:,:);
                end
                fprintf('Ri (fixed)\t%e\n', norm(d.ri.rotations(:)-R(:)) );
                
                fprintf('shifts\t%e\n', norm(d.shifts.shifts(:)-data.others.est_shifts(:)) );
                
                fprintf('\n');
                
                v1 = ReadMRC(sprintf(...
                    '../../aspire_merge/output/80s_89_complete/N%d_nn%02d_g1_3N_method%d_jscores%d_itr1/vol.mrc',...
                    N,nn,M(m+1),j));
                v10 = ReadMRC(sprintf(...
                    '../../aspire_merge/output/80s_89_complete/N%d_nn%02d_g1_3N_method%d_jscores%d_itr1/vol_ref-aligned.mrc',...
                    N,nn,M(m+1),j));
                v2 = ReadMRC(sprintf(...
                    '../../aspire_merge/output/80s_89/N%d_nn%02d_sw%d_jw%d_g1/vol.mrc',...
                    N,nn,m,j));
                v20 = ReadMRC(sprintf(...
                    '../../aspire_merge/output/80s_89/N%d_nn%02d_sw%d_jw%d_g1/vol_ref-aligned.mrc',...
                    N,nn,m,j));
                v0 = ReadMRC('../../aspire_merge/output/vol_80s_89.mrc');
                
                [r1,f1] = fsc(v10,v0);
                [r2,f2] = fsc(v20,v0);
                [~,f12_indirect] = fsc(v10,v20);
                [~,~,v21,~] = cryo_align_densities(v1,v2);
                [~,f12] = fsc(v1,v21);
                h=figure;
                plot_fsc(f1,1.34*359/89,0.143,h); hold on;
                plot_fsc(f2,1.34*359/89,0.143,h); hold on;
                plot_fsc(f12_indirect,1.34*359/89,0.143,h); hold on;
                plot_fsc(f12,1.34*359/89,0.143,h); hold on;
                title(sprintf('N=%d, nn%02d, Sw=%d, Jw=%d\n', N, nn, m, j));
                fprintf('high freq corr\t%f\n', f12(end));
                
            end
        end
    end
end
