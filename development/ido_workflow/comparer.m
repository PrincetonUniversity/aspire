

M = [1,5];

for N = 200
    for nn = 5
        for m = 0:1
            for j = 0:1
                fprintf('N=%d, nn%02d, Sw=%d, Jw=%d\n', N, nn, m, j);

                d = get_reconstruction_data(89,'',N,nn,m,j,1);
                
                load(sprintf(...
                    '../../aspire_merge/output/80s_89_complete/N%d_nn%02d_g1_3N_method%d_jscores%d_itr1/data_final_%d.mat',...
                    N,nn,M(m+1),j,N));
                
                fprintf('projs\t%f\n', norm(d.projs.projs(:)-data.xxN1.projs(:)) );
                fprintf('cls\t%f\n', norm(d.cls.clstack(:)-data.N1N1.clstack(:)) );
                fprintf('Rij0\t%f\n', norm(d.rij.Rij0(:)-data.xxN2.Rij0(:)) );
                fprintf('jsync\t%f\n', norm(d.jsync.J_sync(:)-data.N2.J_sync(:)) );
                fprintf('Rij\t%f\n', norm(d.jsync.Rij(:)-data.xxN2.Rij(:)) );
                
                if m; fprintf('Pij\t%f\n', norm(d.w.Pij(:)-data.N2.Pij(:)) ); end
                
                fprintf('Ri\t%f\n', norm(d.ri.rotations(:)-data.xxN1.R(:)) );
                R = data.xxN1.R;
                s = sign(d.ri.rotations(:,1,1).*data.xxN1.R(:,1,1));
                for i = 1:3
                    R(i,:,:) = s(i)*R(i,:,:);
                end
                fprintf('Ri (fixed)\t%f\n', norm(d.ri.rotations(:)-R(:)) );
                
                fprintf('shifts\t%f\n', norm(d.shifts.shifts(:)-data.others.est_shifts(:)) );
                
                fprintf('\n');
                
            end
        end
    end
end
