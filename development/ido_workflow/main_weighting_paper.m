
% main_wrapper(n, DATA_SET, Ns, NN, S_WEIGHTS, J_WEIGHTS, GROUPS)

%main_wrapper(179, '', 1000, [9,7,5,4,3,2], 1, 0, 1);
%main_wrapper(179, '', 1000, [2], 1, 0, 1);
main_wrapper(179, '', 1000, [9,5], 1, 0, 2);
main_wrapper(89, '', 3000, [2], 1, 0, 1:2);

% gather resolutions
resolutions = [];
tmp = [0,0,0];
N=1000; n=179;
for nn=[9,5]
    for g=1:2
        base_path=sprintf('/home/idog/matlab/results/weighting_paper/80s_%d/N%d_nn%02d_sw1_jw0_g%d',n,N,nn,g);
        x=load(sprintf('%s/resolution.mat',base_path));
        tmp(g+1) = x.res;
        % TODO: for g=1 also use x.fc
    end
    v1=ReadMRC(sprintf('/home/idog/matlab/results/weighting_paper/80s_%d/N%d_nn%02d_sw1_jw0_g%d/vol.mrc',n,N,nn,1));
    v2=ReadMRC(sprintf('/home/idog/matlab/results/weighting_paper/80s_%d/N%d_nn%02d_sw1_jw0_g%d/vol.mrc',n,N,nn,2));
    [~,~,v2,~] = cryo_align_densities(v1,v2);
    [res,fc] = fsc(v1, v2);
    tmp(1)=res;
    resolutions = [resolutions;tmp];
end
N=3000; n=89; nn=2;
for g=1:2
    base_path=sprintf('/home/idog/matlab/results/weighting_paper/80s_%d/N%d_nn%02d_sw1_jw0_g%d',n,N,nn,g);
    x=load(sprintf('%s/resolution.mat',base_path));
    tmp(g+1) = x.res;
    % TODO: for g=1 also use x.fc
end
v1=ReadMRC(sprintf('/home/idog/matlab/results/weighting_paper/80s_%d/N%d_nn%02d_sw1_jw0_g%d/vol.mrc',n,N,nn,1));
v2=ReadMRC(sprintf('/home/idog/matlab/results/weighting_paper/80s_%d/N%d_nn%02d_sw1_jw0_g%d/vol.mrc',n,N,nn,2));
[~,~,v2,~] = cryo_align_densities(v1,v2);
[res,fc] = fsc(v1, v2);
tmp(1)=res;
resolutions = [resolutions;tmp];
disp(resolutions);


% scores hists
a = 2.03;
peak2sigma = 2.31e-2;
hist_intervals = 100;
h = 1/hist_intervals;
hist_x = ((h/2):h:(1-h/2))';
N=1000;
NN=[9,7,5,4,3,2];
fprintf('N\tK\tP\tsigma [deg]\tR2\n');
for i = 1:numel(NN)
    w=load(sprintf(...
        '/home/idog/matlab/results/weighting_paper/80s_179/N1000_nn%02d_sw1_jw0_g1/weights.mat',NN(i)));
    % Initialize variables
    A = (N*(N-1)*(N-2)/2)/hist_intervals*(a+1); % normalization factor of one component of the histogram
    start_values = [0, 0.6^3, 2.5, 0.9]; % B, P, b, x0
    B0 = start_values(2)*(N*(N-1)*(N-2)/2) /... % normalization of 2nd component: B = P*N_delta/sum(f), where f is the component formula
        sum(((1-hist_x).^start_values(3)).*exp(-start_values(3)/(1-start_values(4)).*(1-hist_x)));
    start_values(1) = B0;
    % Fit the distribution
    conf = fitoptions('Method','NonlinearLeastSquares',...
        'Robust','LAR',...
        'Lower',[B0/100, 0, 2, 0.5],... % B, P, b, x0
        'Upper',[2*B0, 0.5, 3, 1],...
        'Startpoint',start_values);
    ft = fittype(['(1-P)*' num2str(A) '*(1-x)^' num2str(a) ' + P*B*((1-x)^b)*exp(-b/(1-x0)*(1-x))'],'options',conf);
    [c,gof] = fit(hist_x, w.scores_hist.hist, ft);
    % Derive P and sigma
    Rsquare = gof.rsquare;
    P = c.P^(1/3);
    peak = c.x0;
    sigma = (1-peak)/peak2sigma;
    Bf = c.P*(N*(N-1)*(N-2)/2) / sum(((1-hist_x).^c.b).*exp(-c.b/(1-c.x0).*(1-hist_x)));
    y = w.scores_hist.hist;
    Rsquare_man = 1 - sum((y-c(hist_x)).^2)/sum((y-mean(y)).^2);
    fprintf('%d\t%d\t%.0f%%\t%.1f\t%.2f\t%.2f\t%.2f\n', ...
        N,NN(i)+1,100*P,sigma,Bf/c.B,Rsquare_man,Rsquare);
    % plot hist
    hfig = figure;
    %title(sprintf('P=%.0f%%, \sigma=%.1f',100*P,sigma), 'fontsize',14);
    hold on;
    bar(hist_x,w.scores_hist.hist,'facecolor','c', 'edgecolor','none');
    hplot = plot(c);
    set(hplot, 'linewidth',1.5);
    xlim([0,1]);
    xlabel('Triplet score');
    ylabel('');
    legend('hide');
    save_figure(hfig, sprintf('/home/idog/matlab/results/weighting_paper/scores_hist_nn%02d',NN(i)));
end
fprintf('\n');
