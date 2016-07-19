function [B,bestfitval]=cryo_estimate_Bfactor(projs,pixA,verbose,prefix,outputdir)
%CRYO_ESTIMATE_BFACTOR Estimate the B factor
%
% [B,bestfitval]=cryo_estimate_Bfactor(projs,pixA)
%   Estimate the B factor of the given set of projections. All projections
%   are assumed to correspond to the same CTF parameters and same B factor.
%   Returns the estimated B factor and quality of the fit using the
%   estimated B factor.
%
% [B,bestfitval]=cryo_estimate_Bfactor(projs,pixA,verbose)
%   verbose=0   No messages are printed (silent).
%   verbose=1   Print debug log messages.
%   verbose=2   Print debug log messages and save debug figures.
%
%   If verbose=2, all figures are saved with prefix 'prefix' (default is
%   bfactor). 
%   If verbose=2, outputdir must be specified.
%
% [B,bestfitval]=cryo_estimate_Bfactor(projs,pixA,verbose,prefix,outputdir)
%   If verbose=2, use the give prefix to all figure names. Figures are
%   saved into outputdir. If verbose~=2, prefix and outputdir are ignored.
%
% Example:
%   See test_cryo_estimate_bfactor
%
% Yoe Shkolnisky, February 2016.

% Verbose=0 silent
% verbose=1 Only printout messages
% verbose=2 save debug figures

if (nargin<2) || (isempty(pixA)) || (pixA<0)
    pixA=1;
    log_message('Pixel size not given. Using dummy size pixA=1.');
end

if nargin<3
    verbose=0;
end

if nargin<4
    prefix='bfactor';
end

if verbose==2
    if nargin<5
        error('Output directory must be specified.');
    end
end

fprefix=fullfile(outputdir,prefix);

% Currently, we assume that the images are odd-sized (assumed when
% estimating the center of each image).
if mod(size(projs,1),2)==0
    projs=projs(1:end-1,1:end-1,:);
end

% Compute the power spectrum of the projections.
log_message('Estimating power spectrum of the projections');
psd_signal_2d=cryo_noise_estimation(projs,0);
orig2d=ceil((size(psd_signal_2d+1)/2)); % Center of the 2D rdially averaged PSD of the projections.
psd_signal_1d=psd_signal_2d(orig2d(1):end,orig2d(2)); % Central ray of the PSD.
psd_signal_1d=psd_signal_1d./max(psd_signal_1d); % Normalize ray to maximum 1.

% Find minima of in the 1d radially averaged power spectrum.
log_message('Detecting minima of 1d radially averaged power spectrum.');
[minpeaklocs] = peakfinder(psd_signal_1d,1.0e-3,100,-1);
log_message('Detected %d minima',numel(minpeaklocs));

% Discard first few minima and make sure we have at most Nminpeaks minima.
% Current Nminpeaks equals the number of detected peaks, but a smaller
% number can be used.
Nminpeaks=numel(minpeaklocs);
Ndiscardpeaks=2;
minpeaklocs=minpeaklocs(1:min(numel(minpeaklocs),Nminpeaks+Ndiscardpeaks));
minpeaklocs=minpeaklocs(Ndiscardpeaks+1:end);

% Compute physical frequencies present in the power spectrum
fmax=1/(2*pixA);
n=numel(psd_signal_1d);
freqs=(0:n-1)/n*fmax;
freqs=freqs(:);

% Interpolate backgorund signal thorugh minima
%backgroundfit = polyfit(freqs(minpeaklocs),psd_signal_1d(minpeaklocs),5);
%background_signal=polyval(backgroundfit,freqs);
background_signal = spline(freqs(minpeaklocs),psd_signal_1d(minpeaklocs),freqs);

if verbose>=2
    % Plot the 1d radially averages power spectrum and fitted background.   
    clf;
    h1a=plot(freqs,psd_signal_1d);
    freqticks=freqs(2:round(n/7):end);
    set(gca,'XTick',freqticks);
    pixres=round(10*(1./freqticks))/10;
    set(gca,'XTickLabel',num2cell(pixres));
    hold on;
    % Plot the peaks we are using. 
    h1b=scatter(freqs(minpeaklocs),psd_signal_1d(minpeaklocs),20,'r','filled');        
    h1c=plot(freqs,background_signal,'r'); % Plot the fitted background signal.
    title('power spectrum and fitted background');
    legend([h1a,h1c,h1b],{'power spectrum','fitted background','minina'});
    fname=sprintf('%s_%s',fprefix,'psd_1d.eps');
    print(fname,'-depsc');
end

% Subtract the background single from the estimted psd.
corrected_psd=psd_signal_1d-background_signal;

if verbose>=2
    clf;
    % Plot the background subtracted power spectrum.
    plot(freqs,corrected_psd);
    freqticks=freqs(2:round(n/7):end);
    set(gca,'XTick',freqticks);
    pixres=round(10*(1./freqticks))/10;
    set(gca,'XTickLabel',num2cell(pixres));
    hold on;
    scatter(freqs(minpeaklocs),corrected_psd(minpeaklocs),20,'r','filled');
    title('Background subtracted power spectrum');
    fname=sprintf('%s_%s',fprefix,'psd_1d_corrected.eps');
    print(fname,'-depsc');
end

% Keep only frequencies after the first minima
cropped_psd=corrected_psd(minpeaklocs(1):minpeaklocs(numel(minpeaklocs)));
cropped_freqs=freqs(minpeaklocs(1):minpeaklocs(numel(minpeaklocs)));

% Find maxima of the corrected power spectrum.
[maxpeaklocs] = peakfinder(cropped_psd,1.0e-4,-1,1);
Nfitmaxpeaks=4; % Number of maxima to use for the fitting.
maxpeaklocs_sav=maxpeaklocs;


bestfitval=1000;  % Value of the objective function of the best fit so far.
bestfitB=-1;      % Estimated B from the best fit so far.
bestfitidx=-1;

for Nskipmaxpeaks=1:3  % Peak to start the fitting (in some cases the first
                       % peaks do not correspond to the exponential decay
                       % part and need to be disacrded.
                       
    maxpeaklocs=maxpeaklocs_sav(Nskipmaxpeaks+1:min(numel(maxpeaklocs_sav),Nfitmaxpeaks+Nskipmaxpeaks));
    if numel(maxpeaklocs)<3
        break % Not enough peaks for fitting
    end
    
    params0=initialBfit(cropped_freqs(maxpeaklocs(1:4)),cropped_psd(maxpeaklocs(1:4)));    
    if ~isempty(params0)
        log_message('Initially estimated parameters (a,B)=(%5.3e,%7.2f)',params0(1),params0(2));
    else
        log_message('Failed to estimate initial parameters (a,B)');
    end
        
        
    params1=refineBfit(double(cropped_freqs(maxpeaklocs)), double(cropped_psd(maxpeaklocs)),params0);
    log_message('Refined parameters (a,B)=(%5.3e,%7.2f)',params1(1),params1(2));
        
    [params2,fitval]=optimizeBfit(cropped_freqs(maxpeaklocs),cropped_psd(maxpeaklocs),params1);
    log_message('Optimized parameters (a,B)=(%5.3e,%7.2f)',params2(1),params2(2));
    log_message('Fitting score fitval=%e',fitval);

    log_message('Current best parameters fitval=%e, B=%7.2f',bestfitval,bestfitB);
    
    if fitval<bestfitval
        bestfitval=fitval;
        bestfitidx=Nskipmaxpeaks;
        bestfitB=params2(2);
        log_message('Better fit detected when skipping %d maxima. Updating parameters to fitval=%e, B=%7.2f',bestfitidx,bestfitval,bestfitB);
    end
    
    
    if verbose>=2
        clf;
        % Plot a zoom in from the first zero of the power spectrum.
        h3a=plot(cropped_freqs,cropped_psd);
        freqticks=cropped_freqs(2:round(numel(cropped_freqs)/7):end);
        set(gca,'XTick',freqticks);
        xlim([min(freqticks) max(freqticks)]);
        pixres=round(10*(1./freqticks))/10;
        set(gca,'XTickLabel',num2cell(pixres));        
        hold on;
        
        % Plot the peaks used for fitting.
        scatter(cropped_freqs(maxpeaklocs_sav),cropped_psd(maxpeaklocs_sav),20,'g','s','filled');
        h3b=scatter(cropped_freqs(maxpeaklocs),cropped_psd(maxpeaklocs),20,'r','filled');
        
        if ~isempty(params0)
            % Plot the initial fit thorught the maxima.
            initial_fit=params0(1)*exp(-params0(2).*cropped_freqs.^2);
            initial_fit(initial_fit>1)=0;
            h3c=plot(cropped_freqs,initial_fit,'g');
        end
        
        % Plot refined fit
%         refined_fit=params1(1)*exp(-params1(2).*cropped_freqs.^2);
%         refined_fit(refined_fit>1)=0;
%         h3d=plot(cropped_freqs,refined_fit,'r');
        
        % Plot optimized fit
        optimizedfit=params2(1)*exp(-params2(2).*cropped_freqs.^2);
        h3e=plot(cropped_freqs,optimizedfit,'m');
        
        if ~isempty(params0)
            % legend([h3a,h3b,h3c,h3d,h3e],{'corrected psd','fitting peaks',...
            %   'initial fit','refined fit','optimized fit'});
            legend([h3a,h3b,h3c,h3e],{'corrected psd','fitting peaks',...
                'initial fit','optimized fit'});
            
        else
            % legend([h3a,h3b,h3d,h3e],{'corrected psd','fitting peaks',...
            %   'refined fit','optimized fit'});
            legend([h3a,h3b,h3e],{'corrected psd','fitting peaks',...
                'optimized fit'});
            
        end
        
        titlestr=sprintf('Curr: B=%4.1f, F=%3.1e, Best:- B=%4.1f, F=%3.1e',params2(2),fitval,bestfitB,bestfitval);
        title(titlestr);
        
        fname=sprintf('%s_fitting_%d.eps',fprefix,Nskipmaxpeaks);
        print(fname,'-depsc');
    end
    
end
B=bestfitB;
log_message('Final estimate B=%7.2f when skipping %d maxima',B,bestfitidx);

end

function params0=initialBfit(x,y)
% Create initial guess for the fitting Gaussian using least sqaures
x=x(y>1.0e-4); % Due to log, values cannot be too small
y=y(y>1.0e-4);

if numel(y)<3
    params0=[]; % Not enough point to estimate an initial guess.
else
    b=log(y);
    A=[ones(numel(x(:)),1) x(:).^2];
    tmp=A\b;
    a0=exp(tmp(1));
    b0=-tmp(2);
    params0=[a0;b0];
    
    if b0<0
        params0=[]; % Something went wrong with the estimation of the
        % initial guess. No initial guess.
    end
end
end

function params1=refineBfit(x,y,params0)
% Refine the initial fit. This function essentially generates a starting
% point for the subsequent optimization procedure, as it returns a fit even
% if the initital fit failed.
[xData, yData] = prepareCurveData(x,y);
ft = fittype( 'a*exp(-b*x^2)', 'independent', 'x', 'dependent', 'y' ); % Set up fittype and options.
opts = fitoptions( 'Method', 'NonlinearLeastSquares','Robust','LAR', 'StartPoint', params0 );
opts.Display = 'Off';
[fitresult,~] = fit( xData, yData, ft, opts ); % Fit model to data.
params1=[fitresult.a;fitresult.b];
end


function [params2,fitval]=optimizeBfit(x,y,initialpoint)
xData=double(x);
yData=double(y);
fun=@(params) (sum((params(1).*exp(-params(2).*xData.^2)-yData).^2));
lb=[0,0];   % Lower bound for (a,B)
ub=[1,800]; % Upper bound for (a,B)
A = []; b = []; Aeq = []; beq = []; % No linear constraints
nonlcon = @(params)Bfitconstraints(params,xData,yData);
options = optimoptions('fmincon','Display','off','Algorithm','active-set');
[params2,fitval] = fmincon(fun,double(initialpoint),A,b,Aeq,beq,lb,ub,nonlcon,options);
end