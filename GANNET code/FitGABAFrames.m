function [MRS_struct] = FitGABAFrames(MRS_struct,ii)
%
% Keep input and output structure names the same to add output data to
% the exisiting MRS data structure.
%Inputs:
% MRS_struct = structure with data loaded from MRSLoadPfiles
% ii = index of spectrum to analyse in MRS_struct

% From GannetFit 120510

% hello from daws

FIT_LSQCURV = 0;
FIT_NLINFIT = 1;

fit_method = FIT_NLINFIT; %FIT_NLINFIT;
waterfit_method = FIT_NLINFIT;

% Needs to be run on GABAerror branch of Gannet - i.e. needs
%    diffSubSpectrum array.
%GABAData=MRS_struct.gabaspec;
GABAData = MRS_struct.diffSubSpectra(:,:,ii);

freq=MRS_struct.freq;


% if strcmp(MRS_struct.Reference_compound,'H2O')
%     WaterData=MRS_struct.waterspec;
% end
MRS_struct.versionfit = '120222a';
disp(['GABA Fit Version is ' MRS_struct.versionfit ]);
fitwater=1;
numscans=size(GABAData);
numscans=numscans(1);

%110624
epsdirname = [ 'MRSfit_' datestr(clock,'yymmdd') ];
numsubspec = MRS_struct.Navg(ii)/4; % Navg / (Nphasecycles * N_ON_OFF)

for jj = 1:numsubspec
    
    
    MRS_struct.pfile{ii};
    %%%%%%%%%%%% GABA fit (Gaussian) %%%%%%%%%%%%%%%
    % Originally (RE version)  gabadata was rescaled to the max across all spectra
    % in the dataset.  Now I normalised the data with a const
    % ... work in progress
    % ...from GaussModel;
    % x(1) = gaussian amplitude
    % x(2) = 1/(2*sigma^2)
    % x(3) = centre freq of peak
    % x(4) = amplitude of linear baseline
    % x(5) = constant amplitude offset

    %Need to set the bounds from the data parameters...
    %MRS_struct.freq(17342)
    %MRS_struct.freq(18000)
    %Hard code it to fit from 2.75 ppm to 3.55 ppm
    z=abs(MRS_struct.freq-3.55);
    lowerbound=find(min(z)==z)
    z=abs(MRS_struct.freq-2.79);%2.75
    upperbound=find(min(z)==z)
    %lowerbound=17342;
    %upperbound=17961;
    %upperbound=18000;
    freqbounds=lowerbound:upperbound;
    plotbounds=(lowerbound-150):(upperbound+150);

    %maxinGABA=max(real(GABAData(freqbounds)));
    maxinGABA=max(real(GABAData(freqbounds)));
    %maxinGABA=1;

    % smarter estimation of baseline params, Krish's idea (taken from Johns
    % code; NAP 121211
    grad_points = (real(GABAData(jj,upperbound)) - real(GABAData(jj,lowerbound))) ./ ...
        (upperbound - lowerbound); %in points
    LinearInit = grad_points ./ (MRS_struct.freq(1) - MRS_struct.freq(2)); %in ppm
    constInit = (real(GABAData(jj,upperbound)) + real(GABAData(jj,lowerbound))) ./2;
    xval = [ 1:(upperbound-lowerbound+1) ];
    linearmodel = grad_points .* xval + GABAData(jj,lowerbound);
    %End copy code

    resnorm=zeros([numscans size(freqbounds,2)]);
    %residuals=resnorm;
    size(resnorm);


    %  GaussModelInit = [10*maxinGABA -90 3.026 0 0];
    % NP; but taken from Johns code, now initialise with parameters declared above
    GaussModelInit = [10*maxinGABA -90 3.026 LinearInit constInit]; %default in 110624

    %OLD INITS; WHY NOT CUT THESE OUT (NP)
    %GaussModelInit = [4.1314 -140.0000 3.0005 -0.8776 0.5684]; %from MINLSQ
    %GaussModelInit = [1 -140.0000 3.0005 -0.8776 0.5684]; %works
    %GaussModelInit = [maxinGABA -90 3.026 0 0]; %works
    %GaussModelInit = [ 1 -90 3.026 0 0]; %fails
    %GaussModelInit = [ 1 -90 3.026 -0.8776 0.5684]; %works

    lb = [0 -200 2.87 -40*maxinGABA -2000*maxinGABA]; %NP; our bounds are 0.03 less due to creatine shift
    ub = [4000*maxinGABA -40 3.12 40*maxinGABA 1000*maxinGABA];



    %NP; From Johns code, initialising steps, just copied over
    options = optimset('lsqcurvefit');
    options = optimset(options,'Display','off','TolFun',1e-10,'Tolx',1e-10,'MaxIter',1e5);
    nlinopts = statset('nlinfit');
    nlinopts = statset(nlinopts, 'MaxIter', 1e5);


    [GaussModelParam(jj,:),resnorm, residg] = lsqcurvefit(@(xdummy,ydummy) GaussModel_area(xdummy,ydummy), ...
        GaussModelInit, ...
        freq(freqbounds),...
        real(GABAData(jj,freqbounds)), ...
        lb,ub,options);
    residg = -residg;

    if(fit_method == FIT_NLINFIT)
        %else  % it's FIT_NLINFIT
        GaussModelInit = GaussModelParam(jj,:);
        % 1111013 restart the optimisation, to ensure convergence

        for fit_iter = 1:100
            [GaussModelParam(jj,:), residg, J, COVB, MSE] = nlinfit(freq(freqbounds), real(GABAData(jj,freqbounds)), ... % J, COBV, MSE edited in
                @(xdummy,ydummy) GaussModel_area(xdummy,ydummy), ...
                GaussModelInit, ...
                nlinopts);
            %MRS_struct.fitparams_iter(fit_iter,:,ii) = GaussModelParam(ii,:);
            GaussModelInit = GaussModelParam(jj,:);
            ci = nlparci(GaussModelParam(jj,:), residg,'covar',COVB); %copied over
        end
    end
    %NP end copy Johns code


    GABAheight = GaussModelParam(jj,1);
    % FitSTD reports the standard deviation of the residuals / gaba HEIGHT
    %MRS_struct.GABAFitError(ii)  =  (ci(1,2)-ci(1,1))/sum(ci(1,:),2)*100;
    MRS_struct.GABAFitError(ii)  =  100*std(residg)/GABAheight;
    % x(2) = 1/(2*sigma^2).  Convert from points to Hz

    % area under gaussian is a * (2 * pi * sigma^2)^(1/2), and
    % GaussModelParam(:,2) = 1 / (2 * sigma^2)
    % This sets GabaArea as the area under the curve.
    MRS_struct.gabaArea(ii)=GaussModelParam(ii,1)./sqrt(-GaussModelParam(ii,2))*sqrt(pi);

    sigma = ( 1 / (2 * (abs(GaussModelParam(ii,2)))) ).^(1/2);
    MRS_struct.GABAFWHM(ii) =  abs( (2* 42.576*3) * sigma);
    MRS_struct.GABAModelFit(ii,:)=GaussModelParam(ii,:);

    figure(22)
    plot(freq(freqbounds),GaussModel_area(GaussModelParam(jj,:),freq(freqbounds)),'r',...
    freq(plotbounds),real(GABAData(jj,plotbounds)), 'b' )
    %, ...
    
    %freq(freqbounds),residg,'b');
    
    legendtxt = regexprep(MRS_struct.pfile{ii}, '_','-');
    title(legendtxt);
    set(gca,'XDir','reverse');
    %set(gca,'YTick',[], 'Xgrid', 'on');
    oldaxis = axis;
    axis( [2.6 3.6 oldaxis(3) oldaxis(4) ] )    
    
    end

% end of MRSGABAfit

%%%%%%%%%%%%%%%%%%%%%%%% GAUSS MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = GaussModel_area(x,freq)

% x(1) = gaussian amplitude
% x(2) = 1/(2*sigma^2)
% x(3) = centre freq of peak
% x(4) = amplitude of linear baseline
% x(5) = constant amplitude offset

%F = x(1)*sqrt(-x(2)/pi)*exp(x(2)*(freq-x(3)).*(freq-x(3)))+x(4)*(freq-x(3))+x(5);
F = x(1)*exp(x(2)*(freq-x(3)).*(freq-x(3)))+x(4)*(freq-x(3))+x(5);



%%%%%%%%%%%%%%%%  OLD LORENTZGAUSSMODEL %%%%%%%%%%%%%%%%%%%%
%function F = LorentzGaussModel(x,freq)
%Lorentzian Model multiplied by a Gaussian.  gaussian width determined by
%x(6). x(7) determines phase.
%F = ((ones(size(freq))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1)*x(1))*cos(x(7))+(ones(size(freq))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1)*x(2).*(freq-x(3)))*sin(x(7))).*(exp(x(6)*(freq-x(3)).*(freq-x(3))))+x(4)*(freq-x(3))+x(5);


%%%%%%%%%%%%%%%%  NEW LORENTZGAUSSMODEL %%%%%%%%%%%%%%%%%%%%
function F = LorentzGaussModel(x,freq)
% CJE 24Nov10 - removed phase term from fit - this is now dealt with
% by the phasing of the water ref scans in MRSLoadPfiles
%Lorentzian Model multiplied by a Gaussian.
% x(1) = Amplitude of (scaled) Lorentzian
% x(2) = 1 / hwhm of Lorentzian (hwhm = half width at half max)
% x(3) = centre freq of Lorentzian
% x(4) = constant baseline amplitude
% x(5) =  -1 / 2 * sigma^2  of gaussian

% 110624; remove linear baseline term
% OBSOLETE:
% x(4) = linear baseline amplitude
% x(5) = constant baseline amplitude
% x(6) =  -1 / 2 * sigma^2  of gaussian

% Lorentzian  = (1/pi) * (hwhm) / (deltaf^2 + hwhm^2)
% Peak height of Lorentzian = 4 / (pi*hwhm)
% F is a normalised Lorentzian - height independent of hwhm
%   = Lorentzian / Peak

%F =((ones(size(freq))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1)*x(1))*cos(x(7))+(ones(size(freq))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1)*x(2).*(freq-x(3)))*sin(x(7))).*(exp(x(6)*(freq-x(3)).*(freq-x(3))))+x(4)*(freq-x(3))+x(5);
% remove phasing
F = (x(1)*ones(size(freq))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1))  ...
    .* (exp(x(6)*(freq-x(3)).*(freq-x(3)))) ... % gaussian
    + x(4)*(freq-x(3)) ... % linear baseline
    +x(5); % constant baseline

% 110624 removed:
% + x(4)*(freq-x(3)) ... % linear baseline
% 111214 replaced:

%%%%%%%%%%%%%%% BASELINE %%%%%%%%%%%%%%%%%%%%%%%
function F = BaselineModel(x,freq)
F = x(2)*(freq-x(1))+x(3);


%%%%%%%%%%%%%%%%%%% INST UNITS CALC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MRS_struct] = MRSGABAinstunits(MRS_struct,ii)
% function [MRS_struct] = MRSGABAinstunits(MRS_struct)
% Convert GABA and Water amplitudes to institutional units
% (pseudo-concentration in mmol per litre).
% March 10: use MRS_struct.

PureWaterConc = 55000; % mmol/litre
WaterVisibility = 0.65; % This is approx the value from Ernst, Kreis, Ross
EditingEfficiency = 0.5;
T1_GABA = 0.80 ; % "empirically determined"...! Gives same values as RE's spreadsheet
% ... and consistent with Cr-CH2 T1 of 0.8 (Traber, 2004)
%Not yet putting in measured GABA T1, but it is in the pipeline - 1.35ish

T2_GABA = 0.13; % from occipital Cr-CH2, Traber 2004
T2_GABA = 0.088; % from JMRI paper 2011 Eden et al.

T1_Water = 1.100; % average of WM and GM, estimated from Wansapura 1999
T2_Water = 0.095; % average of WM and GM, estimated from Wansapura 1999
SiteFactor=0.779; % site specific scaling factor, CUBRIC=0.779
MM=0.45;  % MM correction: fraction of GABA in GABA+ peak. (In TrypDep, 30 subjects: 55% of GABA+ was MM)
%This fraction is platform and implementation dependent, base on length and
%shape of editing pulses and ifis Henry method. 
%
TR=1.8;
TE=0.068;
N_H_GABA=2;
N_H_Water=2;
Nspectra = length(MRS_struct.pfile);
%Nwateravg=8;

T1_factor = (1-exp(-TR./T1_Water)) ./ (1-exp(-TR./T1_GABA));
T2_factor = exp(-TE./T2_Water) ./ exp(-TE./T2_GABA);

MRS_struct.gabaiu(ii) = (MRS_struct.gabaArea(ii)  ./  MRS_struct.waterArea(ii))  ...
    * PureWaterConc*WaterVisibility*T1_factor*T2_factor*(N_H_Water./N_H_GABA) ...
    * MM * SiteFactor * (MRS_struct.Nwateravg ./ MRS_struct.Navg(ii)) ./ EditingEfficiency;

%%%%%%%%%%%%%%% INSET FIGURE %%%%%%%%%%%%%%%%%%%%%%%
function [h_main, h_inset]=inset(main_handle, inset_handle,inset_size)

% The function plotting figure inside figure (main and inset) from 2 existing figures.
% inset_size is the fraction of inset-figure size, default value is 0.35
% The outputs are the axes-handles of both.
%
% An examle can found in the file: inset_example.m
%
% Moshe Lindner, August 2010 (C).

if nargin==2
    inset_size=0.35;
end

inset_size=inset_size*.5;
%figure
new_fig=gcf;
main_fig = findobj(main_handle,'Type','axes');
h_main = copyobj(main_fig,new_fig);
set(h_main,'Position',get(main_fig,'Position'))
inset_fig = findobj(inset_handle,'Type','axes');
h_inset = copyobj(inset_fig,new_fig);
ax=get(main_fig,'Position');
set(h_inset,'Position', [1.3*ax(1)+ax(3)-inset_size 1.001*ax(2)+ax(4)-inset_size inset_size*0.7 inset_size*0.9])







