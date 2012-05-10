function [MRS_struct] = FitGABAFrames(MRS_struct,ii)
%
% Keep input and output structure names the same to add output data to
% the exisiting MRS data structure.
%Inputs:
% MRS_struct = structure with data loaded from MRSLoadPfiles
% ii = index of spectrum to analyse in MRS_struct

% From GannetFit 120510
FIT_LSQCURV = 0;
FIT_NLINFIT = 1;

fit_method = FIT_NLINFIT; %FIT_NLINFIT;


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
numscans=size(GABAData);
numscans=numscans(1);

size(GABAData(:,1))
% remove blank lines
for kk = numscans:-1:1
    if (GABAData(kk,1) == 0)
        GABAData(kk,:) = [];
    end
end

numsubspec = length(GABAData(:,1));

% For SNR, average subspectra in blocks of 8
subspecblock = 8;
npoints = length(GABAData(1,:));
nblocks = fix( numsubspec / subspecblock );
nremain = mod( numsubspec , subspecblock );


blockGABAdata = GABAData(1:(end-nremain),:); % discard remainder, until I think of a better way... 
blockGABAdata = blockGABAdata'; % to [points, subspectra]
size(blockGABAdata)

blockGABAdata = reshape(blockGABAdata, [npoints, subspecblock, nblocks]);
% sum across blocks
blockGABAdata = sum(blockGABAdata,2);
blockGABAdata = reshape(blockGABAdata, [npoints nblocks ]);

GABAData = blockGABAdata';% back to [subspectra, points]

%plot(real(GABAData(:,17000:18500))');

for jj = 1:nblocks
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
    lowerbound=find(min(z)==z);
    z=abs(MRS_struct.freq-2.79);%2.75
    upperbound=find(min(z)==z);
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
    FitError(jj)  =  100*std(residg)/GABAheight;
    % x(2) = 1/(2*sigma^2).  Convert from points to Hz

    % area under gaussian is a * (2 * pi * sigma^2)^(1/2), and
    % GaussModelParam(:,2) = 1 / (2 * sigma^2)
    % This sets GabaArea as the area under the curve.
    MRS_struct.GABA_area_subspec(ii,jj)=GaussModelParam(jj,1)./sqrt(-GaussModelParam(jj,2))*sqrt(pi);


    %figure(22)
    %plot(freq(freqbounds),GaussModel_area(GaussModelParam(jj,:),freq(freqbounds)),'r',...
    %freq(plotbounds),real(GABAData(jj,plotbounds)), 'b' )

%     legendtxt = regexprep(MRS_struct.pfile{ii}, '_','-');
%     title(legendtxt);
%     set(gca,'XDir','reverse');
%     %set(gca,'YTick',[], 'Xgrid', 'on');
%     oldaxis = axis;
%     axis( [2.6 3.6 oldaxis(3) oldaxis(4) ] )

end

% work out all of the fits


plotstart = 17000; plotend =17999;


for jj=1:nblocks
    % calculate fit, for plotting later
    GABAblockfit(jj,:) = GaussModel_area(GaussModelParam(jj,:),freq(plotstart:plotend));
end


startoffset=0;
offsetspec = mean(max(real(GABAData(1,plotstart:plotend))));

endoffset = startoffset + ...
    (nblocks-1) * offsetspec;
offsetvals = [startoffset:offsetspec:endoffset]'

offsetvals = repmat(offsetvals, [1 length([plotstart:plotend])]);
size(offsetvals)
size(GABAData(:,plotstart:plotend))

plotdiffspec = real(GABAData(:,plotstart:plotend)) + offsetvals;
GABAblockfit = GABAblockfit + offsetvals;

figure(23)
plot([plotstart:plotend],plotdiffspec','k', [plotstart:plotend], GABAblockfit, 'r')


%plot(plotdiffspec')


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



