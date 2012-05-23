function MRSplotstack(MRS_struct, verticaloffset)
%function MRSplotstack(MRS_struct, verticaloffset)
%  MRS_struct: struct loaded by Gannet
%  specoffset: vertical offset of spectra when plotted - as a fraction of GABA peak
%              height
% Plots MRS data loaded by MRSLoadPfiles
% 110214:  Scale spectra by the peak _height_ of water
%          Plot multiple spectra as a stack - baselines offset
%            by mean height of GABA

% water~16300, glx ~17150, gaba~17700,  naa~18600, mm09~19500

numspec = length(MRS_struct.gabaspec(:,1));

% Find Water amplitude max, across all Pfiles
waterheight = abs(max(MRS_struct.waterspec,[],2));
heightrescale = repmat((1./waterheight), [1 length(MRS_struct.gabaspec(1,:))]);
SpectraToPlot = MRS_struct.gabaspec .* heightrescale;

% Estimate baseline from between Glx and GABA
specbaseline = mean(real(SpectraToPlot(:,17250:17650)),2);

% averaged gaba height across all scans - to estimate stack spacing
gabaheight = abs(max(SpectraToPlot(:,17250:18000),[],2));
gabaheight = mean(gabaheight);

% fix the offset between spectra here...
specoffset = gabaheight * verticaloffset;
%specoffset = 0; %overlay


plotstackoffset = [ 0 : (numspec-1) ]';
plotstackoffset = plotstackoffset * specoffset;
plotstackoffset = plotstackoffset - specbaseline;

SpectraToPlot = SpectraToPlot + ...
    repmat(plotstackoffset, [ 1  length(MRS_struct.gabaspec(1,:))]);

%figure(99)
plot(MRS_struct.freq, real(SpectraToPlot),'k');
legendtxt = regexprep(MRS_struct.pfile, '_','-');

%legend(legendtxt)
set(gca,'XDir','reverse');
oldaxis = axis;
% yaxis max = top spec baseline + 2*meangabaheight
yaxismax = 1.5*gabaheight + numspec *specoffset; % top spec + 2* height of gaba
yaxismin =  -gabaheight; % extend gaba height below zero
%yaxismin =  - 2* gabaheight; % extend 2* gaba heights below zero


%axis([0 5  yaxismin yaxismax])
axis([2.2 4.0  yaxismin yaxismax])

set(gca, 'YTick', [], 'XTick', [2.2 2.6 3.0 3.4 3.8], 'FontSize', 16); 


%, 'XTickLabel', ...
