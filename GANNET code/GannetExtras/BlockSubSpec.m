function [spec_blocks] = BlockSubSpec(GannetStruct, Nsubspec)
% Take pre-alignment OddFrames and Even Frames

odd = real(GannetStruct.OddFrames);
even = real(GannetStruct.EvenFrames);
diffspec = odd-even;

siz=size(odd);
Npts =siz(1);
Nspec=siz(2);
Nblocks = floor(Nspec/Nsubspec)
Nspecremain=mod(Nspec,Nsubspec);
Nspecused=Nspec-Nspecremain

diffspecused = diffspec(:,(1:Nspecused));

tmp = reshape(diffspecused, [ Npts Nsubspec Nblocks ]);
tmp = mean(tmp,2);
spec_blocks.gabaspec = squeeze(tmp)';
spec_blocks.Reference_compound = GannetStruct.Reference_compound;
spec_blocks.freq = GannetStruct.freq;
spec_blocks.vendor = GannetStruct.vendor;
spec_blocks.waterspec= GannetStruct.waterspec;
for ii=1:Nblocks
    spec_blocks.pfile{ii,:}=GannetStruct.pfile{1};
    spec_blocks.waterspec(ii,:)= GannetStruct.waterspec(1,:) .* (Nsubspec/Nspec);
end


MRSplotstack(spec_blocks)
title(['Subspectra per block = ' num2str(Nsubspec)]);