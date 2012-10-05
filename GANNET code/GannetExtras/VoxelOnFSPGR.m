function VoxelOnFSPGR(pfile, ImagingDataDirectory, MRS_Rot_folder, FSPGR_dcm_folder, FSPGR_nifti)
% function VoxelOnFSPGR(pfile, ...
%                   ImagingDataDirectory,  ...
%                   MRS_Rot_folder, ...
%                   FSPGR_dcm_folder, ...
%                   FSPGR_nifti)
%   pfile and ImagingDataDirectory require full paths 
%   MRS_Rot_folder, FSPGR_dcm_folder, FSPGR_nifti are paths from ImagingDataDirectory
%
% e.g. MRSVoxelOnFSPGR('~/mrs/somato_motor/data/3143_6_P12800.7', ...
%       '/gpfs/home/sapje1/mrs/somato_motor/3143_081007-1', ...
%       'Series_00005', ...
%       'Series_00004_bak', ...
%       'fspgr.nii.gz')
% This is based on Richard version from his visit back to CUBRIC in Dec 09

% Dec 09: turn into a function, get voxel centre information from pfile header generated by print_raw_headerss
% 15 Jan 10: tidy up paths, voxel dimensions from pfile header,use pfile name in Mask filename output
%           Does voxel dim calc work for all rotator orientations?
% v100910: Sept 10: Fix crash when voxel is near edge of FSPGR

if nargin==0
    FSPGR_dcm_folder='/gpfs/home/sapje1/mrs/somato_motor/3143_081008-1/Series_00004_bak';
    MRS_Rot_folder='/gpfs/home/sapje1/mrs/somato_motor/3143_081008-1/Series_00005';
    FSPGR_nifti='~/mrs/somato_motor/3144_081007-1/3144fspgr.nii.gz'; % 15 Dec: use old one or generate new one? not sure yet
    pfile='~/mrs/somato_motor/data/3143_6_P12800.7'; 
    ImagingDataDirectory = '/gpfs/home/sapje1/mrs/somato_motor/3143_081007-1'
end

pfile;
MRS_Rot_folder=[ImagingDataDirectory '/' MRS_Rot_folder];
FSPGR_dcm_folder=[ImagingDataDirectory '/' FSPGR_dcm_folder];
FSPGR_nifti= [ImagingDataDirectory '/' FSPGR_nifti];

GenerateNifti = 0; % from the dcm files

%get from P file header!
%VoxelSize=[25 40 30];

%Tactical Plan
%1. Get voxel centre from MRS prescription
%2. Get FSPGR location (and rotation) from FSPGR dicom header
%3. Build Cubic mask in the coordinate frame of the FSPGR (unrotated).
%4. Get voxel rotation from MRS_Rot dicom header.
%5. Rotate Cubic mask in FSPGR space.


%1. Get voxel centre from MRS prescription
% CJE get this information from the P file header (*.7.hdr) generated by
% copyspectro, rather than from the .shf sage header file
ptr=fopen([pfile '.hdr' ]);
MRSHead=fscanf(ptr,'%c');

[m n]=size(MRSHead); 
fclose(ptr);

k = findstr(' user11:', MRSHead);
b=str2num(MRSHead(k+11:k+18));
k = findstr(' user12:', MRSHead);
c=str2num(MRSHead(k+11:k+18));
k = findstr(' user13:', MRSHead);
a=str2num(MRSHead(k+11:k+18));
Voxel_RAS_coord=[ -b -c a];
k = findstr(' user8:', MRSHead);
d=str2num(MRSHead(k+10:k+17));
k = findstr(' user9:', MRSHead);
e=str2num(MRSHead(k+10:k+17));
k = findstr(' user10:', MRSHead);
f=str2num(MRSHead(k+11:k+18));
VoxelSize = [f e d ]; % works for ob-axial rotator 

%2. Get FSPGR location (and rotation) from FSPGR header
% Get this from the dicom header prior to conversion
dcmfile_first = dir([FSPGR_dcm_folder '/*.1']); % get first dcm filename
FSPGRHead=dicominfo([FSPGR_dcm_folder '/' dcmfile_first.name]);
FSPGR_Coord1=FSPGRHead.ImagePositionPatient;
FSPGR_Rot=FSPGRHead.ImageOrientationPatient;
FSPGR_Rot=reshape(FSPGR_Rot',[3 2]);
FSPGR_Rot(:,3)=cross(FSPGR_Rot(:,1),FSPGR_Rot(:,2));
FSPGR_Rot(:,1:2)=FSPGR_Rot(:,1:2)/2;
FSPGR_Rot(1:2,:)=-FSPGR_Rot(1:2,:);

% get number of files (...so number of images)
dcmfiles = dir([FSPGR_dcm_folder '/*']);
numslices = size(dcmfiles);
numslices = numslices(1) - 2; %don't include  . and .. as files!
dcmfile_last = dir ([ FSPGR_dcm_folder '/*.' num2str(numslices) ]);

FSPGRHead=dicominfo([FSPGR_dcm_folder '/' dcmfile_last.name]);
FSPGR_Coord2=FSPGRHead.ImagePositionPatient;

FSPGR_SideCentre=((FSPGR_Coord1+FSPGR_Coord2)/2);
FSPGR_VolumeCentre=FSPGR_SideCentre+FSPGRHead.ImageOrientationPatient(1:3)*128+FSPGRHead.ImageOrientationPatient(4:6)*128;

sysout = ['mkdir ' ImagingDataDirectory '/MRSvoxelReg' ];
system(sysout);

savename = [ ImagingDataDirectory '/MRSvoxelReg/' FSPGRHead.StudyID '_FSPGR_Coord.mat'];
%savename=mat2str(savename)
save (savename, 'FSPGR_Coord1', 'FSPGR_Coord2', 'FSPGR_Rot', 'FSPGR_SideCentre', 'FSPGR_VolumeCentre');

if ( GenerateNifti )
    % CJE convert to nifti
    cd(FSPGR_dcm_folder);
    sysout = ['geprepfunct ' num2str(numslices) ' 1 ' ImagingDataDirectory '/MRSvoxelReg/' FSPGRHead.StudyID '_FSPGR']
    system(sysout)
    FSPGR_nifti = [ ImagingDataDirectory '/MRSvoxelReg/' FSPGRHead.StudyID '_FSPGR.nii.gz' ]
end

FSPGR_to_Vox=Voxel_RAS_coord-FSPGR_VolumeCentre';
FSPGR_to_Vox_FSL=FSPGR_to_Vox';

%January fiddle step.
FSPGR_to_Vox_FSL(2)=-FSPGR_to_Vox_FSL(2);

Voxel_Centre_FSL=FSPGR_to_Vox_FSL+[128 128 (numslices/2) ]'; 

%SKIP STEP 3 AND GO TO STEP 4 - VOXEL ROTATION.
files = dir(MRS_Rot_folder);
MRSRotHead=dicominfo([MRS_Rot_folder '/' files(3).name]); 
MRS_Rot=MRSRotHead.ImageOrientationPatient;
MRS_Rot=reshape(MRS_Rot',[3 2]);
%MRS_Rot(:,2)=-MRS_Rot(:,2);
%MRS_Rot(1:2,1:2)=MRS_Rot(1:2,1:2);
MRS_Rot(:,3)=cross(MRS_Rot(:,1),MRS_Rot(:,2))
MRS_Rot(1,:)=-MRS_Rot(1,:);
%MRS_Rot=MRS_Rot'
MRS_Rot=-MRS_Rot;
%MRS_Rot(2,1)=-MRS_Rot(2,1);
%MRS_Rot(1,2)=-MRS_Rot(1,2);
%MRS_Rot(3,:)=-MRS_Rot(3,:);
%MRS_Rot(:,3)=-MRS_Rot(:,3);


%Build a mask
Mask=zeros(256,256,numslices);
AxesV=zeros(256,256,256);

% CJE sept 10: need to make sure voxel doesn't poke out of the FSPGR...

ii_min=fix(Voxel_Centre_FSL(1)-50);
if(ii_min<1) 
   ii_min=1;
end
ii_max=fix(Voxel_Centre_FSL(1)+50);
if(ii_max>256)
   ii_max=256;
end

jj_min=fix(Voxel_Centre_FSL(2)-50);
if(jj_min<1) 
   jj_min=1;
end
jj_max=fix(Voxel_Centre_FSL(2)+50);
if(jj_max>256)
   jj_max=256;
end

kk_min=fix(Voxel_Centre_FSL(3)-25);
if(kk_min<1) 
   kk_min=1;
end
kk_max=fix(Voxel_Centre_FSL(3)+25);
if(kk_max>256)
   kk_max=256;
end

for ii=ii_min:ii_max % 50~ 30*sqrt(3)
    for jj=jj_min:jj_max
        for kk=kk_min:kk_max
            Coord_from_Vox_Center=[ii jj kk]'-Voxel_Centre_FSL;
            %Need to take this into the Coordinate frame of the Voxel
            %New_Coord=FSPGR_Rot*Coord_from_Vox_Center;
            New_Coord=[[-1 0 0];[0 -1 0];[0 0 1]]*Coord_from_Vox_Center;
            
            New_Coord=inv(MRS_Rot)*New_Coord;
            Mask(ii,jj,kk)=(New_Coord(1)<=VoxelSize(1)/2)*(New_Coord(2)<=VoxelSize(2)/2)*(New_Coord(3)<=VoxelSize(3)/2)*(New_Coord(1)>-VoxelSize(1)/2)*(New_Coord(2)>-VoxelSize(2)/2)*(New_Coord(3)>-VoxelSize(3)/2);
            %Mask(ii,jj,kk)=Mask(ii,jj,kk)*(New_Coord(1)>-VoxelSize(1)/2)*(New_Coord(2)>-VoxelSize(2)/2)*(New_Coord(3)>-VoxelSize(3)/2);
            %AxesV(ii,jj,kk)=Mask(ii,jj,kk)*((abs(New_Coord(1))<0.5)*(abs(New_Coord(2))<0.5)+((abs(New_Coord(2))<0.5)*(abs(New_Coord(3))<0.5))+((abs(New_Coord(3))<0.5)*(abs(New_Coord(1))<0.5)));
        end
    end
end


FSPGR_nifti
FSPGR=load_nifti(FSPGR_nifti);
ViewVolume=circshift(Mask,[0 0 0])*1000+FSPGR.vol;
 
MaskOutput=FSPGR;

% work out exam & series from last '/' and '_' in pfilename
slsh=findstr('/', pfile);
if isempty(slsh) == 1
   slsh=0;
end
slsh=slsh(end);
undersc=findstr('_', pfile);
undersc=undersc(end);
exam_series_ID=pfile((slsh+1):(undersc-1));


%Arbitrary shift 'up'
MaskOutput.vol=circshift(Mask,[0 0 0]);
save_nifti(MaskOutput,[ImagingDataDirectory '/MRSvoxelReg/' exam_series_ID '_Mask.nii.gz'] );
%save_nifti(MaskOutput,[ImagingDataDirectory '/MRSvoxelReg/Mask.nii.gz'] );

%figure(1)
%mimage(squeeze(ViewVolume(:,fix(Voxel_Centre_FSL(2)),end:-1:1))');
% mimage(squeeze(ViewVolume(:,279,end:-1:1))');   %Choose the relevant slice from FSLview...
% set(gca,'PlotBoxAspectRatio',[256 120 1]);
% set(gca,'DataAspectRatio',[2 1 1]);
% figure(2);
% %mimage(squeeze(ViewVolume(fix(Voxel_Centre_FSL(1)),:,end:-1:1))');
% mimage(squeeze(ViewVolume(163,:,end:-1:1))');
% set(gca,'PlotBoxAspectRatio',[256 120 1]);
% set(gca,'DataAspectRatio',[2 1 1]);
% figure(3)
% %mimage(squeeze(ViewVolume(:,end:-1:1,fix(Voxel_Centre_FSL(3))))');
% mimage(squeeze(ViewVolume(:,end:-1:1,97))');


% %Grey Matter, White Matter, CSF
% 
% % Grey=load_nifti(['../' ExamNo '/Segment_prob_1.nii.gz']);
% % White=load_nifti(['../' ExamNo '/Segment_prob_2.nii.gz']);
% % CSF=load_nifti(['../' ExamNo '/Segment_prob_0.nii.gz']);
% % 
% % GreyVox=Grey.vol.*Mask;
% % WhiteVox=White.vol.*Mask;
% % CSFVox=CSF.vol.*Mask;
% % 
% % Voxel_Segmentation = [sum(CSFVox(:))/sum(Mask(:)) sum(GreyVox(:))/sum(Mask(:)) sum(WhiteVox(:))/sum(Mask(:))]