function [Summary,SleepStage]=SleepSEEG(FileList, patID)
% function [Summary,SleepStage]=SleepSEEG(FileList,ExtraFiles)
    %changing path stuff to work in ubuntu
    addpath /mnt/ernie_main/shared_toolboxes/MatlabProgressBar/
    addpath /mnt/ernie_main/shared_toolboxes/shared_toolboxes/Derek_functions/
    
    % set the patient ID so that we can adaptively get the channels used
    % spreadsheet
    % patID = "Spat36";
    spreadsheet_root_dir = "/mnt/ernie_main/000_Data/SEEG/SEEG_EyesClosed_RestingState/data/";
    subj_resting_state_dir = get_files_in_dir_linux(spreadsheet_root_dir, patID+"*", 1);
    
    if length(subj_resting_state_dir) ~= 1
        error("Unable to find resting state directory");
    else
        spreadsheet_fpath = sprintf("%s/%s_Channels_Used.xlsx", subj_resting_state_dir, patID);
    end

 % [Summary,SleepStage]=SleepSEEG(FileList,ExtraFiles)
 %
 % SleepSEEG performs automatic sleep scoring of intracranial EEG files.
 % The file SleepSEEG_models.mat must be in the same folder as this file.
 % It requires a complete night of iEEG data.
 % 
 % Input (optional): File name or cell array of file names of one complete 
 % night. If excluded, the file selection will be done interactively. It is
 % possible to select one or many files, they do not need to be in any 
 % particular order. Select only files for 8-12 hours of 1 complete night. 
 % ExtraFiles: File name or cell array of file names of extra segments to
 % be scored (e.g. short naps). IF empty they can be chosen interactively.
 % As many files as desired can be selected, or no extra files.
 % 
 % The iEEG files format must be EDF or EDF+C (continuous). If the format
 % is EDF+D (discontinuous), you must use other software to convert it to
 % EDF+C (e.g., EDFbrowser). Accepted sampling rates are 200, 256, 500,
 % 512, 1000, 1024, 2000, or 2048 samples per second. Different files can
 % have different sampling rate from this set. If the sampling rate of the
 % input files is not in this list, a sampling rate conversion must be done
 % with another software. All the different files must have the same 
 % channels, with the same name and in the same order. It is assumed that 
 % the channels are the recording channels in a referential montage. 
 % SleepSEEG needs a bipolar montage and will attemt to form bipolar
 % channels, assuming that the channel names end in a number indicating the
 % contact. The user has to select from a list of channels all the ones
 % that are bipolar iEEG channels, and exclude the ones that are not.
 %
 % Output variables:
 % Summary: A cell array with the list of transitions bewtween different
 % sleep stages. First column is tha file name, second column the date,
 % thrid column the time of the beginning of the sleep stage, fourth column
 % the sleep stage, and fifth column the number of consecutive epochs that
 % this sleep stage is present.
 % SleepStage: A list of the stage for each epoch. First column the file
 % index; second column the date and time in numeric (matlab) format; third
 % column the sleep stage (1 to 5, Stages R, w, N1, N2, N3); fourth column
 % the confidence (the posterior probability of the stage with highest 
 % posterior probability, higher than 0.5 means high confidence).
 %
 % SleepSEEG v1.0
 % Nicolás von Ellenrieder - 2021-12-22
 % nicolas.vonellenrieder@mcgill.ca

% Select version: 
% BEA -blind to epileptic activity, if selected, all channels can be
% incldued (default)
% EEA -exccluding epileptic activity, only channels with IED rate lower
% than 1 per minute should be included (use Janca's automatic detector  
% doi: 10.1007/s10548-014-0379-1 to compute the IED rate).
version='BEA';

% Select night files
if nargin==0
    [FileList,FilePath]=uigetfile('*.SIG;*.EDF;*.REC','Select files of one complete night','MultiSelect','on');
    if iscell(FileList)
        for ii=1:length(FileList), FileList{ii}=[FilePath FileList{ii}]; end
    else
        FileList=[FilePath FileList];
    end
end
if iscell(FileList), file=FileList; else, file=cell(1); file{1}=FileList; end

% DJD Mod -> I don't want to select extra files
ExtraFiles = 0;
ExtraPath = 0;
% % Select extra files
% if nargin<2
%     [ExtraFiles,ExtraPath]=uigetfile('*.SIG;*.EDF;*.REC','Select all extra files to process, e.g. naps','MultiSelect','on');
%     if ~all(ExtraFiles==0)
%         if iscell(ExtraFiles)
%             for ii=1:length(ExtraFiles), ExtraFiles{ii}=[ExtraPath ExtraFiles{ii}]; end
%         else
%             ExtraFiles=[ExtraPath ExtraFiles];
%         end
%     end
% end

sf=size(file);
[nnf,sf]=max(sf);
if ~all(ExtraFiles==0)
    if iscell(ExtraFiles)
        file=cat(sf,file,ExtraFiles);
    else
        file{nff+1}=ExtraFiles;
    end
end
% Read first file and get number of channel
FileName=file{1};
fp=fopen(FileName,'r','ieee-le');
header_ini=char(fread(fp,256,'uchar')');
Nch=str2double(header_ini(253:256));
channel_info.name=char(fread(fp,[16,Nch],'char')');
fclose(fp);
% Verify files are not EDF+D and have the same channels
for nf=1:length(file)
    FileName=file{nf};
    fp=fopen(FileName,'r','ieee-le');
    header_ini=char(fread(fp,256,'uchar')');
    if header_ini(197)=='D'
        disp([FileName ' is an EDF+D file (discontinuous), which is not supported by IntraSleep.']);
        disp('It is possible to convert it to EDF+C with e.g., EDF Browser free software.');
        SleepStage=[]; Summary=[];         
        return;
    end
    num_channels=str2double(header_ini(253:256));
    nam=char(fread(fp,[16,num_channels],'char')');
    if Nch~=num_channels
        disp([FileName ' does not have the same number of channels as ' file{1} '. Unable to continue.']);
        SleepStage=[]; Summary=[];
        return;
    elseif any(any(nam~=channel_info.name))
        disp([FileName ' does not have the same exact channels as ' file{1} '. Unable to continue.']);
        SleepStage=[]; Summary=[];
        return;
    end
end
% Select channels
nmax=0; ele=cell(Nch,1); num=zeros(Nch,1);
for ii=1:Nch
    chan=deblank(channel_info.name(ii,:));
    sp=strfind(chan,' ');
    if ~isempty(sp), chan=chan(sp+1:end); end
    if isnan(str2double(chan(end)))
        ele{ii}=chan; num(ii)=NaN;
    elseif isnan(str2double(chan(end-1)))
        ele{ii}=chan(1:end-1); num(ii)=str2double(chan(end)); 
    else
        ele{ii}=chan(1:end-2); num(ii)=str2double(chan(end-1:end));
    end
    nmax=max(nmax,num(ii));
end
eles=unique(ele,'stable');
% Create bipolars
MM=zeros(Nch); nc=0; ChL=cell(Nch,1);
for ii=1:length(eles)
    for jj=2:nmax
        k=find(strcmp(eles{ii},ele)&num==jj-1,1);
        ke=find(strcmp(eles{ii},ele)&num==jj,1);
        if ~isempty(k)&&~isempty(ke)
            nc=nc+1; MM(nc,ke)=-1; MM(nc,k)=1; 
            ChL{nc}=[eles{ii} num2str(jj-1) '-' eles{ii} num2str(jj)];
        end
    end
end
ChL=ChL(1:nc); MM=MM(1:nc,:);
% old code with GUI interface
% [Channels,l]=listdlg('ListString', ChL,'PromptString','Select channels');

% DJD MOD
% I know what the bipolar montages are, so I'm going to automate this part
% and have it pull directly from the channels used excel sheet
% read in and parse the excel sheet

ch_used_excel = readtable(spreadsheet_fpath, 'ReadVariableNames', false);
ch_used_excel = string(ch_used_excel{:,1});
ch_used_excel = strrep(ch_used_excel, " - ", "-");

% convert all the channels from the edf file into strings
ChL_new = string(ChL);
Channels = [];

% Figure out which bipoles to keep
for pp = 1:length(ChL_new)
    if sum(ChL_new(pp) == ch_used_excel) > 0 
        Channels(end+1) = pp;
    end
end

% set l=1 to denote a success
if sum(Channels) > 0
    l = 1;
end
% DJD MOD

if l==0; return;  end
MM=MM(Channels,:);
Nch=nnz(Channels);
Nfeat=24;    
ne=0;
feature=zeros(Nfeat,Nch,1500);
SleepStage=zeros(1500,5);
night=zeros(1500,1)==0;
for nf=1:length(file)
    FileName=file{nf};
    fp=fopen(FileName,'r','ieee-le');
    header_ini=char(fread(fp,256,'uchar')');
    len=str2double(header_ini(185:192));
    num_records=str2double(header_ini(237:244));
    record_duration=str2double(header_ini(245:252));
    num_channels=str2double(header_ini(253:256));
    channel_info.name=char(fread(fp,[16,num_channels],'char')');
    channel_info.sensor=char(fread(fp,[80,num_channels],'char')');
    channel_info.unit=char(fread(fp,[8,num_channels],'char')');
    channel_info.physmin =str2num(char(fread(fp,[8,num_channels],'char')')); %#ok<*ST2NM>
    channel_info.physmax =str2num(char(fread(fp,[8,num_channels],'char')'));
    channel_info.digimin =str2num(char(fread(fp,[8,num_channels],'char')'));
    channel_info.digimax=str2num(char(fread(fp,[8,num_channels],'char')'));
    channel_info.preproc=char(fread(fp,[80,num_channels],'char')');
    channel_info.samples_per_record=str2num(char(fread(fp,[8,num_channels],'char')'));
    samp=channel_info.samples_per_record(1);
    fs=round(samp/record_duration);
    a=(channel_info.physmax-channel_info.physmin)./(channel_info.digimax-channel_info.digimin);
    b=(channel_info.physmin-a.*channel_info.digimin)';
    b(a<0)=0; a(a<0)=1;
    NumSamps=num_records*samp;
    StartTime=[header_ini(99:109) ' ' header_ini(177:178) ':' header_ini(180:181) ':' header_ini(183:184)];
    Time=datevec(StartTime);
    % Set start time of first epoch at :29 or :59 for compatibility with Stellate
    if Time(6)<27.75
        sta=fs*(29-Time(6));
        Time(6)=29;
    elseif Time(6)<57.75
        sta=fs*(59-Time(6));
        Time(6)=59;
    else
        sta=fs*(60-Time(6)+29);
        Time(6)=29;
    end
    sta=round(sta);
    timenum=datenum(Time);
    Ne=(NumSamps-sta)/fs/30;
    if rem(Ne,1)<1.25/30, Ne=Ne-.05; end
    Ne=floor(Ne);
    SleepStage(ne+(1:Ne),1:2)=[ones(Ne,1)*nf timenum+(0:Ne-1)'/2880];
    % Deal with possible channels with different number of samples per record
    keep_channels=channel_info.samples_per_record==mode(channel_info.samples_per_record);
    num_ch=nnz(keep_channels);
    keep=[];
    for ii=1:num_channels
        if keep_channels(ii)
            keep=cat(1,keep,ones(channel_info.samples_per_record(ii),1)==1);
        else
            keep=cat(1,keep,ones(channel_info.samples_per_record(ii),1)<0);
        end
    end
    disp(['File: ' FileName]);
    disp([num2str(Ne) ' epochs']);
    disp('Computing features...');
    % read buffer
    buffer=ceil((sta+1.25*fs)/samp); 
    fseek(fp,len,-1);
    X=fread(fp,sum(channel_info.samples_per_record)*buffer,'int16');
    X=reshape(X(logical(repmat(keep,buffer,1))),[samp num_ch buffer]);
    X=permute(X,[1 3 2]);
    X=reshape(X,samp*buffer,num_ch);
    X=X(sta-1.25*fs+1:end,:);
    M=diag(a(keep_channels))*MM(:,keep_channels)';
    Xb=(X-repmat(b(keep_channels),size(X,1),1))*M;
    % Compute features
    for ii = progress(1:Ne)
        dur=ceil((fs*30-size(Xb,1)+2.5*fs)/samp);
        X=fread(fp,sum(channel_info.samples_per_record)*dur,'int16');
        X=reshape(X(logical(repmat(keep,dur,1))),[samp num_ch dur]);
        X=permute(X,[1 3 2]);
        X=reshape(X,samp*dur,num_ch);
        X=(X-repmat(b(keep_channels),size(X,1),1))*M;
        X=cat(1,Xb,X);
        Xb=X(30*fs+1:end,:);
        X=X(1:32.5*fs,:);
        X=change_sampling_rate(X,fs);
        X(1:64,:)=[]; X(8193:end,:)=[];
        X=X-repmat(mean(X),8192,1);
        ne=ne+1;
        for nc=1:Nch
            feature(:,nc,ne)=compute_features(X(:,nc));
        end
        night(ne)=nf<=nnf;
        SleepStage(ne,5)=sta+30*fs*(ii-1);
        % disp(Ne-ii);
    end
    fclose(fp);
end
feature=feature(:,:,1:ne);
night=night(1:ne);
feature=permute(feature,[2 3 1]);
SleepStage=SleepStage(1:ne,:);
% Preprocess features
[x,y]=find(sum(isnan(feature)|isinf(feature),3)>0);
for ii=1:length(x), feature(x(ii),y(ii),:)=NaN; end
featfeat=zeros(Nfeat,Nch);
for nch=1:Nch
    for nf=1:Nfeat
        % outlier detection
        f=feature(nch,:,nf);
        in=isnan(f);
        f(in)=[];
        f=f-movmean(f,10);
        fs=sort(f);
        m=fs(round([.25 .75]*size(fs,2)));
        del=isnan(feature(nch,:,nf));
        del(~in)=f<m(1)-2.5*(m(2)-m(1))|f>m(2)+2.5*(m(2)-m(1))|isnan(f);
        feature(nch,del,nf)=NaN;
        % Smoothing
        f=feature(nch,:,nf);
        ff=f;
        in=isnan(ff);
        ff(in)=[];
        ff=movmean(ff,3);
        f(~in)=ff;
        % normalizing
        if any(f), f=(f-nanmean(f(night)))./nanstd(f(night)); end
        feature(nch,:,nf)=f;
        % get features' coordinates
        f=feature(nch,night,nf);
        in=isnan(f);
        f(in)=[];
        fm=movmean(f,10);
        f=norm(f);
        if f, featfeat(nf,nch)=norm(fm)/norm(f); end
    end
end
featfeat=featfeat';

load SleepSEEG_models Model_EEA Model_BEA GC_EEA GC_BEA;
if strcmp(version,'EEA'), Model=Model_EEA; GC=GC_EEA; else, Model=Model_BEA; GC=GC_BEA; end

% Determine cluster for each channel
ch_gr=zeros(Nch,1);
Ng=size(Model,1); Np=size(feature,2);
for nc=1:Nch
    [~,ch_gr(nc)]=min(sum((repmat(featfeat(nc,:),Ng,1)-GC).^2,2));
end
% Score channels
postprob=zeros([Nch Np 7]);
for gc=1:Ng
    ik=ch_gr==gc;
    if nnz(ik)>0
        fe=reshape(feature(ik,:,:),nnz(ik)*Np,Nfeat);
        del=sum(isnan(fe),2)==Nfeat;
        prob=zeros(nnz(ik)*Np,7)+NaN;
        [~,p]=predict(Model{gc,1},fe(~del,:));
        prob(~del,1:4)=p;
        [~,p]=predict(Model{gc,2},fe(~del,:));
        prob(~del,5)=p(:,1);
        [~,p]=predict(Model{gc,3},fe(~del,:));
        prob(~del,6)=p(:,1);
        [~,p]=predict(Model{gc,4},fe(~del,:));
        prob(~del,7)=p(:,1);
        postprob(ik,:,:)=reshape(prob,[nnz(ik),Np,7]);
    end
end
% Combine channels
prop=squeeze(nanmean(postprob,1));
prop(:,1:4)=prop(:,1:4)./repmat(sum(prop(:,1:4),2),[1 4]);
confidence=prop;

% Define stage
[mm,sa]=max(confidence(:,1:4),[],2);
sa(isnan(mm))=2;
mm(isnan(mm))=0;
sa(sa>2)=sa(sa>2)+1;
sa(sa==1&confidence(:,5)<.5)=3;
sa(sa==2&confidence(:,6)<.5)=3;
sa(sa==4&confidence(:,7)<.5)=3;

% Output variables
SleepStage(:,3:4)=[sa mm];
[~,ix]=sort(SleepStage(:,2));
SleepStage=SleepStage(ix,:);
stagename={'R','W','N1','N2','N3'};
ic=[true; SleepStage(1:end-1,3)~=SleepStage(2:end,3)|SleepStage(1:end-1,1)~=SleepStage(2:end,1)];
Sl=SleepStage(ic,:);
icc=find(ic); icc(end+1)=size(ic,1)+1; icc=icc(2:end)-icc(1:end-1); 
Summary=cell(size(icc,1)+1,5);
a=datestr(Sl(1,2));
Summary(1,:)={'File','Date','Time','Sleep stage','# of epochs'};
Summary(2,:)={file{Sl(1,1)},a(1:11),a(13:20),stagename{Sl(1,3)},icc(1)};
for ii=2:size(icc,1)
    a=datestr(Sl(ii,2));
    Summary(ii+1,:)={file{Sl(ii,1)},a(1:11),a(13:20),stagename{Sl(ii,3)},icc(ii)};
end
for ii=size(icc,1)+1:-1:3
    if strcmp(Summary{ii,1},Summary{ii-1,1}), Summary{ii,1}=' '; end
    if strcmp(Summary{ii,2},Summary{ii-1,2}), Summary{ii,2}=' '; end
end

end


function X=change_sampling_rate(X,fs)

od5=54; od2=26; ou2=34;
% hd5=firpm(od5,[0 1/8 1/4 1],[1 1 0 0]); % Defined explicitly below for speed.
% hd2=firpm(od2,[0 3/8 5/8 1],[1 1 0 0]); % Defined explicitly below for speed.
% hu2=firpm(ou2,[0 1/5 2/5 1],[1 1 0 0]); % Defined explicitly below for speed.
hd5 = [-0.000413312132792   0.000384910656353   0.000895384486596   0.001426584098180   0.001572675788393...
        0.000956099017099  -0.000559378457343  -0.002678217568221  -0.004629975982837  -0.005358589238386...
       -0.003933117464092  -0.000059710059922   0.005521319363883   0.010983495478404   0.013840996082966...
        0.011817315106321   0.003905283425021  -0.008768844009700  -0.022682212400564  -0.032498023687148...
       -0.032456772047175  -0.018225658085891   0.011386634156651   0.053456542440034   0.101168250947271...
        0.145263694388270   0.176384224234024   0.187607302744229   0.176384224234024   0.145263694388270...
        0.101168250947271   0.053456542440034   0.011386634156651  -0.018225658085891  -0.032456772047175...
       -0.032498023687148  -0.022682212400564  -0.008768844009700   0.003905283425021   0.011817315106321...
        0.013840996082966   0.010983495478404   0.005521319363883  -0.000059710059922  -0.003933117464092...
       -0.005358589238386  -0.004629975982837  -0.002678217568221  -0.000559378457343   0.000956099017099...
        0.001572675788393   0.001426584098180   0.000895384486596   0.000384910656353  -0.000413312132792];
hd2 = [ 0.001819877350937   0.000000000000000  -0.005222562671417  -0.000000000000000   0.012064004143824...
        0.000000000000000  -0.024375517671448  -0.000000000000001   0.046728257943321   0.000000000000001...
       -0.095109546093322  -0.000000000000001   0.314496714630999   0.500000000000001   0.314496714630999...
       -0.000000000000001  -0.095109546093322   0.000000000000001   0.046728257943321  -0.000000000000001...
       -0.024375517671448   0.000000000000000   0.012064004143824  -0.000000000000000  -0.005222562671417...
        0.000000000000000   0.001819877350937];
hu2 = [-0.000436344318144   0.000871220095338   0.002323148230146   0.002161775863027  -0.001491947315250...
       -0.006921578214988  -0.008225251275241  -0.000175549247298   0.014346107990165   0.022307211309135...
        0.009627222866552  -0.023067228777693  -0.052048618204333  -0.041386934088038   0.030232588120746...
        0.146540456844067   0.255673093905358   0.300317753050023   0.255673093905358   0.146540456844067...
        0.030232588120746  -0.041386934088038  -0.052048618204333  -0.023067228777693   0.009627222866552...
        0.022307211309135   0.014346107990165  -0.000175549247298  -0.008225251275241  -0.006921578214988...
       -0.001491947315250   0.002161775863027   0.002323148230146   0.000871220095338  -0.000436344318144];
fs=round(fs);
switch fs
    case 200
        z=zeros(ou2/2,size(X,2));
        X=2*reshape(permute(cat(3,X,0*X),[3 1 2]),2*size(X,1),size(X,2)); % upsample x2
        X=filter(hu2,1,[X;z]); X(1:ou2/2,:)=[];
        z=zeros(od5/2,size(X,2));
        X=4*reshape(permute(cat(3,X,0*X,0*X,0*X),[3 1 2]),4*size(X,1),size(X,2)); % upsample x4
        X=filter(hd5,1,[X;z]); X(1:od5/2,:)=[];
        X=X(1:5:end,:); % downsample x5
        X=4*reshape(permute(cat(3,X,0*X,0*X,0*X),[3 1 2]),4*size(X,1),size(X,2)); % upsample x4
        X=filter(hd5,1,[X;z]); X(1:od5/2,:)=[];
        X=X(1:5:end,:); % downsample x5
    case 256
        % No change necessary.
    case 500
        z=zeros(od5/2,size(X,2));
        X=4*reshape(permute(cat(3,X,0*X,0*X,0*X),[3 1 2]),4*size(X,1),size(X,2)); % upsample x4
        X=filter(hd5,1,[X;z]); X(1:od5/2,:)=[];
        X=X(1:5:end,:); % downsample x5
        X=4*reshape(permute(cat(3,X,0*X,0*X,0*X),[3 1 2]),4*size(X,1),size(X,2)); % upsample x4
        X=filter(hd5,1,[X;z]); X(1:od5/2,:)=[];
        X=X(1:5:end,:); % downsample x5
        X=4*reshape(permute(cat(3,X,0*X,0*X,0*X),[3 1 2]),4*size(X,1),size(X,2)); % upsample x4
        X=filter(hd5,1,[X;z]); X(1:od5/2,:)=[];
        X=X(1:5:end,:); % downsample x5
    case 512
        z=zeros(od2/2,size(X,2));
        X=filter(hd2,1,[X;z]); X(1:od2/2,:)=[];        
        X=X(1:2:end,:); % downsample x2
    case 1000        
        z=zeros(od5/2,size(X,2));
        X=2*reshape(permute(cat(3,X,0*X),[3 1 2]),2*size(X,1),size(X,2)); % upsample x2
        X=filter(hd5,1,[X;z]); X(1:od5/2,:)=[];
        X=X(1:5:end,:); % downsample x5
        X=4*reshape(permute(cat(3,X,0*X,0*X,0*X),[3 1 2]),4*size(X,1),size(X,2)); % upsample x4
        X=filter(hd5,1,[X;z]); X(1:od5/2,:)=[];
        X=X(1:5:end,:); % downsample x5
        X=4*reshape(permute(cat(3,X,0*X,0*X,0*X),[3 1 2]),4*size(X,1),size(X,2)); % upsample x4
        X=filter(hd5,1,[X;z]); X(1:od5/2,:)=[];
        X=X(1:5:end,:); % downsample x5
    case 1024
        z=zeros(od2/2,size(X,2));
        X=filter(hd2,1,[X;z]); X(1:od2/2,:)=[];        
        X=X(1:2:end,:); % downsample x2
        X=filter(hd2,1,[X;z]); X(1:od2/2,:)=[];        
        X=X(1:2:end,:); % downsample x2
    case 2000      
        z=zeros(od5/2,size(X,2));
        X=filter(hd5,1,[X;z]); X(1:od5/2,:)=[];
        X=X(1:5:end,:); % downsample x5
        X=4*reshape(permute(cat(3,X,0*X,0*X,0*X),[3 1 2]),4*size(X,1),size(X,2)); % upsample x4
        X=filter(hd5,1,[X;z]); X(1:od5/2,:)=[];
        X=X(1:5:end,:); % downsample x5
        X=4*reshape(permute(cat(3,X,0*X,0*X,0*X),[3 1 2]),4*size(X,1),size(X,2)); % upsample x4
        X=filter(hd5,1,[X;z]); X(1:od5/2,:)=[];
        X=X(1:5:end,:); % downsample x5
    case 2048
        z=zeros(od2/2,size(X,2));
        X=filter(hd2,1,[X;z]); X(1:od2/2,:)=[];        
        X=X(1:2:end,:); % downsample x2
        X=filter(hd2,1,[X;z]); X(1:od2/2,:)=[];        
        X=X(1:2:end,:); % downsample x2
        X=filter(hd2,1,[X;z]); X(1:od2/2,:)=[];        
        X=X(1:2:end,:); % downsample x2        
    otherwise
        disp(['Error: Input sampling rate of ' num2str(fs) ' samples per second not supported.']);
        disp('Supported sampling rates are 200,256,500,512,1000,1024,2000, or 2048 samples per second.');
        clear X;
end
end


function f=compute_features(x)

J=8; % scales, so that the widest is 0.5-1 Hz (since the sampling rate is 256).
nd=5;
% 1 are the scale function cofficients
l=[.22641898 .85394354 1.02432694 .19576696 -.34265671 -.04560113 .10970265 -.0088268 -.01779187 .00471742793]';
h=flipud(l).*(-1).^(0:2*nd-1)';
x=x(:);
x=x(1:2^J*floor(size(x,1)*2^-J));
laa=zeros(size(x,1)/2,1); 
C=zeros(3,J);
b=ones(J,1);
mwc=zeros(3*J,1);
for jj=1:J
    % Compute the wavelet leaders
    d=filter(h,1,x); lea=d(1:2:end);
    x=filter(l,1,x); x=x(1:2:end);
    mm=abs(lea(nd+ceil(256/2^jj):end-max(ceil(256/2^jj)-nd,0)));
    mwc(jj)=log(mean(mm)); 
    mwc(jj+J)=sum(mm.^2);
    mwc(jj+2*J)=log(std(mm));
    lea=abs([0; lea; 0]); lea=max(max(lea(1:end-2),lea(2:end-1)),lea(3:end));
    lea=max(lea,laa); laa=max(lea(1:2:end),lea(2:2:end));
    % Get cumulants of ln leaders
    lea=lea(nd+ceil(256/2^jj):end-max(ceil(256/2^jj)-nd,0)); % Here the transients are discarded 
    le=log(lea);
    u1=mean(le); u2=mean(le.^2); u3=mean(le.^3);
    C(:,jj)=[u1; u2-u1^2; u3-3*u1*u2+2*u1^3];
    nj=length(lea);
    b(jj)=nj; % If we want weighted fit for c_p
end
sc=2:6; % Selected scales
C=C(:,sc); b=b(sc);
V0=sum(b); V1=sc*b; V2=sc.^2*b;
w=b.*((V0*sc'-V1)./(V0*V2-V1^2));
f=[log2(exp(1))*C*w; mwc];
f([4 4+J 4+2*J])=[]; % Exclude the 64-128 Hz scale
J=J-1;
f(J+4:2*J+3)=log10(f(J+4:2*J+3)/sum(f(J+4:2*J+3)));
end
