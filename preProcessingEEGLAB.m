% In this Script we will load the EEG data and will do some processing in
% EEGlab.
% This code is runnable for healty subjects now for pateinets we should just
% chnage the name of variables that is commented in the code (MS).
% Before run this script You should run eeglab in your Matlab one time
% Also, you should have the channel location to run the
% code(Copy_of_32channellocation.ced)
% MS -2021
clc
clear
close all
SR=512;
% location of Stimuli for all the three pathways
location_Magno_trigger=[1,3,5,8,9,18,19,22,29,30,32,35,36,40,41];
location_Konio_trigger=[2,4,6,10,11,15,16,21,23,24,27,31,33,38,39,43,47];
location_Parvo_trigger=[7,12,13,14,17,20,25,26,28,34,37,42,44,45,46];
number_of_sample=612;
number_of_channel=9;
number_of_visualpathway=3;
%%
% Merging and preprocessing data for Healthy subjects and People living
% with MS


for subjectsNum=1:16
    
for session='1':'2'
lowband=0.4;
highband=40;

%str=['HC' num2str(subjectsNum) '_' session '.mat'];
str=['MS' num2str(subjectsNum) '_' session '.mat'];
SR=512;
EEG.etc.eeglabvers = '2021.0'; % this tracks which version of EEGLAB is being used, you may ignore it
%   importdata
EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',str,'srate',SR,'pnts',0,'xmin',0);
% import event
EEG_without_event = pop_chanevent(EEG, 35,'edge','leading','edgelen',0);
% remove time channel
EEG_removed_extra_channels = pop_select( EEG_without_event, 'nochannel',1);
% remove Reference channel
EEG_removedChannel = pop_select( EEG_removed_extra_channels, 'nochannel',33);


% load channel location
EEG_withLocation=pop_chanedit(EEG_removedChannel, 'lookup','/Users/masoud/OneDrive - Rutgers University/document/MS journal/Data Sharing/eeglab2021.0/plugins/dipfit/standard_BEM/elec/standard_1005.elc','changefield',{1,'labels',''},'load',{'/Users/masoud/OneDrive - Rutgers University/document/MS journal/Data Sharing/eeglab2021.0/Copy_of_32channellocation.ced','filetype','autodetect'});

% frequency bandpass filter (0.4Hz- 40 Hz)
EEG_frequency_filtered = pop_eegfiltnew(EEG_withLocation, 'locutoff',lowband,'hicutoff',highband,'plotfreqz',0);
% remove faltline
EEG_cleanFlatline=clean_flatlines(EEG_frequency_filtered);
% seperating data to differents epoches and removing baseline for each epoch
EEG_epoched = pop_epoch( EEG_cleanFlatline, {  }, [-1  2], 'epochinfo', 'yes');
EEG_removedBaseline = pop_rmbase(EEG_epoched, [-1000 0] ,[]);
% run ICA
EEG_runICA = pop_runica(EEG_removedBaseline, 'icatype', 'runica', 'extended',1,'interrupt','on');
% giving Labels to the components(Brain, eye, muscle ...)
EEG_icaLabeled = pop_iclabel(EEG_runICA, 'default');
% Giving 70 % thereshold to potentioal noisy components (eye,muscle, ..
EEG_markedComponents  = pop_icflag(EEG_icaLabeled, [NaN NaN;0.7 1;0.7 1;0.7 1;0.7 1;0.7 1;0.7 1]);
% Visualize marked components
% EEG_markedComponents=pop_selectcomps(EEG_markedComponents );
% removing bad or unrelated components 
EEG_removedComponents = pop_subcomp(EEG_markedComponents);
% extract event type
EEG_event=struct2cell(EEG_removedComponents.event);
EEG_event=squeeze(EEG_event);
eventType=EEG_event(2,2);
% extract Epoch for 100 samples before event and 512 samples after that and removing baseline for each epoch
[OUTEEG, indices] = pop_epoch( EEG_removedComponents, eventType, [-100/512 1]);
EEG_removedBaseline2 = pop_rmbase( OUTEEG, [(-100/512*1000) 0] ,[]);
% Finding noisy Trials
[Iin, Iout, newsignal, elec] = eegthresh(  EEG_removedBaseline2.data, 612,24:32, -80, 80, [0 120], 0,120);
EEG_data_preprocessed=OUTEEG.data;
EEG_eventIn=Iin;
EEG_eventOut=Iout;
% Separating events for each pathway
magnoPosition=location_Magno_trigger;
konioPosition=location_Konio_trigger;
parvoPosition=location_Parvo_trigger;
% separating trigger for each visual pathway and removing noisy trials
for i=1:length(EEG_eventOut)
index=find(magnoPosition==EEG_eventOut(i));
magnoPosition(index)=[];
index=find(konioPosition==EEG_eventOut(i));
konioPosition(index)=[];
index=find(parvoPosition==EEG_eventOut(i));
parvoPosition(index)=[];
end
% extracting remain trials in selected channels
Magno=EEG_data_preprocessed(24:32,:,magnoPosition);
Konio=EEG_data_preprocessed(24:32,:,konioPosition);
Parvo=EEG_data_preprocessed(24:32,:,parvoPosition);

    
VEP{1}=Magno;
VEP{2}=Konio;
VEP{3}=Parvo;
% saving the first session of each subject
if session=='1'
subjectsTrials1=VEP;
VEP=[];
end
% saving the second session of each subject
if session=='2'
subjectsTrials2=VEP;
VEP=[];
end

end
for i=1:number_of_visualpathway
tempSubjectsTrials1=cell2mat(subjectsTrials1(i));
tempSubjectsTrials2=cell2mat(subjectsTrials2(i));
tempSubjectsTrials=cat(3,tempSubjectsTrials1,tempSubjectsTrials2);% merging data set of two sessions
subjectsTrials{i}=tempSubjectsTrials;
tempSubjectsTrials=[];
 tempSubjectsTrials2=[];
 tempSubjectsTrials1=[];
end


subjectTrialsFinal(subjectsNum,:)=subjectsTrials;%

end