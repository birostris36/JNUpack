clc; clear

addpath('C:/Users/shjo/OneDrive/Sources/Tools/HHT-Tutorial-master/HuangEMD/')
pth='C:/Users/shjo/OneDrive/mySO/SCP_EEMD/OHC700_1980_2020_230E250E_60S50S/';
myDATA=dir([pth,'*.mat']);

for i=1:length(myDATA)
    disp(['!!! ', myDATA(i).name ,' !!!'])
    mySig=load([myDATA(i).folder,'\',myDATA(i).name]).sgnl;
    myDATA(i).EEMD=eemd_H(mySig,0.2,1000);
    tmp=myDATA(i).EEMD;
    save(['C:/Users/shjo/OneDrive/mySO/SCP_EEMD/OHC700_1980_2020_230E250E_60S50S/EEMD_results/EEMD_',myDATA(i).name],'tmp')
end





















