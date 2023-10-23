clc; clear

addpath('C:/Users/shjo/OneDrive/Sources/Tools/HHT-Tutorial-master/HuangEMD/')

pth='C:/Users/shjo/OneDrive/mySO/EEMD_sigs/OHC700_1960_2020_200E250E_60S53S/';
myDATA=dir([pth,'*.mat']);

mkdir([pth,'EEMD_results/'])
for i=1:length(myDATA)
    disp(['!!! ', myDATA(i).name ,' !!!'])
    mySig=load([myDATA(i).folder,'\',myDATA(i).name]).sgnl;
    myDATA(i).EEMD=eemd_H(mySig,0.2,1000);
    tmp=myDATA(i).EEMD;
    save([pth,'EEMD_results/EEMD_',myDATA(i).name],'tmp')
end
