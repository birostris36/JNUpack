clc; clear
% addpath
addpath('C:/Users/shjo/OneDrive/base142/EOF/pcatool')

pth='C:/Users/shjo/OneDrive/mySO/REOF_OHC/OHC700_1993_2020_0E360E_75S30S/';

N=10; % Total EOF number

myDATA=dir([pth,'*.mat']);

mkdir([pth,'REOF_results/'])
for i=1:length(myDATA)
    disp(['!!! ', myDATA(i).name ,' !!!'])
    mySig=load([myDATA(i).folder,'\',myDATA(i).name]).sgnl;
    mask=myMasking(mySig);
    G = map2mat(mask(:,:),mySig); 
    [pcs,e,expvar] = caleof(G',N,2); % Perform rotate EOF G->G'
    eofs = mat2map(mask(:,:),e) ;
    save([pth,'REOF_results/REOF_',myDATA(i).name],'eofs','pcs','expvar')
end

function [mask]=myMasking(myDATA)
    
    [~,AT,ON] = size(myDATA) ;
    var_mask = myDATA ;
    mask = zeros([AT,ON]) ; 
     
     for i2 = 1 : AT
         for j2 = 1 : ON
    
             if isnan( mean( var_mask(:,i2,j2) ))
                mask(i2,j2) = 0 ;
             else 
                 mask(i2,j2) = 1 ;
             end
         end
     end
     
end