clear; clc

addpath('D:/JNUpack/SOHEAT/matlab/Toolbox/ocean/')
% Filedir = 'D:\7.자료및관측\EN4\Data\';
Filedir = 'G:/EN4/';

% Syear = 1960;
Syear = 1993;
Eyear = 2020;
Smon = 1;
Emon = 12;
Mrho = 28.106331; % sigma(0,0,35) = 28.106331 kg/m^3

Filename = ([Filedir,'EN.4.2.1.f.analysis.g10.',num2str(Syear),num2str(Smon,'%2.2i'),'.nc']);
% Filename = ([Filedir,'EN.4.2.1.f.analysis.l09.',num2str(Syear),num2str(Smon,'%2.2i'),'.nc']);
% Filename = ([Filedir,'EN.4.2.2.f.analysis.l09.',num2str(Syear),num2str(Smon,'%2.2i'),'.nc']);
nc = netcdf(Filename);
depth = nc{'depth'}(:);
db = nc{'depth_bnds'}(:);
lat = nc{'lat'}(:);
lon = nc{'lon'}(:);
clear nc

[lonr,latr] = meshgrid(lon,lat);
[M,L] = size(lonr);

cp = 4000;

ohc = zeros((Eyear-Syear+1).*12,M,L);
nm = 0;
for Year = Syear : 2017
    disp(['========= Processing ',num2str(Year),' Year ========='])
    %% 
    
    for Mon = 1 : 12
        
        disp(['========= Processing ',num2str(Mon),' Mon ========='])
        nm = nm + 1;
        
        if Year==2021 && Mon ==7||Year==2021 && Mon ==8||Year==2021 && Mon ==9 || ...
            Year==2021 && Mon ==10||Year==2021 && Mon ==11||Year==2021 && Mon ==12
            
            disp('No data month')
        else
            
            Filename = ([Filedir,'EN.4.2.1.f.analysis.g10.',num2str(Year),num2str(Mon,'%2.2i'),'.nc']);
%             Filename = ([Filedir,'EN.4.2.1.f.analysis.l09.',num2str(Year),num2str(Mon,'%2.2i'),'.nc']);
%             Filename = ([Filedir,'EN.4.2.2.f.analysis.l09.',num2str(Year),num2str(Mon,'%2.2i'),'.nc']);
            
            nc = netcdf(Filename);
            temp = nc{'temperature'}(:);
            salt = nc{'salinity'}(:);
            clear nc
            
            temp(temp == -32768) = nan;
            salt(salt == -32768) = nan;
            
            % ==============================================================
            % Calculate Density (rho)
            % ==============================================================
            
            dens = zeros(size(temp));
%             cp = zeros(size(temp));
            zet = zeros(M,L);
            
            for iy = 1 : M
                for ix = 1 : L
                    zdz = 0;
                    zet(iy,ix) = 0;
                    for iz = 1 : length(depth)
                        
%                         P = zts2p(depth(iz,1),temp(iz,iy,ix),salt(iz,iy,ix),latr(iy,ix));
                        P = depth(iz,1);
                        dens(iz,iy,ix) = sigma(P,temp(iz,iy,ix),salt(iz,iy,ix))+1000;
%                         cp(iz,iy,ix) = cpsw(P,temp(iz,iy,ix),salt(iz,iy,ix));
                        
                        dz = db(iz,2) - db(iz,1);
                        
                        if zdz <= 2000
%                         if zdz <= 700    
%                         if zdz <= 300
                            
                            zet(iy,ix) = zet(iy,ix) + cp*temp(iz,iy,ix)*dens(iz,iy,ix)*dz; % low dens --> sea level higher
                            zdz = zdz + dz;
                            
                        else
%                             disp(zdz)
                        end
                        
                    end
                    
%                     clear P zdz
                end
            end
            disp(nm)
            ohc(nm,:,:) = zet;
        end
        
    end
end

tmp=squeeze(ohc(1,:,:));

% save EN4_4_2_1_OHC_l09_2000m_93_20_20211116.mat lon lat depth db ohc
% save EN4_4_2_1_OHC_l09_2000m_93_21(06)_20220119.mat lon lat depth db ohc
% save EN4_4_2_1_OHC_l09_2000m_60_21(06)_20220119.mat lon lat depth db ohc
save EN4_4_2_2_OHC_l09_2000m_93_20_20220519.mat lon lat depth db ohc
% save EN4_4_2_1_OHC_l09_700m_60_20_20210818.mat lon lat depth db ohc
% save EN4_4_2_1_OHC_l09_700m_93_20_20210816.mat lon lat depth db ohc
% save EN4_4_2_1_OHC_300m_60_17_20201005.mat lon lat depth db ohc
