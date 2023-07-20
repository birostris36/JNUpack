# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 14:01:33 2023

@author: shjo9
"""

import xarray as xr
import os
import numpy as np


OHC=xr.open_dataset('D:/HEAT/EN4_OHC_GLOBAL_1960_2023.nc',decode_times=True)\
    .loc[dict(time=slice('1980-01','2022-12'))]

### Remove Trend ==============================================================
data=OHC.OHC.assign_coords({'TT':('time',range(len(OHC.time)))})
data=data.swap_dims({"time":"TT"})
data_s=data.polyfit(dim='TT',deg=1,skipna=True)
fit = xr.polyval(data.TT, data_s.polyfit_coefficients)
detrend_OHC=data-fit
detrend_OHC=detrend_OHC.swap_dims({"TT":"time"})
### ===========================================================================
DATA=OHC.OHC
DATA=detrend_OHC

PAC=DATA.loc[dict(lon=slice(185,240),lat=slice(-50,-20))].mean(dim=['lat','lon'])
IND=DATA.loc[dict(lon=slice(70,110),lat=slice(-50,-20))].mean(dim=['lat','lon'])

PAC_2Y=PAC.rolling(time=24,center=True).mean(dropna=True)
IND_2Y=IND.rolling(time=24,center=True).mean(dropna=True)

PAC.plot(figsize=(16,4))
IND.plot()
IND_2Y.plot()
PAC_2Y.plot()
Corr_PAC_IND_2Y=np.corrcoef(PAC_2Y[12:-12],IND_2Y[12:-12])[1][0]

# --> Strong seasonality + realtively week Decadal variations
# --> Replace running average to LPF

### FFT =======================================================================
def JNUFFT(sig):
    from scipy import cos, sin
    from scipy.fftpack import fft, fftfreq, ifft
    sig = sig - np.nanmean(sig)
    sig = sig.reshape(-1)
    x= np.arange(len(sig))
    n = len(sig)
    freqs = fftfreq(n)    # 필요한 모든 진동수를 만든다.
    mask = freqs > 0    # 절반의 값을 무시
    nwaves = freqs*n    # 도메인 길이에 따른 파수
    fft_vals = fft(sig)    # FFT 계산
    fft_norm = fft_vals*(1.0/n)    # FFT 계산된 결과를 정규화
    fft_theo = 2.0*abs(fft_norm)    # 푸리에 계수 계산
    # 계산하고싶은 파수의 범위를 지정 (0~50 사이의 숫자를 입력)
    wavenumber = 50  #int(input("input wavenumber (~50) : ",))
    omg = 1 # Anglural vcelocity
    x0  = np.ones(n)
    origin = fft_norm.real[0]*x0    # 상수부분인 푸리에 계수를 a0 더함
    for k in range(1, wavenumber+1):    # 푸리에계수 an, bn을 이용해 IFFT 구현
        origin +=   2 * fft_norm.real[k] * cos(k*omg*x) + \
                  (-2)* fft_norm.imag[k] * sin(k*omg*x)
    return freqs, fft_theo, mask, origin

freqs, fft_theo, mask, origin = JNUFFT()









