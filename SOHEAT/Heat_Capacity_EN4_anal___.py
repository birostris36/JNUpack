# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 17:27:19 2023

@author: shjo9
"""
# %matplotlib auto
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numpy import arange, ones, pi
from scipy import cos, sin
from scipy.fftpack import fft, fftfreq, ifft

# r_path = 'D:/Working_hub/OneDrive/Projects/Kuroshio_paper/EEMD/EMD_Sigs/'

sig = IND[12:-12].values
# sig = PAC[12:-12].values

# sig = sig[:

def JNUFFT(sig):
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
    x0  = ones(n)
    origin = fft_norm.real[0]*x0    # 상수부분인 푸리에 계수를 a0 더함
    for k in range(1, wavenumber+1):    # 푸리에계수 an, bn을 이용해 IFFT 구현
        origin +=   2 * fft_norm.real[k] * cos(k*omg*x) + \
                  (-2)* fft_norm.imag[k] * sin(k*omg*x)
    return freqs, fft_theo, mask, origin
# 1.   
plt.figure(figsize=(16,7))
plt.subplot(211)
plt.plot(sig, color='k', label='Original',zorder=1)
plt.plot(x, origin, color='g', label='IFFT',zorder=0)
plt.title("Original & IFFT Signal")
plt.legend()
# xtick_location = date.tolist()[::12*2]
# xtick_labels = date.tolist()[::12*2]
# plt.xticks(ticks=xtick_location, labels=xtick_labels, rotation=0, fontsize=12, alpha=.7)
# plt.title("Peak and Troughs of Air Passengers Traffic (1949 - 1969)", fontsize=22)
plt.yticks(fontsize=12, alpha=.7)
# Lighten borders
plt.gca().spines["top"].set_alpha(.0)
plt.gca().spines["bottom"].set_alpha(.3)
plt.gca().spines["right"].set_alpha(.0)
plt.gca().spines["left"].set_alpha(.3)
plt.legend(loc='upper right')
plt.grid(axis='y', alpha=.3)

# 2. 푸리에 계수
plt.subplot(212)
plt.bar(freqs[mask]*n, fft_theo[mask], label="true fft values",color='k')
plt.title("True FFT values")
# plt.xlim([0,80])
plt.axhspan(-0, .02, facecolor='0.5', alpha=0.2)
plt.axhline(y=0, linewidth=1, color='r')
plt.legend()
plt.show()



import math

Periods= len(sig)/(freqs[mask]*n)/12
Power = fft_theo[mask]

Periods2=[math.log(i,12) for i in Periods]

plt.figure(figsize=(16,7))
plt.subplot(211)
plt.plot( sig, color='k', label='Original')
# plt.plot(x, origin, color='g', label='IFFT')
plt.title("Original & IFFT Signal")
plt.legend()
# xtick_location = date.tolist()[::12*2]
# xtick_labels = date.tolist()[::12*2]
# plt.xticks(ticks=xtick_location, labels=xtick_labels, rotation=0, fontsize=12, alpha=.7)
# plt.title("Peak and Troughs of Air Passengers Traffic (1949 - 1969)", fontsize=22)
plt.yticks(fontsize=12, alpha=.7)
# Lighten borders
plt.gca().spines["top"].set_alpha(.0)
plt.gca().spines["bottom"].set_alpha(.3)
plt.gca().spines["right"].set_alpha(.0)
plt.gca().spines["left"].set_alpha(.3)
plt.legend(loc='upper right')
plt.grid(axis='y', alpha=.3)

plt.subplot(212)
# plt.vlines(Periods[2:],[0], Power[2:], label="true fft values",color='k')
plt.plot(Periods[2:], Power[2:], label="true fft values",color='k')

plt.title("True FFT values")
# plt.xlim([0,80])
plt.legend()
plt.show()


# =============================================================================
# 
# ============================================================================