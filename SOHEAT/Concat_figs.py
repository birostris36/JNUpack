# -*- coding: utf-8 -*-
"""
Created on Tue May 30 16:32:40 2023

@author: shjo9
"""

# -*- codilist(ng: utf-8 -*-
"""
Created on Tue Feb 14 20:52:54 2023

@author: birostris
@email : birostris36@gmail.com

Name : 
Reference :
Description :
"""

import os
import cv2
import pandas as pd
# from tqdm import tqdm
import numpy as np

save_path='E:/_tmp/RMnT/Concat_GECOHC_ADT/'

image_path01='E:/_tmp/RMnT/GECOHC/'
image_path02='E:/_tmp/RMnT/SLA_dt/'
# image_path03='D:/HEAT/figs/Pair/ppt/'

IMAGES01=list(np.sort([image_path01+i for i in os.listdir(image_path01) if i.endswith('.png')]))
IMAGES02=list(np.sort([image_path02+i for i in os.listdir(image_path02) if i.endswith('.png')]))
# IMAGES03=list(np.sort([image_path03+i for i in os.listdir(image_path03) if i.endswith('.png')]))

# IMAGES01=IMAGES01[1:]
# DATES=pd.date_range('1980-02','2015-12',freq='m').strftime('%Y_%m')
DATES=range(1,len(IMAGES02)+1)

for i,j,d in zip(IMAGES01,IMAGES02,DATES):
    # Open the first image
    img1 = cv2.imread(i)
    # Open the second image
    img2 = cv2.imread(j)
    # img3 = cv2.imread(k)
    # Concatenate the images horizontally
    result = cv2.hconcat([img1,img2])
    #concatenated_image = np.vstack((img1, img2))
    # Display the concatenated image
    cv2.imwrite(save_path+f'{d:02d}'+'.png', result)

















