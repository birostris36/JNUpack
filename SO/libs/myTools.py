import time
import os
import sys

def myInfo(loc,wpth):
    f = open(wpth+"/myInfo.txt", 'w')
    info='Ran time:\t'+time.strftime('%c')+'\nScript loc:\t'+os.getcwd()+\
        '\nScript name:\t'+loc+'\nwpth:\t\t'+ wpth
    f.write(info)
    f.close()

    
    
    
    
    
    
    