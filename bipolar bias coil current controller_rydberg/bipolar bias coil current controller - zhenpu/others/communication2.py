# -*- coding: utf-8 -*-
"""
Created on Sat Sep 30 15:11:18 2017

@author: wjf1998
"""

import serial
import time 

port='COM10'

#send a frequency in multiples of 1000Hz(i.e. 1kHz).This is so that we need to send only
#2 bytes and gives a maximum frequency of 65 MHz. If multipiled by 32, it gives 2 GHz. 
ser=serial.Serial(port,9600)
time.sleep(1)
 
# Unit: kHz
#processed_list=[1679697]      
processed_list=[160500]      
send_length=len(processed_list)
ser.write(chr(send_length))

for i in range(send_length):
    ser.write(chr((int(processed_list[i])     & 0b111111110000000000000000)>>16))  #高8位 Unit:256^2
    ser.write(chr((int(processed_list[i])<<8  & 0b111111110000000000000000)>>16))  #中8位 Unit:256
    ser.write(chr((int(processed_list[i])<<16 & 0b111111110000000000000000)>>16))  #低8位 Unit:1
    #print teh frequncy in MHz
    print "The output frequnecy is: " + str(float(processed_list[i])/1000) + " MHz."  

#for i in range(send_length):
#    print(ser.readline()
    
    
ser.close()

