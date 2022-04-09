#!/usr/bin/env python3
# -*- coding: utf-8 -*-2
"""
Created on Thu Apr  7 17:04:24 2022

@author: elise
"""
import csv
import pandas as pd
import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt

#reads in .tsv and conversts to .csv in an array
tsv_file='inner.tsv'
csv_table=pd.read_table(tsv_file,sep='\t')
csv_table.to_csv("inner2.csv",index=False)

tsv_file2='outer.tsv'
csv_table2=pd.read_table(tsv_file2,sep='\t')
csv_table2.to_csv("outer2.csv",index=False)

tsv_file3='dnf.tsv'
csv_table3=pd.read_table(tsv_file3,sep='\t')
csv_table3.to_csv("dnf2.csv",index=False)

#--!!check passingvariables!

array_rawsig = np.loadtxt("inner2.csv", delimiter = ",", usecols = (0)) #loads in data file to float array
print(array_rawsig)

array_noisesig = np.loadtxt("outer2.csv", delimiter = ",", usecols = (0)) #loads in data file to float array
print(array_noisesig)

array_dnfsig = np.loadtxt("dnf2.csv", delimiter = ",", usecols = (0)) #loads in data file to float array
print(array_dnfsig)


#calculating powers of raw sig array
n =len(array_rawsig)
square = 0
Araw_sig = 0
mean = 0

for i in range(0,n):
    square += (array_rawsig[i]**2)
    
mean = (square/(float)(n))
Araw_sig = np.sqrt(mean)
print(Araw_sig)    

#calculating powers of noise sig array
n =len(array_rawsig) #same sample length for all
square2 = 0
Anoise_sig = 0
mean2 = 0

for i in range(0,n):
    square2 += (array_noisesig[i]**2)
    
mean2 = (square2/(float)(n))
Anoise_sig = np.sqrt(mean2)
print(Anoise_sig)    

#calculating powers of dnf sig array
n =len(array_rawsig) #same sample length for all
square3 = 0
Adnf_sig = 0
mean3 = 0

for i in range(0,n):
    square3 += (array_dnfsig[i]**2)
    
mean3 = (square3/(float)(n))
Adnf_sig = np.sqrt(mean3)
print(Adnf_sig)    


#SNR of raw sig and noise 
SNR= np.log10(np.square(Araw_sig/Anoise_sig))*10#snr calculation
print(SNR, "dB")

#SNR of dnf sig and noise 
SNR2= np.log10(np.square(Adnf_sig/Anoise_sig))*10#snr calculation
print(SNR2, "dB")


#array to test plotting 
n =len(array_rawsig)
plot = np.abs(array_dnfsig)**2
meanie = plot/(float)(n)
squaret = np.sqrt(meanie)
print(squaret)

snrr= np.log10(np.square(squaret/3))*10#snr calculation, using 3 where noise will be for plotting purposes
print (snrr)

x = snrr
y = np.sort(x)

plt.figure()
plt.plot(x)


#--wav conversion (not working yet)
amplitude = 5 #any val to test

outsig = amplitude*array_rawsig #amplifies signal

sf.write("rawsig.wav",outsig,44100) #converts to .wav (does this too quickly)
#-------------














