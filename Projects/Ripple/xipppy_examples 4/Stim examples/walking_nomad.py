#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 16:10:57 2022

Walking Nomad demo

OVERVIEW

This will take input from the first two channels on an EMG front end,
then stimulate on any stimulating front end based on EMG amplitude.

The EMG and stim FE can be attached to any port. An error will trip if
an EMG and stim FE are not connected.

The system was designed to be used with an LED board that has rows of 8 LEDs.

Maximum EMG magnitude illumates all 8 LEDs. Minimum illuminates 1. 


CALIBRATION

The device is auto-calibrated. It starts with a maximum of 0.

As the muscle activity is recorded by the EMG front end, the maximum value
is reset to the most recent maximum. So, to illumate all 8 LEDs, you must 
provide 7/8 of the most recent maximum reading. 

EXAMPLE USAGE

Connect the two electrodes from channel 1 to either side of a muscle.

Connect EMG to electrode leads, a Nano+stim to the LED board, and to a 
processor.

Start the code. 

All LEDs will initially be on. Flex the muscle and a new maximum value will
be set. You can then flex or extend with varying LED output based on effort.

@author: kyleloizos
"""

import xipppy as xp
from time import sleep

if __name__ == '__main__':
    # Connect to the processor (using TCP protocol)
    
    with xp.xipppy_open(use_tcp=True):
        # Set parameters for calibration
        max_emg = 0
        max_emg2 = 0
        
        # Verify stim FE is connected
        stim_channels = xp.list_elec('stim')
        if not stim_channels:
            print("ERROR - Stim FE not detected")

        # Define recroding channels (EMG front end)            
        rec_channels = xp.list_elec('EMG')
        if not rec_channels:
            print("ERROR - EMG FE not detected")
        
        # Set recording electrodes
        rec_elec1 = rec_channels[0]
        rec_elec2 = rec_channels[1]
        
        # Set filters
        xp.filter_set(rec_elec1, 'hires', 4)
        
        # Set stimulation resolution and enable
        stim_res = 4
        seq = []
        xp.stim_enable_set(False)
        sleep(0.01)
        
        for channel in stim_channels:
            xp.stim_set_res(channel, stim_res)
        
        xp.stim_enable_set(True)
        
        elec = stim_channels[0]
        
        # Design stimulation waveform
        pseg = xp.StimSegment(6,90,-1)
        ipi = xp.StimSegment(1,0,1, enable=False)
        nseg = xp.StimSegment(6,90,1)

        # Measure EMG (hi-res stream)
        while 1:
            # Reset 
            seq = []
            
            # Measure EMG
            emg_data, timestamp = xp.cont_hires(1, [rec_elec1])
            emg_data2, timestamp2 = xp.cont_hires(1, [rec_elec2])
            sleep(0.0001)
            emg_datapoint = emg_data[0]
            emg_datapoint2 = emg_data2[0]
            sleep(0.0001)
            
            # Set new maximum
            emg_datapoint = abs(emg_datapoint)
            emg_datapoint2 = abs(emg_datapoint2)
            
            if emg_datapoint > max_emg:
                max_emg = emg_datapoint
            
            if emg_datapoint2 > max_emg2:
                max_emg2 = emg_datapoint2
            
            # Light first LED in each row 
            for led_row in range(0,4):
                seq0 = xp.StimSeq(led_row*8,50,500,pseg,ipi,nseg)
                seq.append(seq0)
            
            # Light LEDs based on EMG input
            for channel in range(1,8):
                #print (emg_datapoint)
                #print ("max", channel*max_emg/8)
                
                if emg_datapoint > (channel*max_emg/8):
                    seq0 = xp.StimSeq(channel,50,500,pseg,ipi,nseg)
                    seq1 = xp.StimSeq(channel+8,50,500,pseg,ipi,nseg)
                    seq.append(seq0)
                    seq.append(seq1)
                if emg_datapoint2 > (channel*max_emg2/8):
                    seq0 = xp.StimSeq(channel+16,50,500,pseg,ipi,nseg)
                    seq1 = xp.StimSeq(channel+24,50,500,pseg,ipi,nseg)
                    seq.append(seq0)
                    seq.append(seq1)
                
            xp.StimSeq.send_stim_seqs(seq)
