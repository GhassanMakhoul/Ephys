"""
This script will take TTL input and trigger a stimulation.
It will also note the time delay between the sensed TTL pulse
and the stimulation.

The stimulation is a charge-balanced symmetric biphasic pulse.

The script will run until stopped by the user.

Created on Mon Feb  7 17:21:58 2022

@author: kloizos
"""

# Import xipppy API to use in python (with shortened name 'xp' for ease of use)
import xipppy as xp
from time import sleep

if __name__ == '__main__':
    # Connect to the processor (using TCP protocol)
    # For UDP protocol, use "with xp.xipppy_open()"
    with xp.xipppy_open(use_tcp=True):
                
        # Set stimulation parameters (electrode 1 on port A, 5uA/step)
        elec = 0
        res = 3
        xp.stim_enable_set(False)
        sleep(0.001)
        xp.stim_set_res(elec,res)   

        # Enable stimulation on NIP           
        xp.stim_enable_set(True)
        
        # Design stimulation waveform
        # Biphasic cathodic-first (single pulse)
        # amplitude=50 steps, pulse-width=200us (6 clock cycles), interphase interval=100us (3 clock cycles)
        # sampling freq = 30kHz, clock cycle 33.33us
        pseg = xp.StimSegment(6,50,-1)
        ipi = xp.StimSegment(3,0,1, enable=False) 
        nseg = xp.StimSegment(6,50,1) 
        seq0 = xp.StimSeq(elec, 1000,1,pseg,ipi,nseg) 
        
        # Scan input for TTL (in an infinite loop)
        while True:
            # Stim at 1Hz #########
            xp.StimSeq.send(seq0)
            sleep(1)
            #######################
            
            # check digital input for an event and note timestamp
            sleep(0.001)
            [n,event] = xp.digin(max_events=1)
            
            # if at least one event is found, trigger stimulation
            if (n > 0):
                # Send sequence to NIP
                xp.StimSeq.send(seq0)
                