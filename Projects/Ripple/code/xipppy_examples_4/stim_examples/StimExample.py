#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is an example script for sending a train of stimulation pulses.

As an example, the stimulation waveform is defined as a biphasic pulse, 
with properties that can be set by the user. However, the waveform can
be modified as the user sees fit.

Created on Wed Nov 11 2022

@author: kyleloizos, GhassanMakhoul (extended)
"""

CLOCK_CYCLE = 33.33 #us (microseconds)

import xipppy as xp
from time import sleep
from loguru import logger 

logger.add("../logs/ripple_test.log")

def connectToProcessor():
    # Close connection
    xp._close()
    sleep(0.001)
    
    # Connect to processor (Try UDP then TCP). 
    try:
        xp._open()
    except:
        try:
            xp._open(use_tcp=True)
        except:
            print("Failed to connect to processor.")
        else:
            print("Connected over TCP")
    else:
        print("Connected over UDP")
    sleep(0.001)
    

def get_cycles_from_t(t, tol=1e-3):
    """Given a time length, returns how many clock cycles it corresponds with
    If there is a larger difference between expected time and clock_cycles capability
    will raise a warning

    Args:
        t (_type_): length in micro-seconds of desired time
    
    Returns:
        c_t (int) : clock cycles corresponding to amount of time
    """
    c_t = int(t/CLOCK_CYCLE)
    if c_t*CLOCK_CYCLE < t - tol:
        logger.warning(f"Desired time {t} but actual {c_t*CLOCK_CYCLE}")
    return c_t

def get_ma(ma_desired, stim_res=4):
    """Given a target amplitdue, returns the proper amount of sim_mag_steps
    For example, imagine you want to program 3mA of stim:
        - 3mA is equivalent to 3000 ua (microamps)
        - If using stim_res 4 then each step is 10uA
        - thus 3(mA) * 1000/1(uA/mA) *1/10(steps/uA) = 300 steps
    """  

    if stim_res == 1:
        #1uA/Step
        return ma_desired*1000
    if stim_res == 2:
        #2uA/step
        return ma_desired*500
    if stim_res == 3:
        # 5uA/step
        return ma_desired*200
    if stim_res == 4:
        # 10uA/step
        return ma_desired*100
    if stim_res == 5:
        #20uA/step
        return ma_desired*50
    

def stimWaveform(stim_channel, pulse_width, stim_mag_steps, stim_res,stim_timing="bipolar", pole='cathode'):
    """creates a stim waveworm through the xp.StimSegment function.
    This will progressively build up the waveform by defining both phases
    the interphase intervale, the amplitude and the pulsewidth. 

    Args:
        stim_channel (list[int]): should be a list of channel numbers, 
        pulse_width (int): pulse width in clock cycles
        stim_mag_steps (int): magnitude steps of stim
        stim_res (int): index of resolution of stim desired     
            #    (e.g.for nano, 1=1uA/step, 2=2uA/step, 3=5uA/step, 4=10uA/step, 5=20uA/step)
        stim_timing (defult: 'bipolar'): setting this value will determine if stim_channel list is to be used as a 
        bipolar stim specifier (cathode, anode) or if 'simultaneous' will stimulate all specified electrodes as per 'pole'
        pole (str, defulat: 'cathode'): specifies how to construct first in stim_channel (cathode, anode) vs (anode, cathode). If 'stim_timing is
        selected to be simulataneous then all channels will be simultaneously stimulated with setting of 'pole'
    Returns:
        xp.StimSeq: stimulation sequence
    """


    xp.stim_enable_set(False)
    sleep(0.001)
    xp.stim_set_res(stim_channel,stim_res)   

    # Enable stimulation on NIP           
    xp.stim_enable_set(True)
    
    # Design stimulation waveform
    # Biphasic cathodic-first (single pulse)
    pseg = xp.StimSegment(pulse_width,stim_mag_steps,-1)
    ipi = xp.StimSegment(round(pulse_width/2),0,1, enable=False) 
    nseg = xp.StimSegment(pulse_width,stim_mag_steps,1) 
    seq0 = xp.StimSeq(stim_channel, 1000,1,pseg,ipi,nseg)
    #TODO return stim length in real time to determine proper freq and duty cycle
    return seq0

def sendStim(stim_params, max_stim_count):
    
    # Define stimulation waveform
    pw = get_cycles_from_t(stim_params['pw'])
    stim_mag = get_ma(stim_params['amplitude_ma'])
    stim_res = stim_params['stim_res']
    stim_waveform = stimWaveform(stim_channel=stim_params['ch'],pulse_width=pw,stim_mag_steps=stim_mag,stim_res=4)
    stim_count = 0
    
    # Stimulate every second until maximum stim count is reached
    while stim_count < max_stim_count:
        sleep(0.1) #TODO compute the actual needed time to sleep to achieve 1hz
        xp.StimSeq.send(stim_waveform)
        stim_count = stim_count + 1
        print("Spike count: ", stim_count, "out of ", max_stim_count)
        sleep(1);
    
if __name__ == '__main__':
    connectToProcessor()
    stim_params = {}
    stim_params['stim_res'] = 4
    stim_params['amplitude_ma'] = 1
    stim_params['pw'] = 300
    stim_params['ch'] = 3
    sendStim(stim_params, 10)
    xp._close()