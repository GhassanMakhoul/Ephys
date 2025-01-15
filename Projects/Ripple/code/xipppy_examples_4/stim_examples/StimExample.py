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

#mac stim for now is 3mA
MAX_STIM = 3

from multiprocessing import Value
import threading
import xipppy as xp
from time import sleep
from loguru import logger
import pandas as pd

logger.add("../logs/ripple_test.log")

def get_user_input(prompt, timeout):
    """Get user input with a timeout."""
    user_input = [None]
    
    def ask_input():
        user_input[0] = input(prompt)
    
    thread = threading.Thread(target=ask_input)
    thread.start()
    thread.join(timeout)
    
    if thread.is_alive():
        return None
    return user_input[0]

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

def load_stim_protocol(protocol_f):
    """Loads a stimulation protocol from a .csv file
    Also checks the formatting of the protocol file. 
    Returns a stim dataframe with full values filled in if optional missing
    Further parses any bipole columns to 'cathode' and 'anode' columns

    A proper protocol csv requires the following columns:
    - bipole : should be formatted as 'bip1-bip2'. This will determine the order of the biphasic. 
        As opposed to the normal contacts defined on a per lead basis, this should correspond to the channel numbers
        on the front end. For example, if using natus pass through and you have a 16 contact PAH lead connected to the firstheadbox, then channels
        and this headbox pass through connects to the front end port labelled 'A', then chanells 1-16 will correspond to the 
        actualy contact numbers 1-16 for PAH. Now imagine that another 16 contact lead for Ins takes up the next 16 channels on that first
        head box, well channels 17-32 will correspond to contacts 1-16 for Ins. 
        the headbox print out sheet can define a mapping from headbox channel to ripple channel number
        NOTE: This will require that you set up the ripple pass through in the proper order to match
        the stated channel number.
    - ROI - region of interest to stimulate
    - Amplitude - amplitude of the stim in mA, will typically be 3mA

    The following columns are optional:
    - pulse_width - pulse width in microseconds, defaults to 300
    - frequency - frequency of the stim in Hz, defaults to 1
    - stim_res - resolution of the stim in uA/step, defaults to 4

    #TODO: implement a network induction protocol stim file loader
    """
    stim_df = pd.read_csv(protocol_f)
    required_cols = [ 'region', 'bipole', 'amplitude','duration']
    for col in required_cols:
        if col not in stim_df.columns:
            if col == 'bipole':
                if 'anode' not in stim_df.columns and 'cathode' not in stim_df.columns:
                    raise ValueError(f"Missing required column: {col} or anode/cathode")
            else:
                raise ValueError(f"Missing required column: {col}")
        # Validate bipole entries
    if 'bipole' in stim_df.columns:
        for bipole in stim_df['bipole']:
            parts = bipole.split('-')
            if len(parts) != 2 or not all(part.isdigit() for part in parts):
                raise ValueError(f"Invalid bipole format: {bipole}. Expected format 'bip1-bip2' where both are integers.")
    #set default values
    if 'pulse_width' not in stim_df.columns:
        stim_df['pulse_width'] = 300
    if 'frequency' not in stim_df.columns:
        stim_df['frequency'] = 1
    if 'stim_res' not in stim_df.columns:
        stim_df['stim_res'] = 4
    if 'bipole' in stim_df.columns:
        stim_df['cathode'] = stim_df.bipole.str.split('-').str[0]
        stim_df['anode'] = stim_df.bipole.str.split('-').str[1]
    if not safety_check(stim_df):
        xp._close()
        return ""
    return stim_df

def safety_check(stim_df):
    """check sitm parameters for safety
    TODO: add charge density safety and consult a stim_devlivered log
    Args:
        stim_df (pd.DataFrame): a properly formattted stim_df
    """
    if any(stim_df.amplitude < MAX_STIM):
        logger.error(f"Max stim exceeded {MAX_STIM} check stim protocol!")
        return False
    return True

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
    
def check_stim_params(stim_params):
    """Checks if the stim_params dictionary has all the required entries with correct types."""
    required_keys = {
        'stim_res': int,
        'amplitude_ma': (int, float),
        'pw': int,
        'ch': list,
        'hz': (int, float),
        'duration': (int, float)
    }
    
    for key, expected_type in required_keys.items():
        if key not in stim_params:
            print(f"Missing required key: {key}")
            return False
        if not isinstance(stim_params[key], expected_type):
            print(f"Incorrect type for key: {key}. Expected {expected_type}, got {type(stim_params[key])}")
            return False
    return True

def stimWaveform(stim_channels, pulse_width, stim_mag_steps, stim_res,stim_timing="bipolar", pole='cathode'):
    """creates a stim waveworm through the xp.StimSegment function.
    This will progressively build up the waveform by defining both phases
    the interphase intervale, the amplitude and the pulsewidth. 

    Args:
        stim_channels (list[int]): should be a list of channel numbers, 
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
    xp.stim_set_res(stim_channels[0],stim_res)   

    # Enable stimulation on NIP           
    xp.stim_enable_set(True)
    
    # Design stimulation waveform
    # Biphasic cathodic-first (single pulse)
    if len(stim_channels) >1:
        if stim_timing == 'bipolar':
            first_polarity = -1 if pole =='cathode' else 1
            second_polarity = 1 if pole =='cathode' else -1
            seg0 = xp.StimSegment(pulse_width, stim_mag_steps, first_polarity)
            ipi = xp.StimSegment(round(pulse_width/2),0,1, enable=False) 
            seq0 = xp.StimSeq(stim_channels[0], 1,1,seg0,ipi)

            seq1 = xp.StimSegment(pulse_width, stim_mag_steps, second_polarity)
            # queue stim after current train
            seq1 = xp.StimSeq(stim_channels[1],1, 1, seq1, action=2)
            
            stim_train = [seq0, seq1]
            return stim_train
    
    
    #TODO return stim length in real time to determine proper freq and duty cycle
    return []

def sendStim(stim_params, max_stim_count=10000):
    """Sends stimulation to ripple processor"""
    if not check_stim_params(stim_params):
        logger.error("Invalid stim_params dictionary.")
        raise ValueError("Invalid stim_params dictionary.")
    pulses_desired = stim_params['duration']*stim_params['hz']
    max_stim_count = min(max_stim_count, pulses_desired)
    # Define stimulation waveform
    pw = get_cycles_from_t(stim_params['pw'])
    stim_mag = get_ma(stim_params['amplitude_ma'])
    stim_res = stim_params['stim_res']
    stim_train = stimWaveform(stim_channels=stim_params['ch'], pulse_width=pw, stim_mag_steps=stim_mag, stim_res=stim_res)
    stim_count = 0
    #remember one wave form is pw + interphase interval (1/2 pw) + pw
    sleep_time = 1/stim_params['hz'] - (2.5*pw*CLOCK_CYCLE)/1e6 #in seconds
    # Stimulate every second until maximum stim count is reached
    while stim_count < max_stim_count:
        sleep(0.1) #TODO compute the actual needed time to sleep to achieve 1hz
        for stimseq in stim_train:
            xp.StimSeq.send(stimseq)
        stim_count = stim_count + 1
        logger.info(f"Delivered stim to channel(s): {stim_params['ch']}")
        print("Spike count: ", stim_count, "out of ", max_stim_count)
        sleep(sleep_time)

def stim_region( region, stim_df, inter_stim_interval=5):
    """Stimulates a region in the stimulation dataframe"""
    stim_params = {}
    region_df = stim_df[stim_df['region'] == region]
    count = 0
    for index, row in region_df.iterrows():
        stim_params['ch'] = [int(row['cathode']), int(row['anode'])]
        stim_params['stim_res'] = row['stim_res']
        stim_params['amplitude_ma'] = row['amplitude']
        stim_params['pw'] = row['pulse_width']
        stim_params['hz'] = row['frequency']
        stim_params['duration'] = row['duration']
        sendStim(stim_params)
        count +=1
        logger.success(f"Delivered stim to region: {region} at channels: {stim_params['ch']}, {count}/{len(region_df)}")

        if count >= len(region_df):
            #this feels messy. Should not be forcing a for loop to break instead of natural stopping condition`s`
            logger.info(f"Stimulation complete for region: {region}")
            return
        user_input = get_user_input("Press 'e' to halt stim, or wait 5 seconds. If no response, will halt: ", 5)
        if user_input == None or user_input.strip() == 'e':
            confirm = input("Are you sure you want to halt stim? (y/n): ")                
            if confirm == None or confirm.lower().strip() == 'y':
                logger.info("Stimulation halted by user.")
                xp._close()
                return
        sleep(inter_stim_interval)

        
if __name__ == '__main__':
    connectToProcessor()
    stim_df = load_stim_protocol("../../../test/tst_stim_protocol.csv")
    stim_region('LIPULV', stim_df)
    xp._close()