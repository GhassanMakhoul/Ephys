"""
This is an example script for running closed-loop stimulation.

It will send a stimulation to a specified channel whenever a spike is detected on a specified channel.

As an example, the stimulation waveform is defined as a biphasic pulse, with properties that can be set by the user.
However, the waveform can be modified as the user sees fit.

Created on Fri Sep 30 10:20:27 2022

@author: kyleloizos
"""

import xipppy as xp
from time import sleep


def connect_to_processor():
    """
    Connects to the Neural Interface Processor. If UDP connection fails, 
    attempts TCP connection.
    """
    xp._close()
    sleep(0.001)

    try:
        xp._open()
        print("Connected over UDP")
    except:
        try:
            xp._open(use_tcp=True)
            print("Connected over TCP")
        except:
            print("Failed to connect to processor.")

    sleep(0.001)


def stim_waveform(stim_channel, pulse_width, stim_mag_steps, stim_res):
    """
    Defines stimulation waveform as a biphasic pulse with user-defined properties.

    Parameters:
        stim_channel: Channel to send stimulus to
        pulse_width: Width of each phase (units of clock cycles)
        stim_mag_steps: Magnitude of stim pulse (units of steps based on resolution)
        stim_res: Index of resolution of stim desired 
                  (e.g., for nano, 1=1uA/step, 2=2uA/step, 3=5uA/step, 4=10uA/step, 5=20uA/step)

    Returns:
        stim_seq: xp.StimSeq - Stimulation sequence object with defined waveform
    """
    xp.stim_enable_set(False)
    sleep(0.001)
    xp.stim_set_res(stim_channel, stim_res)

    xp.stim_enable_set(True)

    cathodic_phase = xp.StimSegment(pulse_width, stim_mag_steps, -1)
    interphase_interval = xp.StimSegment(round(pulse_width / 2), 0, 1, enable=False)
    anodic_phase = xp.StimSegment(pulse_width, stim_mag_steps, 1)
    stim_seq = xp.StimSeq(stim_channel, 1000, 1, cathodic_phase, interphase_interval, anodic_phase)

    return stim_seq


def stim_when_spike_detected(max_stim_count):
    """
    Triggers stimulation whenever a spike is detected until the maximum stimulation count is reached.

    Parameters:
        max_stim_count: Maximum number of stimulations to perform
    """
    # Stimulation channel and parameters may be adjusted
    stim_waveform_seq = stim_waveform(stim_channel=0, pulse_width=200, stim_mag_steps=50, stim_res=3)
    stim_count = 0

    while stim_count < max_stim_count:
        sleep(0.1)

        # Collect spike data, spk_data parameter indicates channel for recording spikes
        spike_count, spike_data = xp.spk_data(1)

        # If a spike was detected, send stimulation
        if spike_count:
            xp.StimSeq.send(stim_waveform_seq)
            stim_count += 1
            print(f"Spike count: {stim_count} out of {max_stim_count}")
            sleep(1)


if __name__ == '__main__':
    connect_to_processor()
    stim_when_spike_detected(200)
    xp._close()