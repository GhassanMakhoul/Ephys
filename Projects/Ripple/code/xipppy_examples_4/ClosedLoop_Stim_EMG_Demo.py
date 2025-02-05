"""
Closed Loop Stimulation EMG demo

OVERVIEW:

This script captures input from the first two channels on an EMG front
end and stimulates on any stimulating front end based on EMG amplitude.

The EMG and stim FE can be attached to any port. An error will trip if
an EMG and stim FE are not connected.

The system was designed to be used with an LED board that has rows of
8 LEDs. Maximum EMG magnitude illuminates all 8 LEDs. Minimum
illuminates 1.

CALIBRATION:

The device is auto-calibrated. It starts with a maximum of 0.

As muscle activity is recorded by the EMG front end, the maximum value
is reset to the most recent maximum. To illuminate all 8 LEDs, you must
provide 7/8 of the most recent maximum reading.

EXAMPLE USAGE:

Connect the two electrodes from channel 1 to either side of a muscle.

Connect EMG to electrode leads, a Nano+stim to the LED board, and to a
processor.

Start the code. All LEDs will initially be on. Flex the muscle and a new
maximum value will be set. You can then flex or extend with varying LED
output based on effort.

Created on Mon May 23 16:10:57 2022

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


def verify_connections():
    """
    Verifies that both EMG and stimulation front ends are connected.
    """
    stim_channels = xp.list_elec('stim')
    if not stim_channels:
        raise RuntimeError("ERROR - Stim FE not detected")

    rec_channels = xp.list_elec('EMG')
    if not rec_channels:
        raise RuntimeError("ERROR - EMG FE not detected")

    return stim_channels, rec_channels


def setup_stim(stim_channels, stim_res=4):
    """
    Sets up stimulation with specified resolution on each channel.

    Parameters:
        stim_channels: Stim channels available
        stim_res: Resolution for stimulation (default is 4)
    """
    xp.stim_enable_set(False)
    sleep(0.01)

    for channel in stim_channels:
        xp.stim_set_res(channel, stim_res)

    xp.stim_enable_set(True)


def create_stim_waveform():
    """
    Creates and returns a biphasic stimulation waveform.

    Returns:
        Positive segment, inter-pulse interval, negative segment of the waveform
    """
    cathodic_phase = xp.StimSegment(6, 90, -1)
    interphase_interval = xp.StimSegment(1, 0, 1, enable=False)
    anodic_phase = xp.StimSegment(6, 90, 1)
    return cathodic_phase, interphase_interval, anodic_phase


def measure_emg(channel):
    """
    Measures the EMG data from the specified channel.

    Parameters:
        channel: Channel to measure EMG data from

    Returns:
        emg_datapoint: Absolute value of the EMG data point
    """
    emg_data, _ = xp.cont_hires(1, [channel])
    sleep(0.0001)
    emg_datapoint = abs(emg_data[0])
    sleep(0.0001)
    return emg_datapoint


def calibrate_emg(emg_datapoint, max_emg):
    """
    Calibrates the EMG by updating the maximum value if the current reading
    is higher.

    Parameters:
        emg_datapoint: Current EMG data point
        max_emg: Current maximum EMG value

    Returns:
        Updated maximum EMG value
    """
    return max(emg_datapoint, max_emg)


def generate_stim_sequence(emg_datapoint, max_emg, cathodic_phase, interphase_interval, anodic_phase, start_chan, next_rows):
    """
    Generates the stimulation sequence based on EMG input.

    Parameters:
        emg_datapoint: Current EMG data point
        max_emg: Maximum EMG value
        pseg, ipi, nseg: Segments defining the stimulation waveform
        start_chan: First stimulation channel
        next_rows: Flag to stimulate next row of LEDs

    Returns:
        seq: List of stimulation sequences
    """
    seq = []

    for channel in range(1, 8):
        if emg_datapoint > (channel * max_emg / 8):
            for row in [0, 8]:
                seq0 = xp.StimSeq(
                    start_chan + channel + row + next_rows * 16,
                    50, 500, cathodic_phase, interphase_interval, anodic_phase
                )
                seq.append(seq0)

    return seq


if __name__ == '__main__':
    connect_to_processor()

    try:
        stim_channels, rec_channels = verify_connections()
        setup_stim(stim_channels)

        rec_elec1 = rec_channels[0]
        rec_elec2 = rec_channels[1]

        # Set filters
        xp.filter_set(rec_elec1, 'hires', 4)

        max_emg1 = 0
        max_emg2 = 0
        cathodic_phase, interphase_interval, anodic_phase = create_stim_waveform()

        while True:
            emg_datapoint1 = measure_emg(rec_elec1)
            emg_datapoint2 = measure_emg(rec_elec2)

            # Set new maximum
            max_emg1 = calibrate_emg(emg_datapoint1, max_emg1)
            max_emg2 = calibrate_emg(emg_datapoint2, max_emg2)

            # Light up first LED in each row
            stim_sequence = []
            for led_row in range(4):
                seq0 = xp.StimSeq(
                    stim_channels[0] + led_row * 8, 50, 500,
                    cathodic_phase, interphase_interval, anodic_phase
                )
                stim_sequence.append(seq0)

            # Light LEDs based on EMG input
            stim_sequence += generate_stim_sequence(
                emg_datapoint1, max_emg1, cathodic_phase, interphase_interval,
                anodic_phase, stim_channels[0], 0
            )
            stim_sequence += generate_stim_sequence(
                emg_datapoint2, max_emg2, cathodic_phase, interphase_interval,
                anodic_phase, stim_channels[0], 1
            )
            xp.StimSeq.send_stim_seqs(stim_sequence)

    except RuntimeError as e:
        print(e)

    finally:
        print("Closing connection.")
        xp._close()