"""
This is a script that displays a livestream of acquired data from
any selected datastream and channel, as well as a live periodogram
(frequency power spectrum) of the last 5 seconds of live acquired
data. In addition, the computer terminal will print out the calculated
Voltage RMS of various frequency bands (Full: 0.3 Hz - 7500 Hz, LFP: 1 Hz 
- 250 Hz, Spike: 250 Hz - 7500 Hz) every 3.333 seconds.

This is a script that takes live streams of SEEG data and runs inference
on a brain state embedding model on the SEEG data.
1. First uses digital notch filters for SEEG data
2. Feeds data into histogram normalization pipeline
3. Runs hist norm data through GPU and performs inference on model
4. Generates an embedding then feeds embedding of last 5s of data into pacmap-> generates an embedding
5. plots this embedding relative to a preplotted 100 point plot 

Created on November 2024 

@author: richard liu
@author: Ghassan Makhoul (build on initial live_periodogram example code)
"""

import ipdb
import time
from time import sleep
import numpy as np
import matplotlib.pyplot as plt
import xipppy as xp
import sys
from scipy.fft import fft, next_fast_len
from scipy.integrate import trapezoid
import matplotlib
from matplotlib.colors import LinearSegmentedColormap


# Explicitly set the backend to TkAgg, avoid conflict with PyQT6
matplotlib.use('TkAgg')
colors = [[0,0,1,0],[0,0,1,0.5],[0,0.2,0.4,1]]
cmap = LinearSegmentedColormap.from_list("", colors)

class LiveStreamPeriodogram:
    """
    This class sets up and displays a livestream of acquired data and its live periodogram.
    """

    def __init__(self, display_s=5, window_ms=5000, stream_ch=0, stream_ty='hi-res'):
        """
        Initialization function for setting up the LiveStreamPeriodogram class.

        Parameters:
            display_s: Number of seconds of data to display on live plots
            window_ms: Number of milliseconds in FFT sliding window
            stream_ch: Channel number to record from
            stream_ty: Type of data stream to record from ('raw', 'lfp', 'hires', 'hifreq')
        """
        self.display_s = display_s
        self.window_ms = window_ms
        self.stream_ch = stream_ch
        self.stream_ty = stream_ty

        # Setup frequency range (x-axis limits) and decibel range (y-axis limits) for periodogram
        self.f_range = [0, 2000]
        self.db_range = [-50, 50]

        # Initialize connection to Neural Interface Processor
        self.connect_to_processor()

        # Pause for display duration to fill signal buffer
        sleep(self.display_s)

        # Setup parameters for performing FFT calculation
        self.samp_freq, self.window_samp, self.N, self.n_sig, self.t_sig, \
            self.f_psd_total, self.f_ind, self.f_psd = self.setup_stream_fft_settings()

        # Setup parameters for performing RMS noise calculations
        self.f_ind_whole, self.f_ind_lfp, self.f_ind_spike = self.setup_rms_noise_settings()

        # Setup plots for live signal and live periodogram subplots
        self.psd_plot, self.ax_sig, self.ax_psd, \
            self.h_sig, self.h_psd, self.intensity = self.setup_plots()

    def connect_to_processor(self):
        """
        Connect to the Neural Interface Processor. If UDP connection fails, 
        attempt TCP connection.
        """
        # Close Connection
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

    def setup_stream_fft_settings(self):
        """
        This function will enable the selected stream type for the Neural Interface Processor
        and calculate the necessary FFT settings for the live periodogram.

        Returns:
            samp_freq: Sampling frequency in Hz of the selected datastream
            window_samp: Number of samples in the FFT window
            N: Most efficient size of the N-Point FFT
            n_sig: Number of signal points in the display window
            t_sig: Time series of the signal points in the display window
            f_psd_total: Frequency series of the FFT calculation
            f_ind: Frequency indices for the selected frequency range
            f_psd: Frequency series for Periodogram display
        """
        if self.stream_ty == 'raw':
            samp_freq = 30000
        elif self.stream_ty == 'lfp':
            samp_freq = 1000
        elif self.stream_ty == 'hi-res':
            samp_freq = 2000
        elif self.stream_ty == 'hifreq':
            samp_freq = 7500
        else:
            sys.exit(f'{self.stream_ty} is an invalid stream type.\n')

        # Enable stream type on Neural Interface Processor if not enabled
        if not xp.signal(self.stream_ch, self.stream_ty):
            xp.signal_set(self.stream_ch, self.stream_ty, True)

        window_samp = int(np.floor(self.window_ms * samp_freq / 1e3))  # FFT window size in samples
        N = next_fast_len(window_samp // 2)  # Size of PSD
        n_sig = int(samp_freq / 1e3 * round(self.display_s * 1e3))  # Number of signal data points in display window
        t_sig = np.linspace(-self.display_s, 0, n_sig)  # Time for Signal
        f_psd_total = np.linspace(0, samp_freq, 2 * N)  # Frequency for Power Spectral Density
        f_ind = np.where((f_psd_total >= self.f_range[0]) & 
                         (f_psd_total <= self.f_range[1]))[0]  # Frequency indices within selected range
        f_psd = f_psd_total[f_ind]  # Frequency for Periodogram display

        return samp_freq, window_samp, N, n_sig, t_sig, f_psd_total, f_ind, f_psd

    def setup_rms_noise_settings(self):
        """
        This function will calculate the required indices for calculating the RMS Noise
        of various common frequency bands based on application (whole spectrum:
        0.3 - 7.5k Hz, LFP band: 1 - 250 Hz, Spike band: 250 - 7.5k Hz). The user may adjust
        these frequency ranges.

        Returns:
            f_ind_whole: Frequency indices for 0.3 - 7.5k Hz range
            f_ind_lfp: Frequency indices for 1 - 250 Hz range
            f_ind_spike: Frequency indices for 250 - 7.5k Hz range
        """
        f_range_whole = [0.3, 7500]
        f_ind_whole = np.where((self.f_psd_total >= f_range_whole[0]) & 
                               (self.f_psd_total <= f_range_whole[1]))[0]

        f_range_lfp = [1, 250]
        f_ind_lfp = np.where((self.f_psd_total >= f_range_lfp[0]) & 
                             (self.f_psd_total <= f_range_lfp[1]))[0]

        f_range_spike = [250, 7500]
        f_ind_spike = np.where((self.f_psd_total >= f_range_spike[0]) & 
                               (self.f_psd_total <= f_range_spike[1]))[0]

        return f_ind_whole, f_ind_lfp, f_ind_spike

    def setup_plots(self):
        """
        This function sets up the subplots for displaying a live signal and the live
        periodogram.

        Returns:
            psd_plot: Figure object for the plots
            ax_sig: Axes object for the signal plot
            ax_psd: Axes object for the periodogram plot
            h_sig: Line2D object for the signal plot data
            h_psd: Line2D object for the periodogram plot data
        """
        # plot hour of data ahead of time and buffer period
        # then update only the last hour 
        # can PacMAP be meaningful enough to see us go into the preictal funnel
        # perfect world: ghoast of last hour, so the most recent points are darkest
        # flow of data: filter -> hist_eq -> run through model -> PacMAP model -> update in FIFO manner with last 10 minutes-hour as 
        # static plot
        # saving all data to a cold storage archive 
        # code should update buffer and a masterfile
        # don't need to replot the same background points 
        # generate a numpy random array (-5 +5 on all axes, and incorporate time stamps for seizures appropriately)
        x_sig = np.zeros(self.n_sig)  # Signal Buffer for initialization
        x_psd = np.zeros(len(self.f_psd))  # Periodogram Buffer for initialization
        
        psd_plot, (ax_sig, ax_psd) = plt.subplots(2, 1)
        # ipdb.set_trace()
        h_sig = ax_sig.scatter(self.t_sig, x_sig, c=[], cmap=cmap, vmin=0, vmax=1)  # Signal plot
        ax_sig.set_xlim(-2,2)
        ax_sig.set_ylim(-2,2)
        ax_sig.set_title(f'Channel #{self.stream_ch}, {self.stream_ty}')
        ax_sig.set_xlabel('Time (s)')
        ax_sig.set_ylabel('μV')

        h_psd, = ax_psd.plot(self.f_psd, x_psd, 'b', linewidth=0.5)  # Periodogram plot
        ax_psd.set_ylim(self.db_range)
        ax_psd.set_xlim(self.f_range)
        ax_psd.set_title('Periodogram')
        ax_psd.set_xlabel('Frequency (Hz)')
        ax_psd.set_ylabel('Power Spectrum (dB/Hz)')

        plt.subplots_adjust(hspace=0.4)
        intensity = []

        return psd_plot, ax_sig, ax_psd, h_sig, h_psd, intensity

    def live_data_loop(self):
        """
        This function performs the live collection of data and updates the plot
        in a loop while the figure is open.
        """
        # Get current Neural Interface Processor time
        t0 = xp.time()
        t2 = xp.time()

        x_og= np.random.random_sample(size=int(100))
        y_og= np.random.random_sample(size=int(100))

        # While plot is open, continue live plotting data
        while plt.fignum_exists(self.psd_plot.number):
            plt.ion()
            t1 = xp.time()
            # If time since last loop is greater than 1/30 seconds, update loop
            # (Frame rate is capped at 30 FPS)
            if (t1 - t0) >= (3e4 / 30):
                # Data collection based on stream type
                if self.stream_ty == 'raw':
                    x_sig, ts = xp.cont_raw(round(self.display_s * self.samp_freq), [self.stream_ch])
                elif self.stream_ty == 'lfp':
                    x_sig, ts = xp.cont_lfp(round(self.display_s * self.samp_freq), [self.stream_ch])
                elif self.stream_ty == 'hi-res':
                    x_sig, ts = xp.cont_hires(round(self.display_s * self.samp_freq), [self.stream_ch])
                elif self.stream_ty == 'hifreq':
                    x_sig, ts = xp.cont_hifreq(round(self.display_s * self.samp_freq), [self.stream_ch])
            
                # Calculation of FFT and periodogram data
                sig_sample = x_sig[-self.window_samp:]
                psd_sample = (1 / (self.samp_freq * len(sig_sample))) * \
                             np.square(np.abs(fft(sig_sample, 2 * self.N)))  # PSD calculation
                x_psd = 10 * np.log10(psd_sample[self.f_ind])
                
                # Plot updated data and rescale axes as needed
                # ipdb.set_trace()
                x= np.random.random_sample(size=int(len(x_sig)/10000))
                y= np.random.random_sample(size=int(len(x_sig)/10000))

                x_og = np.append(x_og, x)
                y_og = np.append(y_og, y)
                self.h_sig.set_offsets(np.c_[x_og,y_og])
                self.intensity = np.concatenate((np.array(self.intensity)*0.96, np.ones(len(x))))
                self.h_sig.set_array(self.intensity)
                # self.ax_sig.relim()
                self.ax_sig.autoscale_view()
                self.h_psd.set_ydata(x_psd)
                plt.pause(0.1)  # Small pause to update the plot

                t0 = t1

            # If time since last loop is greater than 3.333 seconds, update
            if (t1 - t2) >= (3e4 * 3.333):
                # Noise RMS Calculations
                std_dev = np.std(sig_sample)
                x_rms_whole = np.sqrt(trapezoid(psd_sample[self.f_ind_whole], self.f_psd_total[self.f_ind_whole]))
                x_rms_lfp = np.sqrt(trapezoid(psd_sample[self.f_ind_lfp], self.f_psd_total[self.f_ind_lfp]))
                x_rms_spike = np.sqrt(trapezoid(psd_sample[self.f_ind_spike], self.f_psd_total[self.f_ind_spike]))

                # Print statements for noise RMS values
                print(f'Neural Interface Processor time: {xp.time()}')
                print(f'Std Dev: {std_dev:.5f}')
                print(f'0.3 to 7.5k Hz uV RMS noise: {x_rms_whole:.5f}')
                print(f'1 to 250 Hz uV RMS noise: {x_rms_lfp:.5f}')
                print(f'250 to 7.5k Hz uV RMS noise: {x_rms_spike:.5f}\n')

                t2 = t1
                plt.pause(.001)
            # plt.waitforbuttonpress()


if __name__ == '__main__':
    processor = LiveStreamPeriodogram()
    processor.live_data_loop()
    plt.show()
    xp._close()