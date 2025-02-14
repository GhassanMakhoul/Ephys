"""
This is a script that takes live streams of SEEG data and runs inference
on a brain state embedding model on the SEEG data.
1. First uses digital notch filters for SEEG data
2. Feeds data into histogram normalization pipeline
3. Runs hist norm data through GPU and performs inference on model
4. Generates an embedding then feeds embedding of last 5s of data into pacmap-> generates an embedding
5. plots this embedding relative to a preplotted 100 point plot 

Created on Janurary 2025

@author: richard
@author: ghassan makhoul, extended original code and adapted for tranformers
"""

import ipdb
import time
import sys
sys.path.append("/home/ghassanmakhoul/Documents/Tornadoes_v1/")
sys.path.append("/home/ghassanmakhoul/Documents/Tornadoes_v1/models")
sys.path.append("/home/ghassanmakhoul/Documents/Tornadoes_v1/train")
from time import sleep
import loguru
import numpy as np
import matplotlib.pyplot as plt
import xipppy as xp
import sys

import matplotlib
from loguru import logger

# Explicitly set the backend to TkAgg, avoid conflict with PyQT6
matplotlib.use('TkAgg')
from matplotlib.colors import LinearSegmentedColormap



class LiveStreamBrainStateEmbedding:
    """
    This class sets up and displays a livestream of acquired data and its live periodogram.
    """

    def __init__(self,subject, model_file, pre_proc_pkl, display_s=5, window_ms=5000, stream_ch=[a for a in range(128)], stream_ty='hi-res'):
        """
        Initialization function for setting up the LiveStreamPeriodogram class.

        Parameters:
            subject : subject name (should follow 'pat' convention) Important for streaming and saving of data
            model_file : location of the saved model weights for running inference, TODO: update to reflect how /where model loaded
            pre_proc_pkl : location of a pickle file containing the histogram normalization code
            display_s: Number of seconds of data to display on live plots
            window_ms: Number of milliseconds in FFT sliding window
            stream_ch: Channel number to record from
            stream_ty: Type of data stream to record from ('raw', 'lfp', 'hires', 'hifreq')
        """
        self.subject = subject
        self.display_s = display_s
        self.window_ms = window_ms
        self.stream_ch = stream_ch
        self.stream_ty = stream_ty

        # Setup frequency range (x-axis limits) and decibel range (y-axis limits) for periodogram
        self.f_range = [0, 2000]
        self.db_range = [-50, 100]

        # Initialize connection to Neural Interface Processor
        self.connect_to_processor()

        # Pause for display duration to fill signal buffer
        sleep(self.display_s)

        # Load histogram norm pkl for normalizing data
        self.hist_norm = self.setup_hist_norm(pre_proc_pkl)
        #Load model for inference 
        self.model = self.setup_inference(model_file)

        # Setup parameters for performing live pacmap calculation
        self.samp_freq, self.window_samp, self.n_sig, self.t_sig, = self.setup_stream_pacmap_settings()

        # Setup plots for live signal and live periodogram subplots
        self.pac_plot, self.ax_pre, self.ax_live, \
            self.h_pre, self.h_live = self.setup_plots()

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


    def setup_hist_norm(self, pre_proc_pkl:str):
        """set up histogram normalization runner

        Args:
            pre_proc_pkl (str): Location of the pickle file for the histogram normalization code
        """
        #TODO: implement

        self.hist_norm = []
        return NotImplemented
    

    def setup_inference(self, model_loc: str):
        """loads the pytorch model to perform live inference

        Args:
            model_loc (str): model location
        """
        #TODO implement
        return NotImplemented
    

    def setup_stream_pacmap_settings(self):
        """
        This function will enable the selected stream type for the Neural Interface Processor
        and calculate the necessary Pacmap settings

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
        
        n_sig = int(samp_freq / 1e3 * round(self.display_s * 1e3))  # Number of signal data points in display window
        t_sig = np.linspace(-self.display_s, 0, n_sig)  # Time for Signal
        

        return samp_freq, window_samp, n_sig, t_sig

    

    def setup_plots(self):
        """
        This function sets up the subplots for displaying a live pacmap and a previous pacmap from BSE.

        Returns:
            pac_plot: Figure object for the plots
            ax_pre Axes object for the pretrained pacmap projection
            ax_live: Axes object for the live rendered plot
            h_pre: Line2D object for the signal plot
            h_live: Line2D object for the periodogram plot data
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
        
        pac_plot, (ax_pre, ax_live) = plt.subplots(2, 1)
        # Pretrained pacmap values here 
        h_pre = ax_pre.scatter(self.t_sig, x_sig)  # Add colors
        # TODO find out the dimensions of the PACMAP
        ax_pre.set_xlim(-2, 2)
        ax_pre.set_ylim(-2,2)
        ax_pre.set_title(f'Pretrained Brain State Embedding for {self.subject}')
        ax_pre.set_xlabel('Dimension 1')
        ax_pre.set_ylabel('Dimension 2')
        # Add pacmap
        h_live = ax_live.scatter(self.f_psd, x_psd)  # Periodogram plot
        ax_live.set_ylim(-2,2)
        ax_live.set_xlim(-2,2)
        ax_live.set_title(f'Live Brain State Embedding for {self.subject}')
        ax_live.set_xlabel('Frequency (Hz)')
        ax_live.set_ylabel('Power Spectrum (dB/Hz)')

        plt.subplots_adjust(hspace=0.4)

        return pac_plot, ax_pre, ax_live, h_pre, h_live
    
    def reshape_stream(self, stream_data):
        """Reshape stream data from a [n_samp*n_ch,1] array to a 
        an array of shape [n_ch, n_ch]

        Args:
            stream_data (np.array): Raw stream
        """
        stream_data = np.array(stream_data)
        stream_data.reshape[self.stream_ch,-1]
        return stream_data
    
    def norm_data(self, data: np.array):
        """Runs the hisogram normalization on the data stream and returns

        Args:
            data (np.array): live, cleaned data
        """
        assert self.setup_hist_norm
        return NotImplemented

    def query_model(self, data: np.array):
        """Queries a model to run inference on, 
        should take advantage of GPU and pretrained model
        
        ASSUMES that the model has been instantiated/loaded

        Args:
            data (np.array): N_ch x M_samples array of live ecog data
        """
        if self.model == None: raise ValueError("Cannot run inference if model has not been loaded")

    def live_data_loop(self):
        """
        This function performs the live collection of data and updates the plot
        in a loop while the figure is open.
        """
        # Get current Neural Interface Processor time
        t0 = xp.time()

        # While plot is open, continue live plotting data
        while plt.fignum_exists(self.pac_plot.number):
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
                #TODO double check reshape order check element to element
                sig_reshape = self.reshape_stream(sig_sample)
                # Filter 60hz line noise and resample to 512Hz
                # do the exact preprocessing filtering
                clean_sig = self.clean_sig(sig_reshape)
                # Can apply norm parameters, open pickle and use logic for the right channel
                # The setup goes here 
                norm_sig = self.norm_data(clean_sig)
                # Make sure to load random model and feed in dummy data 
                # TODO - load pangolin, initialize, load weights, torch.rand(with correct shape)
                #           - pull out the mu, check the pacmap forloop where it calles pacmap then fine tunes
                #           - use that inference code run_epoch in training object
                #           - if get_latent_only (diff forward pass)
                #           -  forward passes a second of data, get the latent prediction at 2-samples of stride (248 predictions)
                #           - then take the average of those predictions, and get average mu
                #           - so I'll feed in a second and get 1024 dim x 1_sec of data
                #           - loop time needs to be faster than a second, or batch into 5 seconds, and just run over 5s and get 5 latent vectors
                #           - if loop time takes 3seconds then do a 5s stride so we get time to run
                mu = self.query_model(norm_sig)
                projections  = self.project_data(mu)
                
                
                # Plot updated data and rescale axes as needed
                # assume embeddings have shape N_samp x 2_dim
                self.h_live.set_offsets(np.c_[projections[0,:], projections[1,:]])
                # self.ax_sig.relim()
                self.ax_sig.autoscale_view()
                # self.h_psd.set_offsets(np.c_[x_sig, x_sig])
                plt.pause(0.001)  # Small pause to update the plot
                plt.waitforbuttonpress()
                t0 = t1


if __name__ == '__main__':
    n_ch = [a for a in range(129)] # num of channels for Spat 113
    processor = LiveStreamBrainStateEmbedding(n_ch)
    processor.live_data_loop()
    plt.show()
    xp._close()