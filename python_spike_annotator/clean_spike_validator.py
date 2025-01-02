#!/usr/bin/env python3

"""
clean_spike_validator.py

Instructions:
1. Install required Python packages:
    - pyEDFlib (or mne if you prefer)
    - numpy
    - scipy
    - pandas
    - matplotlib

   Example (in a conda env):
       conda install -c conda-forge pyedflib numpy scipy pandas matplotlib

2. Update EDF_DIR and OUT_FILE below to match your local paths.

3. Run:
    python clean_spike_validator.py

Key bindings:
- up arrow:  decrease gain (channels get closer together vertically)
- down arrow: increase gain (channels spread out more)
- b: switch to bipolar montage
- c: switch to CAR montage
- l / r / m / n: choose "left/right/midline/no spike"
- Enter: confirm choice (must choose l/r/m/n first)

Overwrite settings:
- overwrite = 0: skip previously completed files
- overwrite = 1: re-do all files (overwrite results)
- overwrite = 2: show all files but do NOT save results
"""

import os
import numpy as np
import pandas as pd
import pyedflib
import matplotlib.pyplot as plt
from scipy.signal import iirfilter, filtfilt

###############################################################################
#                              USER PARAMETERS
###############################################################################
EDF_DIR = "/Users/erinconrad/Desktop/research/scalp_eeg_toolbox/results/example_edfs/"   # <- EDIT ME
OUT_FILE = "selections.csv"

# Overwrite settings:
#   0 = skip previously completed files
#   1 = re-write selections of all files
#   2 = show all files but do not save selections
overwrite = 0

###############################################################################
#                     HELPER FUNCTIONS
###############################################################################

def clean_temple(ch_labels):
    """
    Equivalent to your MATLAB 'clean_temple' function.
    Removes or replaces some substrings in the channel names.
    """
    new_labels = []
    for lbl in ch_labels:
        lbl = lbl.replace("EEG", "")
        lbl = lbl.replace("_REF", "")
        lbl = lbl.replace("FP", "Fp")
        lbl = lbl.replace("Z", "z")
        new_labels.append(lbl)
    return new_labels

def read_in_edf_clean(filepath):
    """
    Similar to your MATLAB read_in_edf_clean function.
    Load EDF using pyEDFlib, return (data, labels, fs).
    """
    f = pyedflib.EdfReader(filepath)
    n_channels = f.signals_in_file

    ch_labels = [f.getLabel(i) for i in range(n_channels)]
    fs = f.getSampleFrequency(0)  # assume same sampling rate for all channels

    # Read data into a (samples x channels) numpy array
    signals = []
    for i in range(n_channels):
        signals.append(f.readSignal(i))
    f.close()

    data = np.column_stack(signals)
    return data, ch_labels, fs

def scalp_bipolar_clean(values, ch_labels):
    """
    Creates a standard bipolar montage. This replicates your MATLAB version.
    Returns (bipolar_values, bipolar_labels).
    """
    # Example pairs from your code
    montage_pairs = [
        ('Fp1','F7'),
        ('F7','T3'),
        ('T3','T5'),
        ('T5','O1'),
        ('',''),  # blank
        ('Fp2','F8'),
        ('F8','T4'),
        ('T4','T6'),
        ('T6','O2'),
        ('',''),  # blank
        ('Fp1','F3'),
        ('F3','C3'),
        ('C3','P3'),
        ('P3','O1'),
        ('',''),
        ('Fp2','F4'),
        ('F4','C4'),
        ('C4','P4'),
        ('P4','O2'),
        ('',''),
        ('Fz','Cz'),
        ('',''),
        ('',''),
        ('EKG1','')  # example from your code
    ]

    # Rename ECGL->EKG1, ECGR->EKG2
    ch_labels = [lbl.replace("ECGL","EKG1").replace("ECGR","EKG2") for lbl in ch_labels]

    n_samples = values.shape[0]
    bipolar_values = []
    bipolar_labels = []

    for pair in montage_pairs:
        ch1, ch2 = pair
        if ch1 == '' and ch2 == '':
            # blank row, produce NaNs
            new_signal = np.full(n_samples, np.nan)
            new_label = ""
        else:
            # Try to find the channel indices in ch_labels
            try:
                idx1 = ch_labels.index(ch1)
            except ValueError:
                idx1 = None
            try:
                idx2 = ch_labels.index(ch2)
            except ValueError:
                idx2 = None

            if idx1 is None and idx2 is None:
                # Both channels missing
                new_signal = np.full(n_samples, np.nan)
                new_label = f"{ch1}-{ch2}"
            elif idx1 is not None and idx2 is None and "EKG" in ch1:
                # EKG1 only
                new_signal = values[:, idx1]
                new_label = f"{ch1}"
            else:
                # Standard difference
                val1 = values[:, idx1] if idx1 is not None else np.zeros(n_samples)
                val2 = values[:, idx2] if idx2 is not None else np.zeros(n_samples)
                new_signal = val1 - val2
                new_label = f"{ch1}-{ch2}"

        bipolar_values.append(new_signal)
        bipolar_labels.append(new_label)

    bipolar_values = np.column_stack(bipolar_values)
    return bipolar_values, bipolar_labels

def car_montage_clean(values, ch_labels):
    """
    Common Average Reference montage from your MATLAB code.
    Returns (car_values, car_labels).
    """
    car_labels = ch_labels[:]  # copy
    brain_channels = {
        'Fp1','F7','T3','T5','O1',
        'Fp2','F8','T4','T6','O2',
        'F3','C3','P3','F4','C4',
        'P4','Fz','Cz'
    }

    # Indices of channels that are considered "brain channels"
    is_brain = [lbl in brain_channels for lbl in ch_labels]
    is_brain = np.array(is_brain)

    # Compute average of those channels
    avg = np.nanmean(values[:, is_brain], axis=1)

    # Subtract that average from the channels
    car_values = values.copy()
    car_values[:, is_brain] = car_values[:, is_brain] - avg[:, None]

    # Append '-CAR' to brain channel labels
    for i, lbl in enumerate(ch_labels):
        if is_brain[i]:
            car_labels[i] = lbl + "-CAR"

    return car_values, car_labels

def highpass_filter(data, fs, cutoff=1.0, order=4):
    """
    Implements a simple Butterworth highpass filter. 
    Equivalent to MATLAB's `highpass(data, 1, fs)`.
    """
    from scipy.signal import butter
    nyquist = 0.5 * fs
    low = cutoff / nyquist
    b, a = butter(order, low, btype='high')
    filtered = filtfilt(b, a, data, axis=0)
    return filtered

def bandstop_filter(data, fs, band=(58,62), order=4):
    """
    Implements a simple Butterworth bandstop filter 
    around the specified band. Equivalent to MATLAB's 
    `bandstop(data, [58 62], fs)`.
    """
    from scipy.signal import butter
    nyquist = 0.5 * fs
    low = band[0] / nyquist
    high = band[1] / nyquist
    b, a = butter(order, [low, high], btype='bandstop')
    filtered = filtfilt(b, a, data, axis=0)
    return filtered

def plot_scalp_eeg_clean(ax, data, fs, labels, gain):
    """
    Equivalent to your MATLAB function plot_scalp_eeg_clean:
    Plots each channel offset by some amount (gain).
    """
    ax.clear()
    n_samples, n_chs = data.shape
    duration = n_samples / fs

    offset_accum = 0
    t = np.linspace(0, duration, n_samples)

    for ich in range(n_chs):
        # Shift channel by offset_accum
        signal_offset = data[:, ich] - offset_accum
        ax.plot(t, signal_offset, 'k')
        # Put a channel label out to the right
        ax.text(duration + 0.05, -offset_accum, labels[ich], 
                fontsize=10, verticalalignment='center')

        offset_accum += gain  # Next channel offset
    
    # Draw a vertical line in the middle (similar to your code)
    ax.axvline(x=duration/2, color='r', linestyle='--')
    ax.set_xlabel("Time (s)")
    ax.set_yticks([])
    ax.set_xlim([0, duration + 1])
    # Let the y-limits be from -offset_accum - gain to some small positive margin
    ax.set_ylim([-offset_accum - gain, gain])
    ax.set_title("EEG Plot", fontsize=12)

###############################################################################
#                     CLASS FOR HANDLING INTERACTIVE LOGIC
###############################################################################

class EEGValidatorApp:
    def __init__(self, edf_files, responses, overwrite):
        """
        edf_files: list of full file paths
        responses: current or empty DataFrame with columns ['Filename','Response']
        overwrite: 0, 1, or 2
        """
        self.edf_files = edf_files
        self.responses = responses
        self.overwrite = overwrite

        # Default montage is 'bipolar' (like MATLAB code)
        self.montage_type = 'bipolar'

        # Default gain
        self.gain = 50

        # User input (spike type)
        self.user_input = None

        # A flag to indicate if user pressed Enter to confirm
        self.confirmed = False

        # We'll store references to the data for re-plotting
        self.bipolar_values = None
        self.bipolar_labels = None
        self.car_values = None
        self.car_labels = None
        self.fs = None

        # We'll set up a figure and axes once
        self.fig, self.ax = plt.subplots(figsize=(14, 8))
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)

    def run(self):
        """
        Main loop over all EDF files.
        """
        for file in self.edf_files:
            edf_name = os.path.basename(file)

            # Check if we already have a response for this file
            if self.overwrite in [0, 2]:
                if edf_name in self.responses['Filename'].values:
                    print(f"\nAlready did {edf_name}, skipping.")
                    continue

            # Otherwise, process the file
            print(f"\nDisplaying {edf_name}")
            self.process_file(file)

            # Save the user input if not in overwrite=2
            if self.overwrite != 2 and self.user_input in ['l','r','m','n']:
                # Append or replace row in self.responses
                new_row = pd.DataFrame([[edf_name, self.user_input]], 
                                       columns=['Filename','Response'])
                self.responses = pd.concat([self.responses, new_row], ignore_index=True)
                # Save partial results
                self.responses.to_csv(os.path.join(EDF_DIR, OUT_FILE), index=False)
            
            # Reset user input for next file
            self.user_input = None

        print("All done!")

    def process_file(self, filepath):
        """
        Reads the EDF, does filtering, sets up montage data, and plots.
        Then waits for user input to confirm a response.
        """
        # Read EDF
        values, ch_labels, fs = read_in_edf_clean(filepath)
        ch_labels = clean_temple(ch_labels)

        # Replace NaNs with channel mean
        for col in range(values.shape[1]):
            col_nans = np.isnan(values[:, col])
            if np.any(col_nans):
                mean_val = np.nanmean(values[:, col])
                values[col_nans, col] = mean_val

        # Demean each channel
        values = values - np.nanmean(values, axis=0)

        # Channels that are all nan -> make zero
        all_nan_chs = np.all(np.isnan(values), axis=0)
        values[:, all_nan_chs] = 0

        # Filter
        values = highpass_filter(values, fs, cutoff=1.0, order=4)
        values = bandstop_filter(values, fs, band=(58, 62), order=4)

        # Get bipolar montage
        bipolar_values, bipolar_labels = scalp_bipolar_clean(values, ch_labels)
        # Get CAR montage
        car_values, car_labels = car_montage_clean(values, ch_labels)

        # Store for re-plotting on key events
        self.bipolar_values = bipolar_values
        self.bipolar_labels = bipolar_labels
        self.car_values = car_values
        self.car_labels = car_labels
        self.fs = fs

        # Reset for this file
        self.montage_type = 'bipolar'
        self.gain = 50
        self.user_input = None
        self.confirmed = False

        # Initial plot
        self.update_plot()

        # Show the figure without blocking so we can handle key events
        plt.show(block=False)

        # Block until the user selects a valid input (l/r/m/n) AND presses Enter
        while not self.confirmed:
            plt.pause(0.1)

        # Once confirmed == True, we know they pressed Enter with a valid choice
        # Clear the figure before next file
        self.ax.cla()

    def update_plot(self):
        """
        Redraw the plot based on current montage, gain, and user_input.
        """
        # Which data to plot?
        if self.montage_type == 'bipolar':
            data = self.bipolar_values
            labels = self.bipolar_labels
        else:
            data = self.car_values
            labels = self.car_labels

        plot_scalp_eeg_clean(self.ax, data, self.fs, labels, self.gain)

        # Title or super-title prompt
        prompt_str = "Select spike: l (left), r (right), m (midline), n (no spike). Then press Enter."
        self.fig.suptitle(prompt_str, fontsize=14, color='blue')

        # Display current selection in red text
        if self.user_input is None:
            selection_str = "Selected: None"
        else:
            selection_str = f"Selected: {self.user_input}"
        # We'll place it in the upper-left corner
        self.ax.text(0.01, 0.98, selection_str, transform=self.ax.transAxes,
                     fontsize=14, color='red', verticalalignment='top')

        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def on_key_press(self, event):
        """
        Key press callback. 
        - up/down: gain
        - b/c: montage type
        - l/r/m/n: user choice (must press Enter to confirm)
        - Enter: confirm choice
        """
        if event.key in ['l','r','m','n']:
            # Store the user input but wait for Enter
            self.user_input = event.key
            self.update_plot()
            print(f"You pressed '{self.user_input}'. Press Enter to confirm.")

        elif event.key == 'enter':
            # Only confirm if we have a valid input
            if self.user_input in ['l','r','m','n']:
                print(f"User input confirmed: {self.user_input}")
                self.confirmed = True  # This allows the main loop to exit
            else:
                print("Please press one of [l, r, m, n] before pressing Enter.")

        elif event.key == 'up':
            self.gain = max(0, self.gain - 10)  # Decrease gain
            self.update_plot()

        elif event.key == 'down':
            self.gain += 10  # Increase gain
            self.update_plot()

        elif event.key == 'b':
            self.montage_type = 'bipolar'
            self.update_plot()

        elif event.key == 'c':
            self.montage_type = 'car'
            self.update_plot()

        else:
            print(f"Other key pressed: {event.key}")

###############################################################################
#                     MAIN SCRIPT LOGIC
###############################################################################

def main():
    # 1) Load existing CSV if it exists (unless overwrite=1)
    if overwrite == 0:
        # If file exists, load it
        csv_path = os.path.join(EDF_DIR, OUT_FILE)
        if os.path.isfile(csv_path):
            responses = pd.read_csv(csv_path)
        else:
            responses = pd.DataFrame(columns=['Filename','Response'])
    elif overwrite == 1:
        # Start fresh
        responses = pd.DataFrame(columns=['Filename','Response'])
    else:
        # overwrite == 2 (do not save)
        responses = pd.DataFrame(columns=['Filename','Response'])

    # 2) Get list of .edf files
    edf_files = [os.path.join(EDF_DIR, f) 
                 for f in os.listdir(EDF_DIR) 
                 if f.lower().endswith('.edf')]

    # 3) Initialize App
    app = EEGValidatorApp(edf_files, responses, overwrite)

    # 4) Run
    app.run()

    # If overwrite != 2, results are saved as we go
    print("Done. Final responses:")
    print(app.responses)

if __name__ == "__main__":
    main()
