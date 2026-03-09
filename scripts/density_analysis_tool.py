#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Set the backend to 'Agg' for non-interactive plotting
import matplotlib.pyplot as plt
from matplotlib import gridspec
import panedr
from tabulate import tabulate
from pprint import pprint
import os
np.set_printoptions(threshold = np.inf)
from matplotlib.ticker import LogFormatter, LogLocator, LogFormatterExponent, FuncFormatter

def read_density_data(file_path):
    df = panedr.edr_to_df(file_path, verbose=True)
    # Ensure Time is correctly set if it's not present in the dataframe
    if 'Time' not in df.columns:
        df['Time'] = df.index
    return df[['Time', 'Density']].reset_index(drop=True)

def calculate_rolling_stats(df, window_size):
    df['Extended Rolling Mean'] = df['Density'].rolling(window=window_size, min_periods=1).mean()
    df['Extended Rolling Std'] = df['Density'].rolling(window=window_size, min_periods=1).std()
    return df

def find_stationary_point(df, window_size, mean_threshold, min_consecutive_points, distance):
    # Calculate the change in the extended rolling mean at a specified distance
    df['Mean Change Distance'] = df['Extended Rolling Mean'].diff(periods=distance).abs()
    
    df['EMA of Mean Change Distance'] = df['Mean Change Distance'].ewm(alpha=0.005).mean()
    
    # Find indices where the change over the distance is below the threshold
    below_threshold = df['EMA of Mean Change Distance'] < mean_threshold
    
    # Group consecutive points below threshold to identify sequences long enough
    streaks = below_threshold.cumsum() - below_threshold.cumsum().where(~below_threshold).ffill().fillna(0)
    
    # Check if any streak of sufficient length exists
    if (streaks.value_counts().max() < min_consecutive_points):
        raise ValueError("No sufficient consecutive stationary points found.")

    # Find the first index where the required minimum consecutive points are met
    first_stationary = streaks[streaks >= min_consecutive_points].first_valid_index()

    if first_stationary is None:
        return 0
    
    return first_stationary - min_consecutive_points + 1  # Adjust to get the start of the sequence

def calculate_average_density(df, stationary_start):
    stationary_data = df.loc[stationary_start:]
    return stationary_data['Density'].mean()

def identify_points_within_range(df, average_density, tolerance=0.01):
    print(f"Average observable: {average_density}")
    lower_bound = average_density - tolerance * abs(average_density)
    upper_bound = average_density + tolerance * abs(average_density)

    print(f"Lower bound: {lower_bound}")
    print(f"Upper bound: {upper_bound}")
    filtered_df = df[(df['Density'] >= lower_bound) & (df['Density'] <= upper_bound)]

    return filtered_df  # Reset index to avoid issues in downstream operations

def count_consecutive_series(df, window_size, min_consecutive_points, distance):
    threshold_range = np.linspace(0.000001, 0.00001, 20)
    results = {}
    for mean_threshold in threshold_range:
        # Calculate the change in the extended rolling mean at a specified distance
        df['Mean Change Distance'] = df['Extended Rolling Mean'].diff(periods=distance).abs()

        # Find indices where the change over the distance is below the threshold
        below_threshold = df['Mean Change Distance'] < mean_threshold
        # Group consecutive points below threshold
        streaks = below_threshold.cumsum() - below_threshold.cumsum().where(~below_threshold).ffill().fillna(0)
        all_lengths = streaks.value_counts()
        long_streaks = all_lengths[all_lengths >= min_consecutive_points]
        try:
            results[mean_threshold] = long_streaks.iloc[1]
        except IndexError as e:
            results[mean_threshold] = 0

    return results

def plot_density(df, stationary_start, average_density, within_range_indices, best_points_indices, file_dir, all_lengths, mean_threshold):
    # Set up a figure with custom GridSpec for layout
    fig = plt.figure(figsize=(15, 10))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[2, 1])

    # Top subplot for density spanning full width
    ax1 = plt.subplot(gs[0, :])
    ax1.plot(df['Time'], df['Density'], label='Density', color='blue', marker='o', linestyle='-', markersize=5, zorder=1)
    ax1.plot(df['Time'], df['Extended Rolling Mean'], color='orange', label='Extended Rolling Mean', zorder=9)
    ax1.axvline(df['Time'].loc[stationary_start], color='green', linestyle='--', label=f'Stationary Start @ {df["Time"].loc[stationary_start]} ps', zorder=4)
    ax1.axhline(y=average_density, color='purple', linestyle='-', label=f'Average Density: {average_density:.2f}', zorder=5)
    
    valid_within_range = df.loc[within_range_indices]
    ax1.scatter(valid_within_range['Time'], valid_within_range['Density'], color='cyan', marker='o', label='Within Tolerance Range', zorder=6)
    
    if best_points_indices is not None and len(best_points_indices) > 0:
        valid_best = df.loc[best_points_indices]
        ax1.scatter(valid_best['Time'], valid_best['Density'], color='red', edgecolors='black', marker='*', s=150, label='Extracted Replicas', zorder=10)

    ax1.fill_between(df['Time'], df['Density'].min(), df['Density'], where=(df.index < stationary_start), color='red', alpha=0.3, label='Discarded Data', zorder=2)
    ax1.set_ylabel('Density (kg/m^3)')
    ax1.set_title('Density Analysis with Stationary Point')
    ax1.legend()
    ax1.grid(True)

    # Bottom left subplot for rolling standard deviation
    ax2 = plt.subplot(gs[1, 0])
    ax2.set_yscale("log")
    ax2.plot(df['Time'], df['Mean Change Distance'], color='blue', label='Mean Change Distance', zorder=3)
    ax2.plot(df['Time'], df['EMA of Mean Change Distance'], color='orange', label='EMA', zorder=4)
    ax2.axhline(y=mean_threshold, color="gray", linestyle='--', label=f'Threshold value: {mean_threshold}', zorder=4)
    ax2.set_xlabel('Time (ps)')
    ax2.set_ylabel('EMA of Mean change Distance (kg/m^3)')
    ax2.set_title('Mean Change Distance')
    ax2.legend()
    ax2.grid(True)
    
    # Bottom right subplot for the histogram
    ax3 = plt.subplot(gs[1, 1])
    bins = np.logspace(np.log10(df['EMA of Mean Change Distance'].min()), np.log10(df['EMA of Mean Change Distance'].max()), 30)
    ax3.hist(df['EMA of Mean Change Distance'], bins=bins, edgecolor='black', rwidth=0.7)
    ax3.set_xlabel('EMA Mean Change Distance')
    ax3.set_title('Histogram of EMA')
    if df['Mean Change Distance'].max() / df['Mean Change Distance'].min() > 100:
        ax3.set_xscale('log')

    plt.tight_layout()
    plt.savefig(os.path.join(file_dir, 'density_analysis.png'))
    plt.close(fig)

def parse_range(arg_range):
    try:
        start, end = map(float, arg_range.split(':'))
        return start, end
    except ValueError:
        raise argparse.ArgumentTypeError("Range must be in the format 'start:end'.")

def parse_index_range(arg_range):
    try:
        start, end = map(int, arg_range.split(':'))
        return start, end
    except ValueError:
        raise argparse.ArgumentTypeError("Index range must be in the format 'start:end'.")

def main(args):
    # Load Data
    df = read_density_data(args.file)
    file_dir = os.path.dirname(os.path.abspath(args.file))
    base_name = os.path.splitext(os.path.basename(args.file))[0]
    
    # Define GROMACS file paths early so they are available for the loop
    topol_file = os.path.join(file_dir, base_name + '.tpr')
    traj_file = os.path.join(file_dir, base_name + '.xtc')
    if not os.path.exists(traj_file):
        traj_file = os.path.join(file_dir, base_name + '.trr')

    # Apply Cropping
    if args.time_crop:
        start_time, end_time = args.time_crop
        print(f"Cropping data from time {start_time} to {end_time}.")
        if end_time > 0:
            df = df[(df['Time'] >= start_time) & (df['Time'] <= end_time)]
        else:
            df = df[df['Time'] >= start_time]

    if args.index_crop:
        start_index, end_index = args.index_crop
        print(f"Cropping data from index {start_index} to {end_index}.")
        df = df.iloc[start_index:end_index]
        
    df = df.reset_index(drop=True)
    df = calculate_rolling_stats(df, args.window_size)
    
    # Identify Stationarity
    try:
        stationary_start = find_stationary_point(df, args.window_size, args.mean_threshold, args.min_consecutive_points, args.distance)
    except ValueError as e:
        print(e)
        return

    average_density = calculate_average_density(df, stationary_start)
    within_range = identify_points_within_range(df.loc[stationary_start:], average_density, args.tolerance)
    
    if not within_range.empty:
        within_range = within_range.copy()
        within_range['Delta to Average'] = abs(within_range['Density'] - average_density)

        # Determine best points
        num_to_extract = args.num_replicas
        best_points = within_range.nsmallest(10, 'Delta to Average') # Keep top 10 for printing
        points_to_extract = within_range.nsmallest(num_to_extract, 'Delta to Average').sort_values(by='Time')

        # Fancy Print Output (Original Logic)
        print(f"Points within {args.tolerance*100}% of the average density:")
        points_to_show = within_range.nsmallest(args.show_points, 'Delta to Average').sort_values(by='Time')
        
        for idx, row in points_to_show.iterrows():
            if idx in best_points.index:
                rank = list(best_points.index).index(idx) + 1
                red = int(255 * rank / 10); green = int(255 - (255 * rank / 10))
                color = f"\033[38;2;{red};{green};0m"
                delta_display = f"{row['Delta to Average']:.4f} (BEST: #{rank})"
            else:
                color = "\033[0m"; delta_display = f"{row['Delta to Average']:.4f}"
            print(f"{color}Index: {idx}, Time: {row['Time']}, Density: {row['Density']:.4f}, Delta: {delta_display}\033[0m")

        # --- AUTOMATIC EXTRACTION ---
        print(f"\n--- Automatically extracting top {num_to_extract} frames ---")
        
        for i, (idx, row) in enumerate(points_to_extract.iterrows()):
            output_file = os.path.join(file_dir, f"start_replica_{i+1}.gro")
            selected_time = row['Time']
            
            # 1. Use the full binary name: gmx_mpi
            # 2. Use mpirun -np 1 to prevent the PMI2 crash
            # 3. Add -quiet to reduce log clutter
            command = (
                f"echo 0 | srun --ntasks=1 --cpus-per-task=1 gmx_mpi trjconv "
                f"-s {topol_file} -f {traj_file} -o {output_file} "
                f"-dump {selected_time} -quiet"
            )
            
            print(f"Extracting Replica {i+1} at {selected_time} ps...")
            
            # We use a trick here: setting an environment variable just for this command
            # to prevent MPI from trying to initialize complex networking
            os.system(f"export I_MPI_SHM_LMT=shm; {command}")
        
        # Plotting
        all_lengths = count_consecutive_series(df, args.window_size, args.min_consecutive_points, args.distance)
        plot_density(df, stationary_start, average_density, within_range.index, points_to_extract.index, file_dir, all_lengths, args.mean_threshold)
        
    else:
        print("No points found within the specified range.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze density data to find stationary points.')
    parser.add_argument('file', type=str, help='Path to the EDR file.')
    parser.add_argument('--num_replicas', type=int, default=3, help='Number of top frames to extract for replicas.')
    parser.add_argument('--window_size', type=int, default=10, help='Window size for rolling calculations.')
    parser.add_argument('--mean_threshold', type=float, default=0.1, help='Threshold for changes in the rolling mean.')
    parser.add_argument('--tolerance', type=float, default=0.01, help='Tolerance for average density range.')
    parser.add_argument('--show_points', type=int, default=20, help='Number of points to show in terminal.')
    parser.add_argument('--distance', type=int, default=1, help='Distance between points to check for mean change.')
    parser.add_argument('--min_consecutive_points', type=int, default=1, help='Min consecutive points for convergence.')
    parser.add_argument('--time_crop', type=parse_range, help="Format 'start:end'.")
    parser.add_argument('--index_crop', type=parse_index_range, help="Format 'start:end'.")
    args = parser.parse_args()
    main(args)
