#!/usr/bin/env python3

# Name: analyze_docprofs.py
# Description: This is a python script will take in the document
#              array profiles in csv format and analyze the compressibility
#              of the profiles.
# Date: December 8th, 2022

import os 
import sys
import argparse
import random
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

def parse_doc_profiles(file_path):
    """ Take *.csv file and return a list of lists where each position correponds to a profile """
    profiles = []
    with open(file_path, "r") as in_fd:
        for line in in_fd:
            curr_profile = [int(x) for x in line.strip().split(",")]
            profiles.append(curr_profile)
            if len(profiles) == 10000:
                break
    return profiles

def convert_to_increases(profile_list, starting_end):
    """ Focus on monotonic increases and convert profiles, and from a certain end """
    new_profiles = []
    for profile in profile_list:
        # reverse it if we are going right to left ..
        if starting_end == "right":
            profile = profile[::-1]
        new_profile = [min(1000,profile[0])]
        for i in range(1, len(profile)):
            val = profile[i] if profile[i] > new_profile[i-1] else new_profile[i-1]
            val = min(1000, val)
            new_profile.append(val)
        if starting_end == "right":
            new_profile = new_profile[::-1]
        new_profiles.append(new_profile)
    return new_profiles

def compute_increases_per_profile(profile_list, direction):
    """ Return the number of increases per profile """
    increase_list = []
    for profile in profile_list:
        increases = 0
        if direction == "right":
            curr_profile = profile[::-1]
        else:
            curr_profile = profile.copy()
        for i in range(1, len(curr_profile)):
            if curr_profile[i] > curr_profile[i-1]:
                increases+=1
        increase_list.append(increases)
    return increase_list   

def compute_sig_lcps_per_profile(all_profiles):
    """ Find the number of significant lcps (>15) in each profile """
    sig_list = []
    for profile in all_profiles:
        count = 0
        for val in profile:
            if val >= 15:
                count += 1
        sig_list.append(count)
    return sig_list

def main(args):
    """ analyze the document array profiles """
    
    # Step 1: Parse the start and end profiles from *.csv
    start_profiles = parse_doc_profiles(args.input_prefix + ".sdap.csv")
    end_profiles = parse_doc_profiles(args.input_prefix + ".edap.csv")
    print("[log] parsed both the start and end doc profiles.")

    # Step 2: Convert to the the monotonic increase version of the profile from left or right
    start_mono_profiles_left = convert_to_increases(start_profiles, "left")
    start_mono_profiles_right = convert_to_increases(start_profiles, "right")

    end_mono_profiles_left = convert_to_increases(end_profiles, "left")
    end_mono_profiles_right = convert_to_increases(end_profiles, "right")
    print("[log] analyzed the profiles from both left and right directions.")

    # Step 3: Analyze the number of monotonic increases from each direction 
    start_profiles_left_increases = compute_increases_per_profile(start_mono_profiles_left, "left")
    start_profiles_right_increases = compute_increases_per_profile(start_mono_profiles_right, "right")

    end_profiles_left_increases = compute_increases_per_profile(end_mono_profiles_left, "left")
    end_profiles_right_increases = compute_increases_per_profile(end_mono_profiles_right, "right")
    print("[log] computed the number of monotonic increases for each profile.")

    # Step 4: Analyze the number of significant lcps at each profile
    start_profiles_sig_lcps = compute_sig_lcps_per_profile(start_profiles)
    end_profiles_sig_lcps = compute_sig_lcps_per_profile(end_profiles)

    # Step 5: Plot the data ...

    # Plot 1: Look at the monotonic increases
    plt.figure(figsize=(10, 10))

    ## Sub-plot 1:
    plt.subplot(2, 2, 1)
    plt.xlabel("Monotonic LCP")
    plt.ylabel("Document Number")
    plt.title("Start of Run Profiles - Left to Right")

    count = 0
    for index in random.sample(range(1, len(start_mono_profiles_left)), 50):
        profile = start_mono_profiles_left[index]
        if profile[-1] == 1000:
            plt.plot([x for x in range(1, len(profile)+1)], profile, label=f"profile_{index}")
            count += 1
            if count >= 10:
                break

    ## Sub-plot 2:
    plt.subplot(2, 2, 2)
    plt.xlabel("Monotonic LCP")
    plt.ylabel("Document Number")
    plt.title("Start of Run Profiles - Right to Left")

    count = 0
    for index in random.sample(range(1, len(start_mono_profiles_right)), 50):
        profile = start_mono_profiles_right[index]
        if profile[0] == 1000:
            plt.plot([x for x in range(1, len(profile)+1)], profile, label=f"profile_{index}")
            count += 1
            if count >= 10:
                break

    ## Sub-plot 3:
    plt.subplot(2, 2, 3)
    plt.xlabel("Monotonic LCP")
    plt.ylabel("Document Number")
    plt.title("End of Run Profiles - Left to Right")

    count = 0
    for index in random.sample(range(1, len(end_mono_profiles_left)), 50):
        profile = end_mono_profiles_left[index]
        if profile[-1] == 1000:
            plt.plot([x for x in range(1, len(profile)+1)], profile, label=f"profile_{index}")
            count += 1
            if count >= 10:
                break

    ## Sub-plot 4:
    plt.subplot(2, 2, 4)
    plt.xlabel("Monotonic LCP")
    plt.ylabel("Document Number")
    plt.title("End of Run Profiles - Right to Left")

    count = 0
    for index in random.sample(range(1, len(end_mono_profiles_right)), 50):
        profile = end_mono_profiles_right[index]
        if profile[0] == 1000:
            plt.plot([x for x in range(1, len(profile)+1)], profile, label=f"profile_{index}")
            count += 1
            if count >= 10:
                break

    plt.savefig(args.input_prefix+".monotonic_profiles_plot.png", dpi=800, bbox_inches="tight")
    print("[log] plotted the monotonic increases of 50 random profiles.")

    # Plot 2: Looking at the number of increases from each direction
    sns.set_style("darkgrid")
    plt.figure(figsize=(10, 10))

    ## Sub-plot 1: Start of run profiles, going from left to right
    plt.subplot(2, 2, 1)
    plt.xlabel("Number of Increases in LCP")
    plt.ylabel("Density")
    plt.title("Start of Run Profiles - Left to Right")

    values, counts = np.unique(start_profiles_left_increases, return_counts=True)
    sns.barplot(x=values, y=(counts/np.sum(counts)))

    ## Sub-plot 2: Start of run profiles, going from right to left
    plt.subplot(2, 2, 2)
    plt.xlabel("Number of Increases in LCP")
    plt.ylabel("Density")
    plt.title("Start of Run Profiles - Right to Left")

    values, counts = np.unique(start_profiles_right_increases, return_counts=True)
    sns.barplot(x=values, y=(counts/np.sum(counts)))

    ## Sub-plot 3: End of run profiles, going from left to right
    plt.subplot(2, 2, 3)
    plt.xlabel("Number of Increases in LCP")
    plt.ylabel("Density")
    plt.title("End of Run Profiles - Left to Right")

    values, counts = np.unique(end_profiles_left_increases, return_counts=True)
    sns.barplot(x=values, y=(counts/np.sum(counts)))

    ## Sub-plot 4: Start of run profiles, going from right to left
    plt.subplot(2, 2, 4)
    plt.xlabel("Number of Increases in LCP")
    plt.ylabel("Density")
    plt.title("End of Run Profiles - Right to Left")

    values, counts = np.unique(end_profiles_right_increases, return_counts=True)
    sns.barplot(x=values, y=(counts/np.sum(counts)))

    plt.savefig(args.input_prefix+".monotonic_increases_plot.png", dpi=800, bbox_inches="tight")
    print("[log] plotted the monotonic increases distributions of 1000 profiles.")

    # Plot 3: Look at the amount of significant lcps
    plt.figure(figsize=(11, 5))

    ## Sub-plot 1: Start of run profiles
    plt.subplot(1, 2, 1)
    plt.xlabel("Number of Significant LCPs (>15)")
    plt.ylabel("Density")
    plt.title("Start of Run Profiles")
    plt.xticks(np.arange(0, max(end_profiles_sig_lcps), step=250))

    sns.kdeplot(data=np.array(start_profiles_sig_lcps), color="blue", fill=True)

    ## Sub-plot 2: End of run profiles
    plt.subplot(1, 2, 2)
    plt.xlabel("Number of Significant LCPs (>15)")
    plt.ylabel("Density")
    plt.title("End of Run Profiles")
    plt.xticks(np.arange(0, max(end_profiles_sig_lcps), step=250))

    sns.kdeplot(data=np.array(end_profiles_sig_lcps), color="blue", fill=True)

    plt.savefig(args.input_prefix+".sig_lcps_plot.png", dpi=800, bbox_inches="tight")
    print("[log] plotted the amount of significant lcps in each profile.\n")

def parse_arguments():
    """ Defines the command-line argument parser, and return arguments """
    parser = argparse.ArgumentParser(description="This script allows users to analyze the compressiblity of the document array profiles.")
    parser.add_argument("-i", "--input", dest="input_prefix", required=True, help="input profile prefix")
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Validate the command-line arguments """
    if not os.path.isfile(args.input_prefix + ".edap.csv"):
        print("Error: the *.edap.csv file is not present.")
        exit(1)
    if not os.path.isfile(args.input_prefix + ".edap.csv"):
        print("Error: the *.edap.csv file is not present.")
        exit(1)


if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)