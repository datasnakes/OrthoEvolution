#!/bin/bash
# For Debian/Ubuntu users
sudo apt-get update; sudo apt-get upgrade

# Install the packages
sudo apt-get install clustalo paml phyml ncbi-blast+ phylip
echo "All packages installed."

# Set the PATHs