#!/bin/bash
# Assumes this script is in the new project folder

# Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt

