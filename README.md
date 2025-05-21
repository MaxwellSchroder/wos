# Wos Solver

This project solves a 2D heat diffusion problem using a C++ solver and provides Python tools for data generation and visualization.

## Prerequisites

* Python 3.8
* C++ compiler (e.g., g++)

## Installation

1.  **Clone the repository:**

    ```bash
    git clone https://github.com/MaxwellSchroder/wos.git
    ```

2.  **Create and activate a Python virtual environment:**

    ```bash
    python3 -m venv venv
    source venv/bin/activate
    pip install -r requirements.txt
    ```

3.  **Install Python dependencies:**

    ```bash
    pip install -r requirements.txt
    ```

## Usage

### 1. Generate a 2D shape (CSV):

1.  Navigate to the `wos-solver` directory:

    ```bash
    cd wos-solver
    ```

2.  Modify the `point_generator.py` script to define your desired 2D shape. The shape is hardcoded within the script.
3.  Run the point generator:

    ```bash
    python3 point_generator.py
    ```

    This will generate a CSV file (e.g., `combined_coordinates.csv`) containing the shape's coordinates.

### 2. Solve the heat diffusion problem (C++):

1.  Compile the C++ solver:

    ```bash
    g++ -std=c++17 -O3 -pedantic -Wall -I./include wos_fileread.cpp -o wos_solver
    ```

2.  Run the solver, providing the generated CSV file as input:

    ```bash
    ./wos_solver combined_coordinates.csv
    ```

    This will generate a CSV file (e.g., `output.csv`) containing the solved temperature values.

### 3. Visualize the results (Python):

1.  Run the visualization script, providing the solver's output CSV file as input:

    ```bash
    python3 visualise.py output.csv
    ```

    This will display a visualization of the solved temperature distribution.

## Project Structure

* `venv/`: Python virtual environment.
* `requirements.txt`: Python dependencies.
* `wos-solver/`: Contains the core solver.
    * `point_generator.py`: Generates 2D shape CSV files.
    * `WoS-fileread.cpp`: C++ solver.
    * `visualise.py`: Python visualization script.
* `README.md`: This file.

## Notes

* Ensure you have a C++ compiler installed.
* The `point_generator.py` script requires manual modification to define the desired shape.
* The `visualise.py` script provides a basic visualization. You can modify it for more advanced visualization.
* The C++ solver accepts two command line arguments. The input filename, and the output filename.
