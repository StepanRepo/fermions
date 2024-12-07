


# fermions

This program calculates statistical properties of fermions in equilibrium.  It's designed to compute macroscopic parameters such as number density (*n*), energy density (*ρ*), and pressure (*p*) as functions of temperature (*T*) and chemical potential (*μ*) for a given fermion mass (*m*) and spin (*s*).  The calculations are performed for different temperature regimes (lowT, normal, highT), allowing for the exploration of different analytical approximation methods.

## Installation and Usage

### Requirements

**Essential Software:**
-   **Fortran Compiler:**  A Fortran compiler supporting Fortran 2008 or later (e.g., gfortran, ifort, pgfortran). Preferable compiler is **gfortran**
-   **CMake:**  The CMake build system is required for building the project.


**Additional Requirements (for plotting and optional performance):**

To generate plots of the results using the provided  `plot.py`  script and to optionally improve performance, you will need the following:

-   **Python 3.x:**  The Python 3 programming language.
-   **NumPy:**  The fundamental package for scientific computing in Python.
-   **Matplotlib:**  A comprehensive library for creating static, interactive, and animated visualizations in Python.
-   **LAPACK (Recommended):**  The Linear Algebra PACKage (LAPACK) library is highly recommended for optimal performance. While not strictly mandatory, its inclusion significantly enhances the speed of certain calculations.

**Installation Instructions (Linux):**

The commands to install the required software packages will vary depending on your Linux distribution.  **Note: You will likely need  `sudo`  rights (administrator privileges) to install system-level packages.**  Here are examples for some common distributions:

**Essential Packages:**

**Debian/Ubuntu (apt):**

```bash
sudo apt update
sudo apt install gfortran cmake
```

**Fedora/CentOS/RHEL (dnf/yum):**

```bash
sudo dnf update  #(use 'sudo yum update' for older RHEL/CentOS versions)
sudo dnf install gfortran cmake
```

**Arch Linux (pacman):**

```bash
sudo pacman -Syu
sudo pacman -S gfortran cmake
```

**Additional Packages (Python and LAPACK):**

After installing the essential packages above, install the remaining requirements:

```bash
pip3 install numpy matplotlib
sudo apt install liblapack-dev  # Debian/Ubuntu; use appropriate command for your distribution (see below)

```

**Installing LAPACK (Distribution-Specific):**

-   **Debian/Ubuntu:**  `sudo apt install liblapack-dev`
-   **Fedora/CentOS/RHEL:**  `sudo dnf install lapack-devel`  (or  `sudo yum install lapack-devel`  for older versions)
-   **Arch Linux:**  `sudo pacman -S lapack`


### Step 1: Clone the Repository

```bash
git clone https://github.com/StepanRepo/fermions.git
cd fermions
```

### Step 2: Build the Project

1.  Ensure you have CMake and a Fortran compiler installed.
2.  Create a build directory and compile the project:
    
```bash
cd build
cmake ..
make
```
    
3.  The compiled executable  `fermions`  will be placed in the root directory of the repository.

### Step 3: Configure the Program

The `fermions` program uses `default.cfg` for customization. This file allows users to set physical and numerical parameters for calculations. The file follows the Fortran `namelist` format and must begin with `&config_nml` and end with `/`.

**Configuration Parameters:**

General Format:
```cfg
&config_nml
parameter_name = value
/
```

**Parameters:**

- **`output_dir`**: Directory for output files.  
  - Example: `output_dir = "results"`

- **`regime`**: Temperature regime.  
  - Options: `"lowT"`, `"normal"`, `"highT"`.  
  - Example: `regime = "highT"`

- **`Tmin`** and **`Tmax`**: Temperature range (in Kelvin).  
  - Example: `Tmin = 0.5`, `Tmax = 100.0`

- **`discr`**: Number of points in the temperature range.  
  - Example: `discr = 200`

- **`mu`**: Chemical potential (in energy units).  
  - Example: `mu = 1e-10`

- **`m`**: Particle mass (in grams).  
  - Example: `m = 1.6726219e-24` (proton mass).

- **`s`**: Particle spin.  
  - Example: `s = 1.5`

- **`nu`**: De Broglie frequency (in Hz). Only used when `m = 0`.  
  - Example: `nu = 1e16`

**Example Configuration**

```cfg
&config_nml
output_dir = "results"
regime = "highT"
Tmin = 0.5
Tmax = 100.0
discr = 200
mu = 1e-10
m = 1.6726219e-24
s = 0.5
nu = 1e16
/
```

**Notes**

- **Comments**: Add comments with `!`. Example: `Tmin = 0.1  ! Minimum temperature`.  
- **Output Directory**: Ensure the directory exists or create it manually.  
- **Massless Particles**: Set `m = 0` and define `nu`.  
- **Discretization**: Higher `discr` increases resolution but may slow computation.  


## How to Use the Results

The  `fermions`  program generates output files containing calculated physical properties of fermions based on the configuration provided in  `default.cfg`. These results can be used for analysis, visualization, and further study.

### Output Files

#### File Format

-   The program produces  `.dat`  files in the directory specified by the  `output_dir`  parameter in  `default.cfg`.
-   Each  `.dat`  file contains tabulated data with the following columns:
    1.  **Temperature (`T`)**: The temperature values (in Kelvin).
    2.  **Chemical Potential (`μ`)**: The chemical potential values (in energy units).
    3.  **Number Density (`n`)**: The number density of fermions (in cm⁻³).
    4.  **Mass Density (`ρ`)**: The mass density of fermions (in g·cm⁻³).
    5.  **Pressure (`p`)**: The pressure of fermions (in dyn·cm⁻²).


-   Example:
    
    ```txt
    T          μ              n              ρ            p
    0.1        1.00e-15       1.23e+23       1.98e-12     1.45e+10
    0.2        1.00e-15       2.34e+23       3.21e-12     2.76e+10
    ```
    

### Plotting the Results

-   Use the provided Python script (`plot.py`) to generate plots for quick visualization.
-   The script processes all  `.dat`  files in the output directory and produces  `.pdf`  files with the following plots:
    1.  **Number Density (`n`) vs Temperature (`T`)**.
    2.  **Energy Density (`ρ`) vs Temperature (`T`)**.
    3.  **Pressure (`p`) vs Temperature (`T`)**.

#### Running the Python Script

```bash
python3 plot.py
```

-   The  `.pdf`  files will be saved in the same directory as the  `.dat`  files.
-   Example output file:  `results/data.pdf`  (if  `output_dir`  is set to  `results`).

The generated plots will include:

-   X-axis: Temperature (`T`, in Kelvin).
-   Y-axis: The respective physical property (`n`,  `ρ`, or  `p`).
