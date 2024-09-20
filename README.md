# Calculation of Theoretical CWSI in Real-Time
This repository provides a comprehensive set of functions in both C++ and Python to accurately estimate theoretical Crop Water Stress Index (CWSI), actual transpiration (Ta), and potential transpiration (Tp) for apple trees. These calculations are based on microclimate data and surface canopy temperatures. While initially calibrated for apple trees, the functions can be adapted and recalibrated to suit other fruit trees or crops, provided plant physiological characteristics are considered.

The underlying equations are derived from the following research papers:

1) Osroosh et al., 2015. Estimating actual transpiration of apple trees based on infrared thermometry. J. Irri Drain Eng., 141(8): 04014084. [Link to article](https://ccda4084-320d-4212-92c7-12721a347519.filesusr.com/ugd/b17667_e0d580341d4745d6b61b51f5c1cf90da.pdf)
2) Osroosh et al., 2015. Estimating potential transpiration of apple trees using theoretical non-water-stressed baselines. J. Irri Drain Eng., 141(9): 04015009. [Link to article](https://ccda4084-320d-4212-92c7-12721a347519.filesusr.com/ugd/b17667_1b59dd5c23fc476384621b855555d9c7.pdf)
3) Osroosh et al., 2016. Daylight crop water stress index for continuous monitoring of water status in apple trees. Irri Science, 34(3): 209–219.[Link to article](https://ccda4084-320d-4212-92c7-12721a347519.filesusr.com/ugd/b17667_10f416167c4048a2b1d59b3f3d8cd75a.pdf)

For related research articles and blog posts, please visit EnviTronics Lab website:
[EnviTronics Lab](https://envitronicslab.com)

## How to run the Python code

The Python code in the `cwsi_fucntions.py` module demonstrates the functionality of the CWSI model implemented in the `CWSIModel` class. You can run the `cwsi_test.py` module to test the model and see example outputs.

**Requirements:**

* Python 3 (tested with version 3.10)

**Instructions:**

1. Make sure you have Python installed on your system.
2. Download or clone the repository containing the code.
3. Open a terminal window and navigate to the directory containing the script (`cwsi_test.py`).
4. Run the following command to execute the script:

```
python cwsi_test.py
```

**Example Output:**

The script will perform two tests: calculating actual and potential transpiration, followed by calculating CWSI for a specific canopy temperature and a range of temperatures. Here's an example output:

```
*********** Transpiration Models Test **********
Actual Transpiration = 2.4 (mm/day)
Potential Transpiration = 5.2 (mm/day)

**************** CWSI Model Test *****************
CWSI (cal) = 0.216767
CWSI (exp) = 0.216767

Test passed!

Num        Tc        CWSI
0         24.1       0.22
1         24.2       0.23
2         24.3       0.24
... (output continues for multiple temperatures)

******************  End  ********************
```

This is just an example, and the actual values you see may differ depending on `Tc` values.

## How to compile/build and run the C++ program (Linux only)
The CWSI_Function has only been tested on Linux (Ubuntu 22.04.2 LTS). A fully functional application on any other OS or OS version is not guaranteed. To find the OS name and version in Linux, run any of the following commands:
```
cat /etc/os-release
lsb_release -a
hostnamectl
```

To test the program, open a terminal window and type the follwoing commands:
```
cd CWSI_Function
chmod u+x run.sh
./run.sh
```

You should see the following output in the terminal:
```
*********** CWSI Model Test **********
CWSI (cal) = 0.216767
CWSI (exp) = 0.216767

Test passed!

Num    Tc     CWSI
0      24.2   0.24
1      24.3   0.27
2      24.4   0.29
... (output continues for multiple temperatures)

************** End **************
```