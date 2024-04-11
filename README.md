# Calculation of Theoretical CWSI in Real-Time
This repository provides required functions to calculate theoretical valaues for CWSI, actual transpiration (Ta) and potential transpiration (Tp) using microclimate data and surface canopy temperatures. The functions are currently calibrated for apple trees, however, they can be modified and/or re-calibrated to work with other fruit trees or crops. Obviously, any change to the equations needs to take the plant physiology into consideration. 

The equations are extracted from the following research papers:
1) Osroosh et al., 2015. Estimating actual transpiration of apple trees based on infrared thermometry. J. Irri Drain Eng., 141(8): 04014084. [Link to article](https://ccda4084-320d-4212-92c7-12721a347519.filesusr.com/ugd/b17667_e0d580341d4745d6b61b51f5c1cf90da.pdf)
2) Osroosh et al., 2015. Estimating potential transpiration of apple trees using theoretical non-water-stressed baselines. J. Irri Drain Eng., 141(9): 04015009. [Link to article](https://ccda4084-320d-4212-92c7-12721a347519.filesusr.com/ugd/b17667_1b59dd5c23fc476384621b855555d9c7.pdf)
3) Osroosh et al., 2016. Daylight crop water stress index for continuous monitoring of water status in apple trees. Irri Science, 34(3): 209–219.[Link to article](https://ccda4084-320d-4212-92c7-12721a347519.filesusr.com/ugd/b17667_10f416167c4048a2b1d59b3f3d8cd75a.pdf)

For related research articles and blog posts, please visit EnviTronics Lab website:
[EnviTronics Lab](https://envitronicslab.com)


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
3      24.5   0.31
4      24.6   0.34
5      24.7   0.36
6      24.8   0.39
7      24.9   0.41
8      25.0   0.43
9      25.1   0.46
10     25.2   0.48
11     25.3   0.51
12     25.4   0.53
13     25.5   0.55
14     25.6   0.58
15     25.7   0.60
16     25.8   0.63
17     25.9   0.65
18     26.0   0.68
19     26.1   0.70
************** End **************
```
