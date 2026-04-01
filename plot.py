#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import argparse

def handleUserInput():
    p = argparse.ArgumentParser(
        description = "",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        "-i", "--input",
        default = "output.dat",
        help = "data file to plot"
    )
    
    return p.parse_args()



def main(): 
    args = handleUserInput()
    file = args.input

    with open(file) as f:
        r, v = np.loadtxt( f, unpack = True );

    # plt.scatter( r, v, color = "red", lw = 10.0 )
    plt.scatter( r, v, color = "red", 
                 lw = 0.1, s = 2.5 )
    plt.grid()
    plt.show()

if __name__ == "__main__":
    main()



