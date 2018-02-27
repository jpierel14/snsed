**********************
Core-Collapse Extrapolations
**********************

Getting Started
===============

Read in a lightcurve:


    >>> from __future__ import print_function
    >>> import snsedextend
    >>> sedFile=snsedextend.example_sed
    >>> myLc=snsedextend.load_example_data()
    >>> print(myLC)

Generate a Color Table
======================
Produce a color table from the example data with some assumptions. You can set any parameter that you would like that is used for SNCosmo fitting.
    
    >>> colorTable=snsedextend.curveToColor(myLC,colors=['U-B', 'r-J', 'r-H', 'r-K'], snType='Ic', zpsys='vega', bounds={'hostebv': (-1, 1), 't0': (53787.94, 53797.94)},
    constants={'mwr_v': 3.1, 'mwebv': '0.1267', 'z': '0.033529863', 'hostr_v': 3.1}, dust='CCM89Dust', effect_frames=['rest', 'obs'], effect_names=['host', 'mw'])
    >>> print(colorTable)
    time                U-B          ...        rK_err      
    ------------------- -------------------- ... -------------------
    -4.373621809507313                   -- ... 0.14190728435615696
    -4.373621809507313                   -- ...                  --
    -4.373621809507313                   -- ...                  --
    -4.366291809514223  -0.3018462428016375 ...                  --
    -2.4106118095078273 -0.05909438979689674 ...                  --
    -2.4006218095091754                   -- ...                  --
    -2.4006218095091754                   -- ... 0.11584654528008068
    ...                  ... ...                 ...
    0.5863781904918142                   -- ...                  --
    3.58512819049065  0.46099767884640763 ...                  --
    4.586158190490096  0.40038613488523844 ...                  --
    5.611378190485993                   -- ...                  --
    5.611378190485993                   -- ...                  --
    5.611378190485993                   -- ... 0.10417369509400125
    5.628438190491579  0.29095420901175495 ...                  --
    Length = 23 rows

Color Curve Fitting
===================
Now we can fit this color table and get a best model by minimizing BIC.
This function returns a python dictionary with colors as keys and an astropy Table object
with time and color vectors as values.

    >>>curveDict=snsedextend.fitColorCurve(colorTable)
    >>>print(curveDict.keys())

If we print the column names from the color table generated above, you'll see they line up.

    >>>print(colorTable.colnames)

SED Extrapolation
=================
Now you can provide an SED to be extrapolated, and let it do the work (This is a type Ic). This will return an
sncosmo.Source object and simultanously save the new SED to a file defined by newFileLoc (default current directory):

    >>>newSED=snsedextend.extendCC(colorTable,curveDict,sedlist=[sedFile],zpsys='vega',showplots=True,verbose=True)
    >>>print(newSED)

Plotting from Timeseries
========================
You can directly plot an SED file.

    >>>snsedextend.plotSED(sedFile)
