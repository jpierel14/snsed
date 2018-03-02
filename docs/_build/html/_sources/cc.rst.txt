****************************
Core-Collapse Extrapolations
****************************

Getting Started
===============

Read in a lightcurve:


    >>> from __future__ import print_function
    >>> import snsedextend
    >>> sedFile=snsedextend.example_sed
    >>> myLC=snsedextend.load_example_data()
    >>> print(myLC)

    name  band     time     mag   magerr
    ------ ---- ----------- ------ ------
    2006aj    U 53788.16733 18.039  0.078
    2006aj    U 53790.12301 17.995  0.047
    2006aj    U 53790.13652 17.982  0.057
    2006aj    U 53796.11875 18.599  0.112
    2006aj    U 53797.11978 18.679  0.054
    2006aj    U 53798.16206 18.738  0.058
    2006aj    B 53788.11683 18.438   0.05
    ...  ...         ...    ...    ...
    2006aj    H    53793.12 16.752  0.047
    2006aj    H   53798.145 16.635  0.037
    2006aj    K    53788.16 17.183   0.12
    2006aj    K   53790.133 16.714  0.094
    2006aj    K   53792.123 16.406  0.145
    2006aj    K    53793.12 16.393   0.08
    2006aj    K   53798.145 16.769   0.08
    Length = 105 rows


Generate a Color Table
======================
Produce a color table from the example data with some assumptions. You can set any parameter that you would like that is used for SNCosmo fitting.
    
    >>> colorTable=snsedextend.curveToColor(myLC,colors=['U-B', 'r-J', 'r-H', 'r-K'], snType='Ic', zpsys='vega', bounds={'hostebv': (-1, 1), 't0': (53787.94, 53797.94)},constants={'mwr_v': 3.1, 'mwebv': '0.1267', 'z': '0.033529863', 'hostr_v': 3.1}, dust='CCM89Dust', effect_frames=['rest', 'obs'], effect_names=['host', 'mw'])
   
    Getting best fit for: U-B,r-J,r-H,r-K
    No model provided, running series of models.
    Best fit model is "snana-2004fe", with a Chi-squared of 310.573042
    
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

    >>> curveDict=snsedextend.fitColorCurve(colorTable)
    
    Running BIC for color:
     U-B
     r-J
     r-H
     r-K
     
    >>> print(curveDict.keys())
    
['r-J', 'r-K', 'r-H', 'U-B']
    


If we print the column names from the color table generated above, you'll see they line up.

    >>> print(colorTable.colnames)

['time', 'U-B', 'UB_err', 'r-J', 'rJ_err', 'r-H', 'rH_err', 'r-K', 'rK_err']

SED Extrapolation
=================
Now you can provide an SED to be extrapolated, and let it do the work (This is a type Ic). This will return an
sncosmo.Source object and simultanously save the new SED to a file defined by newFileLoc (default current directory):

    >>> newSED=snsedextend.extendCC(colorTable,curveDict,sedlist=[sedFile],zpsys='vega',showplots=False,verbose=True)
    >>> print(newSED)

Plotting from Timeseries
========================
You can directly plot an SED file.

    >>> snsedextend.plotSED(sedFile)
