J.R. Pierel & S.Rodney 

2017.04.26

__SUMMARY__

Extrapolate SED up to J, H and K bands, allowing user to define V-J, V-H and/or V-K. Extrapolation is improvement upon previous method in that
it chooses a linear slope such that the integration of the J (H,K) band will correspond to the user-defined V-J (V-H,V-K). The script assumes
that the defined color is in ABmag, but Vega can be defined using the --vega flag (see below).

__SETUP__

The script assumes that you have set the SNDATA_ROOT environment variable (automatically set if you
run the normal SNDATA setup process). Then you can place this script inside any directory, and it will
look for your filenames (or .SED files) in the '~/SNDATA_ROOT/snsed/NON1A' directory. The directory you
place this file in will be populated with the extended SED files (same filename). Careful, as this means
that the SED files will be overwritten if you are running the script from within the NON1A directory. An
error.log file is also created, containing error messages from running the script.
	
__SYNOPSIS__

	
	python snsedextend.py -i <file1,file2,...> -p <day# or all> -v <vTrans.dat> -j <jTrans.dat> -h <hTrans.dat> -k <kTrans.dat> --vh <V-H> --vk <V-K> --jh <J-H> --jk <J-K> --vj <V-J> --vega
	

__DESCRIPTION__

	The options are as follows:

	-i  	    This allows user to define a filename (example.SED) or list of filenames
		    (example.SED,example2.SED,...) separated by commas (no spaces). If this flag
		    is not set, then all .SED files in the SNDATA_ROOT/snsed/NON1A directory will
		    be extrapolated.


	-p	    This allows user to flag to plot or not. If you want to plot, add the day (6)
		    or you can plot all epochs sequentially (all)

	-v	    This allows user to define a transmission file to use for the v-band.

	-j	    This allows user to define a transmission file to use for the j-band.

	-h	    This allows user to define a transmission file to use for the h-band.

	-k	    This allows user to define a transmission file to use for the k-band.


	--vh	    User defined V-H


	--vk	    User defined V-K


	--jh	    User defined J-H


	--jk	    User defined J-K

	--vj	    User defined V-J

	--vega	    This allows user to input color in Vega instead of (assumed) AB


__Transmission Files__

You may use your own transmission files to define the filters used in the extrapolation. The default filters are tophat filters for J,H, and K, and Bessell for V.
To use your own transmission file, use the flags described above, and place the files in the appropriate folder (i.e. transmission file for V should go in vBand folder, etc.)
Alternatively, you may change the dictionary 'filters' at the top of the script:


       filters={
		'V':'vBand/bessellv.dat',
    	'J':'jBand/tophatJ.dat',
    	'H':'hBand/tophatH.dat',
		'K':'kBand/tophatK.dat'
	}


This is the default dictionary, simply change the '.dat' filenames to the filenames of your transmission files, and place the files in the appropriate folder. The
format of the transmission files should be two-column, with the first column being wavelength and the second column being transmission.