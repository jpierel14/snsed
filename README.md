J.R. Pierel & S.Rodney 

2017.03.31

Extrapolate SED up to H and K bands, allowing user to define V-H and V-K. Extrapolation is improvement upon previous method in that it chooses a linear slope such that the integration of the H (V) band will correspond to the user-defined V-H (V-K).

SETUP
	The script assumes that you have set the SNDATA_ROOT environment variable (automatically set if you
	run the normal SNDATA setup process). Then you can place this script inside any directory, and it will
	look for your filenames (or .SED files) in the '~/SNDATA_ROOT/snsed/NON1A' directory. The directory you
	place this file in will be populated with the extended SED files (same filename). Careful, as this means
	that the SED files will be overwritten if you are running the script from within the NON1A directory. An
	error.log file is also created, containing error messages from running the script. 
SYNOPSIS
	python snsedextend.py -i <file1,file2,...> -p <day# or all> --vh <V-H> --vk <V-K> --jh <J-H> --jk <J-K>

DESCRIPTION
	The options are as follows:

	-i  	    This allows user to define a filename (example.SED) or list of filenames
		    (example.SED,example2.SED,...) separated by commas (no spaces). If this flag
		    is not set, then all .SED files in the SNDATA_ROOT/snsed/NON1A directory will
		    be extrapolated.


	-p	    This allows user to flag to plot or not. If you want to plot, add the day (6)
		    or you can plot all epochs sequentially (all)


	--vh	    User defined V-H


	--vk	    User defined V-K


	--jh	    User defined J-H


	--jk	    User defined J-K


OTHER INFO
      At the top of the script the following global variables are set. You can change this to affect the
      respective parameter (in angstroms):

      		 Center of V,H,K bands:
			VBAND=5500
			HBAND=15414
			KBAND=22000
		 Width of V,H,K bands:
		        vWidth=3000
			hWidth=3390
			kWidth=4000
