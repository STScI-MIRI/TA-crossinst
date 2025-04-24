# TA-crossinst
code to test TA algorithms and procedures for the JWST instruments

This repository contains:
- jwst_ta.py: a set of functions for centroid calculation using the same method as GENTALOCATE. the code was adapted from an idl code written by P. Goudfrooij for NIRISS
- miri_test.py: a script that shows the calling sequence, used with MIRI data and parameters. can be used for testing.
- nircam_test.py: as above, for NIRCam.

Setup:
- Clone the repository to a local folder using `git clone https://github.com/STScI-MIRI/TA-crossinst.git`
- If necessary, activate the appropriate conda environment.
- `cd` to the folder where you cloned the repository and run `pip install .` (or `pip install -e .` for editable installations).
- Import the code with `from jwst_ta import jwst_ta`


Status (May 10th 2018):

- Code returns the same result as the idl script on the MIRI test data
- Code supports:
    * different checkbox and centroid window sizing
    * using full array or custom-sized region of interest for calculation
	* background subtraction: fractional or absolute value.
	* flat fielding
	* creating a TA image from a ramp or ingesting a pre-made TA image
    * creating a TA image from a custom list of group indices (and calculation of the scale factor)
   
- Not yet implemented:
    * random offset of target within ROI
    

Contributors:
- Bryan Hilbert
- Sarah Kendrew
- Jonathan Aguilar

Comments and contriutions welcome. Please use Issues and Pull requests.



