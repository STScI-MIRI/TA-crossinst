# TA-crossinst
code to test TA algorithms and procedures for the JWST instruments

This repository contains:
- jwst_ta.py: a set of functions for centroid calculation using the same method as GENTALOCATE. the code was adapted from an idl code written by P. Goudfrooij for NIRISS
- miri_test.py: a script that shows the calling sequence, used with MIRI data and parameters.
- nircam_test.py: as above, for NIRCam.

Status (Feb 27th 2018):

- Code returns the same result as the idl script on the MIRI test data
- Code supports:
    * different checkbox and centroid window sizing
    * using full array or custom-sized region of interest for calculation
	* background subtraction: fractional or absolute value.
	* flat fielding
	* creating a TA image from a ramp or ingesting a pre-made TA image
   
- Not yet implemented:
    * multiple iterations
    * random offset of target within ROI
    

Contributors:
- Bryan Hilbert
- Sarah Kendrew

Comments welcome!



