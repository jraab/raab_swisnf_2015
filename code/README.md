This directory contains the underlying code for this paper. If you just want to recreate figures see the two scripts
    runfirst.sh
    createFigs.sh

 in the top directory. All the figure creating scripts are located in `figures/`

 If you are interested in the early processing steps (i.e. alignment, peak calling, etc) see the scripts that start 

     pipeline.<SOMETHING>.sh 

 These were executed on our cluster (Rocks) and the `call.pipline.<SOMETHING>.sh` scripts were used to launch the jobs. 


 The ENCODE files were processed using the `call.calc_enrichment_at_peaks.py`. This script launched jobs that paired each encode file with its correct input so that normalization could be done and calculates the signal over the ARID peaks. These scripts need several changes to run on another machine because I hardcoded the location of our downloads of the ENCODE files. The results of these are stored in `output/encode_coverages/` and I've provided the output of those scripts because they take a while to generate. 

 Please contact open an issue or contact me at jesse.r.raab@gmail.com with any issues. 

