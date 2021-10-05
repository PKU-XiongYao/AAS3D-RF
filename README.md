      ******************************************************************************
      **                                                                          **                                         
      **                              AAS3D-RF                                    **
      ** Random Forest based disease-associated Amino Acid Substitution predictor **  
      **              with 3D structural features incorporated                    **
      **                                                                          **
      **                                                                          **
      **                         README for AAS3D-RF                              **
      **                                                                          **
      ******************************************************************************

The following instructions are for 64-bit GNU/Linux operating systems (e.g., Ubuntu 16.04.3 LTS) 
and Windows is currently not supported. All of the Linux commands (starting with a "$" ) are typed 
in the Bash shell. 

I. INSTALLATION

Caution: 1. Please carefully follow the instructions below to meet the prerequisites to make the package work. 
         2. The location of the "run_AAS3D-RF.py" script is denoted as <PATH_TO_AAS3D-RF>.
         
** PREREQUISITES **
The following softwares and databases need to be presented in the system before attempting to use AAS3D-RF.

*** Anaconda (https://www.anaconda.com/distribution/#download-section)
The AAS3D-RF package is based on Python 2, and a set of related python modules are
required. We recommend to use Anaconda and its environment to ease this process:
    $   cd <PATH_TO_AAS3D-RF>/
    $   conda create --name AAS3DRF python==2.7.13
    $   conda activate AAS3DRF

And then use the following command to install the required modules:
    $   conda install -c conda-forge numpy==1.13.1 scipy==0.19.1 scikit-learn==0.19.0 biopython==1.69 pandas==0.20.3


The Users can also use the Python 2.7 shipped with the operation system and install these modules by using pip.


*** PSI-BLAST && DSSP
A default version of PSI-BLAST (blast-2.2.29+, https://blast.ncbi.nlm.nih.gov/Blast.cgi) and DSSP (DSSP-v2.2.1, https://swift.cmbi.umcn.nl/gv/dssp/) 
have been included in the <PATH_TO_AAS3D-RF>/Tools/.

The copyrights belong to the original authors and developers. 


*** UniRef90 Database
Download UniRef90 cluster database from UniProt, move it to the preset directory in AAS3D-RF, and format it
    $   wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
    $   mv uniref90.fasta.gz  <PATH_TO_AAS3D-RF>/Tools/blast-2.2.29+/UniRef90
    $   cd <PATH_TO_AAS3D-RF>/Tools/blast-2.2.29+/UniRef90
    $   gunzip uniref90.fasta.gz 
    $   <PATH_TO_AAS3D-RF>/Tools/blast-2.2.29+/bin/makeblastdb -in uniref90.fasta -dbtype prot -parse_seqids -out uniref90 -title uniref90

This may need 1-2 hours depending on the performance of the server.

If you would like to use a UniRef90 database at a different location, you need to set the "BLAST_DATABASE" variable 
in the <PATH_TO_AAS3D-RF>/PredMod/userConfig.py file.
  
AAS3D-RF was developed based on UniRef90 Release 2018_06, so we recommend to use UniRef90 Release 2018_06 in
users' application for better consistency.


*** RING package
Please visit http://protein.bio.unipd.it/ring/ website to learn more about RING (Residue Interaction Network Generator).
To download its standalone package, you need to submit the download request at http://protein.bio.unipd.it/download/.
It may take one week to get the download link of the academically licensed RING package.

Download the compressed package
Move it to <PATH_TO_AAS3D-RF>/Tools/
Unpack it
Rename it as "RING", so that "RING" contains the "bin", "data", and "lib"
subdirectories.

Add the following line in your .bashrc file, and source your .bashrc file:
        export VICTOR_ROOT=<PATH_TO_AAS3D-RF>/Tools/RING/     # Note: the trailing '/' is required!

If you have already had RING software at a different location and want to use it, you need to set "RING" variable
    in <PATH_TO_AAS3D-RF>/PredMod/userConfig.py file.
  
AAS3D-RF is developed based on RING version 2.0, so we recommend to use this version for
better consistency.


*** Naccess package
Visit http://www.bioinf.manchester.ac.uk/naccess/ to download, decrypt, and unpack the Naccess V2.1.1.
The users need to write an email to the author for the decryption key.
Move it to <PATH_TO_AAS3D-RF>/Tools/
Rename it to "Naccess"
Compile it (commands provided as below)

Note: 1) We have increased some of the parameters in the original accall.pars file for AAS3D-RF predictor to handle 
         larger proteins, and the modified accall.pars file is at <PATH_TO_AAS3D-RF>/Tools/Naccess_Para_File/accall.pars. We 
         recommend to copy this modified file to <PATH_TO_AAS3D-RF>/Tools/Naccess to
		 replace the original one before you compile it.
      2) You may need to change "write(4,'(a,i)')" to "write(4,*)" in accall.f file (line 255) in order to compile it properly.

Commands:
    $   cd <PATH_TO_AAS3D-RF>/Tools/Naccess
    $   cp <PATH_TO_AAS3D-RF>/Tools/Naccess_Para_File/accall.pars  ./
    $   csh ./install.scr


If you have already had Naccess package at a different location, we recommand to recompile it using the modified accall.pars file in 
    <PATH_TO_AAS3D-RF>/Tools/Naccess_Para_File/ and you need to set "Naccess" variable in <PATH_TO_AAS3D-RF>/PredMod/userConfig.py
    file.
  

*** FoldX package
Visit http://foldxsuite.crg.eu/ to download the FoldX package
Move it to <PATH_TO_AAS3D-RF>/Tools/
Unpack it
Create a directory named "FoldX" in <PATH_TO_AAS3D-RF>/Tools/:
    $   mkdir FoldX
Put the "foldx" executable and also the "rotabase.txt" file into <PATH_TO_AAS3D-RF>/Tools/FoldX/ 


If you have already had the FoldX package in a different location and want to use it, you need to set "FoldX" variable
in <PATH_TO_AAS3D-RF>/PredMod/userConfig.py.



II. AAS3D-RF USAGE

For predicting disease-associated amino acid substitutions, you need to:

1) Create or choose a work path (your work path);

2) Integrate a tab-delimited input file containing 4 columns, representing UniProt_AC, Structure_ID, Chain, and AAS;
   An example file named AAS_example is in <PATH_TO_AAS3D-RF>/Example;

3) Download UniProt_AC.xml file from UniProt database for every UniProt_AC in your input file 
   (e.g. https://www.uniprot.org/uniprot/P20933.xml), and put all of them in the folder named UniProt_XML;

4) Provide a folder named STRUCTURES containing all the Structure_ID.pdb file in your input file;

5) And then run AAS3D-RF:
    $   python <PATH_TO_AAS3D-RF>/run_AAS3D-RF.py --in input_file  --out output_file --workPath work_path --num 16
    
    --in        Input amino acid subsitution file
    --out       Output file name
    --workPath  Working Path (please use absolute path !!!)
    --num       Number of threads (CPU cores) to use, default=1


Note: 1) Each structure file contains only one chain without expression tags, insertion codes, and non-standard residues
      2) The sequence fragment in structure file must exactly match the sequence from UniProt_AC.xml file with consistent 
         residue numbers
      3) Please put tab-delimited input file, UniProt_XML folder, and STRUCTURES folder in your work path
      4) Please run AAS3D-RF in your work path
      5) The output file will be in your work path 



Before starting your own work, firstly, we recommend you to run example given in <PATH_TO_AAS3D-RF>/Example to ensure 
all softwares have been installed properly and the AAS3D-RF works normally.

Example:
    $   cd <PATH_TO_AAS3D-RF>/Example
    $   python <PATH_TO_AAS3D-RF>/run_AAS3D-RF.py --in AAS_example --out example_out --workPath <PATH_TO_AAS3D-RF>/Example --num 32



III. LICENSE

AAS3D-RF is free to academic users. AAS3D-RF relies on several prerequisite packages, whose copyrights
belong to their authors or institutions.


IV. CITATION

Manuscript is under review.
