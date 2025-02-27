# Data Processing, Training and Inference famework 


## Basic Setup 

1. Log into LXPLUS server (CERN computers)

```bash
ssh username@lxplus.cern.ch
```

2. Set a recent CMSSW version 

```bash
cmsrel CMSSW_13_3_0
cd CMSSW_13_3_0/src
cmsenv
```

3. Clone the repository  and compile 

```bash
git clone https://github.com/castaned/ML-integration-CMSSW DeepNTuples
scram b -j 4
```


## DATA Processing

### Filter nanoAOD orignal files


1. Set up the proxy (to use samples stored in the GRID):

If you dont have a valid grid certificate follow the instructions here:

[Grid certificate](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookStartingGrid#ObtainingCert)


Then generate the grid certificate file and store in .globus directory (this is needed to access files stored in grid servers)

```bash
voms-proxy-init --voms cms --valid 192:00 --out $HOME/.globus/x509up_u$(id -u)
```

To check that the certificate was generate correctly and that the file was stored in the .globus directory:

```bash
voms-proxy-info --all
```

2. Move to the directory to submit jobs to process nanoAOD orginal files

```bash
cd  cd MyNanoAODTools/scripts/
```

3.

- Check that the datasets to be processed are in the datasets.yaml file, the format should be consistent with the one found in the DAS ( https://cmsweb.cern.ch/das/)  (e.g. /WprimeToWZToWlepZlep_narrow_M1000_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM)

- Check that the list of branches to keep are updated in the branchsel.txt (according to the nanoAOD version the name of the branches may change, so it is always better to check the list of branches in the original file  https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/home)


4.

- update submit_condor.py   to replace the location where the files will be saved (e.g. change /eos/user/c/castaned by /eos/user/u/username)
- update run_filter.sh or.py

   - to replace the location where the code is located (e.g. change /afs/cern.ch/work/c/castaned/CMSSW_13_3_0/src by /afs/cern.ch/user/u/username)
   - to replace the location where files will be saved (e.g. change EOS_DIR="/eos/user/c/castaned/NanoAOD_Filtered/${DATASET_FOLDER}"  by  EOS_DIR="/eos/user/u/username/NanoAOD_Filtered/${DATASET_FOLDER}" )


5. Submit the condor jobs

```bash
condor_submit condor_submit.jdl
```

6. Check the progress of the jobs 


```bash
condor_q
```

7. After the jobs finished you should look at the EOS directory to verify the skimmed samples were created 


### Merge directories (randomly) and produce h5 files

1. Samples to merge are located in datasets directory, use mergeSamples script to merge into single root files

```bash
mergeSamples.py [events per output file] [output dir] [path to the filelist produced in step 1]
```
e.g.,
```bash
cd DeepNTuples
export OUTDIR=$PWD/datasets 
export MERGEDIR=$PWD/output
mergeSamples.py 200000 ${MERGEDIR} ${OUTDIR}/signal.txt ${OUTDIR}/bkg.txt
```

2. Split into training and testing samples (e.g. separate from 10 files, 7 for training and the rest for test)

```bash
export TRAINDIR=${MERGEDIR}/train
export TESTDIR=${MERGEDIR}/test
mkdir -p $TRAINDIR $TESTDIR
mv ${MERGEDIR}/ntuple_merged_[.0-7.].root ${TRAINDIR}/
mv ${MERGEDIR}/ntuple_merged_*.root ${TESTDIR}/
```


Then you can run


```bash
convert-uproot-opendata.py [input file (.root)] [output file (.h5)]
```
e.g.,
```
convert-uproot-opendata.py ${TRAINDIR}/ntuple_merged_5.root ${TRAINDIR}/ntuple_merged_5.h5
```
which produces `HDF5` files with different arrays for each output variable.



## Training





## Inference

 Ensure to have the requied packages

- onnxruntime: For running the ONNX model.
- uproot: To read NanoAOD files in pure Python.
- numpy: For handling input arrays.

```bash
pip install onnxruntime uproot numpy
```

 Execute script 

```bash
python test_onnx_nanoaod.py
```

The functionalities are: 
- Opens a NanoAOD ROOT file using uproot.
- Extracts electron kinematic variables (pt, eta, phi, mass).
- Loops over each event and selects the first electron.
- Formats the data as an input array for the ONNX model.
- Runs the ONNX model on each event and prints the output.



