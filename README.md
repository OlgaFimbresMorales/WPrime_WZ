### Filtering, Mixing of datasets (nanoAOD) and converting from ROOT to HD5 (arrays) format for use in ML training

1. Log into LXPLUS server


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

4. Set up the proxy (to use samples stored in the GRID):

If you dont have a valid grid certificate follow the instructions here:

[Grid certificate](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookStartingGrid#ObtainingCert)


Then try the following command

```bash
voms-proxy-init -voms cms -valid 192:00
```




5. Test the filter script locally

```bash
reduce_nanoaod.py  root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/WprimeToWZToWlepZlep_narrow_M1000_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/270000/6AB66CC6-59BA-FF43-89AA-FF10B6E170C7.root output2.root
```


6. Samples to merge are located in datasets directory, use mergeSamples script to merge into single root files

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

7. Split into training and testing samples (e.g. separate from 10 files, 7 for training and the rest for test)

```bash
export TRAINDIR=${MERGEDIR}/train
export TESTDIR=${MERGEDIR}/test
mkdir -p $TRAINDIR $TESTDIR
mv ${MERGEDIR}/ntuple_merged_[.0-7.].root ${TRAINDIR}/
mv ${MERGEDIR}/ntuple_merged_*.root ${TESTDIR}/
```


### Convert `ROOT` files to `HDF5` files using `uproot`  (up to maximum 100 particles) 


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

### Ensure to have the requied packages

- onnxruntime: For running the ONNX model.
- uproot: To read NanoAOD files in pure Python.
- numpy: For handling input arrays.

```bash
pip install onnxruntime uproot numpy
```

### Execute script 

```bash
python test_onnx_nanoaod.py
```

The functionalities are: 
- Opens a NanoAOD ROOT file using uproot.
- Extracts electron kinematic variables (pt, eta, phi, mass).
- Loops over each event and selects the first electron.
- Formats the data as an input array for the ONNX model.
- Runs the ONNX model on each event and prints the output.



