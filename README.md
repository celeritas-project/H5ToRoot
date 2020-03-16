README
======


# About

The `H5ToRoot` class takes a tally group from a Shift HDF5 input file and creates
a ROOT file containing a TTree with the information from the tally group.

An executable with 3 different examples on how to use the class is included in
the package.



# Dependencies

- ROOT
- HDF5
- SCALE
  - Nemesis
- MPI



# How to install

From the H5ToRoot folder:

```
$ ./build.sh
```



# Example

Check the `main.cpp` and run the `./h5ToRoot` example.



# Class Implementation


## Public variables
```
TString                h5InputFilename;   
TString                rootFilename;   
TString                tallyGroupName;   
std::shared_ptr<TFile> rootFile;   
std::shared_ptr<TTree> tallyTree;   
std::vector<TString>   datasetList;   
```


## Constructors

`H5ToRoot(TString   h5InputFilename,
          TString   tallyGroupName);`   
Case the user wants the .root filename to be the same as the .h5 file
(apart from the file extension).


`H5ToRoot(TString   h5InputFilename,
          TString   tallyGroupName,
          TString   rootFilename);`   
Case of a user-specified root filename.


`H5ToRoot(TString   h5InputFilename,
          TString   tallyGroupName,
          std::shared_ptr<TFile> rootFile);`   
Case the user wants to point to an existing root file. This adds the possibility
to write several tallies to a single root file.


## Public class methods

`TString CreateRootFilename(TString &h5InputFilename);`   
Assigns to the class public variable TString rootFilename the same name of the
h5 input filename but replaces the .h5 or .hdf5 extension by the .root extension.


`void SelectDataset(TString dataset);`   
Adds a user-specified dataset to the dataset list. This list is written to the
ROOT file by invoking the H5ToRoot::AddBranches() method.


`void RemoveDataset(TString dataset);`   
Removes a user-specified dataset to be written to the ROOT file.


`void PrintSelectedDatasets();`   
Prints the current selections of datasets.   


`void SelectAllDatasets();`   
Selects all possible datasets existent in any tally group. This list is written
into the ROOT file by invoking H5ToRoot::AddBranches().   
Any dataset included in the list that is not found in the .h5 file is
automatically skipped.


`void AddBranch(TString &dataset);`   
Adds a single dataset to the TTree as a branch variable. Mind that if the
dataset does not exists in the .h5 file, this function will simply skip it and
will not return an error.


`void AddBranches(std::vector<TString> &datasetList);`   
Add all datasets specified in an external std::vector<TString> as branches in
the ROOT TTree.


`void AddBranches();`   
Add all datasets specified in the public class member
std::vector<TString> datasetList.


`void End();`   
Fills the ROOT Tree, writes the ROOT file to disk, and closes the TFile object.
Both the TFile and TTree objects are public members of the class and, thus,
can be manipulated at will.



___
**Stefano Tognini**  
**Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.**
