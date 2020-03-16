//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   include/H5ToRoot.hh
 * \author Stefano Tognini
 * \brief  Class description for reading a tally group from a .h5 file and
 *         producing a .root file.
 * \note   Copyright (c) 2019 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//


#ifndef h5ToRoot_h
#define h5ToRoot_h


// C++
#include <string>
#include <vector>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

// Nemesis
#include "Nemesis/comm/global.hh"

#include "Nemesis/harness/HDF5_Diagnostics.hh"

#include "Nemesis/hdf5/core/File.hh"
#include "Nemesis/hdf5/core/Reader.hh"
#include "Nemesis/hdf5/core/Definitions.hh"
#include "Nemesis/hdf5/io/Arithmetic.hh"
#include "Nemesis/hdf5/io/Vector.hh"
#include "Nemesis/hdf5/io/View_Field.hh"
#include "Nemesis/hdf5/io/Vector_Bool.hh"
#include "Nemesis/hdf5/io/Bool.hh"
#include "Nemesis/hdf5/io/String.hh"
#include "Nemesis/hdf5/io/Vector_String.hh"
#include "Nemesis/hdf5/io/Hyperslab_Vector.hh"
#include "Nemesis/hdf5/dt/String_Datatype.hh"
#include "Nemesis/hdf5/dt/Vlen_Datatype.hh"

#include "Nemesis/utils/Range.hh"



namespace h5toroot
{

//===========================================================================//
/*!
 * \brief Class implementation for reading a tally group from a .h5 file and
 *        producing a .root file.
 */
//===========================================================================//
class H5ToRoot
{
    
  //--------------------------- Public class members ------------------------//
  public:
    TString                h5InputFilename;
    TString                rootFilename;
    TString                tallyGroupName;
    std::shared_ptr<TFile> rootFile;
    std::shared_ptr<TTree> tallyTree;
    std::vector<TString>   datasetList;
    
    
    //----------------------------- Constructors ----------------------------//
    /*!
     * Constructor for the case of a user-specified root filename.
     */
    H5ToRoot(TString   h5InputFilename,
             TString   tallyGroupName,
             TString   rootFilename);
    
    /*!
     * Constructor for the case when the user wants the .root filename to be
     * the same as the .h5 file. (Apart from the file extension.)
     */
    H5ToRoot(TString   h5InputFilename,
             TString   tallyGroupName);

    /*!
     * Constructor for when user wants to point to an existing root file.
     * This adds the possibility to write several tallies to a single root file.
     */
    H5ToRoot(TString   h5InputFilename,
             TString   tallyGroupName,
             std::shared_ptr<TFile> rootFile);

    
    
    //----------------------------- Destructor ------------------------------//
    ~H5ToRoot();
    
    
    //------------------------- Public class methods ------------------------//
    /*!
     * Assigns to the class public variable TString rootFilename the same name
     * of the h5 input filename but replaces the .h5 or .hdf5 extension by the
     * .root extension.
     */
    TString CreateRootFilename(TString &h5InputFilename);
    
    /*!
     * Adds a user-specified dataset to dataset list. This list is written to
     * the ROOT file by invoking the H5ToRoot::AddBranches() method.
     */
    void SelectDataset(TString dataset);
    
    /*!
     * Removes a user-specified dataset to be written to the ROOT file.
     */
    void RemoveDataset(TString dataset);
    
    /*!
     * Prints the current selections of datasets.
     */
    void PrintSelectedDatasets();
    
    /*!
     * Selects all possible datasets existent in any tally group. This list
     * is written into the ROOT file by invoking H5ToRoot::AddBranches().
     * Any dataset included in the list that is not found in the .h5 file is
     * automatically skipped.
     */
    void SelectAllDatasets();

    /*!
     * Fills the ROOT Tree, writes the ROOT file to disk, and closes the
     * TFile object. Both the TFile and TTree objects are public members
     * of the class and, thus, can be manipulated at will.
     */
    void End();
    
    /*!
     * Adds a single dataset to the TTree as a branch variable. Mind that
     * if the dataset does not exists in the .h5 file, this function will
     * simply skip it and will not return an error.
     */
    void AddBranch(TString &dataset);
    
    /*!
     * Add all datasets specified in an external std::vector<TString> as
     * branches in the ROOT TTree.
     */
    void AddBranches(std::vector<TString> &datasetList);
    
    /*!
     * Add all datasets specified in the public class member
     * std::vector<TString> datasetList.
     */
    void AddBranches();
    
    
    
  //-------------------------- Private class members ------------------------//
  private:
    std::string                          type;
    std::vector<std::string>             dimension_labels_binned;
    nemesis::Hyperslab_Vector<double, 4> hyperslabBinned_4;
    nemesis::Hyperslab_Vector<double, 6> hyperslabBinned_6;
    std::vector<unsigned long>           binnedDims;
    std::string                          description;
    std::vector<double>                  group_bounds_n;
    int                                  max_encountered_bins;
    std::vector<std::string>             mesh_cell;
    std::vector<std::string>             mesh_stat;
    std::vector<double>                  mesh_x;
    std::vector<double>                  mesh_y;
    std::vector<double>                  mesh_z;
    std::vector<std::string>             multiplier_descs;
    std::vector<std::string>             multiplier_names;
    double                               normalization;
    double                               num_histories;
    std::vector<std::string>             dimension_labels_total;
    nemesis::Hyperslab_Vector<double, 3> hyperslabTotal_3;
    nemesis::Hyperslab_Vector<double, 5> hyperslabTotal_5;
    std::vector<unsigned long>           totalDims;
    std::vector<unsigned int>            union_cellids;
    std::vector<unsigned int>            union_lengths;
    std::vector<std::string>             dimension_labels_volumes;
    nemesis::Hyperslab_Vector<double, 1> hyperslabVolumes_1;
    nemesis::Hyperslab_Vector<double, 2> hyperslabVolumes_2;
    nemesis::Hyperslab_Vector<double, 3> hyperslabVolumes_3;
    std::vector<unsigned long>           volumeDims;
    
    
    //------------------------ Private class methods ------------------------//
    /*!
     * Private class member functions for reading individual datasets and
     * attributes.
     *
     * To add a new function, three steps need to be taken to ensure that the
     * class will work as intended:
     *
     * 1) Add a new function
     *    AddBranch_datasetName(nemesis::hdf5::Reader &reader) { ... }.
     * 2) Include a new if statement calling the aforementioned function in
     *    WhichDatasetToAdd(...).
     * 3) Include the new datasetName in public member function
     *    SelectAllDatasets().
     */
    
    
    inline void AddBranch_type(nemesis::hdf5::Reader &reader);
    
    inline void AddBranch_binned(nemesis::hdf5::Reader &reader);
    
    inline void AddBranch_description(nemesis::hdf5::Reader &reader);
    
    inline void AddBranch_group_bounds_n(nemesis::hdf5::Reader &reader);
    
    inline void AddBranch_max_encountered_bins(nemesis::hdf5::Reader &reader);
    
    inline void AddBranch_mesh_stat(nemesis::hdf5::Reader &reader);
    
    inline void AddBranch_mesh_cell(nemesis::hdf5::Reader &reader);
    
    inline void AddBranch_mesh_x(nemesis::hdf5::Reader &reader);
    
    inline void AddBranch_mesh_y(nemesis::hdf5::Reader &reader);
    
    inline void AddBranch_mesh_z(nemesis::hdf5::Reader &reader);
    
    inline void AddBranch_multiplier_descs(nemesis::hdf5::Reader &reader);
    
    inline void AddBranch_multiplier_names(nemesis::hdf5::Reader &reader);
    
    inline void AddBranch_normalization(nemesis::hdf5::Reader &reader);
    
    inline void AddBranch_num_histories(nemesis::hdf5::Reader &reader);
    
    inline void AddBranch_total(nemesis::hdf5::Reader &reader);
    
    inline void AddBranch_union_cellids(nemesis::hdf5::Reader &reader);
    
    inline void AddBranch_union_lenghts(nemesis::hdf5::Reader &reader);
    
    inline void AddBranch_volumes(nemesis::hdf5::Reader &reader);

    inline void WhichDatasetToAdd(nemesis::hdf5::Reader &reader,
                                  TString &dataset);
    
}; /* end H5ToRoot class */

} /* end namespace h5toroot */

#endif
