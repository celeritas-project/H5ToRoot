//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   src/H5ToRoot.cpp
 * \author Stefano Tognini
 * \brief  Class implementation for reading a tally group from a .h5 file and
 *         producing a .root file.
 * \note   Copyright (c) 2019 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//



#include "H5ToRoot.hh"



//------------------------------ Constructors -------------------------------//
h5toroot::H5ToRoot::H5ToRoot(TString h5InputFilename,
                             TString tallyGroupName,
                             TString rootFilename)
: h5InputFilename(h5InputFilename), tallyGroupName(tallyGroupName),
  rootFilename(rootFilename)
{
    this->rootFile = std::make_shared<TFile>(rootFilename.Data(), "RECREATE");
    this->tallyTree = std::make_shared<TTree>(tallyGroupName.Data(),
                                              tallyGroupName.Data());
}


h5toroot::H5ToRoot::H5ToRoot(TString h5InputFilename,
                             TString tallyGroupName)
: h5InputFilename(h5InputFilename), tallyGroupName(tallyGroupName)
{
    this->rootFilename = CreateRootFilename(h5InputFilename);
    this->rootFile = std::make_shared<TFile>(rootFilename.Data(), "RECREATE");
    this->tallyTree = std::make_shared<TTree>(tallyGroupName.Data(),
                                              tallyGroupName.Data());
}


h5toroot::H5ToRoot::H5ToRoot(TString   h5InputFilename,
                             TString   tallyGroupName,
                             std::shared_ptr<TFile> rootFile)
: h5InputFilename(h5InputFilename),  tallyGroupName(tallyGroupName),
  rootFile(rootFile)
{
    this->tallyTree = std::make_shared<TTree>(tallyGroupName.Data(),
                                              tallyGroupName.Data());
}



//------------------------------ Destructor ---------------------------------//
h5toroot::H5ToRoot::~H5ToRoot()
{
}



//-------------------------- Public class methods ---------------------------//
TString h5toroot::H5ToRoot::CreateRootFilename(TString &h5InputFilename)
{
    // the filename is the same as the hdf5, except for the extension
    TString rootOutputFilename = h5InputFilename;
    
    // Adding the ROOT filename extension
    if (rootOutputFilename.Contains(".h5"))
    {
        rootOutputFilename.ReplaceAll(".h5", ".root");
    }
    
    if (rootOutputFilename.Contains(".hdf5"))
    {
        rootOutputFilename.ReplaceAll(".hdf5", ".root");
    }
    
    return rootOutputFilename;
}

void h5toroot::H5ToRoot::SelectDataset(TString dataset)
{
    this->datasetList.push_back(dataset);
}

void h5toroot::H5ToRoot::RemoveDataset(TString dataset)
{
    for (int i = 0; i < this->datasetList.size(); i++)
    {
        TString thisDataset = this->datasetList.at(i);
        
        if (thisDataset == dataset)
        {
            datasetList.erase(this->datasetList.begin() + i);
        }
    }
}

void h5toroot::H5ToRoot::PrintSelectedDatasets()
{
    for (auto dataset : this->datasetList)
    {
        std::cout << dataset << std::endl;
    }
}

void h5toroot::H5ToRoot::SelectAllDatasets()
{
    this->datasetList.push_back("type");
    this->datasetList.push_back("binned");
    this->datasetList.push_back("description");
    this->datasetList.push_back("group_bounds_n");
    this->datasetList.push_back("max_encountered_bins");
    this->datasetList.push_back("mesh_cell");
    this->datasetList.push_back("mesh_stat");
    this->datasetList.push_back("mesh_x");
    this->datasetList.push_back("mesh_y");
    this->datasetList.push_back("mesh_z");
    this->datasetList.push_back("multiplier_descs");
    this->datasetList.push_back("multiplier_names");
    this->datasetList.push_back("normalization");
    this->datasetList.push_back("num_histories");
    this->datasetList.push_back("total");
    this->datasetList.push_back("union_cellids");
    this->datasetList.push_back("union_lengths");
    this->datasetList.push_back("volumes");
}

void h5toroot::H5ToRoot::End()
{
    this->tallyTree->Fill();
    this->rootFile->Write();
    this->rootFile->Close();
}

//---------------------------------------------------------------------------//
//                       AddBranch() / AddBranches()
//---------------------------------------------------------------------------//
void h5toroot::H5ToRoot::AddBranch(TString &dataset)
{
    // Opening the HDF5 file
    auto h5File = std::make_shared<nemesis::hdf5::File>
    (h5InputFilename.Data(),
     nemesis::hdf5::File::READ,
     nemesis::hdf5::File::SERIAL);
    
    // Creating the reader
    nemesis::hdf5::Reader reader(h5File);
    
    // Opening tally group
    reader.open_group("tally");
    reader.open_group(tallyGroupName.Data());
    
    // Updating the TFile opening option
    rootFile->ReOpen("UPDATE");
    tallyTree->GetEntry(0);
    
    // Special case. Type is an attribute, so the other if statement fails
    if (dataset == "type")
    {
        WhichDatasetToAdd(reader, dataset);
    }
    
    if (!reader.location().exists(dataset.Data()))
    {
        reader.close_group();
        reader.close_group();
        return;
    }
    
    WhichDatasetToAdd(reader, dataset);
    
    reader.close_group();
    reader.close_group();
}

void h5toroot::H5ToRoot::AddBranches(std::vector<TString> &datasetList)
{
    // Opening the HDF5 file
    auto h5File = std::make_shared<nemesis::hdf5::File>
    (h5InputFilename.Data(),
     nemesis::hdf5::File::READ,
     nemesis::hdf5::File::SERIAL);
    
    // Creating the reader
    nemesis::hdf5::Reader reader(h5File);
    
    // Opening tally group
    reader.open_group("tally");
    reader.open_group(tallyGroupName.Data());
    
    // Updating the TFile opening option
    rootFile->ReOpen("UPDATE");
    tallyTree->GetEntry(0);
    
    for (auto dataset : datasetList)
    {
        // Special case. Type is an attribute, so the other if statement fails
        if (dataset == "type")
        {
            WhichDatasetToAdd(reader, dataset);
        }
        
        if (reader.location().exists(dataset.Data()))
        {
            WhichDatasetToAdd(reader, dataset);
        }
    }
    
    reader.close_group();
    reader.close_group();
}

void h5toroot::H5ToRoot::AddBranches()
{
    // Opening the HDF5 file
    auto h5File = std::make_shared<nemesis::hdf5::File>
    (h5InputFilename.Data(),
     nemesis::hdf5::File::READ,
     nemesis::hdf5::File::SERIAL);
    
    // Creating the reader
    nemesis::hdf5::Reader reader(h5File);
    
    // Opening tally group
    reader.open_group("tally");
    reader.open_group(tallyGroupName.Data());
    
    // Updating the TFile opening option
    rootFile->ReOpen("UPDATE");
    tallyTree->GetEntry(0);
    
    for (auto dataset : this->datasetList)
    {
        // Special case. Type is an attribute, so the other if statement fails
        if (dataset == "type")
        {
            WhichDatasetToAdd(reader, dataset);
        }
        
        if (reader.location().exists(dataset.Data()))
        {
            WhichDatasetToAdd(reader, dataset);
        }
    }
    
    reader.close_group();
    reader.close_group();
}



//------------------------- Private class methods ---------------------------//
inline void h5toroot::H5ToRoot::AddBranch_type(nemesis::hdf5::Reader &reader)
{
    reader.open_attr();
    read(reader, "type", type);
    reader.close_attr();
    
    TBranch *thisBranch = tallyTree->Branch("type", &type);
    thisBranch->Fill();
}

inline void h5toroot::H5ToRoot::AddBranch_binned(nemesis::hdf5::Reader &reader)
{
    reader.open_attr("binned");
    read(reader, "DIMENSION_LABELS", dimension_labels_binned);
    reader.close_attr();
    
    if (dimension_labels_binned.size() == 6)
    {
        
        nemesis::hdf5::read(reader,
                            nemesis::hdf5::Arithmetic_Datatype<double>(),
                            "binned",
                            hyperslabBinned_6);
        
        auto dim = hyperslabBinned_6.dims();
        
        for (int i = 0; i < dim.size(); i++)
            binnedDims.push_back(dim[i]);
        
        double binned[binnedDims.at(0)][binnedDims.at(1)]
        [binnedDims.at(2)][binnedDims.at(3)]
        [binnedDims.at(4)][binnedDims.at(5)];
        
        for (int i = 0; i < binnedDims[0]; i++)
            for (int j = 0; j < binnedDims[1]; j++)
                for (int k = 0; k < binnedDims[2]; k++)
                    for (int l = 0; l < binnedDims[3]; l++)
                        for (int m = 0; m < binnedDims[4]; m++)
                            for (int n = 0; n < binnedDims[5]; n++)
                                binned[i][j][k][l][m][n] =
                                hyperslabBinned_6[i][j][k][l][m][n];
        
        TBranch *thisBranch =
        tallyTree->Branch("binned",
                          binned,
                          Form("binned[%lu][%lu][%lu][%lu][%lu][%lu]/D",
                               binnedDims[0],
                               binnedDims[1],
                               binnedDims[2],
                               binnedDims[3],
                               binnedDims[4],
                               binnedDims[5])
                          );
        
        thisBranch->Fill();
        
    }
    
    else if (dimension_labels_binned.size() == 4)
    {
        static nemesis::Hyperslab_Vector<double, 4> hyperslabBinned_4;
        static std::vector<unsigned long> binnedDims;
        
        nemesis::hdf5::read(reader,
                            nemesis::hdf5::Arithmetic_Datatype<double>(),
                            "binned",
                            hyperslabBinned_4);
        
        auto dim = hyperslabBinned_4.dims();
        
        for (int i = 0; i < dim.size(); i++)
            binnedDims.push_back(dim[i]);
        
        double binned[binnedDims.at(0)][binnedDims.at(1)]
        [binnedDims.at(2)][binnedDims.at(3)];
        
        for (int i = 0; i < binnedDims[0]; i++)
            for (int j = 0; j < binnedDims[1]; j++)
                for (int k = 0; k < binnedDims[2]; k++)
                    for (int l = 0; l < binnedDims[3]; l++)
                        binned[i][j][k][l] =
                        hyperslabBinned_4[i][j][k][l];
        
        TBranch *thisBranch =
        tallyTree->Branch("binned",
                          binned,
                          Form("binned[%lu][%lu][%lu][%lu]/D",
                               binnedDims[0],
                               binnedDims[1],
                               binnedDims[2],
                               binnedDims[3])
                          );
        
        thisBranch->Fill();
    }
    
    TBranch *thisBranch = tallyTree->Branch("dimension_labels_binned",
                                            &dimension_labels_binned);
    thisBranch->Fill();
}

inline void h5toroot::H5ToRoot::AddBranch_description(nemesis::hdf5::Reader &reader)
{
    read(reader, "description", description);
    
    TBranch *thisBranch = tallyTree->Branch("description", &description);
    thisBranch->Fill();
}

inline void h5toroot::H5ToRoot::AddBranch_group_bounds_n(nemesis::hdf5::Reader
                                                         &reader)
{
    read(reader, "group_bounds_n", group_bounds_n);
    
    TBranch *thisBranch = tallyTree->Branch("group_bounds_n", &group_bounds_n);
    thisBranch->Fill();
}

inline void h5toroot::H5ToRoot::AddBranch_max_encountered_bins(nemesis::hdf5::Reader
                                                               &reader)
{
    read(reader, "max_encountered_bins", max_encountered_bins);
    
    TBranch *thisBranch = tallyTree->Branch("max_encountered_bins",
                                            &max_encountered_bins);
    thisBranch->Fill();
}

inline void h5toroot::H5ToRoot::AddBranch_mesh_cell(nemesis::hdf5::Reader &reader)
{
    read(reader, "mesh_cell", mesh_cell);
    
    TBranch *thisBranch = tallyTree->Branch("mesh_cell", &mesh_cell);
    thisBranch->Fill();
}

inline void h5toroot::H5ToRoot::AddBranch_mesh_stat(nemesis::hdf5::Reader &reader)
{
    read(reader, "mesh_stat", mesh_stat);
    
    TBranch *thisBranch = tallyTree->Branch("mesh_stat", &mesh_stat);
    thisBranch->Fill();
}

inline void h5toroot::H5ToRoot::AddBranch_mesh_x(nemesis::hdf5::Reader &reader)
{
    read(reader, "mesh_x", mesh_x);
    
    TBranch *thisBranch = tallyTree->Branch("mesh_x", &mesh_x);
    thisBranch->Fill();
}

inline void h5toroot::H5ToRoot::AddBranch_mesh_y(nemesis::hdf5::Reader &reader)
{
    read(reader, "mesh_y", mesh_y);
    
    TBranch *thisBranch = tallyTree->Branch("mesh_y", &mesh_y);
    thisBranch->Fill();
}

inline void h5toroot::H5ToRoot::AddBranch_mesh_z(nemesis::hdf5::Reader &reader)
{
    read(reader, "mesh_z", mesh_z);
    
    TBranch *thisBranch = tallyTree->Branch("mesh_z", &mesh_z);
    thisBranch->Fill();
}

inline void h5toroot::H5ToRoot::AddBranch_multiplier_descs(nemesis::hdf5::Reader
                                                           &reader)
{
    read(reader, "multiplier_descs", multiplier_descs);
    
    TBranch *thisBranch = tallyTree->Branch("multiplier_descs",
                                            &multiplier_descs);
    thisBranch->Fill();
}

inline void h5toroot::H5ToRoot::AddBranch_multiplier_names(nemesis::hdf5::Reader
                                                           &reader)
{
    read(reader, "multiplier_names", multiplier_names);
    
    TBranch *thisBranch = tallyTree->Branch("multiplier_names",
                                            &multiplier_names);
    thisBranch->Fill();
}

inline void h5toroot::H5ToRoot::AddBranch_normalization(nemesis::hdf5::Reader
                                                        &reader)
{
    read(reader, "normalization", normalization);
    
    TBranch *thisBranch = tallyTree->Branch("normalization",
                                            &normalization,
                                            "normalization/D");
    thisBranch->Fill();
}

inline void h5toroot::H5ToRoot::AddBranch_num_histories(nemesis::hdf5::Reader
                                                        &reader)
{
    read(reader, "num_histories", num_histories);
    
    TBranch *thisBranch = tallyTree->Branch("num_histories",
                                            &num_histories,
                                            "num_histories/D");
    thisBranch->Fill();
}

inline void h5toroot::H5ToRoot::AddBranch_total(nemesis::hdf5::Reader &reader)
{
    reader.open_attr("total");
    read(reader, "DIMENSION_LABELS", dimension_labels_total);
    reader.close_attr();
    
    if (dimension_labels_total.size() == 5)
    {
        nemesis::hdf5::read(reader,
                            nemesis::hdf5::Arithmetic_Datatype<double>(),
                            "total",
                            hyperslabTotal_5);
        
        auto dims = hyperslabTotal_5.dims();
        
        for (int i = 0; i < dims.size(); i++)
            totalDims.push_back(dims[i]);
        
        double total[totalDims.at(0)][totalDims.at(1)]
        [totalDims.at(2)][totalDims.at(3)][totalDims.at(4)];
        
        for (int i = 0; i < totalDims[0]; i++)
            for (int j = 0; j < totalDims[1]; j++)
                for (int k = 0; k < totalDims[2]; k++)
                    for (int l = 0; l < totalDims[3]; l++)
                        for (int m = 0; m < totalDims[4]; m++)
                            total[i][j][k][l][m] =
                            hyperslabTotal_5[i][j][k][l][m];
        
        TBranch *thisBranch =
        tallyTree->Branch("total",
                          total,
                          Form("total[%lu][%lu][%lu][%lu][%lu]/D",
                               totalDims[0],
                               totalDims[1],
                               totalDims[2],
                               totalDims[3],
                               totalDims[4])
                          );
        
        thisBranch->Fill();
        
    }
    
    else if (dimension_labels_total.size() == 3)
    {
        nemesis::hdf5::read(reader,
                            nemesis::hdf5::Arithmetic_Datatype<double>(),
                            "total",
                            hyperslabTotal_3);
        
        auto dims = hyperslabTotal_3.dims();
        
        for (int i = 0; i < dims.size(); i++)
            totalDims.push_back(dims[i]);
        
        double total[totalDims.at(0)][totalDims.at(1)][totalDims.at(2)];
        
        for (int i = 0; i < totalDims[0]; i++)
            for (int j = 0; j < totalDims[1]; j++)
                for (int k = 0; k < totalDims[2]; k++)
                    total[i][j][k] = hyperslabTotal_3[i][j][k];
        
        TBranch *thisBranch =
        tallyTree->Branch("total",
                          total,
                          Form("total[%lu][%lu][%lu]/D",
                               totalDims[0],
                               totalDims[1],
                               totalDims[2])
                          );
        
        thisBranch->Fill();
    }
    
    TBranch *thisBranch = tallyTree->Branch("dimension_labels_total",
                                            &dimension_labels_total);
    thisBranch->Fill();
}

inline void h5toroot::H5ToRoot::AddBranch_union_cellids(nemesis::hdf5::Reader &reader)
{
    read(reader, "union_cellids", union_cellids);
    
    TBranch *thisBranch = tallyTree->Branch("union_cellids", &union_cellids);
    thisBranch->Fill();
}

inline void h5toroot::H5ToRoot::AddBranch_union_lenghts(nemesis::hdf5::Reader
                                                        &reader)
{
    read(reader, "union_lengths", union_lengths);
    
    TBranch *thisBranch = tallyTree->Branch("union_lengths", &union_lengths);
    thisBranch->Fill();
}

inline void h5toroot::H5ToRoot::AddBranch_volumes(nemesis::hdf5::Reader &reader)
{
    reader.open_attr("volumes");
    read(reader, "DIMENSION_LABELS", dimension_labels_volumes);
    reader.close_attr();
    
    if (dimension_labels_volumes.size() == 1)
    {
        nemesis::hdf5::read(reader,
                            nemesis::hdf5::Arithmetic_Datatype<double>(),
                            "volumes",
                            hyperslabVolumes_1);
        
        auto dims = hyperslabVolumes_1.dims();
        
        for (int i = 0; i < dims.size(); i++)
            volumeDims.push_back(dims[i]);
        
        double volumes[volumeDims.at(0)];
        
        for (int i = 0; i < volumeDims[0]; i++)
            volumes[i] = hyperslabVolumes_1[i];
        
        TBranch *thisBranch =
        tallyTree->Branch("volumes",
                          volumes,
                          Form("volumes[%lu]/D", volumeDims[0]));
        
        thisBranch->Fill();
    }
    
    else if (dimension_labels_volumes.size() == 2)
    {
        nemesis::hdf5::read(reader,
                            nemesis::hdf5::Arithmetic_Datatype<double>(),
                            "volumes",
                            hyperslabVolumes_2);
        
        auto dims = hyperslabVolumes_2.dims();
        
        for (int i = 0; i < dims.size(); i++)
            volumeDims.push_back(dims[i]);
        
        double volumes[volumeDims.at(0)][volumeDims.at(1)];
        
        for (int i = 0; i < volumeDims[0]; i++)
            for (int j = 0; j < volumeDims[1]; j++)
                volumes[i][j] = hyperslabVolumes_2[i][j];
        
        TBranch *thisBranch =
        tallyTree->Branch("volumes",
                          volumes,
                          Form("volumes[%lu][%lu]/D",
                               volumeDims[0],
                               volumeDims[1])
                          );
        
        thisBranch->Fill();
    }
    
    else if (dimension_labels_volumes.size() == 3)
    {
        nemesis::hdf5::read(reader,
                            nemesis::hdf5::Arithmetic_Datatype<double>(),
                            "volumes",
                            hyperslabVolumes_3);
        
        auto dims = hyperslabVolumes_3.dims();
        
        for (int i = 0; i < dims.size(); i++)
            volumeDims.push_back(dims[i]);
        
        double volumes[volumeDims.at(0)][volumeDims.at(1)][volumeDims.at(2)];
        
        for (int i = 0; i < volumeDims[0]; i++)
            for (int j = 0; j < volumeDims[1]; j++)
                for (int k = 0; k < volumeDims[2]; k++)
                    volumes[i][j][k] = hyperslabVolumes_3[i][j][k];
        
        TBranch *thisBranch =
        tallyTree->Branch("volumes",
                          volumes,
                          Form("volumes[%lu][%lu][%lu]/D",
                               volumeDims[0],
                               volumeDims[1],
                               volumeDims[2])
                          );
        
        thisBranch->Fill();
    }
    
    TBranch *thisBranch = tallyTree->Branch("dimension_labels_volumes",
                                            &dimension_labels_volumes);
    thisBranch->Fill();
}



//---------------------------------------------------------------------------//
//                          WhichDatasetToAdd()
//---------------------------------------------------------------------------//
inline void h5toroot::H5ToRoot::WhichDatasetToAdd(nemesis::hdf5::Reader &reader,
                                                  TString &dataset)
{
    if (dataset == "type")
    {
        AddBranch_type(reader);
    }
    
    if (dataset == "binned")
    {
        AddBranch_binned(reader);
    }
    
    if (dataset == "description")
    {
        AddBranch_description(reader);
    }
    
    if (dataset == "group_bounds_n")
    {
        AddBranch_group_bounds_n(reader);
    }
    
    if (dataset == "max_encountered_bins")
    {
        AddBranch_max_encountered_bins(reader);
    }
    
    if (dataset == "mesh_cell")
    {
        AddBranch_mesh_cell(reader);
    }
    
    if (dataset == "mesh_stat")
    {
        AddBranch_mesh_stat(reader);
    }
    
    if (dataset == "mesh_x")
    {
        AddBranch_mesh_x(reader);
    }
    
    if (dataset == "mesh_y")
    {
        AddBranch_mesh_y(reader);
    }
    
    if (dataset == "mesh_z")
    {
        AddBranch_mesh_z(reader);
    }
    
    if (dataset == "multiplier_descs")
    {
        AddBranch_multiplier_descs(reader);
    }
    
    if (dataset == "multiplier_names")
    {
        AddBranch_multiplier_names(reader);
    }
    
    if (dataset == "normalization")
    {
        AddBranch_normalization(reader);
    }
    
    if (dataset == "num_histories")
    {
        AddBranch_num_histories(reader);
    }
    
    if (dataset == "total")
    {
        AddBranch_total(reader);
    }
    
    if (dataset == "union_cellids")
    {
        AddBranch_union_cellids(reader);
    }
    
    if (dataset == "union_lengths")
    {
        AddBranch_union_lenghts(reader);
    }
    
    if (dataset == "volumes")
    {
        AddBranch_volumes(reader);
    }
}
