//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   H5ToRoot/main.cpp
 * \author Stefano Tognini
 * \brief  Reads a tally group from a .h5 file and produces a .root file.
 *         3 examples are presented.
 * \note   Copyright (c) 2019 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//


#include "TString.h"
#include "Nemesis/comm/global.hh"
#include "H5ToRoot.hh"


int main(int argc, char *argv[])
{
    
    nemesis::initialize(argc, argv);
    
    // Printing help
    if (argc < 3 || argc > 5)
    {
        std::cout << std::endl;
        std::cout << "This program presents 3 examples on how the class ";
        std::cout << "H5ToROOT can be implemented." << std::endl;
        std::cout << std::endl;
        std::cout << "  1. Creates a root file with the same name as the h5, ";
        std::cout << "and adds one tally group:" << std::endl;
        std::cout << "     ./h5ToRoot [h5File] [tallyGroup]" << std::endl;
        std::cout << "  2. Creates a root file with a user-specified filename, ";
        std::cout << "and adds one tally group:" << std::endl;
        std::cout << "     ./h5ToRoot [h5File] [tallyGroup] [rootFilename]";
        std::cout << std::endl;
        std::cout << "  3. Creates a root file with a user-specified filename ";
        std::cout << "and adds 2 tally groups, each in a different tree:";
        std::cout << std::endl;
        std::cout << "     ./h5ToRoot [h5File] [tallyGroup] [tallyGroup2] ";
        std::cout << "[rootFilename]" << std::endl;
        std::cout << std::endl;
        
        return 0;
    }
    
    
    //-----------------------------------------------------------------------//
    /*!
     * EXAMPLE 1:
     * ----------
     *
     * Creates a root file with the same name as the h5, and adds one tally
     * group:
     *
     * $ ./h5ToRoot [h5 file] [tally group]
     */
    //-----------------------------------------------------------------------//
    if (argc == 3)
    {
        // Fetching the HDF5 filename
        TString h5InputFilename = argv[1];
        
        // Fetching the tally group name
        TString tallyGroupName = argv[2];
        
        // Creating the root file with the tally tree
        h5toroot::H5ToRoot rootOutput(h5InputFilename, tallyGroupName);
        
        // Selecting all possible datasets
        rootOutput.SelectAllDatasets();
        
        // Adding existing datasets to the root tree. Skipping non-existent ones
        rootOutput.AddBranches();
        
        // Filling the root tree, writing the file to disk, and closing it
        rootOutput.End();
        
        return 0;
    }
    
    
    //-----------------------------------------------------------------------//
    /*!
     * EXAMPLE 2:
     * ----------
     *
     * Creates a root file with a user-specified filename, and adds one tally
     * group:
     *
     * $ ./h5ToRoot [h5 file] [tally group] [root filename]
     */
    //-----------------------------------------------------------------------//
    else if (argc == 4)
    {
        // Fetching the HDF5 filename
        TString h5InputFilename = argv[1];
        
        // Fetching the tally group name
        TString tallyGroupName = argv[2];
        
        // Fetching the root filename
        TString rootFilename = argv[3];
        
        // Creating the root file with the tally tree
        h5toroot::H5ToRoot rootOutput(h5InputFilename,
                                      tallyGroupName,
                                      rootFilename);
        
        // Selecting all possible datasets
        rootOutput.SelectAllDatasets();
        
        // Adding existing datasets to the root tree. Skipping non-existent ones
        rootOutput.AddBranches();
        
        // Filling the root tree, writing the file to disk, and closing it
        rootOutput.End();

        return 0;
    }
    
    
    //-----------------------------------------------------------------------//
    /*!
     * EXAMPLE 3:
     * ----------
     *
     * Creates a root file with a user-specified filename and adds 2 tally
     * groups, each in a different tree:
     *
     * $ ./h5ToRoot [h5 file] [tally group] [tally group 2] [root filename]
     */
    //-----------------------------------------------------------------------//
    else if (argc == 5)
    {
        // Fetching the HDF5 filename
        TString h5InputFilename = argv[1];
        
        // Fetching the first tally group name
        TString tallyGroupName = argv[2];

        // Fetching the second tally group name
        TString tallyGroupName2 = argv[3];
        
        // Fetching the root filename
        TString rootFilename = argv[4];

        // Creating the root file with the tally tree
        h5toroot::H5ToRoot rootOutput(h5InputFilename,
                                      tallyGroupName,
                                      rootFilename);
        
        // Selecting all possible datasets
        rootOutput.SelectAllDatasets();
        
        // Adding existing datasets to the root tree. Skipping non-existent ones
        rootOutput.AddBranches();
        
        // Filling the root tree, writing the file to disk, and closing it
        rootOutput.End();
        
        // Reopening the previous root file
        std::shared_ptr<TFile> rootFile(new TFile(rootFilename.Data(), "update"));
                       
        // Adding the second branch to the reopened root file
        h5toroot::H5ToRoot rootOutput2(h5InputFilename,
                                       tallyGroupName2,
                                       rootFile);
        
        // Selecting all possible datasets
        rootOutput2.SelectAllDatasets();
        
        // Adding existing datasets to the root tree. Skipping non-existent ones
        rootOutput2.AddBranches();
        
        // Filling the root tree, writing the file to disk, and closing it=
        rootOutput2.End();
        
        return 0;
    }
}
