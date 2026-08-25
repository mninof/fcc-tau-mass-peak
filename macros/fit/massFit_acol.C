///////////////////////////////////////////////////////////////////////////////////
///// Future Circular Collider Software - Tau-lepton Mass Peak Reconstruction /////
///// This macro performs the fit of the tau-lepton mass for several bins of  /////
///// the variable "acolinearity". Used to study cuts on this quantity        /////
///////////////////////////////////////////////////////////////////////////////////

#include "TTree.h"
#include "RooRealVar.h"
#include "RooDataSet.h"

#include "RooAddPdf.h"

#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"

#include <iostream>
#include <iomanip>

#include <sstream>

/* This macro can be called once the rootfile from the simulation has been loaded to
   ROOT. The argument, a TTree, corresponds to the name of the tree in the rootfile
   where the events are stored. Typically, this name is just "events".
  
   This macro requires a new branch, called "Acolinearity", which can be produced with
   the files and setting provided in the /further_analyses directory.
  
   The produced result is a text file containing all the information from the bins */
void massFit_acol(TTree* events)
{       
    // 1. Define the observable
    double centre = 1.7768199; double ext = 0.2;
    double xmin = centre - ext; double xmax = centre + ext;

    // 2. Define Acolinearity bins
    const int nBins = 1;
    double loEdge = -1.0;
    double hiEdge = -0.999905;
    double binWidth = (hiEdge - loEdge) / nBins;
    double edges[nBins + 1];
    for (int i = 0; i <= nBins; i++)
        edges[i] = loEdge + i * binWidth;

    // 3. Read all entries once into vectors
    TTreeReader reader(events);
    TTreeReaderValue<double> massReco(reader, "MassReco");
    TTreeReaderValue<double> massRecoCC(reader, "MassReco_CC");
    TTreeReaderValue<double> acolinearity(reader, "Acolinearity");

    // Store as pairs (mass, acolinearity)
    std::vector<std::pair<double,double>> entries;
    while (reader.Next()) {
        double acol = *acolinearity;
        if (*massReco > xmin && *massReco < xmax)
            entries.push_back({*massReco, acol});
        if (*massRecoCC > xmin && *massRecoCC < xmax)
            entries.push_back({*massRecoCC, acol});
    }
    
    std::ofstream outFile("acol_results.txt");
    outFile << std::fixed << std::setprecision(7);
    outFile << "# acol_o    acol_f    entries\n";
    std::cout << std::fixed << std::setprecision(7);
    
    std::vector<double> meanValues, meanErrValues, sigmaValues, sigmaErrValues;

    // 4. Loop over bins
    for (int b = 0; b < nBins; b++) {
        double lo = edges[b], hi = edges[b + 1];
        std::cout << "\n=== Bin " << b + 1 << ": Acolinearity ∈ ["
                  << lo << ", " << hi << ") ===" << std::endl;

        RooRealVar gMass("gMass", "Reconstructed Mass", xmin, xmax);
        RooArgSet vars(gMass);
        RooDataSet gData("data", "dataset from tree", vars);

        for (auto& [mass, acol] : entries) {
            if (acol >= lo && acol < hi) {
                gMass.setVal(mass);
                gData.add(vars);
            }
        }

        std::cout << "Entries in dataset: " << gData.numEntries() << std::endl;
        gMass.setRange("fitRange", xmin, xmax);
        gMass.setRange("plotRange", xmin, xmax);

        RooRealVar mean("mean", "mean", 1.7768199, 1.75, 1.85);
        RooRealVar sigma("sigma", "sigma", 0.007, 1e-5, 0.02);
        RooRealVar alphaL("alphaL", "alphaL", 0.5, 0.1, 10.0);
        RooRealVar nL("nL", "nL", 1.0, 0.1, 10.0);
        RooRealVar alphaR("alphaR", "alphaR", 0.5, 0.1, 10.0);
        RooRealVar nR("nR", "nR", 1.0, 0.1, 10.0);
        RooCrystalBall sig("sig", "Signal function", gMass, mean, sigma, alphaL, nL, alphaR, nR);

        auto fitResult = sig.fitTo(gData, RooFit::Strategy(2),
                                   RooFit::Range("fitRange"), RooFit::Save(), RooFit::PrintLevel(-1));
        fitResult->Print();
        std::cout << "Mean: " << mean.getVal() << "\n";
        std::cout << "Mean error: " << mean.getError() << "\n";
        std::cout << "Sigma: " << sigma.getVal() << "\n";
        std::cout << "Sigma error: " << sigma.getError() << "\n";
        
        outFile << std::setprecision(5) << lo << "     " << hi << "     "
                << std::setprecision(7) << gData.numEntries() << "\n";
        
        meanValues.push_back(mean.getVal()); meanErrValues.push_back(mean.getError());
        sigmaValues.push_back(sigma.getVal()); sigmaErrValues.push_back(sigma.getError());
    }
    
    outFile << "\n";
    outFile << "[";
    for (int i = 0; i < nBins; i++) {
        outFile << meanValues[i] * 1000; // MeV
        if (i != nBins - 1) outFile << ", ";
    }
    outFile << "]\n";

    outFile << "[";
    for (int i = 0; i < nBins; i++) {
        outFile << meanErrValues[i] * 1000;
        if (i != nBins - 1) outFile << ", ";
    }
    outFile << "]\n";
    
        outFile << "[";
    for (int i = 0; i < nBins; i++) {
        outFile << sigmaValues[i] * 1000;
        if (i != nBins - 1) outFile << ", ";
    }
    outFile << "]\n";

    outFile << "[";
    for (int i = 0; i < nBins; i++) {
        outFile << sigmaErrValues[i] * 1000;
        if (i != nBins - 1) outFile << ", ";
    }
    outFile << "]";
    
    outFile.close();
}