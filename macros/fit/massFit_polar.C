///////////////////////////////////////////////////////////////////////////////////
///// Future Circular Collider Software - Tau-lepton Mass Peak Reconstruction /////
///// This macro performs the fit of the tau-lepton mass for several bins of  /////
///// the polar angle. Used to study cuts on this quantity. This was part of  /////
///// an study of polar-angle impact originally requested by Alberto Lusiani  /////
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
  
   This macro requires two new branches, corresponding to the polar angles of the taus,
   which can be produced with the files and setting provided in the /further_analyses directory.
  
   The produced result is a text file containing all the information from the bins */
void massFit_polar(TTree* events)
{       
    // 1. Define the observable
    double centre = 1.7768199; double ext = 0.3;
    double xmin = centre - ext; double xmax = centre + ext;

    // 2. Define cos(theta) bins
    const int nBins = 10;
    double edges[nBins + 1];
    for (int i = 0; i <= nBins; i++)
        edges[i] = i * 0.1;

    // 3. Read all entries once into vectors
    TTreeReader reader(events);
    TTreeReaderValue<double> massReco(reader, "MassReco");
    TTreeReaderValue<double> massRecoCC(reader, "MassReco_CC");
    TTreeReaderValue<double> cosTheta(reader, "CosTheta");
    TTreeReaderValue<double> cosThetaCC(reader, "CosTheta_CC");

    // Store as pairs (mass, |cosTheta|)
    std::vector<std::pair<double,double>> entries;
    while (reader.Next()) {
        double absCos = std::abs(*cosTheta);
        double absCosCC = std::abs(*cosThetaCC);
        if (*massReco > xmin && *massReco < xmax)
            entries.push_back({*massReco, absCos});
        if (*massRecoCC > xmin && *massRecoCC < xmax)
            entries.push_back({*massRecoCC, absCosCC});
    }
    
    std::ofstream outFile("polar_results.txt");
    outFile << std::fixed << std::setprecision(7);
    outFile << "# cosTheta_o  cosTheta_f  mean       mean_error  sigma      sigma_err   entries\n";
    std::cout << std::fixed << std::setprecision(7);

    // 4. Loop over bins
    for (int b = 0; b < nBins; b++) {
        double lo = edges[b], hi = edges[b + 1];
        std::cout << "\n=== Bin " << b << ": |cosθ| ∈ ["
                  << lo << ", " << hi << ") ===" << std::endl;

        RooRealVar gMass("gMass", "Reconstructed Mass", xmin, xmax);
        RooArgSet vars(gMass);
        RooDataSet gData("data", "dataset from tree", vars);

        for (auto& [mass, absCos] : entries) {
            if (absCos >= lo && absCos < hi) {
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
        
        bool correlations = false; // when true, renders fit params. correlation matrices of some bins
        if (correlations == true && (b == 1 || b == 3 || b == 5 || b == 7) ) {
            TCanvas map(Form("map_%d", b), "map", 600, 600);
            map.SetRightMargin(0.18);
            TH2* hcorr = (TH2*)fitResult->correlationHist()->Clone(Form("corr_bin%d", b));
            hcorr->Draw("colz text"); hcorr->SetTitle(""); hcorr->SetStats(0);
            std::string filename = "fit_correlations_bin" + std::to_string(b) + ".pdf";
            map.SaveAs(filename.c_str());
        }
        
        outFile << lo << "     " << hi << "   " << mean.getVal() << "  " << mean.getError()
                << "   " << sigma.getVal() << "  " << sigma.getError() << "  " << gData.numEntries() << "\n";
    }
    
    outFile.close();
}