///////////////////////////////////////////////////////////////////////////////////////
///// Future Circular Collider Software - Tau-lepton Mass Peak Reconstruction     /////
///// This macro calculates the systematic uncertainty in the mass reconstruction /////
///// associated to the fixed (nominal) value used for E_tau. Result: 1.7 keV     /////
///////////////////////////////////////////////////////////////////////////////////////

#include "TTree.h"
#include "RooRealVar.h"
#include "RooDataSet.h"

#include "RooAddPdf.h"

#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"

#include <iostream>

/* The function below should NOT be called from ROOT. The only correct call from ROOT is
   to the function "sysFromEtau.C" */
double myRooFit(TTree* events, const char* branch1, const char* branch2)
{       
    // 1. Define the observable
    double centre = 1.77683; double ext = 0.2;
    double xmin = centre - ext; double xmax = centre + ext;
    RooRealVar gMass("gMass", "Reconstructed Mass", xmin, xmax);
    RooArgSet vars(gMass);

    // 2. Create an empty dataset and fill it manually from both branches
    RooDataSet gData("data", "dataset from tree", vars);

    TTreeReader reader(events);
    TTreeReaderValue<double> massReco(reader, branch1);
    TTreeReaderValue<double> massRecoCC(reader, branch2);

    while (reader.Next()) {
        if (*massReco > xmin && *massReco < xmax) {
            gMass.setVal(*massReco);
            gData.add(vars);
        }
        if (*massRecoCC > xmin && *massRecoCC < xmax) {
            gMass.setVal(*massRecoCC);
            gData.add(vars);
        }
    }

    std::cout << "Entries in dataset: " << gData.numEntries() << std::endl;
    gMass.setRange("fitRange", xmin, xmax);
    gMass.setRange("plotRange", xmin, xmax);

    // 3. Build signal mode (Double-sided Crystal Ball)
    RooRealVar mean("mean","mean", 1.77683, 1.75, 1.85);
    RooRealVar sigma("sigma","sigma", 0.007, 0.001, 0.02);
    RooRealVar alphaL("alphaL","alphaL", 0.5, 0.1, 1.0);
    RooRealVar nL("nL","nL", 2.0, 1.0, 100.0);
    RooRealVar alphaR("alphaR","alphaR", 0.5, 0.1, 1.0);
    RooRealVar nR("nR","nR", 2.0, 1.0, 100.0);

    RooCrystalBall sig("sig","Signal function",gMass, mean, sigma,
                        alphaL, nL, alphaR, nR);
    
    // 4. Perform the Fit
    auto fitmodel = sig.fitTo(gData, RooFit::Strategy(1), RooFit::Range("fitRange"), RooFit::Save());
    //fitmodel->Print();
    
    //std::cout << mean.getVal() << endl;
    return mean.getVal();
}

/* This macro can be called once the rootfile from the simulation has been loaded to
   ROOT. The argument, a TTree, corresponds to the name of the tree in the rootfile
   where the events are stored. Typically, this name is just "events".
  
   This macro requires producing a new rootfile (apart from the original one produced
   with the steering file) which contains four new branches: "MassReco_m",
   "MassReco_CC_m", "MassReco_p" and "MassReco_CC_p". Those are the reco'ed masses
   using the extreme values for the nominal E_tau (i.e. + and minus the reported
   uncertainty on (mass of Z boson)/2). These new branches may be produced with
   the files and setting provided in the /further_analyses directory */
void sysFromEtau(TTree* events)
{
    double M = myRooFit(events, "MassReco", "MassReco_CC");
    double M_m = myRooFit(events, "MassReco_m", "MassReco_CC_m");
    double M_p = myRooFit(events, "MassReco_p", "MassReco_CC_p");
    
    std::cout << "M: " << M << endl;
    std::cout << "M_m: " << M_m << endl;
    std::cout << "M_p: " << M_p << endl;
    
    double M_sys = std::max( std::abs(M_m - M),  std::abs(M_p - M) );
    std::cout << "M_sys: " << M_sys << endl;
}