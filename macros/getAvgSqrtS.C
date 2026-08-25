///////////////////////////////////////////////////////////////////////////////////
///// Future Circular Collider Software - Tau-lepton Mass Peak Reconstruction /////
///// This macro computes the average value of the centre-of-mass energy at   /////
///// the MC level. Prints and returns the value as a double                  /////
///////////////////////////////////////////////////////////////////////////////////

#include <iomanip>

// This macro can be called once the rootfile from the simulation has been loaded to
// ROOT. The argument, a TTree, corresponds to the name of the tree in the rootfile
// where the events are stored. Typically, this name is just "events"
double GetSqrtSMean(TTree* tree) {
    tree->SetAlias("sqrt_s", "TauMC_E + TauMC_E_CC");
    tree->Draw("sqrt_s>>htemp", "", "goff"); // goff: graphics off
    TH1F* htemp = (TH1F*)gDirectory->Get("htemp");
    double mean = htemp->GetMean();
    std::cout << std::setprecision(6) << "Mean of sqrt_s = " << mean << std::endl;

    return mean;
}