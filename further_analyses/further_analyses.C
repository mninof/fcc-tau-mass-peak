///////////////////////////////////////////////////////////////////////////////////////
///// Future Circular Collider Software - Tau-lepton Mass Peak Reconstruction     /////
///// The purpose of this macro is producing new rootfiles with data or variables /////
///// that are not available in the original rootfile produced with the steering  /////
///// file. e.g combinations of the branches into new kinematic variables, or     /////
///// a re-computation of the tau mass with other assumptions. This file should   ///// 
///// be modifed according to the purpose of the new analysis to perform          /////
///////////////////////////////////////////////////////////////////////////////////////

#include <ROOT/RDataFrame.hxx>
#include <TRandom3.h>
#include <TMath.h>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "helperFunctions.C" // This file contains several crucial functions for the new analyses!

/* No argument is required in order to call this macro. It must be modified internally according
   to what new data is required for the new analysis */
void further_analyses()
{
    //std::cout << "Beam energy lower limit: " << E_TAU_M << " GeV" << std::endl;
    //std::cout << "Beam energy upper limit: " << E_TAU_P << " GeV" << std::endl;
    //std::cout << "Tau mass lower limit: " << TAU_MASS_M * 1000 << " MeV" << std::endl;
    //std::cout << "Tau mass lower limit: " << TAU_MASS_P * 1000 << " MeV" << std::endl;
    
    /* This is the path to the original rootfile produced by the steering file. The new analyses
       are typically combinations and/or modifications produced through the data cotained in that
       original rootfile */
    ROOT::RDataFrame df("events", "result/p8_ee_Ztautau_ecm91.root");

    /* Here below are the computation of the new variables for the analysis. The new branches defined
       in this block can later be written into a new rootfile.
    
       All the helper functions imported from "helperFunctions.C" should be called here */
    auto df2 = df
        .Define("x_sm", [](double x) { return smear_vtx(x, kX); }, {"TauMCVertex_x"})
        .Define("y_sm", [](double y) { return smear_vtx(y, kY); }, {"TauMCVertex_y"})
        .Define("z_sm", [](double z) { return smear_vtx(z, kZ); }, {"TauMCVertex_z"})
        .Define("dx","TauMCVertex_x - TauVertex_x")
        .Define("dx_sm","TauMCVertex_x - x_sm")
        .Define("dy","TauMCVertex_y - TauVertex_y")
        .Define("dy_sm","TauMCVertex_y - y_sm")
        .Define("dz","TauMCVertex_z - TauVertex_z")
        .Define("dz_sm","TauMCVertex_z - z_sm")
        .Define("sqrt_s", "TauMC_E + TauMC_E_CC")
        .Define("NewEnergy", Etau, {"Photons_E_tot"})
        .Define("MassMC", MassWithEnergy_MC, {"TauMCVertex_x", "TauMCVertex_y", "TauMCVertex_z",
                                              "TauMCVertex_x_CC", "TauMCVertex_y_CC", "TauMCVertex_z_CC",
                                              "invMass_3pi_MC", "TauMC_E",
                                              "PiMC_px", "PiMC_py", "PiMC_pz"})
        .Define("MassMC_CC", MassWithEnergy_MC, {"TauMCVertex_x", "TauMCVertex_y", "TauMCVertex_z",
                                                 "TauMCVertex_x_CC", "TauMCVertex_y_CC", "TauMCVertex_z_CC",
                                                 "invMass_3pi_MC_CC", "TauMC_E_CC",
                                                 "PiMC_px_CC", "PiMC_py_CC", "PiMC_pz_CC"})
        .Define("MassReco", MassWithEnergy, {"TauVertex_x", "TauVertex_y", "TauVertex_z",
                                             "TauVertex_x_CC", "TauVertex_y_CC", "TauVertex_z_CC",
                                             "invMass_3pi_fit", "NewEnergy",
                                             "Pi_px", "Pi_py", "Pi_pz"})
        .Define("MassReco_CC", MassWithEnergy, {"TauVertex_x", "TauVertex_y", "TauVertex_z",
                                                "TauVertex_x_CC", "TauVertex_y_CC", "TauVertex_z_CC",
                                                "invMass_3pi_fit_CC", "NewEnergy",
                                                "Pi_px_CC", "Pi_py_CC", "Pi_pz_CC"})
        .Define("Acolinearity", DotProductUni, {"Pi_px", "Pi_py", "Pi_pz",
                                                "Pi_px_CC", "Pi_py_CC", "Pi_pz_CC"});
    
    // The line below stores the requested branches of the new dataframe into a new rootfile
    df2.Snapshot("events", "result/mass_corrected_by_photons.root", {"MassReco", "MassReco_CC", "Acolinearity"});

    // The code block below was used to study several combination of truth-level with reco-level kinematic variables
    /*double centre = 1.77693; double range = 0.7; int nBins = 75;
    auto h1 = df2.Histo1D({"h1", ";Mass [GeV];Events", nBins, centre - range, centre + range}, "MCvtx_MCmom_MCm");
    auto h2 = df2.Histo1D({"h2", ";Mass [GeV];Events", nBins, centre - range, centre + range}, "MCvtx_MCmom_Rm");
    auto h3 = df2.Histo1D({"h3", ";Mass [GeV];Events", nBins, centre - range, centre + range}, "MCvtx");
    auto h4 = df2.Histo1D({"h4", ";Mass [GeV];Events", nBins, centre - range, centre + range}, "R");
    auto h5 = df2.Histo1D({"h5", ";Mass [GeV];Events", nBins, centre - range, centre + range}, "MassReco");
    TCanvas c("c", "c", 800, 600);
    h1->Draw(); h1->SetLineWidth(2);
    h2->Draw("SAME"); h2->SetLineWidth(2); h2->SetLineColor(kRed);
    h3->Draw("SAME"); h3->SetLineWidth(2); h3->SetLineColor(kGreen);
    h4->Draw("SAME"); h4->SetLineWidth(2); h4->SetLineColor(kBlack);
    h5->Draw("SAME"); h5->SetLineWidth(2); h5->SetLineColor(kMagenta);
    c.SaveAs("plots/sm.pdf");*/
    
    bool plot_res = false; /* Set true to plot the vertex resolution
                              in this case, this will plot the beam spot resolution in y direction */
    if (plot_res == true) {
        auto h2_ptr = df2.Histo1D({"h2", ";Resolution in y [mm];Events", 100, -1e-4, 1e-4}, "ZMCVertex_y");
        TH1D& h2 = *h2_ptr;
        TCanvas c("c", "c", 800, 800);
        gStyle->SetOptStat(0);
        TPad* pad1 = new TPad("pad1", "", 0, 0, 1, 1);
        pad1->SetBottomMargin(0.12); //pad1->SetTopMargin(0.05);
        pad1->SetLeftMargin(0.18); //pad1->SetRightMargin(0.05);
        pad1->SetTickx(1); pad1->SetTicky(1);
        //pad1->SetLogy(); // log-scale for y axis (useful for sqrt(s) plots)
        pad1->Draw(); pad1->cd();

        h2.SetLineColor(kRed); h2.SetLineWidth(2);
        h2.GetXaxis()->SetTitleOffset(1.2); h2.GetYaxis()->SetTitleOffset(2.0);
        h2.GetXaxis()->SetLabelSize(0.045); h2.GetXaxis()->SetTitleSize(0.05); h2.GetXaxis()->SetNdivisions(505);
        h2.GetYaxis()->SetLabelSize(0.045); h2.GetYaxis()->SetTitleSize(0.05);

        h2.SetMinimum(1); h2.Draw("HIST");
        //TLegend* leg = new TLegend(0.58, 0.75, 0.81, 0.88); // upper-right corner
        TLegend* leg = new TLegend(0.20, 0.75, 0.43, 0.88); // upper-left corner

        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.05);
        leg->AddEntry(&h2, "MC truth", "l");
        //leg->Draw();

        c.SaveAs("plots/BS_res_y.pdf");
    }
}