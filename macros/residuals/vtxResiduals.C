///////////////////////////////////////////////////////////////////////////////////
///// Future Circular Collider Software - Tau-lepton Mass Peak Reconstruction /////
///// This macro produces the plots of the residuals of the 3-momentum of the /////
///// 3-pion system                                                           /////  
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

#include "lhcbStyle.C"
#include <sstream>

/* This macro can be called once the rootfile from the simulation has been loaded to
   ROOT. The argument, a TTree, corresponds to the name of the tree in the rootfile
   where the events are stored. Typically, this name is just "events" */
void vtxResiduals(TTree* events)
{
    TTreeReader reader(events);

    TTreeReaderValue<double> vx_MC(reader, "TauMCVertex_x");
    TTreeReaderValue<float> vx_reco(reader, "TauVertex_x");
    TTreeReaderValue<double> vx_MC_CC(reader, "TauMCVertex_x_CC");
    TTreeReaderValue<float> vx_reco_CC(reader, "TauVertex_x_CC");

    TTreeReaderValue<double> vy_MC(reader, "TauMCVertex_y");
    TTreeReaderValue<float> vy_reco(reader, "TauVertex_y");
    TTreeReaderValue<double> vy_MC_CC(reader, "TauMCVertex_y_CC");
    TTreeReaderValue<float> vy_reco_CC(reader, "TauVertex_y_CC");

    TTreeReaderValue<double> vz_MC(reader, "TauMCVertex_z");
    TTreeReaderValue<float> vz_reco(reader, "TauVertex_z");
    TTreeReaderValue<double> vz_MC_CC(reader, "TauMCVertex_z_CC");
    TTreeReaderValue<float> vz_reco_CC(reader, "TauVertex_z_CC");

    std::vector<double> vres_vx, vres_vy, vres_vz;

    while (reader.Next()) {
        vres_vx.push_back((*vx_reco - *vx_MC) / *vx_MC);
        vres_vx.push_back((*vx_reco_CC - *vx_MC_CC) / *vx_MC_CC);

        vres_vy.push_back((*vy_reco - *vy_MC) / *vy_MC);
        vres_vy.push_back((*vy_reco_CC - *vy_MC_CC) / *vy_MC_CC);

        vres_vz.push_back((*vz_reco - *vz_MC) / *vz_MC);
        vres_vz.push_back((*vz_reco_CC - *vz_MC_CC) / *vz_MC_CC);
    }

    lhcbStyle();
    //gStyle->SetOptStat(1); gStyle->SetOptStat("MR"); // add stats box - M: mean, R: stdev
    // -- might clash with lhcbStyle

    double plot_range = 0.3;
    auto makeHist = [plot_range](const char* name, std::vector<double>& v, int color) {
        TH1D* h = new TH1D(name, "", 100, -plot_range, plot_range);
        for (double val : v) h->Fill(val);
        h->SetLineColor(color);
        return h;
    };

    TH1D* hvx = makeHist("hvx", vres_vx, kRed);
    TH1D* hvy = makeHist("hvy", vres_vy, kRed);
    TH1D* hvz = makeHist("hvz", vres_vz, kRed);

    TH1D* hists[3] = {hvx, hvy, hvz};
    const char* xTitles[3] = {"(v^{reco}_{x} - v^{MC}_{x}) / v^{MC}_{x}",
                               "(v^{reco}_{y} - v^{MC}_{y}) / v^{MC}_{y}",
                               "(v^{reco}_{z} - v^{MC}_{z}) / v^{MC}_{z}"};

    TCanvas c("c", "c", 1800, 600);
    c.Divide(3, 1, 0, 0);

    for (int i = 0; i < 3; i++)
    {
        TVirtualPad* pad = c.cd(i + 1);
        //pad->SetLeftMargin(i == 0 ? 0.25 : 0.15);
        pad->SetLeftMargin(0.25);

        hists[i]->GetXaxis()->SetTitle(xTitles[i]);
        hists[i]->GetYaxis()->SetTitle(i == 0 ? "Entries" : "");
        hists[i]->GetYaxis()->SetTitleOffset(1.8);

        hists[i]->Draw("HIST");

        TLegend* leg = new TLegend(0.55, 0.75, 0.92, 0.88);
        //leg->AddEntry(hists[i], "Reco", "L");

        auto label = LHCbStyle::create_label("FCC-ee");
        //label->SetFillStyle(0);
    }

    gSystem->mkdir("plots", kTRUE);
    c.SaveAs("plots/vtxResiduals.pdf");
}