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
void momResiduals(TTree* events)
{
    TTreeReader reader(events);
    
    TTreeReaderValue<float> px_MC(reader, "PiMC_px");
    TTreeReaderValue<float> px_pre(reader, "Tracks_px");
    TTreeReaderValue<double> px_pos(reader, "Pi_px");
    TTreeReaderValue<float> px_MC_CC(reader, "PiMC_px_CC");
    TTreeReaderValue<float> px_pre_CC(reader, "Tracks_px_CC");
    TTreeReaderValue<double> px_pos_CC(reader, "Pi_px_CC");
    
    TTreeReaderValue<float> py_MC(reader, "PiMC_py");
    TTreeReaderValue<float> py_pre(reader, "Tracks_py");
    TTreeReaderValue<double> py_pos(reader, "Pi_py");
    TTreeReaderValue<float> py_MC_CC(reader, "PiMC_py_CC");
    TTreeReaderValue<float> py_pre_CC(reader, "Tracks_py_CC");
    TTreeReaderValue<double> py_pos_CC(reader, "Pi_py_CC");
    
    TTreeReaderValue<float> pz_MC(reader, "PiMC_pz");
    TTreeReaderValue<float> pz_pre(reader, "Tracks_pz");
    TTreeReaderValue<double> pz_pos(reader, "Pi_pz");
    TTreeReaderValue<float> pz_MC_CC(reader, "PiMC_pz_CC");
    TTreeReaderValue<float> pz_pre_CC(reader, "Tracks_pz_CC");
    TTreeReaderValue<double> pz_pos_CC(reader, "Pi_pz_CC");

    std::vector<double> vres_px_pos, vres_px_pre;
    std::vector<double> vres_py_pos, vres_py_pre;
    std::vector<double> vres_pz_pos, vres_pz_pre;

    while (reader.Next()) {
        vres_px_pos.push_back((*px_pos  - *px_MC)  / *px_MC);
        vres_px_pos.push_back((*px_pos_CC - *px_MC_CC) / *px_MC_CC);
        vres_px_pre.push_back((*px_pre  - *px_MC)  / *px_MC);
        vres_px_pre.push_back((*px_pre_CC - *px_MC_CC) / *px_MC_CC);
        
        vres_py_pos.push_back((*py_pos  - *py_MC)  / *py_MC);
        vres_py_pos.push_back((*py_pos_CC - *py_MC_CC) / *py_MC_CC);
        vres_py_pre.push_back((*py_pre  - *py_MC)  / *py_MC);
        vres_py_pre.push_back((*py_pre_CC - *py_MC_CC) / *py_MC_CC);
        
        vres_pz_pos.push_back((*pz_pos  - *pz_MC)  / *pz_MC);
        vres_pz_pos.push_back((*pz_pos_CC - *pz_MC_CC) / *pz_MC_CC);
        vres_pz_pre.push_back((*pz_pre  - *pz_MC)  / *pz_MC);
        vres_pz_pre.push_back((*pz_pre_CC - *pz_MC_CC) / *pz_MC_CC);
    }
    
    lhcbStyle();
    // Uncomment to render boxes with errors included (caps here below) -- might clash with lhcbStyle
    //gStyle->SetOptStat(1); gStyle->SetOptStat("MR");
    
    auto makeHist = [](const char* name, std::vector<double>& v, int color) {
        TH1D* h = new TH1D(name, "", 100, -0.005, 0.005);
        for (double val : v) h->Fill(val);
        h->SetLineColor(color);
        return h;
    };

    TH1D* hpx_pos = makeHist("hpx_pos", vres_px_pos, kRed);
    TH1D* hpx_pre = makeHist("hpx_pre", vres_px_pre, kBlue+2);
    TH1D* hpy_pos = makeHist("hpy_pos", vres_py_pos, kRed);
    TH1D* hpy_pre = makeHist("hpy_pre", vres_py_pre, kBlue+2);
    TH1D* hpz_pos = makeHist("hpz_pos", vres_pz_pos, kRed);
    TH1D* hpz_pre = makeHist("hpz_pre", vres_pz_pre, kBlue+2);

    TH1D* histsRow1[3] = {hpx_pre, hpy_pre, hpz_pre};
    TH1D* histsRow2[3] = {hpx_pos, hpy_pos, hpz_pos};
    const char* xTitles[3] = {"(p^{reco}_{x} - p^{MC}_{x}) / p^{MC}_{x}",
                              "(p^{reco}_{y} - p^{MC}_{y}) / p^{MC}_{y}",
                              "(p^{reco}_{z} - p^{MC}_{z}) / p^{MC}_{z}"};
    
    TCanvas c("c", "c", 1800, 1200);
    c.Divide(3, 2, 0, 0);

    for (int i = 0; i < 3; i++) 
    {
        TVirtualPad* pad1 = c.cd(i + 1);
        //pad1->SetLeftMargin(i == 0 ? 0.25 : 0.15);
        pad1->SetLeftMargin(0.25);

        histsRow1[i]->GetXaxis()->SetTitle(xTitles[i]);
        histsRow1[i]->GetYaxis()->SetTitle(i == 0 ? "Entries" : "");
        histsRow1[i]->GetYaxis()->SetTitleOffset(1.8);

        histsRow1[i]->Draw("HIST");

        TLegend* leg1 = new TLegend(0.55, 0.75, 0.92, 0.88);
        leg1->AddEntry(histsRow1[i], "Pre-fit", "L");
        //leg1->Draw();

        auto label1 = LHCbStyle::create_label("FCC-ee");
        label1->SetFillStyle(0);
        //label1->Draw();
        
        TVirtualPad* pad2 = c.cd(i + 4);
        //pad2->SetLeftMargin(i == 0 ? 0.25 : 0.15);
        pad2->SetLeftMargin(0.25);

        histsRow2[i]->GetXaxis()->SetTitle(xTitles[i]);
        histsRow2[i]->GetYaxis()->SetTitle(i == 0 ? "Entries" : "");
        histsRow2[i]->GetYaxis()->SetTitleOffset(1.8);

        histsRow2[i]->Draw("HIST");

        TLegend* leg2 = new TLegend(0.55, 0.75, 0.92, 0.88);
        leg2->AddEntry(histsRow2[i], "Post-fit", "L");
        //leg2->Draw();

        auto label2 = LHCbStyle::create_label("FCC-ee");
        label2->SetFillStyle(0);
        //label2->Draw();
    }

    gSystem->mkdir("plots", kTRUE);
    c.SaveAs("plots/momResiduals.pdf");
}