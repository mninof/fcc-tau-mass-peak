///////////////////////////////////////////////////////////////////////////////////
///// Future Circular Collider Software - Tau-lepton Mass Peak Reconstruction /////
///// This macro produces the plots of the residuals of the invariant mass of /////
///// the 3-pion system                                                       /////
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
void invMassResiduals(TTree* events)
{
    TTreeReader reader(events);

    TTreeReaderValue<double> invMass_3pi_MC(reader, "invMass_3pi_MC");
    TTreeReaderValue<double> invMass_3pi_raw(reader, "invMass_3pi_raw");
    TTreeReaderValue<double> invMass_3pi_fit(reader, "invMass_3pi_fit");
    TTreeReaderValue<double> invMass_3pi_MC_CC(reader, "invMass_3pi_MC_CC");
    TTreeReaderValue<double> invMass_3pi_raw_CC(reader, "invMass_3pi_raw_CC");
    TTreeReaderValue<double> invMass_3pi_fit_CC(reader, "invMass_3pi_fit_CC");

    std::vector<double> invMass_fit, invMass_raw;

    while (reader.Next()) {
        invMass_raw.push_back((*invMass_3pi_raw - *invMass_3pi_MC) / *invMass_3pi_MC);
        invMass_raw.push_back((*invMass_3pi_raw_CC - *invMass_3pi_MC_CC) / *invMass_3pi_MC_CC);

        invMass_fit.push_back((*invMass_3pi_fit - *invMass_3pi_MC) / *invMass_3pi_MC);
        invMass_fit.push_back((*invMass_3pi_fit_CC - *invMass_3pi_MC_CC) / *invMass_3pi_MC_CC);
    }

    lhcbStyle();
    //gStyle->SetOptStat(1); gStyle->SetOptStat("MR");

    double plot_range = 0.014;
    auto makeHist = [plot_range](const char* name, std::vector<double>& v, int color) {
        TH1D* h = new TH1D(name, "", 100, -plot_range, plot_range);
        for (double val : v) h->Fill(val);
        h->SetLineColor(color);
        return h;
    };

    TH1D* hraw = makeHist("hraw", invMass_raw, kBlue+2);
    TH1D* hfit = makeHist("hfit", invMass_fit, kRed);

    TH1D* histsRow1[2] = {hraw, hfit};
    const char* xTitles[2] = {"(M_{3#pi}^{pre} - M_{3#pi}^{MC}) / M_{3#pi}^{MC}",
                              "(M_{3#pi}^{pos} - M_{3#pi}^{MC}) / M_{3#pi}^{MC}"};

    TCanvas c("c", "c", 1300, 600);
    c.Divide(2, 1, 0, 0);

    for (int i = 0; i < 2; i++)
    {
        TVirtualPad* pad1 = c.cd(i + 1);
        //pad1->SetLeftMargin(i == 0 ? 0.25 : 0.15); pad1->SetRightMargin(0.13);
        pad1->SetLeftMargin(0.25);
        
        histsRow1[i]->GetXaxis()->SetTitle(xTitles[i]);
        histsRow1[i]->GetYaxis()->SetTitle(i == 0 ? "Entries" : "");
        histsRow1[i]->GetYaxis()->SetTitleOffset(1.8);

        histsRow1[i]->Draw("HIST");

        TLegend* leg1 = new TLegend(0.55, 0.75, 0.92, 0.88);
        leg1->AddEntry(histsRow1[i], i == 0 ? "Pre-fit" : "Post-fit", "L");

        auto label1 = LHCbStyle::create_label("FCC-ee");
        label1->SetFillStyle(0);
    }

    gSystem->mkdir("plots", kTRUE);
    c.SaveAs("plots/invMassResiduals.pdf");
}