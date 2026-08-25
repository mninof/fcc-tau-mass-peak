////////////////////////////////////////////////////////////////////////////////////
///// Future Circular Collider Software - Tau-lepton Mass Peak Reconstruction  /////
///// This macro produces the plots of the residuals of the variable cos(beta) /////
////////////////////////////////////////////////////////////////////////////////////

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
void cosResiduals(TTree* events)
{
    TTreeReader reader(events);
    
    bool pre_and_pos = false; // true - uses pre- and post-fit momentum, false - only post-fit momentum 
    if (pre_and_pos == true) {
        TTreeReaderValue<double> cosBeta_MC(reader, "cosBeta_MC");
        TTreeReaderValue<double> cosBeta_pre(reader, "cosBeta_pre");
        TTreeReaderValue<double> cosBeta_pos(reader, "cosBeta_pos");
        TTreeReaderValue<double> cosBeta_MC_CC(reader, "cosBeta_MC_CC");
        TTreeReaderValue<double> cosBeta_pre_CC(reader, "cosBeta_pre_CC");
        TTreeReaderValue<double> cosBeta_pos_CC(reader, "cosBeta_pos_CC");

        std::vector<double> vres_cos_pos, vres_cos_pre;

        while (reader.Next()) {
            vres_cos_pre.push_back((*cosBeta_MC - *cosBeta_pre) / *cosBeta_MC);
            vres_cos_pre.push_back((*cosBeta_MC_CC - *cosBeta_pre_CC) / *cosBeta_MC_CC);

            vres_cos_pos.push_back((*cosBeta_MC - *cosBeta_pos) / *cosBeta_MC);
            vres_cos_pos.push_back((*cosBeta_MC_CC - *cosBeta_pos_CC) / *cosBeta_MC_CC);
        }

        lhcbStyle();
        //gStyle->SetOptStat(1); gStyle->SetOptStat("MR");

        double plot_range = 2e-5;
        auto makeHist = [plot_range](const char* name, std::vector<double>& v, int color) {
            TH1D* h = new TH1D(name, "", 100, -plot_range, plot_range);
            for (double val : v) h->Fill(val);
            h->SetLineColor(color);
            return h;
        };

        TH1D* hcos_pre = makeHist("hcos_pre", vres_cos_pre, kBlue + 2);
        TH1D* hcos_pos = makeHist("hcos_pos", vres_cos_pos, kRed);

        TH1D* histsRow1[2] = {hcos_pre, hcos_pos};
        const char* xTitles[2] = {"(cos#beta^{MC} - cos#beta^{reco}_{pre}) / cos#beta^{MC}",
                                  "(cos#beta^{MC} - cos#beta^{reco}_{pos}) / cos#beta^{MC}"};

        TCanvas c("c", "c", 1300, 600);
        c.Divide(2, 1, 0, 0);

        for (int i = 0; i < 2; i++)
        {
            TVirtualPad* pad1 = c.cd(i + 1);
            pad1->SetLeftMargin(i == 0 ? 0.25 : 0.15);
            pad1->SetRightMargin(0.13);

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
        c.SaveAs("plots/cosResiduals.pdf");
    }
    
    else
    {
        TTreeReaderValue<double> CosBeta_MC(reader, "CosBeta_MC");
        TTreeReaderValue<double> CosBeta(reader, "CosBeta");
        TTreeReaderValue<double> CosBeta_MC_CC(reader, "CosBeta_MC_CC");
        TTreeReaderValue<double> CosBeta_CC(reader, "CosBeta_CC");

        std::vector<double> vres_cos;

        while (reader.Next()) {
            vres_cos.push_back((*CosBeta - *CosBeta_MC) / *CosBeta_MC);
            vres_cos.push_back((*CosBeta_CC - *CosBeta_MC_CC) / *CosBeta_MC_CC);
        }

        lhcbStyle();
        //gStyle->SetOptStat(1); gStyle->SetOptStat("MR");

        double plot_range = 1e-4;
        auto makeHist = [plot_range](const char* name, std::vector<double>& v, int color) {
            TH1D* h = new TH1D(name, "", 100, -plot_range, plot_range);
            for (double val : v) h->Fill(val);
            h->SetLineColor(color); h->SetLineWidth(4);
            return h;
        };

        TH1D* hcos = makeHist("hcos", vres_cos, kRed);

        TCanvas c("c", "c", 700, 600);
        c.SetLeftMargin(0.25);
        c.SetRightMargin(0.13);
        hcos->GetXaxis()->SetTitle("(cos#beta^{reco} - cos#beta^{MC}) / cos#beta^{MC}");
        hcos->GetYaxis()->SetTitle("Entries");
        hcos->GetYaxis()->SetTitleOffset(1.8);

        hcos->Draw("HIST");
        
        gSystem->mkdir("plots", kTRUE);
        c.SaveAs("plots/cosResiduals.pdf");
    }
}