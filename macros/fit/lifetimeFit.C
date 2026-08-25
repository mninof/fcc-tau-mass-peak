///////////////////////////////////////////////////////////////////////////////////
///// Future Circular Collider Software - Tau-lepton Mass Peak Reconstruction /////
/////       This macro performs the fit of the tau-lepton lifetime            /////
///////////////////////////////////////////////////////////////////////////////////

#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "RooRealVar.h"
#include "RooDataSet.h"

#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"

#include "TPad.h"
#include "TPaveText.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"

#include <iostream>

/* This macro can be called once the rootfile from the simulation has been loaded to
   ROOT. The argument, a TTree, corresponds to the name of the tree in the rootfile
   where the events are stored. Typically, this name is just "events" */
void lifetimeFit(TTree* events)
{       
    // 1. Define the observable
    double xmin = 0; double xmax = 3000;
    RooRealVar t("t", "Reconstructed lifetime", xmin, xmax);
    RooArgSet vars(t);

    // 2. Create an empty dataset
    RooDataSet data("data", "dataset from tree", vars);
    
    /* Variable "mogens": When set to false, fits the lifetime from the two tau lifetimes
                          separated. Requires an estimation of the PV. Uses two branches.
                          When set to true, fits the lifetime from the sum of the two taus'
                          lifetimes. It is preferrable, since no estimation of the PV is needed, 
                          but requires a slight modification of the steering file. Uses one branch */
    bool mogens = false;
    if (mogens == false) {
        TTreeReader reader(events);
        TTreeReaderValue<double> lifetime(reader, "TauLifet");
        TTreeReaderValue<double> lifetimeCC(reader, "TauLifet_CC");

        while (reader.Next()) {
            if (*lifetime > xmin && *lifetime < xmax) {
                t.setVal(*lifetime);
                data.add(vars);
            }
            if (*lifetimeCC > xmin && *lifetimeCC < xmax) {
                t.setVal(*lifetimeCC);
                data.add(vars);
            }
        }

        std::cout << "Entries in dataset: " << data.numEntries() << std::endl;
        t.setRange("fitRange", xmin, xmax);
        t.setRange("plotRange", xmin, xmax);

        // 3. Build signal mode (exponential decay)
        RooRealVar tau("tau", "tau", 290.3, 250.0, 350.0);
        RooGenericPdf sig("pdf", "exp(-t/tau)", RooArgList(t, tau));

        // 4. Perform the Fit
        auto fitmodel = sig.fitTo(data, RooFit::Minos(true), RooFit::Strategy(2), RooFit::Range("fitRange"), RooFit::Save());
        fitmodel->Print();
    }
    else {
        TTreeReader reader(events);
        TTreeReaderValue<double> lifetime(reader, "LifetimeMogens");

        while (reader.Next()) {
            if (*lifetime > xmin && *lifetime < xmax) {
                t.setVal(*lifetime);
                data.add(vars);
            }
        }

        std::cout << "Entries in dataset: " << data.numEntries() << std::endl;
        t.setRange("fitRange", xmin, xmax);
        t.setRange("plotRange", xmin, xmax);

        // 3. Build signal mode (exponential decay)
        RooRealVar gammaShape("gammaShape", "shape", 2.0); gammaShape.setConstant(true); // fixed: sum of 2 exponentials
        RooRealVar mu("mu", "threshold", 0.0); mu.setConstant(true);

        RooRealVar tau("tau", "tau", 290.3, 250.0, 350.0);  // this plays the role of "beta"
        RooGamma sig("pdf", "Gamma pdf", t, gammaShape, tau, mu);

        // 4. Perform the Fit
        auto fitmodel = sig.fitTo(data, RooFit::Minos(true), RooFit::Strategy(2), RooFit::Range("fitRange"), RooFit::Save());
        fitmodel->Print();
        
        // 5. Plotting
        bool pulls = true; // true - produces plot with pull frame, false - plots only exp data    
        if (pulls == true) {
            TCanvas c("c", "c", 700, 600);
            double prop_main = 0.75; // proportion of the whole graph corresponding to main plot
            double fs = 0.06; // fontsize as a fraction of the main frame height
            //double prop_main = 0.99;
            double scale = prop_main / (1 - prop_main);
            TPad* pad1 = new TPad("pad1", "", 0, 1 - prop_main, 1, 1);
            TPad* pad2 = new TPad("pad2", "", 0, 0,   1, 1 - prop_main);
            pad1->SetBottomMargin(0.03); pad1->SetLeftMargin(0.15); pad1->SetRightMargin(0.05);
            pad2->SetTopMargin(0.02); pad2->SetLeftMargin(0.15); pad2->SetRightMargin(0.05);
            pad2->SetBottomMargin(0.4);
            pad1->Draw(); pad2->Draw();

            // Main frame
            RooPlot* frame = t.frame(RooFit::Title(" "));
            frame->GetYaxis()->SetTitleOffset(1.3); // push label away from axis
            frame->GetXaxis()->SetLabelSize(0);
            frame->GetXaxis()->SetTitleSize(0);
            //frame->GetXaxis()->SetTitle("MassReco [GeV]");
            frame->GetYaxis()->SetLabelSize(fs);
            frame->GetYaxis()->SetTitleSize(fs);

            data.plotOn(frame, RooFit::Name("data"));
            sig.plotOn(frame, RooFit::Name("model"), RooFit::LineColor(kRed), RooFit::LineWidth(2));

            ////double chi2ndf = frame->chiSquare(fitmodel->floatParsFinal().getSize());
            TPaveText* stats = new TPaveText(0.65, 0.75, 0.85, 0.85, "NDC");
            stats->SetFillColor(0); stats->SetBorderSize(1);
            stats->SetTextFont(42); stats->SetTextSize(0.04);
            ////stats->AddText(Form("#chi^{2}/ndf = %.3f", chi2ndf));

            RooHist* pullHist = frame->pullHist();
            RooPlot* pullFrame = t.frame(RooFit::Title(" "));
            pullFrame->addPlotable(pullHist, "P");
            pullFrame->GetYaxis()->SetTitle("Pulls");
            pullFrame->GetYaxis()->SetNdivisions(5);

            pullFrame->GetXaxis()->SetTitle("Reconstructed #tau_{#tau} [fs]");
            pullFrame->GetXaxis()->SetLabelSize(fs * scale);
            pullFrame->GetXaxis()->SetTitleSize(fs * scale);
            pullFrame->GetYaxis()->SetLabelSize(fs * scale);
            pullFrame->GetYaxis()->SetTitleSize(fs * scale);
            pullFrame->GetYaxis()->SetTitleOffset(1.3 / scale);

            /*model.plotOn(frame, RooFit::Components("sig"), RooFit::LineStyle(kDashed),
                           RooFit::LineColor(kBlue), RooFit::LineWidth(2));
             model.plotOn(frame, RooFit::Components("bkg"), RooFit::LineStyle(kDashed),
                          RooFit::LineColor(kGreen), RooFit::LineWidth(2));*/

            pad1->cd(); frame->drawBefore("data", "model"); frame->Draw(); ////stats->Draw();
            pad2->cd(); pullFrame->Draw();
            gSystem->mkdir("plots", kTRUE);
            c.SaveAs("plots/lifetime_plot_fit.pdf");
        }
        else {
            TCanvas c("c", "c", 900, 600);

            TPad* pad1 = new TPad("pad1", "", 0, 0, 1, 1);
            pad1->SetBottomMargin(0.15);
            pad1->SetLeftMargin(0.14);
            pad1->Draw();

            RooPlot* frame = t.frame(RooFit::Title(" "));

            frame->GetYaxis()->SetTitleOffset(1.4);
            frame->GetXaxis()->SetTitleOffset(1.4);
            frame->GetXaxis()->SetTitle("Reconstructed #tau_{#tau} [fs]");
            frame->GetXaxis()->SetLabelSize(0.05);
            frame->GetXaxis()->SetTitleSize(0.05);
            frame->GetYaxis()->SetLabelSize(0.05);
            frame->GetYaxis()->SetTitleSize(0.05);

            data.plotOn(frame, RooFit::Name("data"));
            //model.plotOn(frame, RooFit::Name("model"), RooFit::LineColor(kRed), RooFit::LineWidth(2));

            pad1->cd(); frame->drawBefore("data", "model"); frame->Draw();

            gSystem->mkdir("plots", kTRUE);
            c.SaveAs("plots/lifetime_plot_raw.pdf");
        }
    }
}