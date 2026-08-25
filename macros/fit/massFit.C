///////////////////////////////////////////////////////////////////////////////////
///// Future Circular Collider Software - Tau-lepton Mass Peak Reconstruction /////
/////          This macro performs the fit of the tau-lepton mass             /////
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
void massFit(TTree* events)
{       
    // 1. Define the observable
    double centre = 1.7768199; double ext = 0.2;
    double xmin = centre - ext; double xmax = centre + ext;
    RooRealVar gMass("gMass", "Reconstructed Mass", xmin, xmax);
    //RooRealVar gMass("gMass", "Reconstructed Mass", -ext, ext);
    RooArgSet vars(gMass);

    // 2. Create an empty dataset and fill it manually from both branches
    RooDataSet gData("data", "dataset from tree", vars);

    TTreeReader reader(events);
    TTreeReaderValue<double> massReco(reader, "MassReco");
    TTreeReaderValue<double> massRecoCC(reader, "MassReco_CC");
    //TTreeReaderValue<double> cosTheta(reader, "CosTheta");
    //TTreeReaderValue<double> cosThetaCC(reader, "CosTheta_CC");
    double cosLo = 0.5; double cosHi = 0.6;

    while (reader.Next()) {
        //if (*massReco > xmin && *massReco < xmax &&
        //    std::abs(*cosTheta) > cosLo && std::abs(*cosTheta) < cosHi) {
        if (*massReco > xmin && *massReco < xmax) {
            gMass.setVal(*massReco);
            //gMass.setVal(*massReco - centre);
            gData.add(vars);
        }
        //if (*massRecoCC > xmin && *massRecoCC < xmax &&
        //    std::abs(*cosThetaCC) > cosLo && std::abs(*cosThetaCC) < cosHi) {
        if (*massRecoCC > xmin && *massRecoCC < xmax) {
            gMass.setVal(*massRecoCC);
            //gMass.setVal(*massRecoCC - centre);
            gData.add(vars);
        }
    }

    std::cout << "Entries in dataset: " << gData.numEntries() << std::endl;
    gMass.setRange("fitRange", xmin, xmax);
    gMass.setRange("plotRange", xmin, xmax);
    //gMass.setRange("fitRange", -ext, ext);
    //gMass.setRange("plotRange", -ext, ext);

    // 3. Build signal mode (Double-sided Crystal Ball)
    RooRealVar mean("mean","mean", 1.7768199, 1.75, 1.85);
    //RooRealVar mean("mean","mean", 0, -0.03, 0.03);
    RooRealVar sigma("sigma","sigma", 0.007, 1e-5, 0.02);
    RooRealVar alphaL("alphaL","alphaL", 0.5, 0.1, 10.0);
    RooRealVar nL("nL","nL", 1.0, 0.1, 10.0);
    RooRealVar alphaR("alphaR","alphaR", 0.5, 0.1, 10.0);
    RooRealVar nR("nR","nR", 1.0, 0.1, 10.0);

    RooCrystalBall sig("sig","Signal function",gMass, mean, sigma,
                        alphaL, nL, alphaR, nR);
    
    // 4. Perform the Fit
    auto fitmodel = sig.fitTo(gData, RooFit::Strategy(2), RooFit::Minos(false),
                              RooFit::Range("fitRange"), RooFit::Save());
    fitmodel->Print();
    std::cout << std::fixed << std::setprecision(7);
    std::cout << "Mean: " << mean.getVal() << std::endl;
    std::cout << "Stat. uncertainty: " << mean.getError() << std::endl;

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
        pad1->SetBottomMargin(0.03); pad1->SetLeftMargin(0.11); pad1->SetRightMargin(0.05);
        pad2->SetTopMargin(0.02); pad2->SetLeftMargin(0.11); pad2->SetRightMargin(0.05);
        pad2->SetBottomMargin(0.4);
        pad1->Draw(); pad2->Draw();

        // Main frame
        RooPlot* frame = gMass.frame(RooFit::Title(" "));
        frame->GetYaxis()->SetTitleOffset(0.9); // push label away from axis
        frame->GetXaxis()->SetLabelSize(0);
        frame->GetXaxis()->SetTitleSize(0);
        //frame->GetXaxis()->SetTitle("MassReco [GeV]");
        frame->GetYaxis()->SetLabelSize(fs);
        frame->GetYaxis()->SetTitleSize(fs);

        gData.plotOn(frame, RooFit::Name("data"));
        sig.plotOn(frame, RooFit::Name("model"), RooFit::LineColor(kRed), RooFit::LineWidth(2));

        ////double chi2ndf = frame->chiSquare(fitmodel->floatParsFinal().getSize());
        TPaveText* stats = new TPaveText(0.65, 0.75, 0.85, 0.85, "NDC");
        stats->SetFillColor(0); stats->SetBorderSize(1);
        stats->SetTextFont(42); stats->SetTextSize(0.04);
        ////stats->AddText(Form("#chi^{2}/ndf = %.3f", chi2ndf));

        RooHist* pullHist = frame->pullHist();
        RooPlot* pullFrame = gMass.frame(RooFit::Title(" "));
        pullFrame->addPlotable(pullHist, "P");
        pullFrame->GetYaxis()->SetTitle("Pulls");
        pullFrame->GetYaxis()->SetNdivisions(5);

        pullFrame->GetXaxis()->SetTitle("Reconstructed m_{#tau} [GeV/c^{2}]");
        pullFrame->GetXaxis()->SetLabelSize(fs * scale);
        pullFrame->GetXaxis()->SetTitleSize(fs * scale);
        pullFrame->GetYaxis()->SetLabelSize(fs * scale);
        pullFrame->GetYaxis()->SetTitleSize(fs * scale);
        pullFrame->GetYaxis()->SetTitleOffset(0.9 / scale);

        pad1->cd(); frame->drawBefore("data", "model"); frame->Draw(); ////stats->Draw();
        pad2->cd(); pullFrame->Draw();
        
        gSystem->mkdir("plots", kTRUE);
        c.SaveAs("plots/mass_plot_MC_alt.pdf");
    }
    else {
        TCanvas c("c", "c", 900, 600);
        
        TPad* pad1 = new TPad("pad1", "", 0, 0, 1, 1);
        pad1->SetBottomMargin(0.15);
        pad1->SetLeftMargin(0.14);
        pad1->Draw();
        
        RooPlot* frame = gMass.frame(RooFit::Title(" "));

        frame->GetYaxis()->SetTitleOffset(1.4);
        frame->GetXaxis()->SetTitleOffset(1.4);
        frame->GetXaxis()->SetTitle("Reconstructed m_{#tau} [GeV/c^{2}]");
        frame->GetXaxis()->SetLabelSize(0.05);
        frame->GetXaxis()->SetTitleSize(0.05);
        frame->GetYaxis()->SetLabelSize(0.05);
        frame->GetYaxis()->SetTitleSize(0.05);

        gData.plotOn(frame, RooFit::Name("data"));
        //model.plotOn(frame, RooFit::Name("model"), RooFit::LineColor(kRed), RooFit::LineWidth(2));

        pad1->cd(); frame->drawBefore("data", "model"); frame->Draw();

        gSystem->mkdir("plots", kTRUE);
        c.SaveAs("plots/mass_plot_raw_alt.pdf");
    }
    
    bool lhcb_plot_style = false;
    if (lhcb_plot_style == true) {
        lhcbStyle();

        TCanvas c("c", "c", 600, 600);
        RooPlot* frame = gMass.frame(RooFit::Title(""), RooFit::Bins(100));
        double binWidth = (gMass.getMax() - gMass.getMin()) / frame->GetNbinsX();
        std::ostringstream yTitle; yTitle << "Events / (" << binWidth << " GeV/c^{2})";
        frame->GetXaxis()->SetTitle("m_{reco} [GeV/c^{2}]");
        frame->GetYaxis()->SetTitle(yTitle.str().c_str());
        frame->GetYaxis()->SetTitleOffset(1.2);
        c.SetLeftMargin(0.17);
        gData.plotOn(frame, RooFit::Name("data"));
        sig.plotOn(frame,
                   RooFit::Name("signal"),
                   RooFit::LineColor(kRed),
                   RooFit::LineWidth(2));
        frame->drawBefore("data", "signal"); frame->Draw();
        TLegend* leg = new TLegend(0.62, 0.74, 0.92, 0.88);
        leg->AddEntry(frame->findObject("data"), "Simulation", "EP"); // E=error bar, P=marker
        leg->AddEntry(frame->findObject("signal"), "Fit", "L");  // L=line
        leg->Draw();

        // LHCb label — use the function from the macro
        // Options: "LHCb", "LHCb Preliminary", "LHCb Unofficial", "LHCb Simulation"
        auto lhcbLabel = LHCbStyle::create_label("FCC-ee");
        lhcbLabel->Draw();
        lhcbLabel->SetFillStyle(0); // Gets rid off an annoying black rectangle
        
        gSystem->mkdir("plots", kTRUE);
        c.SaveAs("plots/mass_plot_fit_alt.pdf");
    }
}