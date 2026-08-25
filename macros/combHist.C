///////////////////////////////////////////////////////////////////////////////////
///// Future Circular Collider Software - Tau-lepton Mass Peak Reconstruction /////
///// Combined histogram macro - Used to combine the data of a branch and its /////
////  CC branch from one event into a single histogram                        /////
///////////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>

// No argument is required in order to call this macro
void combHist()
{
    // Once the root file with the simulation data has been produced,
    // the path of the resulting file should be indicated below
    TFile f("result/p8_ee_Ztautau_ecm91.root");
    // Name of the branch containing the events data. Typically "events"
    TTreeReader reader("events", &f);

    // Names of the branches to join. Names inside quotes correspond to
    // the branch names from the rootfile
    TTreeReaderValue<double> MassReco(reader, "MassReco");
    TTreeReaderValue<double> MassReco_CC(reader, "MassReco_CC");
    TTreeReaderValue<double> MassReco_truem(reader, "MassReco_truem");
    TTreeReaderValue<double> MassReco_truem_CC(reader, "MassReco_truem_CC");

    double centre = 1.77693; double range = 0.4;
    TH1D h("hMass", "Tau mass distribution;Mass [GeV];Events", 150, centre - range, centre + range);
    TH1D h2("hMass2", "Tau mass distribution;Mass [GeV];Events", 150, centre - range, centre + range);

    while (reader.Next()) {
        h.Fill(*MassReco);
        h.Fill(*MassReco_CC);
        h2.Fill(*MassReco_truem);
        h2.Fill(*MassReco_truem_CC);
    }
    
    TCanvas c("c", "c", 900, 600);
    gStyle->SetOptStat(0);
    TPad* pad1 = new TPad("pad1", "", 0, 0, 1, 1);
    pad1->SetBottomMargin(0.14);
    pad1->SetLeftMargin(0.14);
    pad1->Draw(); pad1->cd();

    h.SetLineColor(kBlack); h.SetLineWidth(2);
    h2.SetLineColor(kRed); h2.SetLineWidth(2);

    h2.SetTitle(" ");

    h2.GetXaxis()->SetTitle("MassReco [GeV]");

    h2.GetXaxis()->SetTitleOffset(1.4); h2.GetYaxis()->SetTitleOffset(1.0);
    h2.GetXaxis()->SetLabelSize(0.05); h2.GetXaxis()->SetTitleSize(0.05);
    h2.GetYaxis()->SetLabelSize(0.05); h2.GetYaxis()->SetTitleSize(0.05);
    
    h2.Draw("HIST"); h.Draw("HIST SAME");
    TLegend* leg = new TLegend(0.58, 0.75, 0.81, 0.88);

    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.05);
    leg->AddEntry(&h,  "Reco. momentum", "l");
    leg->AddEntry(&h2, "Truth momentum", "l");
    leg->Draw();

    gSystem->mkdir("plots", kTRUE);
    c.SaveAs("plots/mass.pdf");
}