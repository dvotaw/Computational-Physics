//  File: roundoff_plotter.C
//
//  This is a macro to plot the results of roundoff.x in ROOT, rather than GNUPlot.
//
//  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
//
//  Revision history:
//      28-Jan-2019  Written from scratch.
//
//  Notes:  
//   * Does not need to be compiled. Run in ROOT, from command line:
//
//     	root -l roundoff_plotter.C
//
//*********************************************************************//

void roundoff_plotter(void)
{
  TCanvas* c = new TCanvas("c", "c");
  c->cd();

  TGraph* g = new TGraph("roundoff_output.txt");
  g->SetTitle("Roundoff error demonstration");
  g->GetXaxis()->SetTitle("log_{10}(N)");
  g->GetXaxis()->CenterTitle();
  g->GetYaxis()->SetTitle("log_{10}(#epsilon_{rel})");
  g->GetYaxis()->CenterTitle();
  g->Draw("AC*");
  g->SetMarkerStyle(20);
  g->SetLineColor(kBlue);

  g->Fit("pol1", "", "", 4, 7);

  c->Print("roundoff.pdf", "pdf");
}
