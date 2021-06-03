void potenziali_di_arresto()
{
//=========Macro generated from canvas: cVlam/Vc0(f)
//=========  (Sat May  1 19:06:47 2021) by ROOT version6.08/06
   TCanvas *cVlam = new TCanvas("cVlam", "Vc0(f)",0,44,782,756);
   cVlam->Range(3.65093e+14,0.380625,9.249715e+14,1.762275);
   cVlam->SetFillColor(0);
   cVlam->SetBorderMode(0);
   cVlam->SetBorderSize(2);
   cVlam->SetGridx();
   cVlam->SetGridy();
   cVlam->SetTickx(1);
   cVlam->SetTicky(1);
   cVlam->SetFrameBorderMode(0);
   cVlam->SetFrameBorderMode(0);
   
   Double_t Graph0_fx1001[5] = {
   4.68135e+14,
   4.950376e+14,
   5.452896e+14,
   6.300582e+14,
   7.041672e+14};
   Double_t Graph0_fy1001[5] = {
   0.6383814,
   0.7729454,
   1.017435,
   1.297934,
   1.457137};
   Double_t Graph0_fex1001[5] = {
   9.908398e+12,
   1.178063e+13,
   3.985348e+13,
   3.473753e+13,
   1.274205e+14};
   Double_t Graph0_fey1001[5] = {
   0.026,
   0.14,
   0.1,
   0.08,
   0.076};
   TGraphErrors *gre = new TGraphErrors(5,Graph0_fx1001,Graph0_fy1001,Graph0_fex1001,Graph0_fey1001);
   gre->SetName("Graph0");
   gre->SetTitle("Graph");
   gre->SetFillStyle(1000);
   
   TH1F *Graph_Graph1001 = new TH1F("Graph_Graph1001","Graph",100,4.210809e+14,8.689836e+14);
   Graph_Graph1001->SetMinimum(0.51879);
   Graph_Graph1001->SetMaximum(1.62411);
   Graph_Graph1001->SetDirectory(0);
   Graph_Graph1001->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph1001->SetLineColor(ci);
   Graph_Graph1001->GetXaxis()->SetTitle("Frequenza [Hz]");
   Graph_Graph1001->GetXaxis()->SetLabelFont(42);
   Graph_Graph1001->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph1001->GetXaxis()->SetTitleSize(0.03);
   Graph_Graph1001->GetXaxis()->SetTitleOffset(1.33);
   Graph_Graph1001->GetXaxis()->SetTitleFont(42);
   Graph_Graph1001->GetYaxis()->SetTitle("Potenziale d'arresto [V]");
   Graph_Graph1001->GetYaxis()->SetLabelFont(42);
   Graph_Graph1001->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph1001->GetYaxis()->SetTitleSize(0.03);
   Graph_Graph1001->GetYaxis()->SetTitleOffset(1.51);
   Graph_Graph1001->GetYaxis()->SetTitleFont(42);
   Graph_Graph1001->GetZaxis()->SetLabelFont(42);
   Graph_Graph1001->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1001->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1001);
   
   
   TF1 *retta1002 = new TF1("retta","[0]*x+[1]",4.210809e+14,8.689836e+14);
   retta1002->SetFillColor(19);
   retta1002->SetFillStyle(0);
   retta1002->SetLineColor(2);
   retta1002->SetLineWidth(2);
   retta1002->SetChisquare(0.2009144);
   retta1002->SetNDF(3);
   retta1002->GetXaxis()->SetLabelFont(42);
   retta1002->GetXaxis()->SetLabelSize(0.035);
   retta1002->GetXaxis()->SetTitleSize(0.035);
   retta1002->GetXaxis()->SetTitleFont(42);
   retta1002->GetYaxis()->SetLabelFont(42);
   retta1002->GetYaxis()->SetLabelSize(0.035);
   retta1002->GetYaxis()->SetTitleSize(0.035);
   retta1002->GetYaxis()->SetTitleFont(42);
   retta1002->SetParameter(0,4.104607e-15);
   retta1002->SetParError(0,9.319914e-16);
   retta1002->SetParLimits(0,0,0);
   retta1002->SetParameter(1,-1.28163);
   retta1002->SetParError(1,0.4554966);
   retta1002->SetParLimits(1,0,0);
   gre->GetListOfFunctions()->Add(retta1002);
   gre->Draw("AP");
   
   TPaveText *pt = new TPaveText(0.4440134,0.9342,0.5485204,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->Draw();
   cVlam->Modified();
   cVlam->cd();
   cVlam->SetSelected(cVlam);
}
