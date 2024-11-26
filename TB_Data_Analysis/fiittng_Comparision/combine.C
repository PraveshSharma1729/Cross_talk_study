void combine()
{

for(int k=1;k<=12;k++){

    // Open the ROOT files
    TFile* inFile1 = TFile::Open("/home/pravesh/Desktop/Cross_Talk/Run1695564190/Run1695564190.root");
    TFile* inFile2 = TFile::Open("/home/pravesh/Desktop/Cross_Talk/My_Cross_talk_simulation_0/100GeV_Pion_without_abs_100K_out0.root");

    // Retrieve the histograms from each file
    TH1F* hist1 = (TH1F*)inFile1->Get(Form("x2_%i_L2",k));
    TH1F* hist2 = (TH1F*)inFile2->Get(Form("x2_%i_L2",k));
    
    // Set statbox options
    gStyle->SetOptStat(1111); // Set to display all stats info

    // Scale the histograms
    hist1->Scale(1.0 / hist1->Integral());
    hist2->Scale(1.0 / hist2->Integral());

    // Set the colors of the histograms
    hist1->SetLineColor(kBlue); // Set hist1 color to blue
    hist2->SetLineColor(kRed);  // Set hist2 color to red

    // Create canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    // Set the x-axis range for both histograms
    hist1->GetXaxis()->SetRangeUser(-0.2, 0.2);
    hist2->GetXaxis()->SetRangeUser(-0.2, 0.2);

    // Draw first histogram without error bars and connected lines
    hist1->Draw("HIST");
    // Set y-axis limits
    hist1->SetMaximum(0.015); // Set the upper limit of the y-axis
    hist1->SetMinimum(0);
    c1->Update();

    // Retrieve and adjust statbox for hist1
    TPaveStats *st1 = (TPaveStats*)hist1->FindObject("stats");
    st1->SetX1NDC(0.1); // New x start position
    st1->SetY1NDC(0.7); // New y start position
    st1->SetX2NDC(0.3); // New x end position
    st1->SetY2NDC(0.9); // New y end position

    // Draw the second histogram on the main canvas without error bars and connected lines
    hist2->Draw("HIST SAME");
    c1->Update();

    // Retrieve and adjust statbox for hist2
    TPaveStats *st2 = (TPaveStats*)hist2->FindObject("stats");
    if (st2) {
        st2->SetX1NDC(0.7); // New x start position
        st2->SetY1NDC(0.7); // New y start position
        st2->SetX2NDC(0.9); // New x end position
        st2->SetY2NDC(0.9); // New y end position
        st2->Draw();  // Ensure the stats box is redrawn
    }

    // Add legend to the canvas
    TLegend *legend = new TLegend(0.5, 0.7, 0.6, 0.8);
    legend->AddEntry(hist1, "Data", "l");
    legend->AddEntry(hist2, "MC", "l");
    legend->Draw();

    // Finalize
    c1->Modified();
    c1->Update();

    // Save the canvas as an image
    c1->SaveAs(Form("/home/pravesh/Desktop/Cross_Talk/My_Cross_talk_simulation_0/images/x2_%i_L2.jpg",k)); // Change the path and filename as needed
}


}
