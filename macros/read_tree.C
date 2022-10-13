void read_tree(const char *fileName = "outFileMCid__myFirstTest.root"){
    TFile *fIn = new TFile(fileName, "READ");

    //float p, eta, phi, start_coord_x, start_coord_y, start_coord_z, E;
    //int PDGID, CLsize, ProcessID, Layer;

    TTreeReader treeReader("MCtree", fIn);
    TTreeReaderValue<Int_t> CLsize(treeReader, "CLsize");
    TTreeReaderValue<Int_t> PDGID(treeReader, "PDGID");
    TTreeReaderValue<Int_t> Layer(treeReader, "Layer");

    TH1F *histCLsizeNuclei = new TH1F("histCLsizeNuclei", "", 500, -0.5, 499.5);
    histCLsizeNuclei -> SetLineColor(kRed);
    TH1F *histCLsizeNonNuclei = new TH1F("histCLsizeNonNuclei", "", 500, -0.5, 499.5);
    histCLsizeNonNuclei -> SetLineColor(kBlue);
    TH1F *histCLsizeAll = new TH1F("histCLsizeAll", "", 500, -0.5, 499.5);
    histCLsizeAll -> SetLineColor(kBlack);
    
    TH1F *histCLsize[10];
    for (int i = 0;i < 10;i++) {
        histCLsize[i] = new TH1F(Form("histCLsizeNuclei_%i", i), "", 500, -0.5, 499.5);
    }

    while (treeReader.Next()) {
        std::cout << *CLsize << " " << *PDGID << " " << *Layer << std::endl;
        histCLsizeAll -> Fill(*CLsize);
        if(TMath::Abs(*PDGID) == 1000020030) {
            histCLsizeNuclei -> Fill(*CLsize);
        } else {
            histCLsizeNonNuclei -> Fill(*CLsize);
        }
        histCLsize[*Layer] -> Fill(*CLsize);
    }

    TH1F *histCLsizeRatio = (TH1F*) histCLsizeNuclei -> Clone("histCLsizeRatio");
    histCLsizeRatio -> Divide(histCLsizeAll);
    histCLsizeRatio -> SetLineColor(kBlack);

    TCanvas *canvas = new TCanvas("canvas", "", 600, 600);
    histCLsizeNonNuclei -> Draw("H");
    histCLsizeNuclei -> Draw("Hsame");

    TCanvas *canvasRatio = new TCanvas("canvasRatio", "", 600, 600);
    histCLsizeRatio -> Draw("H");
    
    TLine *lineUnity = new TLine(-0.5, 1, 499.5, 1);
    lineUnity -> SetLineColor(kGray+1);
    lineUnity -> Draw("same");

    TCanvas *canvasCLperLayer = new TCanvas("canvasCLperLayer", "", 3000, 1200);
    canvasCLperLayer -> Divide(5, 2);
    for (int i = 0;i < 10;i++) {
        canvasCLperLayer -> cd(i+1);
        gPad -> SetLogy(1);
        histCLsize[i] -> Draw("H");
    }

    TFile *fOut = new TFile("AnalysisResults.root", "RECREATE");
    canvas -> Write();
    canvasRatio -> Write();
    canvasCLperLayer -> Write();
    fOut -> Close();
}