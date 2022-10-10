void read_tree(const char *fileName = "outFileMCid__myFirstTest.root"){
    TFile *fIn = new TFile(fileName, "READ");

    //float p, eta, phi, start_coord_x, start_coord_y, start_coord_z, E;
    //int PDGID, CLsize, ProcessID, Layer;

    TTreeReader treeReader("MCtree", fIn);
    TTreeReaderValue<Int_t> CLsize(treeReader, "CLsize");
    TTreeReaderValue<Int_t> PDGID(treeReader, "PDGID");
    TTreeReaderValue<Int_t> Layer(treeReader, "Layer");

    TH1I *histCLsizeNuclei = new TH1I("histCLsizeNuclei", "", 100, -0.5, 99.5);
    histCLsizeNuclei -> SetLineColor(kRed);
    TH1I *histCLsizeNonNuclei = new TH1I("histCLsizeNonNuclei", "", 100, -0.5, 99.5);
    histCLsizeNonNuclei -> SetLineColor(kBlue);
    
    TH1I *histCLsize[10];
    for (int i = 0;i < 10;i++) {
        histCLsize[i] = new TH1I(Form("histCLsizeNuclei_%i", i), "", 100, -0.5, 99.5);
    }

    while (treeReader.Next()) {
        std::cout << *CLsize << " " << *PDGID << " " << *Layer << std::endl;
        if(TMath::Abs(*PDGID) == 1000020030) {
            histCLsizeNuclei -> Fill(*CLsize);
        } else {
            histCLsizeNonNuclei -> Fill(*CLsize);
        }
        histCLsize[*Layer] -> Fill(*CLsize);
    }

    TCanvas *canvas = new TCanvas("canvas", "", 600, 600);
    histCLsizeNonNuclei -> Draw("H");
    histCLsizeNuclei -> Draw("Hsame");

    TCanvas *canvasCLperLayer = new TCanvas("canvasCLperLayer", "", 3000, 1200);
    canvasCLperLayer -> Divide(5, 2);
    for (int i = 0;i < 10;i++) {
        canvasCLperLayer -> cd(i+1);
        gPad -> SetLogy(1);
        histCLsize[i] -> Draw("H");
    }
}