#if !defined(CLING) || defined(ROOTCLING)

#include <iostream>

#include "../utils/ClusterStudyUtils.h"

#include "CommonDataFormat/RangeReference.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCTrack.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "ITSBase/GeometryTGeo.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TrkClusRef.h"

#include "ITStracking/IOUtils.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"

#include "MFTAssessment/MFTAssessment.h"
#include "Framework/InputSpec.h"
#include "MFTBase/GeometryTGeo.h"

#include <gsl/gsl>
#include "TSystemDirectory.h"
#include <TLorentzVector.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "CommonDataFormat/RangeReference.h"
#include "DetectorsVertexing/DCAFitterN.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CCDBTimeStampUtils.h"

#endif

using GIndex = o2::dataformats::VtxTrackIndex;
using V0 = o2::dataformats::V0;
using MCTrack = o2::MCTrack;
using Cascade = o2::dataformats::Cascade;
using RRef = o2::dataformats::RangeReference<int, int>;
using VBracket = o2::math_utils::Bracket<int>;
using namespace o2::itsmft;
using CompClusterExt = o2::itsmft::CompClusterExt;
using ITSCluster = o2::BaseCluster<float>;
using Vec3 = ROOT::Math::SVector<double, 3>;
using MCLabCont = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
//using mMFTMapping = o2::itsmft::ChipMappingMFT;
o2::itsmft::ChipMappingMFT mMFTMapping;
o2::itsmft::ChipMappingMFT mMFTChipMapper;

void mcOrigin(std::string inPath="/home/lmichele/alice/mft_sim_full/myDir/", std::string outLabel="_myFirstTest", bool isOldData=true, unsigned int pix_thr = 2, bool verbose=false)
{
    /*
     - outLabel: label added to final root file 
     - isOldData: bool to adopt old/new dictionary
     - verbose: allow additional print
    */
    // "------------------ GLOBAL info ------------------"
    o2::base::GeometryManager::loadGeometry("../utils/o2_geometry.root");
    auto gman = o2::its::GeometryTGeo::Instance();
    gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G));
    std::vector<int> nStaves{12, 16, 20, 24, 30, 42, 48};
    std::vector<double> deltaEta{0.3, 0.44, 0.6, 0.72, 0.77, 0.72, 0.6, 0.44, 0.3};
    std::vector<int> PDGcodeOutsider;

    auto outFile = TFile(Form("outFileMCid_%s.root", outLabel.data()), "recreate");
    TTree *MCtree = new TTree("MCtree", "MCtree");

    float p, eta, phi, start_coord_x, start_coord_y, start_coord_z, E;
    int PDGID, CLsize, ProcessID, Layer;

    MCtree->Branch("CLsize", &CLsize);
    MCtree->Branch("p", &p);
    MCtree->Branch("phi", &phi);
    MCtree->Branch("eta", &eta);
    MCtree->Branch("X", &start_coord_x);
    MCtree->Branch("Y", &start_coord_y);
    MCtree->Branch("Z", &start_coord_z);
    MCtree->Branch("E", &E);
    MCtree->Branch("PDGID", &PDGID);
    MCtree->Branch("ProcessID", &ProcessID);
    MCtree->Branch("Layer", &Layer);

    TH1D *hCLid = new TH1D("hCLid", ";MC Track ID; entries", 8000, -0.5, 7999.5);
    TH1D *hCLsizeAll = new TH1D("hCLsizeAll", ";CL size; entries", 500, -0.5, 499.5);
    TH1D *hCLsizeTracks = new TH1D("hCLsizeTracks", ";CL size; entries", 500, -0.5, 499.5);
    TH1D *hTrackCLid = new TH1D("hTrackCLid", ";MC Track ID; entries", 20000, -0.5, 19999.5);
    TH1D *hCLsizeNuclTrk = new TH1D("hCLsizeNuclTrk", ";CL size; entries", 500, -0.5, 499.5);
    TH1D *hCLsizePionTrk = new TH1D("hCLsizePionTrk", ";CL size; entries", 500, -0.5, 499.5);
    TH1F *hCLsizeMean = new TH1F("hCLsizeMean", ";CL size; entries", 500, -0.5, 499.5);
    TH1D *hPDGcode = new TH1D("hPDGcode", ";PDG code; entries", 20000, -0.5, 19999.5);
    TH1D *hPDGcodeOut = new TH1D("hPDGcodeOut", ";PDG code; entries", 20000, -0.5, 19999.5);
    TH1D *hPDGmass = new TH1D("hPDGmass", ";PDG mass (); entries", 10000, -0.5, 99.5);
    TH1D *hPDGeloss = new TH1D("hPDGeloss", ";PDG dE/dx (); entries", 10000, -0.5, 99.5);
    TH1D *hPDGp = new TH1D("hPDGp", ";PDG p (); entries", 50000, -0.5, 49.5);
    hCLid->SetDirectory(nullptr);
    hTrackCLid->SetDirectory(nullptr);
    hPDGmass->SetDirectory(nullptr);
    hPDGeloss->SetDirectory(nullptr);
    hPDGp->SetDirectory(nullptr);
    hPDGcode->SetDirectory(nullptr);
    hCLsizeAll->SetDirectory(nullptr);
    hCLsizeNuclTrk->SetDirectory(nullptr);
    hCLsizePionTrk->SetDirectory(nullptr);
    hCLsizeMean->SetDirectory(nullptr);

    LOG(info) << "------------------ LOADING INPUT FILES ------------------";
    // Topology dictionary
    if (verbose){
        LOG(info) << "Loading topology dictionary";
    }
    o2::itsmft::TopologyDictionary mdict;
    o2::itsmft::ChipMappingITS chipMapping;
    if (isOldData)
    {
        LOG(info) << "Loading OLD dictionary: if you are analysing data older than JUNE should be fine";
        mdict.readFromFile(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(o2::detectors::DetID::MFT, "../utils/ITS"));
    }
    else
    {
        LOG(info) << "Loading LATEST dictionary: if you are analysing data older than JUNE check out the dictionary";
        auto f = TFile("../utils/o2_itsmft_TopologyDictionary_1653153873993.root");
        mdict = *(reinterpret_cast<o2::itsmft::TopologyDictionary *>(f.Get("ccdb_object")));
    }

    // Define the PB input file
    if (verbose)
    {
        LOG(info) << "Loading PB data file from " << inPath;
    }
    TSystemDirectory dir("MyDir", inPath.data());
    auto files = dir.GetListOfFiles();
    std::vector<std::string> dirs;
    for (auto fileObj : *files)
    {
        std::string file = ((TSystemFile *)fileObj)->GetName();
        if (verbose)
        {
            LOG(info) << "Keeping " << file;
        }
        dirs.push_back(inPath+file);
    }
    std::sort(dirs.begin(), dirs.end());

    std::vector<std::string> fulldirs;
    for (auto &dir : dirs)
    {
        TSystemDirectory subdir("MyDir2", dir.data());
        auto files = subdir.GetListOfFiles();
        for (auto fileObj : *files)
        {
            std::string file = ((TSystemFile *)fileObj)->GetName();
            if (file.substr(0, 2) == "tf")
            {
                if (verbose)
                {
                    LOG(info) << "Keeping " << file;
                }
                fulldirs.push_back(dir + "/" + file);
            }
        }
    }
    std::sort(fulldirs.begin(), fulldirs.end());
    int counter = 0;
    for (auto &dir : fulldirs)
    {
        if (true)
        {
            //if (counter > 10)
            //{
                //continue;
            //}
            counter ++;
        }
        //counter ++;
        LOG(info) << "Analysing directory: " << dir <<"  counter:"<<counter;
        std::string o2clus_its_file = dir + "/" + "mftclusters.root";
        std::string o2trac_its_file = dir + "/" + "mfttracks.root";
        LOG(info) << dir.substr(44, dir.size()).data();
        //std::string o2kine_file = dir + "/" + Form("sgn_%s_Kine.root", dir.substr(66, dir.size()).data());
        //std::string o2kine_file = dir + "/" + Form("sgn_%s_Kine.root", dir.substr(1, dir.size()).data());
        std::string o2kine_file = dir + "/" + Form("sgn_%s_Kine.root", dir.substr(44, dir.size()).data());
        //std::string o2kine_file = dir + "/" + "sgn_10_Kine.root";

        auto fITSclus = TFile::Open(o2clus_its_file.data());
        auto fITStrac = TFile::Open(o2trac_its_file.data());
        auto fMCTracks = TFile::Open(o2kine_file.data());

        if ( !fITSclus || !fITStrac || !fMCTracks)
        {
            LOG(info) << "SKIPPING: missing file!";
            continue;
        }

        auto treeITSclus = (TTree *)fITSclus->Get("o2sim");
        auto treeITStrac = (TTree *)fITStrac->Get("o2sim");
        auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");

        std::vector<CompClusterExt> *ITSclus = nullptr;
        std::vector<unsigned char> *ITSpatt = nullptr;
        o2::dataformats::MCTruthContainer<o2::MCCompLabel> *clusLabArr = nullptr;
        std::vector<int> *ITSTrackClusIdx = nullptr;
        std::vector<o2::mft::TrackMFT> *ITStracks = nullptr;
        std::vector<o2::MCTrack> *MCtracks = nullptr;
        std::vector<o2::BaseCluster<float>> mMFTClustersGlobal;
        //mMFTClustersGlobal.clear();
        //mMFTClustersGlobal.reserve(mMFTClusters.size());

        treeITSclus->SetBranchAddress("MFTClusterComp", &ITSclus);
        treeITSclus->SetBranchAddress("MFTClusterPatt", &ITSpatt);
        treeITSclus->SetBranchAddress("MFTClusterMCTruth", &clusLabArr);
        treeITStrac->SetBranchAddress("MFTTrack", &ITStracks);
        treeITStrac->SetBranchAddress("MFTTrackClusIdx", &ITSTrackClusIdx);
        treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);

        std::vector<int> LargeCLTrackID;
        std::vector<int> LargeCLEvID;
        std::vector<int> CLsiezes;
        std::vector<int> LayerID;
        for (int frame = 0; frame < treeITSclus->GetEntriesFast(); frame++)
        {   // LOOP OVER FRAMES  
            if (!treeITSclus->GetEvent(frame) || !treeITStrac->GetEvent(frame))
            {
                if(verbose)
                {
                    LOG(info) << "Skipping frame: " << frame;
                }
                continue;
            }

            if(verbose)
            {
                LOG(info) << "Frame: " << frame;
            }
            std::vector<o2::itsmft::ClusterPattern> pattVec;
            getClusterPatterns(pattVec, ITSclus, ITSpatt, mdict, gman);
            for (unsigned int iClus{0}; iClus < ITSclus->size(); iClus++)
            {   // LOOP OVER CLUSTERS
                auto &patt = pattVec[iClus];
                auto &clus = ITSclus->at(iClus);
                auto chipID = clus.getChipID();
                auto clsLayer = mMFTChipMapper.chip2Layer(chipID); //WARNING! With this I get 10 layers, which is the position of each layer? Print it
                int layer, sta, ssta, mod, chipInMod;
                auto pattID = clus.getPatternID();
                layer = gman->getLayer(clus.getSensorID());
                int npix = patt.getNPixels();
                hCLsizeAll->Fill(npix);
                if (npix > pix_thr)
                {
                    auto &labCls = (clusLabArr->getLabels(iClus))[0];
                    int  trackID, evID, srcID;
                    bool fake;
                    labCls.get(trackID, evID, srcID, fake);
                    if (verbose)
                    {
                        LOG(info) << "(NPIX="<<npix<<") Labels info: trackID="<<trackID<<", eventID="<<evID<<", srcID="<<srcID;
                    }
                    if (!labCls.isNoise() && labCls.isValid() && labCls.isCorrect() && !labCls.isFake())
                    {
                        LargeCLTrackID.push_back(trackID);
                        LargeCLEvID.push_back(evID);
                        hCLid->Fill(trackID);
                        CLsiezes.push_back(npix);
                        LayerID.push_back(clsLayer);
                    }
                }

            }
        }
            o2::mft::TrackMFT ITStrack;
            std::vector<o2::itsmft::ClusterPattern> pattVec;
            getClusterPatterns(pattVec, ITSclus, ITSpatt, mdict, gman);
            double meanCLsize = 0;
            bool allCLsOk = true;
            for (unsigned int iTrack{0}; iTrack < ITStracks->size(); iTrack++)
            {   // LOOP OVER TRACKS
                if (iTrack%10 == 0 && verbose)
                {
                    LOG(info) << "iTrack: " << iTrack;
                }

                auto &patt = pattVec[iTrack];
                ITStrack = (*ITStracks)[iTrack];
                //auto firstClus = ITStrack.getFirstClusterEntry();
                auto firstClus = ITStrack.getExternalClusterIndexOffset();
                //auto ncl = ITStrack.getNumberOfClusters();
                auto ncl = ITStrack.getNumberOfPoints();
                //LOG(info) << "CL offset = " << firstClus << "; N. CLs = " << ncl;

                meanCLsize = 0;
                std::vector <int> vectCLsize;
                allCLsOk = true;
                for (int icl = 0; icl < ncl; icl++)
                {   // LOOP OVER CLUSTERS
                    //LOG(info) << "LOOP OVER CLUSTERS";
                    //LOG(info) << "---------> " << (*ITSTrackClusIdx)[firstClus + icl] << " ITSClus size = " << (*ITSclus).size() ;
                    //auto &clus = (*ITSclus)[(*ITSTrackClusIdx)[firstClus + icl]];
                    auto &patt = pattVec[(*ITSTrackClusIdx)[firstClus + icl]];

                    ////////////////////
                    //auto clsEntry = (*ITSTrackClusIdx)[firstClus + icl];
                    //auto globalCluster = mMFTClustersGlobal[clsEntry];
                    //auto layer1 = mMFTMapping.ChipID2Layer[globalCluster.getSensorID()];
                    //LOG(info) << "LAYER1 =" << layer1;
                    //LayerID.push_back(layer1);
                    ////////////////////

                    //if ((*ITSTrackClusIdx)[firstClus + icl] > pattVec.size()) {
                        //LOG(info) << "---------------------------> " << (*ITSTrackClusIdx)[firstClus + icl] << " " << pattVec.size();
                    //}
                    
                    int npix = patt.getNPixels();
                    vectCLsize.push_back(npix);
                    meanCLsize += npix;
                    //LOG(info) << "NUMBER of PIXELS = " << npix;
                    auto &clus = (*ITSclus)[(*ITSTrackClusIdx)[firstClus + icl]];
                    if (npix > pix_thr)
                    {
                        hCLsizeTracks->Fill(npix);

                        auto &labCls = (clusLabArr->getLabels(ITSTrackClusIdx->at(firstClus+icl)))[0];
                        int  trackID, evID, srcID;
                        bool fake;
                        labCls.get(trackID, evID, srcID, fake);
                        if (verbose)
                        {
                            LOG(info) << "Labels info: trackID="<<trackID<<", eventID="<<evID<<", srcID="<<srcID;
                        }
                        hTrackCLid->Fill(trackID);
                    } else {
                        allCLsOk = false;
                    }
                }
                if (ncl > 8 && allCLsOk) {
                    meanCLsize = meanCLsize / ncl;
                    hCLsizeMean->Fill(meanCLsize);
                    for (int i = 0;i < int (vectCLsize.size());i++) {
                        printf("%i ", vectCLsize.at(i));
                    }
                    printf("\n");
                }
                vectCLsize.clear();
            }

        if (LargeCLTrackID.size() == 0)
        {
            if(true)
            {
                LOG(info) <<"Skipping: no large cluster found!";
            }
            continue;
        }

        std::vector<std::vector<o2::MCTrack>> mcTracksMatrix;
        auto nev = treeMCTracks->GetEntriesFast();
        mcTracksMatrix.resize(nev);
        for (int n = 0; n < nev; n++)
        { // loop over MC events
            treeMCTracks->GetEvent(n);
            mcTracksMatrix[n].resize(MCtracks->size());
            if (verbose)
            {
                LOG(info) << "N MC ev.=" << nev <<", N MC tracks="<<MCtracks->size();
            }
            for (unsigned int mcI{0}; mcI < MCtracks->size(); ++mcI)
            {   // LOOP over MC tracks
                mcTracksMatrix[n][mcI] = MCtracks->at(mcI);
            }
        }

        //LOG(info) << "---- GETTING MC tracks info ----";
        for (int i = 0; i < LargeCLEvID.size(); i++)
        {
            
            auto evID = LargeCLEvID.at(i);
            auto trID = LargeCLTrackID.at(i);
            if (verbose)
            {
                LOG(info) <<"evID="<<evID<<", trID="<<trID;
            }
            auto trPDG = mcTracksMatrix[evID][trID].GetPdgCode();

            //if (trPDG == 1000020030) {
                //LOG(info) << "----------> NUCLEO TROVATO";
            //}

            //if (trPDG < 1000000000)
            //{
                //double mass = TDatabasePDG::Instance()->GetParticle(trPDG)->Mass();

            if (trPDG == 1000020030 || trPDG == 211) {
                if (mcTracksMatrix[evID][trID].GetP() > 1) {
                    continue;
                } else {
                    if (TMath::Abs(trPDG) == 1000020030) {hCLsizeNuclTrk->Fill(CLsiezes.at(i));}
                    if (TMath::Abs(trPDG) == 211) {hCLsizePionTrk->Fill(CLsiezes.at(i));}
                    CLsize = CLsiezes.at(i);
                    p = mcTracksMatrix[evID][trID].GetP();
                    phi = mcTracksMatrix[evID][trID].GetPhi();
                    eta = mcTracksMatrix[evID][trID].GetEta();
                    start_coord_x = mcTracksMatrix[evID][trID].GetStartVertexCoordinatesX();
                    start_coord_y = mcTracksMatrix[evID][trID].GetStartVertexCoordinatesY();
                    start_coord_z = mcTracksMatrix[evID][trID].GetStartVertexCoordinatesZ();
                    E = mcTracksMatrix[evID][trID].GetEnergy();
                    PDGID = mcTracksMatrix[evID][trID].GetPdgCode();
                    ProcessID = mcTracksMatrix[evID][trID].getProcess();
                    Layer = LayerID.at(i);
                    MCtree->Fill();
                }
            }

            //if (verbose)


            
            //{
            //    LOG(info) <<"m="<<mass<<", p="<<p<<", PDG="<<trPDG;
            //}
            //hPDGmass->Fill(mass);
            //}
            //else
            //{
            //    PDGcodeOutsider.push_back(trPDG);
            //    LOG(info) << "------------------ PDG OUTSIDER INFO ------------------";
            //    LOG(info) << "Tot PDG OUTSIDER"<<PDGcodeOutsider.size();
            //    LOG(info) << "PDG OUTSIDER"<<PDGcodeOutsider[PDGcodeOutsider.size()];
            //}
            hPDGp->Fill(p);
            hPDGcode->Fill(trPDG);
        }
    }

    LOG(info) << "------------------ SAVING OUTFILE ------------------";
    outFile.cd();
    hCLid->Write();
    hTrackCLid->Write();
    hCLsizeNuclTrk->Write();
    hCLsizePionTrk->Write();
    hCLsizeMean->Write();
    hPDGcode->Write();
    hPDGmass->Write();
    hPDGp->Write();
    hCLsizeAll->Write();
    hCLsizeTracks->Write();
    MCtree->Write();

    //LOG(info) << "------------------ PDG OUTSIDER INFO ------------------";
    //LOG(info) << "Tot PDG OUTSIDER"<<PDGcodeOutsider.size();
    //for (int i = 0; i < PDGcodeOutsider.size(); i++)
    //{
    //    LOG(info) << "PDG OUTSIDER"<<PDGcodeOutsider[i];   
    //}
    outFile.Close();

}
