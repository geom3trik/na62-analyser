#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <TChain.h>
#include <math.h>
#include <list>
#include "spectro2.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"


#include "TRecoGigaTrackerEvent.hh"
#include "TRecoSpectrometerEvent.hh"
using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

int COUNTNUMBER = 0;

spectro2::spectro2( Core::BaseAnalysis *ba ) : Analyzer( ba, "spectro" )
{   //SOMETHING ELSE
    RequestTree( "GigaTracker", new TRecoGigaTrackerEvent );
    RequestTree( "Spectrometer", new TRecoSpectrometerEvent );
    RequestTree( "LKr", new TRecoLKrEvent );
    RequestTree( "MUV3", new TRecoMUV3Event );
    RequestTree( "CEDAR", new TRecoCedarEvent );
    RequestTree( "LAV", new TRecoLAVEvent );
    RequestTree( "RICH", new TRecoRICHEvent );
    RequestTree( "SAC", new TRecoSACEvent );
    RequestTree( "CHANTI", new TRecoCHANTIEvent );
    RequestTree( "CEDAR", new TRecoCedarEvent  );
}

void spectro2::InitOutput()
{

}

void spectro2::CreateHist1D( TString name, TString title, int nbins, double low, double high )
{
    TH1D* h1 = new TH1D( name, title, nbins, low, high );
    BookHisto( h1 );
}

void spectro2::SetHistAxisLabels( TString name, TString xlabel, TString ylabel )
{
    TH1* h = fHisto.GetHisto( name );
    if( h != NULL )
    {
        h -> GetXaxis() -> SetTitle(xlabel);
        h -> GetYaxis() -> SetTitle(ylabel);
    }

}

void spectro2::InitHist()
{
    // Reconstructed momentum histogram
    CreateHist1D( "MomentumHist", "Reconstructed Total Momentum STRAW Before Magnet", NumberOfBins, 0, 80 );
    SetHistAxisLabels( "Reconstructed Total Momentum In STRAW Before Magnet","Momentum GeV/c","Number of Entries");

    // True momentum histogram plotted with same bins as reconstructed momentum
    CreateHist1D( "CompareTrueMomentumHist", " True Total Momentum STRAW Before Magnet", NumberOfBins, 0, 80 );
    SetHistAxisLabels("True Total Momentum STRAW Before Magnet","Momentum GeV/c","Number of Entries");

    // Reconstructed x momentum histogram
    CreateHist1D( "xMomentumHist", "Compare x Momentum", NumberOfBins, -0.3, 0.3 );
    SetHistAxisLabels("xMomentumHist","Momentum GeV","Number of Entries");


    // True x momentum histogram plotted with same bins as reconstructed x momentum
    CreateHist1D( "CompareTruexMomentumHist", "Compare True x Momentum", NumberOfBins, -0.3, 0.3 );
    SetHistAxisLabels("CompareTruexMomentumHist","Momentum GeV","Number of Entries");

    // Reconstructed y momentum histogram
    CreateHist1D( "yMomentumHist", "y Momentum", NumberOfBins, -0.3, 0.3 );
    SetHistAxisLabels("yMomentumHist","Momentum GeV","Number of Entries");

    // True y momentum histogram plotted with same bins as reconstructed y momentum
    CreateHist1D( "CompareTrueyMomentumHist", "Compare True y Momentum", NumberOfBins, -0.3, 0.3 );
    SetHistAxisLabels("CompareTrueyMomentumHist","Momentum GeV","Number of Entries");

    // Reconstructed z momentum histogram
    CreateHist1D( "zMomentumHist", "z Momentum", NumberOfBins, 0, 80 );
    SetHistAxisLabels("zMomentumHist","Momentum GeV","Number of Entries");

    // True z momentum histogram plotted with same bins as reconstructed z momentum
    CreateHist1D( "CompareTruezMomentumHist", "Compare True z Momentum", NumberOfBins, 0, 80 );
    SetHistAxisLabels("CompareTruezMomentumHist","Momentum GeV","Number of Entries");

    // Reconstructed transverse momentum histogram
    CreateHist1D( "TransverseMomentumHist", "Transverse Momentum", NumberOfBins, 0, 300 );
    SetHistAxisLabels("TransverseMomentumHist","Momentum MeV","Number of Entries");

    // True transverse momentum histogram plotted with same bins as reconstructed transverse momentum
    CreateHist1D( "CompareTrueTransverseMomentumHist", " Compare True Transverse Momentum", NumberOfBins, 0, 300 );
    SetHistAxisLabels("CompareTrueTransverseMomentumHist","Momentum MeV","Number of Entries");

    // Reconstructed azimuthal angle histogram
    CreateHist1D( "AzimuthalMomentumHist", "Azimuthal Angle", NumberOfBins, -pi, pi );
    SetHistAxisLabels("AzimuthalMomentumHist","Azimuthal Angle Radians","Number of Entries");

    // True azimuthal angle histogram with same bins as reconstructed azimuthal angle histogram
    CreateHist1D( "CompareTrueAzimuthalMomentumHist", "Compare True Azimuthal Angle", NumberOfBins, -pi, pi );
    SetHistAxisLabels("CompareTrueAzimuthalMomentumHist","Azimuthal Angle Radians","Number of Entries");

    // Reconstructed polar angle histogram
    CreateHist1D( "PolarMomentumHist", "Polar Angle", NumberOfBins, 0, 0.017 );
    SetHistAxisLabels("PolarMomentumHist","Polar Angle Radians","Number of Entries");

    // True polar angle histogram with same bins as reconstructed polar angle histogram
    CreateHist1D( "CompareTruePolarMomentumHist", "Compare True Polar Angle", NumberOfBins, 0, 0.017 );
    SetHistAxisLabels("CompareTruePolarMomentumHist","Polar Angle Radians","Number of Entries");

    TH2D* h7 = new TH2D( "EnergyVsAzimuthal", "Total Momentum and Azimuthal Angle", NumberOfBins, -pi, pi, NumberOfBins, 0, 80 );
    h7 -> GetXaxis() -> SetTitle( "Azimuthal Angle Radians" );
    h7 -> GetYaxis() -> SetTitle( "Momentum GeV" );
    BookHisto( h7 );

    TH2D* h107 = new TH2D( "CompareTrueEnergyVsAzimuthal", "Compare True Total Momentum and Azimuthal Angle", NumberOfBins, -pi, pi, NumberOfBins, 0, 80 );
    h107 -> GetXaxis() -> SetTitle( "Azimuthal Angle Radians" );
    h107 -> GetYaxis() -> SetTitle( "Momentum GeV" );
    BookHisto( h107 );

    TH2D* h8 = new TH2D( "EnergyVsPolar", "Total Momentum and polar Angle", NumberOfBins, 0, 0.017, NumberOfBins, 0, 80 );
    h8 -> GetXaxis() -> SetTitle( "Polar Angle Radians" );
    h8 -> GetYaxis() -> SetTitle( "Momentum GeV" );
    BookHisto( h8 );

    TH2D* h108 = new TH2D( "CompareTrueEnergyVsPolar", "Compare True Total Momentum and Polar Angle", NumberOfBins, 0, 0.017, NumberOfBins, 0, 80 );
    h108 -> GetXaxis() -> SetTitle( "Azimuthal Angle Radians" );
    h108 -> GetYaxis() -> SetTitle( "Momentum GeV" );
    BookHisto( h108 );

    TH2D* h9 = new TH2D( "TranverseEnergyVsAzimuthal", "Tranverse Momentum and Azimuthal Angle", NumberOfBins, -pi, pi, NumberOfBins, 0, 0.3 );
    h9 -> GetXaxis() -> SetTitle( "Azimuthal Angle Radians" );
    h9 -> GetYaxis() -> SetTitle( "Momentum GeV" );
    BookHisto( h9 );

    TH2D* h109 = new TH2D( "CompareTrueTranverseEnergyVsAzimuthal", "Compare Tranverse Momentum and Azimuthal Angle", NumberOfBins, -pi, pi, NumberOfBins, 0, 0.3 );
    h109 -> GetXaxis() -> SetTitle( "Azimuthal Angle Radians" );
    h109 -> GetYaxis() -> SetTitle( "Momentum GeV" );
    BookHisto( h109 );

    CreateHist1D( "TrueMomentumHist", " True Momentum", NumberOfBins, 0, 0 );
    SetHistAxisLabels("TrueMomentumHist","Momentum GeV","Number of Entries");

    CreateHist1D( "TruexMomentumHist", "True x Momentum", NumberOfBins, 0, 0 );
    SetHistAxisLabels("TruexMomentumHist","Momentum GeV","Number of Entries");

    CreateHist1D( "TrueyMomentumHist", "True y Momentum", NumberOfBins, 0, 0 );
    SetHistAxisLabels("TrueyMomentumHist","Momentum GeV","Number of Entries");

    CreateHist1D( "TruezMomentumHist", "True z Momentum", NumberOfBins, 0, 0 );
    SetHistAxisLabels("TruezMomentumHist","Momentum GeV","Number of Entries");

    CreateHist1D( "TrueTransverseMomentumHist", " True Transverse Momentum", NumberOfBins, 0, 0 );
    SetHistAxisLabels("TrueTransverseMomentumHist","Momentum MeV","Number of Entries");

    CreateHist1D( "TrueAzimuthalMomentumHist", " True Azimuthal Angle", NumberOfBins, 0, 0 );
    SetHistAxisLabels("TrueAzimuthalMomentumHist","Azimuthal Angle Radians","Number of Entries");

    CreateHist1D( "TruePolarMomentumHist", " True Polar Angle", NumberOfBins, 0, 0 );
    SetHistAxisLabels("TruePolarMomentumHist","Polar Angle Radians","Number of Entries");

    TH2D* h22 = new TH2D( "TrueEnergyVsAzimuthal", "True Total Momentum and Azimuthal Angle", NumberOfBins, 0, 0, NumberOfBins, 0, 0 );
    h22 -> GetXaxis() -> SetTitle( "Azimuthal Angle Radians" );
    h22 -> GetYaxis() -> SetTitle( "Momentum GeV" );
    BookHisto( h22 );

    TH2D* h23 = new TH2D( "TrueEnergyVsPolar", "True Total Momentum and Polar Angle", NumberOfBins, 0, 0, NumberOfBins, 0, 0 );
    h23 -> GetXaxis() -> SetTitle( "Azimuthal Angle Radians" );
    h23 -> GetYaxis() -> SetTitle( "Momentum GeV" );
    BookHisto( h23 );

    TH2D* h24 = new TH2D( "TrueTranverseEnergyVsAzimuthal", "Tranverse Momentum and Azimuthal Angle", NumberOfBins, 0, 0, NumberOfBins, 0, 0 );
    h24 -> GetXaxis() -> SetTitle( "Azimuthal Angle Radians" );
    h24 -> GetYaxis() -> SetTitle( "Momentum GeV" );
    BookHisto( h24 );

    CreateHist1D( "ClosestPointFromBeamAxis", "Closest Point From Beam Axis", NumberOfBins, 0, 300 );
    SetHistAxisLabels("ClosestPointFromBeamAxis","m","Number of Entries");

    CreateHist1D( "ClosestxPointFromBeamAxis", " Closest x Point From Beam Axis", NumberOfBins, -200, 200 );
    SetHistAxisLabels("ClosestxPointFromBeamAxis","mm","Number of Entries");

    CreateHist1D( "ClosestyPointFromBeamAxis", " Closest y Point From Beam Axis", NumberOfBins, -100, 100 );
    SetHistAxisLabels("ClosestyPointFromBeamAxis","mm","Number of Entries");

    CreateHist1D( "ClosestzPointFromBeamAxis", " Closest z Point From Beam Axis", NumberOfBins, 0, 300 );
    SetHistAxisLabels("ClosestzPointFromBeamAxis","m","Number of Entries");

    CreateHist1D( "MinimumDistanceToBeamAxis", " Minimum Distance From Beam Axis", NumberOfBins, -100, 100 );
    SetHistAxisLabels("MinimumDistanceToBeamAxis","mm","Number of Entries");

    CreateHist1D( "ClosestDistanceToBeamAxis", " Minimum Distance(TWO) From Beam Axis", NumberOfBins, 0, 100 );
    SetHistAxisLabels("ClosestDistanceToBeamAxis","mm","Number of Entries");

    CreateHist1D( "ClosestxDistanceToBeamAxis", " Minimum xDistance From Beam Axis", NumberOfBins, -100, 100 );
    SetHistAxisLabels("ClosestxDistanceToBeamAxis","mm","Number of Entries");

    CreateHist1D( "ClosestyDistanceToBeamAxis", " Minimum yDistance From Beam Axis", NumberOfBins, -100, 100 );
    SetHistAxisLabels("ClosestyDistanceToBeamAxis","mm","Number of Entries");

    CreateHist1D( "ClosestzDistanceToBeamAxis", " Minimum zDistance From Beam Axis", NumberOfBins, -0.1, 0.1 );
    SetHistAxisLabels("ClosestzDistanceToBeamAxis","mm","Number of Entries");

    TH2D* h34 = new TH2D( "DecayPoisition", " Decay Poisition ", NumberOfBins, 0, 270, NumberOfBins, -180 , 180 );
    h34 -> GetXaxis() -> SetTitle( "z metres" );
    h34 -> GetYaxis() -> SetTitle( "x mm" );
    BookHisto( h34 );

    CreateHist1D( "ClosestPointOfMuon", " ClosestPointOfMuon", NumberOfBins*10, 0, 0 );
    SetHistAxisLabels("ClosestPointOfMuon","m","Number of Entries");

    CreateHist1D( "ClosestxPointOfMuon", " ClosestxPointOfMuon", NumberOfBins*10, 0, 0 );
    SetHistAxisLabels("ClosestxPointOfMuon","mm","Number of Entries");

    CreateHist1D( "ClosestyPointOfMuon", " ClosestyPointOfMuon", NumberOfBins*10, 0, 0 );
    SetHistAxisLabels("ClosestyPointOfMuon","mm","Number of Entries");

    CreateHist1D( "ClosestzPointOfMuon", " ClosestzPointOfMuon", NumberOfBins*10, 0, 0 );
    SetHistAxisLabels("ClosestzPointOfMuon","m","Number of Entries");

    CreateHist1D( "ClosestTimeOfMuon", " ClosestTimeOfMuon", NumberOfBins*10, 0, 0 );
    SetHistAxisLabels("ClosestTimeOfMuon","s","Number of Entries");

    CreateHist1D( "ClosestTimeOfKaon", " ClosestTimeOfKaon", NumberOfBins*10, 0, 0 );
    SetHistAxisLabels("ClosestTimeOfKaon","s","Number of Entries");

    CreateHist1D( "ClosestDistanceOfMuonToKaon", " ClosestDistanceOfMuonToKaon", NumberOfBins*10, 0, 0 );
    SetHistAxisLabels("ClosestDistanceOfMuonToKaon","mm","Number of Entries");

    CreateHist1D( "ClosestxDistanceOfMuonToKaon", " ClosestxDistanceOfMuonToKaon", NumberOfBins*10, 0, 0 );
    SetHistAxisLabels("ClosestxDistanceOfMuonToKaon","mm","Number of Entries");

    CreateHist1D( "ClosestyDistanceOfMuonToKaon", " ClosestyDistanceOfMuonToKaon", NumberOfBins*10, 0, 0 );
    SetHistAxisLabels("ClosestyDistanceOfMuonToKaon","mm","Number of Entries");

    CreateHist1D( "ClosestzDistanceOfMuonToKaon", " ClosestzDistanceOfMuonToKaon", NumberOfBins*10, 0, 0 );
    SetHistAxisLabels("ClosestzDistanceOfMuonToKaon","mm","Number of Entries");

    CreateHist1D( "ClosestSpaceTimeInterval", " ClosestSpaceTimeInterval", NumberOfBins*10, 0, 0 );
    SetHistAxisLabels("ClosestSpaceTimeInterval","m^2","Number of Entries");

    TH2D* h46 = new TH2D( "ParticleProductionPosition", "StartingPositionOfJetParticles", NumberOfBins, 0, 0, NumberOfBins, 0, 0 );
    h46 -> GetXaxis() -> SetTitle( "m" );
    h46 -> GetYaxis() -> SetTitle( "mm" );
    BookHisto( h46 );

    TH2D* h47 = new TH2D( "KaonEndingPosition", "EndingPositionOfKaons", NumberOfBins, 0, 0, NumberOfBins, 0, 0 );
    h47 -> GetXaxis() -> SetTitle( "m" );
    h47 -> GetYaxis() -> SetTitle( "mm" );
    BookHisto( h47 );

    CreateHist1D( "TrueProductionPositionX", " MuonProductionPointx", NumberOfBins, 0, 0 );
    SetHistAxisLabels("TrueProductionPositionx","m","Number of Entries");

    CreateHist1D( "TrueProductionPositionY", " MuonProductionPointy", NumberOfBins, 0, 0 );
    SetHistAxisLabels("TrueProductionPositiony","m","Number of Entries");

    CreateHist1D( "TrueProductionPositionZ", " MuonProductionPointz", NumberOfBins, 0, 0 );
    SetHistAxisLabels("TrueProductionPositionZ","m","Number of Entries");

    TH2D* h51 = new TH2D( "TrueProductionPosition", "MuonProductionPoint", NumberOfBins, 0, 0, NumberOfBins, 0, 0 );
    h51 -> GetXaxis() -> SetTitle( "m" );
    h51 -> GetYaxis() -> SetTitle( "mm" );
    BookHisto( h51 );

    CreateHist1D("MissingMass", "Missing Mass Squared", NumberOfBins, -0.2, 0.2);
    SetHistAxisLabels("MissingMass","Missing Mass Squared, GeV /c","Number of Entries");

    CreateHist1D("TrueMissingMass", "True Missing Mass Squared", 1000, 0, 0);
    SetHistAxisLabels("TrueMissingMass","True Missing Mass Squared","Number of Entries");

    CreateHist1D("SpectrometerXMomentumResolution", "Spectrometer x Momentum Resolution", NumberOfBins, 0, 0);
    SetHistAxisLabels("SpectrometerXMomentumResolution","x Momentum MeV","Number of Entries");

    CreateHist1D("SpectrometerYMomentumResolution", "Spectrometer y Momentum Resolution", NumberOfBins, 0, 0);
    SetHistAxisLabels("SpectrometerYMomentumResolution","y Momentum MeV","Number of Entries");

    CreateHist1D("SpectrometerZMomentumResolution", "Spectrometer z Momentum Resolution", NumberOfBins, 0, 0);
    SetHistAxisLabels("SpectrometerZMomentumResolution","z Momentum GeV","Number of Entries");

    CreateHist1D("SpectrometerMomentumResolution", "Spectrometer Momentum Resolution", NumberOfBins, 0, 0);
    SetHistAxisLabels("SpectrometerMomentumResolution","Momentum GeV","Number of Entries");

    CreateHist1D("XZAngle", "XZAngle", NumberOfBins, 0, 0);
    SetHistAxisLabels("XZAngle","Radians","Number of Entries");

    CreateHist1D("YZAngle", "YZAngle", NumberOfBins, 0, 0);
    SetHistAxisLabels("YZAngle","Radians","Number of Entries");

    CreateHist1D("ReconstructedMomentumHist","Momentum",NumberOfBins,0,80);
    SetHistAxisLabels("ReconstructedMomentumHist","Momentum GeV","Number of Entries");

    CreateHist1D("ReconstructedxMomentumHist","xMomentum",NumberOfBins,-0.3,0.3);
    SetHistAxisLabels("ReconstructedxMomentumHist","Momentum GeV","Number of Entries");

    CreateHist1D("ReconstructedyMomentumHist","yMomentum",NumberOfBins,-0.3,0.3);
    SetHistAxisLabels("ReconstructedyMomentumHist","Momentum GeV","Number of Entries");


    CreateHist1D("ReconstructedzMomentumHist","zMomentum",NumberOfBins,0,80);
    SetHistAxisLabels("ReconstructedzMomentumHist","Momentum GeV","Number of Entries");

    CreateHist1D("TrueXZAngle", "TrueXZAngle", NumberOfBins, 0, 0);
    SetHistAxisLabels("TrueXZAngle","Radians","Number of Entries");

    CreateHist1D("TrueYZAngle", "TrueYZAngle", NumberOfBins, 0, 0);
    SetHistAxisLabels("TrueYZAngle","Radians","Number of Entries");

    CreateHist1D("cluster energy", "cluster energy", NumberOfBins, 0, 0);
    SetHistAxisLabels("Cluster Energy GeV","Counts","Number of Entries");

    CreateHist1D("Energyclusterenergy", "Energyclusterenergy", NumberOfBins, 0, 0);
    SetHistAxisLabels("Energy / ClusterEnergy","Counts","Number of Entries");

    CreateHist1D("TotalDecays", "TotalDecays", NumberOfBins, 0, 0);
    SetHistAxisLabels("TotalDecays","Counts","Number of Entries");

    CreateHist1D("ClusterEnergyAndIntersectionLength", "ClusterEnergyAndIntersectionLength", NumberOfBins, 0, 0);
    SetHistAxisLabels("ClusterEnergyAndIntersectionLength","Counts","Number of Entries");

    CreateHist1D("LkrEntryx", "LkrEntryx", NumberOfBins, 0, 0);
    SetHistAxisLabels("LkrEntryx","m","Number of Entries");

    CreateHist1D("LkrEntryy", "LkrEntryy", NumberOfBins, 0, 0);
    SetHistAxisLabels("LkrEntryy","m","Number of Entries");

    CreateHist1D("LkrEntryz", "LkrEntryz", NumberOfBins, 0, 0);
    SetHistAxisLabels("LkrEntryz","m","Number of Entries");

    CreateHist1D("closestxDistanceToMUV3", "closestxDistanceToMUV3", NumberOfBins, 0, 0);
    SetHistAxisLabels("closestxDistanceToMUV3","m","Number of Entries");

    CreateHist1D("closestyDistanceToMUV3", "closestyDistanceToMUV3", NumberOfBins, 0, 0);
    SetHistAxisLabels("closestyDistanceToMUV3","m","Number of Entries");

    CreateHist1D("closestzDistanceToMUV3", "closestzDistanceToMUV3", NumberOfBins, 0, 0);
    SetHistAxisLabels("closestzDistanceToMUV3","m","Number of Entries");

    CreateHist1D("closestDistanceToMUV3", "closestDistanceToMUV3", NumberOfBins, 0, 0);
    SetHistAxisLabels("closestDistanceToMUV3","m","Number of Entries");

    CreateHist1D("closestxPointToMUV3", "closestxPointToMUV3", NumberOfBins, 0, 0);
    SetHistAxisLabels("closestxPointToMUV3","m","Number of Entries");

    CreateHist1D("closestyPointToMUV3", "closestyPointToMUV3", NumberOfBins, 0, 0);
    SetHistAxisLabels("closestyPointToMUV3","m","Number of Entries");

    CreateHist1D("closestzPointToMUV3", "closestzPointToMUV3", NumberOfBins, 0, 0);
    SetHistAxisLabels("closestzPointToMUV3","m","Number of Entries");

    CreateHist1D("MUV3candidates", "MUV3candidates", NumberOfBins, 0, 0);
    SetHistAxisLabels("MUV3candidates","numberofcandidates","Number of Entries");

    CreateHist1D("timedifference", "timedifference", 5000, 0, 0);
    SetHistAxisLabels("timedifference","time","Number of Entries");

    CreateHist1D("spectrotime", "spectrotime", 5000, 0, 0);
    SetHistAxisLabels("spectrotime","time","Number of Entries");

    CreateHist1D("MUV3time", "MUV3time", 5000, 0 , 0 );
    SetHistAxisLabels("MUV3time","time","Number of Entries");


    CreateHist1D("LKrTimeDifference", "LKrTimeDifference", 5000, 0 , 0 );
    SetHistAxisLabels("LKrTimeDifference","time","Number of Entries");

    CreateHist1D("LKrcandidates", "LKrTimeDifference", 5000, 0 , 0 );
    SetHistAxisLabels("LKrTimeDifference","time","Number of Entries");

    CreateHist1D("LKrtime", "LKrtime", 5000, 0 , 0 );
    SetHistAxisLabels("LKrtime","time","Number of Entries");



    TH2D* h52 = new TH2D( "SpectrometerXMomentumResolutionAgainstXMomentum", "x Momentum Resolution Against x Momentum", NumberOfBins, 0, 0, NumberOfBins, 0, 0 );
    h52 -> GetXaxis() -> SetTitle( "MeV" );
    h52 -> GetYaxis() -> SetTitle( "MeV" );
    BookHisto( h52 );

    TH2D* h53 = new TH2D( "SpectrometerYMomentumResolutionAgainstYMomentum", "y Momentum Resolution Against y Momentum", NumberOfBins, 0, 0, NumberOfBins, 0, 0 );
    h53 -> GetXaxis() -> SetTitle( "MeV" );
    h53 -> GetYaxis() -> SetTitle( "MeV" );
    BookHisto( h53 );

    TH2D* h54 = new TH2D( "SpectrometerZMomentumResolutionAgainstZMomentum", "z Momentum Resolution Against z Momentum", NumberOfBins, 0, 0, NumberOfBins, 0, 0 );
    h54 -> GetXaxis() -> SetTitle( "GeV" );
    h54 -> GetYaxis() -> SetTitle( "GeV" );
    BookHisto( h54 );

    TH2D* h55 = new TH2D( "SpectrometerMomentumResolutionAgainstMomentum", "Momentum Resolution Against Momentum", NumberOfBins, 0, 0, NumberOfBins, 0, 0 );
    h55 -> GetXaxis() -> SetTitle( "GeV" );
    h55 -> GetYaxis() -> SetTitle( "GeV" );
    BookHisto( h55 );


    for(int i=0;i<13;i++)
    {

        string intstr = to_string(i);
        TString num = TString(intstr);
        CreateHist1D(TString("ResolutionTemp") + num, "Title", 500,0,0);
        CreateHist1D(TString("ResolutionTempX") + num, "Title", 500,0,0);
        CreateHist1D(TString("ResolutionTempY") + num, "Title", 500,0,0);
        CreateHist1D(TString("ResolutionTempZ") + num, "Title", 500,0,0);
        CreateHist1D(TString("ResolutionTempXZ") + num, "Title", 500,0,0);
        CreateHist1D(TString("ResolutionTempYZ") + num, "Title", 500,0,0);
    }

    TH2D* h56 = new TH2D( "MissingMassVsZPosition", "Missing Mass Vs Z position", NumberOfBins, 0, 0, NumberOfBins, 0, 0 );
    h56 -> GetXaxis() -> SetTitle( "GeV Squared" );
    h56 -> GetYaxis() -> SetTitle( "m" );
    BookHisto( h56 );

    TH2D* h57 = new TH2D( "MUV3Position", "MUV3Position", NumberOfBins, 0, 0, NumberOfBins, 0, 0 );
    h57 -> GetXaxis() -> SetTitle( "m" );
    h57 -> GetYaxis() -> SetTitle( "m" );
    BookHisto( h57 );

    TH2D* h66 = new TH2D("ClusterPosition","Cluster Position",NumberOfBins,0,0,NumberOfBins,0,0);
    h66 -> GetXaxis() -> SetTitle( "mm" );
    h66 -> GetYaxis() -> SetTitle( "mm" );
    BookHisto( h66 );

    TH2D* h67 = new TH2D("TrackPosition","Track Position",NumberOfBins,0,0,NumberOfBins,0,0);
    h67 -> GetXaxis() -> SetTitle( "mm" );
    h67 -> GetYaxis() -> SetTitle( "mm" );
    BookHisto( h67 );

    CreateHist1D("ClusterTrackDistanceX", "Distance Between Cluster X & Track X", NumberOfBins, 0 , 0 );
    SetHistAxisLabels("ClusterTrackDistance","mm","Number of Entries");

    CreateHist1D("ClusterTrackDistanceY", "Distance Between Cluster Y & Track Y", NumberOfBins, 0 , 0 );
    SetHistAxisLabels("ClusterTrackDistanceY","mm","Number of Entries");


    /*
    TH1D* h14 = new TH1D( "TrueMuonxMomentumHist", "True Muon x Momentum", NumberOfBins, 0, 0 );
    h14 -> GetXaxis() -> SetTitle( "Momentum MeV" );
    h14 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h14 );
   TH1D* h15 = new TH1D( "TrueMuonyMomentumHist", "True Muon y Momentum", NumberOfBins, 0, 0 );
    h15 -> GetXaxis() -> SetTitle( "Momentum MeV" );
    h15 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h15 );
   TH1D* h16 = new TH1D( "TrueMuonzMomentumHist", "True Muon z Momentum", NumberOfBins, 0, 0 );
    h16 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h16 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h16 );
    TH1D* h17 = new TH1D( "TrueMuonEnergyHist", "True Muon z Momentum", NumberOfBins, 0, 0 );
    h17 -> GetXaxis() -> SetTitle( "Energy GeV" );
    h17 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h17 );
*/
}

void spectro2::DefineMCSimple()
{

}

void spectro2::StartOfRunUser()
{

}

void spectro2::StartOfBurstUser()
{

}

void spectro2::SaveAllPlotsPDF()
{
    // Get iterator for CanvasOrganizer map
    NA62Analysis::NA62Map<TString,CanvasOrganizer*>::type::iterator ptr;
    // Get the CanvasOrganizer map
    NA62Analysis::NA62Map<TString,CanvasOrganizer*>::type canvases = GetCanvases();
    // Loop through all CanvasOrganizer map
    for( ptr = canvases.begin(); ptr != canvases.end(); ptr++ )
    {
        // Retrieve ROOT canvas and save as .pdf file
        ptr->second->GetCanvas()->SaveAs( TString( ptr->first + ".pdf" ) );
    }
}


TVector3 ClosestPointOnVectorToOtherVector( TVector3 Vector1Position, TVector3 Vector1Direction, TVector3 Vector2Position, TVector3 Vector2Direction ) //14th November
{
    double a,b,c,d,e, VectorScale1, VectorScale2;
    a = Vector1Direction.Dot( Vector1Direction );
    b = Vector1Direction.Dot( Vector2Direction );
    c = Vector2Direction.Dot( Vector2Direction );
    d = Vector1Direction.Dot( ( Vector1Position - Vector2Position ) );
    e = Vector2Direction.Dot( ( Vector1Position - Vector2Position ) );
    //f = ( Vector1Position - Vector2Position ).Dot( ( Vector1Position - Vector2Position ) );
    VectorScale1 = ( b * e - c * d ) / ( a * c - pow( b, 2 ) );
    VectorScale2 = ( a * e - b * d ) / ( a * c - pow( b, 2 ) );

    TVector3 ClosestPointOnVector1 = Vector1Position + VectorScale1 * Vector1Direction;
    TVector3 ClosestPointOnVector2 = Vector2Position + VectorScale2 * Vector2Direction;
    return ClosestPointOnVector1;
}

TVector3 ClosestPointOnVectorToPoint( TVector3 VectorPosition, TVector3 VectorDirection, TVector3 Point ) //30th jan
{
    double t;
    t = ( VectorDirection.Dot( Point ) - VectorDirection.Dot( VectorPosition ) ) / VectorDirection.Dot( VectorDirection );
    TVector3 ClosestPointOnVector = VectorPosition + VectorDirection * t;
    return ClosestPointOnVector;
}


TLorentzVector ClosestSpaceTimePointOnVectorToOtherVector( TVector3 Vector1Position, TVector3 Vector1Momentum, double Mass1, double Time1, TVector3 Vector2Position, TVector3 Vector2Momentum, double Mass2, double Time2 ) //14th November
{
    double b,d, VectorScale1, VectorScale2, ClosestTime1, ClosestTime2;
    TVector3  A, B, C, D;
    //MAY NEED TO SCALE TIME1 AND TIME2 TO GET THEM IN SECONDS. DON'T KNOW WHAT THEY OUTPUT.
    A = Vector1Position * pow( 1000., -1 );
    B = Vector1Momentum.Unit();
    C = Vector2Position * pow( 1000., -1 );
    D = Vector2Momentum.Unit();
    b = Vector1Momentum.Mag() * pow ( 10 , 6 ) * 5.3442883 * pow( 10, -28 );
    d = Vector2Momentum.Mag() * pow ( 10 , 6 ) * 5.3442883 * pow( 10, -28 );
    //Time1 = Time1;
    //Time2 = Time2;
    VectorScale1 = ( ( pow( d, 2 ) - pow( SpeedOfLight, 2 ) ) * ( pow( SpeedOfLight, 2 ) * ( Time2 - Time1 ) - ( C - A ).Dot( b * B ) ) + ( ( d * D ).Dot( b * B ) - pow( SpeedOfLight, 2 ) ) * ( ( C - A ).Dot( d * D ) - pow( SpeedOfLight, 2 ) * ( Time2 - Time1 ) ) ) / ( pow( SpeedOfLight, 2 ) * pow ( ( d - b ), 2 ) );
    VectorScale2 = ( ( pow( b, 2 ) - pow( SpeedOfLight, 2 ) ) * VectorScale1 - ( C - A ).Dot( b * B ) + pow( SpeedOfLight , 2 ) * ( Time2 - Time1 ) ) / ( (d * D ).Dot( b * B ) - pow( SpeedOfLight, 2 ) );
    TVector3 ClosestPointOnVector1 = A + VectorScale1 * B * b * pow( sqrt( pow( Mass1, 2 ) + pow( b, 2 ) / pow ( SpeedOfLight , 2 ) ), -1 );
    TVector3 ClosestPointOnVector2 = C + VectorScale2 * D * d * pow( sqrt( pow( Mass2, 2 ) + pow( d, 2 ) / pow ( SpeedOfLight , 2 ) ), -1 );
    ClosestTime1 = Time1 + VectorScale1;
    ClosestTime2 = Time2 + VectorScale2;
    TLorentzVector ClosestSpaceTimePointOnVector1, ClosestSpaceTimePointOnVector2;
    ClosestSpaceTimePointOnVector1.SetVect( ClosestPointOnVector1 );
    ClosestSpaceTimePointOnVector1(3) = ClosestTime1;
    ClosestSpaceTimePointOnVector2.SetVect( ClosestPointOnVector2 );
    ClosestSpaceTimePointOnVector2(3) = ClosestTime2;
    return ClosestSpaceTimePointOnVector1;
}

bool intersection( double bxmin, double bxmax, double bymin, double bymax, double bzmin, double bzmax, TVector3 VectorOrigin, TVector3 VectorDirection )  //30th jan

{
    double tmin = -INFINITY, tmax = INFINITY;

    if ( VectorDirection[0] != 0.0 )
    {
        double tx1 = ( bxmin - VectorOrigin[0] ) / VectorDirection[0];
        double tx2 = ( bxmax - VectorOrigin[0] ) / VectorDirection[0];

        tmin = max( tmin, min( tx1, tx2 ) );
        tmax = min( tmax, max( tx1, tx2 ) );
    }

    if ( VectorDirection[1] != 0.0 )
    {
        double ty1 = ( bymin - VectorOrigin[1] ) / VectorDirection[1];
        double ty2 = ( bymax - VectorOrigin[1] ) / VectorDirection[1];

        tmin = max( tmin, min( ty1, ty2 ) );
        tmax = min( tmax, max( ty1, ty2 ) );
    }

    if ( VectorDirection[2] != 0.0 )
    {
        double tz1 = ( bzmin - VectorOrigin[2] ) / VectorDirection[2];
        double tz2 = ( bzmax - VectorOrigin[2] ) / VectorDirection[2];

        tmin = max( tmin, min( tz1, tz2 ) );
        tmax = min( tmax, max( tz1, tz2 ) );
    }
    return tmax >= tmin;
}

double intersection_length( double bxmin, double bxmax, double bymin, double bymax, double bzmin, double bzmax, TVector3 VectorOrigin, TVector3 VectorDirection ) //30th jan

 {
    double tmin = -INFINITY, tmax = INFINITY;
    double tx1,tx2,ty1,ty2,tz1,tz2;
    if ( VectorDirection[0] != 0.0 )
    {
        tx1 = ( bxmin - VectorOrigin[0] ) / VectorDirection[0];
        tx2 = ( bxmax - VectorOrigin[0] ) / VectorDirection[0];

        tmin = max( tmin, min( tx1, tx2 ) );
        tmax = min( tmax, max( tx1, tx2 ) );
    }

    if ( VectorDirection[1] != 0.0 )
    {
        ty1 = ( bymin - VectorOrigin[1] ) / VectorDirection[1];
        ty2 = ( bymax - VectorOrigin[1] ) / VectorDirection[1];

        tmin = max( tmin, min( ty1, ty2 ) );
        tmax = min( tmax, max( ty1, ty2 ) );
    }

    if ( VectorDirection[2] != 0.0 )
    {
        tz1 = ( bzmin - VectorOrigin[2] ) / VectorDirection[2];
        tz2 = ( bzmax - VectorOrigin[2] ) / VectorDirection[2];

        tmin = max( tmin, min( tz1, tz2 ) );
        tmax = min( tmax, max( tz1, tz2 ) );
    }
    double length = sqrt( pow( ( tx2 - tx1 ), 2 ) + pow( ( ty2 - ty1 ), 2 ) + pow( ( tz2 - tz1 ), 2 ) );
    return length;
}




void spectro2::Process( int iEvent )
{
    //if( fMCSimple.fStatus == MCSimple::kMissing ){printIncompleteMCWarning( iEvent );return;}
    //if( fMCSimple.fStatus == MCSimple::kEmpty ){printNoMCWarning();return;}

    // Get the events from each detector
    TRecoSpectrometerEvent *SpectrometerEvent = ( TRecoSpectrometerEvent* )GetEvent( "Spectrometer" );
    TRecoLKrEvent *LKrEvent = ( TRecoLKrEvent* )GetEvent( "LKr" );
    TRecoMUV3Event *MUV3Event = ( TRecoMUV3Event* )GetEvent( "MUV3" );
    TRecoCedarEvent *CEDAREvent = ( TRecoCedarEvent* )GetEvent( "CEDAR" );
    TRecoLAVEvent *LAVEvent = ( TRecoLAVEvent* )GetEvent( "LAV" );
    TRecoRICHEvent *RICHEvent = ( TRecoRICHEvent*)GetEvent( "RICH" );
    TRecoSACEvent *SACEvent = ( TRecoSACEvent*)GetEvent( "SAC" );
    TRecoCHANTIEvent *CHANTIEvent = ( TRecoCHANTIEvent*)GetEvent( "CHANTI" );

    //Get truth event
    Event *MCTruthEvent = GetMCEvent();
    bool cluster_energetic_enough = false;
    bool ring_correct_size = false;
    int spectro_charge = 0;
    int decay_area = 0; //Initialise variable to determine where the particle decayed from. 0 is not in the beam, 1 is before the fiducial, 2 is in the fiducial.
    int LKrCandidates = LKrEvent-> GetNCandidates();
    double cluster_energy[ 100 ] = { 0 };
    double ring_radius = { 0 };
    double LKr_time[ 100 ] = { 0 };
    double LKr_time_difference[ 100 ] = { 0 };

    double MUV3_time[ 20 ] = { 0 };
    double MUV3_time_difference[ 20 ] = { 0 };
    double Spectro_time = 0;
    TVector3 true_muon_velocity_squared;
    bool been_detected = false;
    bool Intersect_LKr = false;
    bool detection_in_LKr = false;
    bool detection_in_rich = false;
    bool detection_in_CEDAR = false;
    bool detection_in_CHANTI = false;
    bool CEDAR_enough_sectors = false;
    bool detection_in_MUV3 = false;
    int sectors [ 20 ] = { 0 };
    TVector3 MUV3_position[ 20 ];
    TVector3 closest_point_MUV3[ 20 ];
    TVector3 distance_to_MUV3[ 20 ];
    TVector3 position_start;
    TVector3 position_after;
    TVector3 momentum;
    TVector3 momentum_after;
    TLorentzVector muon_momentum;
    TVector3 muon_velocity_squared;

    double min_dist_to_baxis_before_fiducial;
    TVector3 closest_point_from_baxis_before_fiducial;
    TVector3 closest_point_of_beam_approach_before_fiducial;
    TVector3 dist_to_baxis_before_fiducial;

    double min_dist_to_baxis_after_fiducial;
    TVector3 closest_point_from_baxis_after_fiducial;
    TVector3 closest_point_of_beam_approach_after_fiducial;
    TVector3 dist_to_baxis_after_fiducial;
    double reco_xz_angle;
    double reco_yz_angle;
    TLorentzVector missing_mass;
    //Check to see if an event was detected in the spectrometer


    if( SpectrometerEvent->GetNCandidates() >= 1 )
    {
        //Loop through each detected particle in the spectrometer
        for ( int k = 0; k < SpectrometerEvent->GetNCandidates(); k++ )
        {
            //Get the candidate from the spectrometer event
            TRecoSpectrometerCandidate *SpectroCandidate = ( TRecoSpectrometerCandidate* )SpectrometerEvent->GetCandidate( k );
            SpectroCandidate->SetEvent( SpectrometerEvent );

            spectro_charge = 1;
            //Get the three momentum of the candidate from the spectrometer before the magnet
            momentum = SpectroCandidate->GetThreeMomentumBeforeMagnet();
            momentum_after = SpectroCandidate->GetThreeMomentumAfterMagnet();

            //Get the three position of the candidate from the spectrometer before the magnet
            position_start = SpectroCandidate->GetPositionBeforeMagnet();
            position_after = SpectroCandidate->GetPositionAfterMagnet();

            Spectro_time = SpectroCandidate->GetTime();
            //Lots of calculations
            min_dist_to_baxis_before_fiducial = ( (b.fiducial_entry - position_start ).Dot( b.beam_axis.Cross( momentum ) ) ) / b.beam_axis.Cross( momentum ).Mag();
            closest_point_from_baxis_before_fiducial = ClosestPointOnVectorToOtherVector( position_start, momentum, b.fiducial_entry, b.beam_axis );
            closest_point_of_beam_approach_before_fiducial = ClosestPointOnVectorToOtherVector( b.fiducial_entry, b.beam_axis, position_start, momentum );
            dist_to_baxis_before_fiducial = - closest_point_of_beam_approach_before_fiducial + closest_point_from_baxis_before_fiducial;

            min_dist_to_baxis_after_fiducial = ( (b.fiducial_entry - position_start ).Dot( b.beam_axis_rotated.Cross( momentum ) ) ) / ( b.beam_axis_rotated.Cross( momentum ) ).Mag();
            closest_point_from_baxis_after_fiducial = ClosestPointOnVectorToOtherVector( position_start, momentum, b.fiducial_entry, b.beam_axis_rotated );
            closest_point_of_beam_approach_after_fiducial = ClosestPointOnVectorToOtherVector( b.fiducial_entry, b.beam_axis_rotated, position_start, momentum );
            dist_to_baxis_after_fiducial = - closest_point_of_beam_approach_after_fiducial + closest_point_from_baxis_after_fiducial;
            /////////////////////
            // Geometric Stuff //
            /////////////////////
            //Finds whether the muon comes from before or after the first magnet
            if( closest_point_from_baxis_before_fiducial[2] <= 104000
                &&( closest_point_from_baxis_after_fiducial[2] < 104000 || abs ( dist_to_baxis_before_fiducial.Mag() ) < abs( dist_to_baxis_after_fiducial.Mag() ) ) )
            {
                decay_area = 1;
            }
            else if ( closest_point_from_baxis_after_fiducial[2] >= 140000 /*104000*/ && closest_point_from_baxis_after_fiducial[2] <= 166000 )
            {
                decay_area = 2;
            }
            reco_xz_angle = SpectroCandidate->GetSlopeXBeforeMagnet();
            reco_yz_angle = SpectroCandidate->GetSlopeYBeforeMagnet();
            ////////////////////////
            // Missing mass stuff //
            ////////////////////////
            //Calculate energy of kaon
            double kaon_energy = sqrt( ( 5625000000 ) + ( kaon_mass_2 ) );
            TLorentzVector kaon_momentum;
            kaon_momentum[2] = 75e3; //Set z component to beam energy
            kaon_momentum[3] = kaon_energy; //Set time component to kaon energy
            kaon_momentum.RotateY( -BeamAngleFromZAxis ); //Transform to beam frame
            //Calculate energy of muon
            double muon_energy = sqrt( ( momentum.Mag2() ) + ( muon_mass_2 ) );
            muon_momentum.SetVect( momentum ); //Set xyz components to muon 3-momentum
            muon_momentum[3] = muon_energy; //Set time component to muon energy
            muon_velocity_squared = muon_momentum.Vect() * ( speed_of_light_2 / muon_energy );

            //Calculate missing mass squared
            missing_mass = kaon_momentum - muon_momentum;
            if  (  decay_area == 2 && been_detected == false )
            {
                FillHisto( "TotalDecays", 1);
                been_detected == true;
            }
        }
    }

    if ( LKrEvent->GetNCandidates() >= 1 )
    {
        detection_in_LKr = true;
        for ( int iLKr = 0; iLKr < LKrEvent->GetNCandidates(); iLKr++ )
        {
            TRecoLKrCandidate *LKrCandidate = ( TRecoLKrCandidate* )LKrEvent->GetCandidate( iLKr );
            LKrCandidate->SetEvent( LKrEvent );

            cluster_energy[iLKr] = LKrCandidate->GetClusterEnergy();
            LKr_time[iLKr] = LKrCandidate->GetTime();
            LKr_time_difference[iLKr] = LKr_time[iLKr] - Spectro_time;
            if ( cluster_energy[iLKr] <= muon_momentum[3] * 0.001 && cluster_energy[iLKr] >= muon_momentum[3] * 0.000007142857142857142857 )
            {
                cluster_energetic_enough = true;
            }
        }
        TRecoLKrCandidate *LKrCandidate = ( TRecoLKrCandidate* )LKrEvent->GetCandidate( 0 );
        LKrCandidate->SetEvent( LKrEvent );


        //Project ray onto lkr front face
        TVector3 ray = position_after + (240000-position_after[2])*(momentum_after*(1/momentum_after[2]));


        //calculate distance between projection point and lkr centre
        TVector3 LKr_centre(-7.1,8.975,240000);
        TVector3 rad = ray - LKr_centre;
        double dist = rad.Mag();
        //make sure is within lkr face
        if( dist >= 150 && dist <= 1200 )
        {
            double lkr_x = LKrCandidate->GetClusterX()*10;
            double lkr_y = LKrCandidate->GetClusterY()*10;
            FillHisto("TrackPosition",ray[0],ray[1]);
            FillHisto("ClusterPosition",lkr_x,lkr_y);

            double clust_dist2 = pow( ( lkr_x-ray[0] ), 2 ) + pow( ( lkr_y-ray[1] ), 2 );
            double clust_dist = sqrt( clust_dist2 );
            FillHisto("ClusterTrackDistanceX",lkr_x-ray[0]);
            FillHisto("ClusterTrackDistanceY",lkr_y-ray[1]);
            //std::cout << clust_dist << std::endl;
            if(clust_dist2 <= 2500)
            {

                Intersect_LKr = true;
            }
            else Intersect_LKr = false;
        }
        else
        {
            Intersect_LKr = false;
        }

        //calculate distance between lkr cluster and
        /*
        if  (   intersection( LKr1minx, Lkr1maxx, Lkr1miny, LKr1maxy, LKr1minz, Lkr1maxz, position_after, momentum_after ) == true
                || intersection( LKr2minx, Lkr2maxx, Lkr2miny, LKr2maxy, LKr2minz, Lkr2maxz, position_after, momentum_after ) == true )
        {
            Intersect_LKr = true; //If ray intersects cuboid of LKr, set intersects.
            if  (   intersection( LKrhole1minx, Lkrhole1maxx, Lkrhole1miny, LKrhole1maxy, LKrhole1minz, Lkrhole1maxz, position_after, momentum_after ) == true
                    && intersection( LKrhole2minx, Lkrhole2maxx, Lkrhole2miny, LKrhole2maxy, LKrhole2minz, Lkrhole2maxz, position_after, momentum_after ) == true )
            {
                Intersect_LKr = false; //If ray intersects the entrance to the hole in the LKr AND the exit of the hole, it did not hit the LKr.
            }
        }
        */
    }

    if ( RICHEvent->GetNRingCandidates() == 1 )
    {
        detection_in_rich = true;
        TRecoRICHCandidate *RICHCandidate = (TRecoRICHCandidate*)RICHEvent->GetCandidate(0);
        RICHCandidate->SetEvent(RICHEvent);
        ring_radius = RICHCandidate->GetRingRadius();
        double velocity = sqrt( ( pow( ( ring_radius * RICH_focal_length_inverse ), 2 ) + 1 ) ) * refractive_index_inverse;
        double RICH_velocity_over_muon_assumption_velocity = ( sqrt( muon_velocity_squared.Mag() ) * speed_of_light_inverse ) / velocity;
        if ( RICH_velocity_over_muon_assumption_velocity <= 1.000002 && RICH_velocity_over_muon_assumption_velocity >= 0.999998 )
        {
            ring_correct_size = true;
        }
    }

    if ( MUV3Event->GetNCandidates() >= 1 )
    {
        detection_in_MUV3 = true;
        for ( int iMUV3 = 0; iMUV3 < MUV3Event->GetNCandidates(); iMUV3++ )
        {
            TRecoMUV3Candidate *MUV3Candidate = ( TRecoMUV3Candidate* )MUV3Event->GetCandidate( iMUV3 );
            MUV3Candidate->SetEvent( MUV3Event );
            MUV3_position[iMUV3] = MUV3Candidate->GetPosition();
            closest_point_MUV3[iMUV3] = ClosestPointOnVectorToPoint( position_after, momentum_after, MUV3_position[iMUV3] );
            distance_to_MUV3[iMUV3] = -closest_point_MUV3[iMUV3] + MUV3_position[iMUV3];
            MUV3_time[iMUV3] = MUV3Candidate-> GetTime() ;
            MUV3_time_difference[iMUV3] = MUV3_time[iMUV3] - Spectro_time;
        }
    }

    if ( CEDAREvent->GetNCandidates() >= 1 )
    {
        detection_in_CEDAR = true;
        for ( int iCEDAR = 0; iCEDAR<CEDAREvent->GetNCandidates(); iCEDAR++ )
        {
            TRecoCedarCandidate *CEDARCandidate = (TRecoCedarCandidate*)CEDAREvent->GetCandidate(iCEDAR);
            CEDARCandidate->SetEvent(CEDAREvent);
            sectors[iCEDAR] = CEDARCandidate->GetNSectors();
            if ( sectors[iCEDAR] >= 5 )
            {
                CEDAR_enough_sectors = true;
            }
        }
    }

    if ( CHANTIEvent->GetNCandidates() >= 1 )
    {
        detection_in_CHANTI = true;
    }
    ///////////////////////////////////
    // Plot Histos with Restrictions //
    ///////////////////////////////////
    //Attempts to select only the muon
    int accepted_number_of_LKr_candidates = 1, accepted_number_of_spectrometer_candidates = 1, accepted_number_of_MUV3_candidates = 1;

    if  (   spectro_charge == 1  //Positive Charge
            && SpectrometerEvent->GetNCandidates() == accepted_number_of_spectrometer_candidates //Single Detection In Spectrometer
            && decay_area == 2  //Decay in Fiducial Region //34728 particles when running on 100,000 munu.
            && (LKrEvent->GetNCandidates() == 1 && Intersect_LKr) || (LKrEvent->GetNCandidates() == 0 && !Intersect_LKr)
            //&& ( LKrEvent->GetNCandidates() <= accepted_number_of_LKr_candidates /*|| Intersect_LKr == false*/ )  // Be detected in LKr, if the muon should hit LKr //27429 with just getncandidates, 31099 with or intersect false.
            && ( cluster_energetic_enough == true || LKrEvent->GetNCandidates() == 0 /*|| Intersect_LKr == false*/ ) //IF a cluster is in the LKr, it should have energy greater than 1/1000 of muon energy and less than 1/140000 //27273, 30943 with or intersect false
            && ( ring_correct_size == true || detection_in_rich == false ) //IF a ring is in the RICH, it should be consistent with a muon. //27273, 30943 with or intersect false for LKr
            && LAVEvent->GetNCandidates() == 0 //No muon should be detected in the LAV. //26112, 27557 with or intersect false for LKr
            && MUV3Event->GetNCandidates() >= accepted_number_of_MUV3_candidates // Be detected in MUV3 IF DETECTED PARTICLE WOULD HIT IT //26076 , 26086 with or intersect false for LKr
            && dist_to_baxis_after_fiducial.Mag() <= 40 && dist_to_baxis_after_fiducial[0] <= 40 && dist_to_baxis_after_fiducial[1] <= 40 && dist_to_baxis_after_fiducial[2] <= 40 //Be linked closely to the beam. // 26069 with or intersect false for LKr
            && abs ( MUV3_time_difference[0] ) <= 15
            && abs ( LKr_time_difference[0] ) <= 15
            && SACEvent->GetNCandidates() == 0
            && abs(distance_to_MUV3[0][0]) <= 200
            && abs(distance_to_MUV3[0][1]) <= 500
            && abs(distance_to_MUV3[0][2]) <= 10
            && abs(distance_to_MUV3[0].Mag()) <= 500
            //&& ( momentum.Mag() <= 35000 && momentum.Mag() >= 15000 && detection_in_rich == true ) //RICH can only distinguish in this range. 8211 with or intersect false for LKr
            //&& missing_mass.Mag2() / pow( 1000, 2 ) < 2000 //Have Missing Mass Correct for Decay
            //&& CEDAREvent->GetNCandidates() >=1 // Be detected in CEDAR
            //&& CEDAREvent->GetNHits() >=10 //Have atleast 10 hits in the CEDAR
            //&& CEDAR_enough_sectors == true
            && detection_in_CHANTI == false
            //&& detection_in_rich == true
        )
    {
        /*
        //true stuff here

        if ( MCTruthEvent -> GetNKineParts() >= 1 )
        {

            //Get the first decay product associated with the decay event (hopefully a muon)
            KinePart *TrueCandidate = ( KinePart* )MCTruthEvent -> GetKineParts() -> At( 1 );
            //Get the kaon associated with the decay event
            KinePart *KaonCandidate = ( KinePart* )MCTruthEvent -> GetKineParts() -> At( 0 );

            TLorentzVector true_kaon_momentum = KaonCandidate->GetFinal4Momentum();

            TVector3 true_momentum = TrueCandidate->GetMomAtCheckPoint( 2 ).Vect();
            TVector3 true_position_start = TrueCandidate->GetProdPos().Vect();
            TVector3 true_position_end = TrueCandidate->GetEndPos().Vect();
            TVector3 true_LKr_entry = TrueCandidate ->GetPosLKrEntry();
            double true_xz_angle = atan( true_momentum[0] / true_momentum[2] );
            double true_yz_angle = atan( true_momentum[1] / true_momentum[2] );

            //Make sure the first true event decay product is a muon
            if ( true_momentum.Mag() != 0 && TrueCandidate -> GetPDGcode() == -13 && abs( true_momentum.Theta() ) > 0 )
            {
                true_momentum.RotateY( BeamAngleFromZAxis );
                momentum.RotateY( BeamAngleFromZAxis );
                for( int i = 0; i < 13; i++ )
                {
                    if( true_momentum.Mag() * 0.001 < 5 * ( i + 3 ) && true_momentum.Mag() * 0.001 >= 5 * ( i + 2 ) )
                    {
                        string intstr = to_string( i );
                        TString num = TString( intstr );
                        FillHisto( TString( "ResolutionTemp" ) + num, ( true_momentum.Mag() - momentum.Mag() ) * 0.001 );
                    }
                    if( true_momentum[0] * 0.001 < 0.026666666666666666 * ( i + 1 ) - 0.2 && true_momentum[0] * 0.001 >= 0.026666666666666666 * ( i ) - 0.2 )
                    {
                        string intstr = to_string( i );
                        TString num = TString( intstr );
                        FillHisto( TString( "ResolutionTempX" ) + num, ( true_momentum[0] - momentum[0] ) * 0.001 );
                    }
                    if( true_momentum[1] * 0.001 < 0.026666666666666666 * ( i + 1 ) - 0.2 && true_momentum[1] * 0.001 >= 0.026666666666666666 * ( i ) - 0.2 )
                    {
                        string intstr = to_string( i );
                        TString num = TString( intstr );
                        FillHisto( TString( "ResolutionTempY" ) + num, ( true_momentum[1] - momentum[1] ) * 0.001 );
                    }
                    if( true_momentum[2] * 0.001 < 5 * ( i + 2 ) && true_momentum[2] * 0.001 >= 5 * ( i + 1 ) )
                    {
                        string intstr = to_string( i );
                        TString num = TString( intstr );
                        FillHisto( TString( "ResolutionTempZ" ) + num, ( true_momentum[2] - momentum[2] ) * 0.001 );
                    }
                }
                true_momentum.RotateY(-BeamAngleFromZAxis);
                momentum.RotateY(-BeamAngleFromZAxis);
                for( int i = 0; i < 13; i++ )
                {
                    if( true_xz_angle < 0.0015333 * ( i + 1 ) - 0.008 && true_xz_angle >= 0.0015333 * ( i ) - 0.008 )
                    {
                        string intstr = to_string( i );
                        TString num = TString( intstr );
                        FillHisto( TString( "ResolutionTempXZ" ) + num, ( true_xz_angle - reco_xz_angle ) );
                    }
                    if( true_yz_angle < ( 0.0016 ) * ( i + 1 ) - 0.012 && true_yz_angle >= ( 0.0016 ) * ( i ) - 0.012 )
                    {
                        string intstr = to_string( i );
                        TString num = TString( intstr );
                        FillHisto( TString( "ResolutionTempYZ" ) + num, ( true_yz_angle - reco_yz_angle ) );
                    }
                }
                ////////////////////////////////////
                // Calculate Momentum Resolutions //
                ////////////////////////////////////
                true_momentum.RotateY( BeamAngleFromZAxis );
                momentum.RotateY( BeamAngleFromZAxis );
                TVector3 ResolutionTemp = true_momentum - momentum;
                FillHisto( "SpectrometerXMomentumResolution",  ResolutionTemp[0] );
                FillHisto( "SpectrometerYMomentumResolution", ResolutionTemp[1] );
                FillHisto( "SpectrometerZMomentumResolution", ResolutionTemp[2] * 0.001 );
                FillHisto( "SpectrometerMomentumResolution", ResolutionTemp.Mag() * 0.001 );

                FillHisto( "SpectrometerXMomentumResolutionAgainstXMomentum", true_momentum[0], ResolutionTemp[0] );
                FillHisto( "SpectrometerYMomentumResolutionAgainstYMomentum", true_momentum[1], ResolutionTemp[1] );
                FillHisto( "SpectrometerZMomentumResolutionAgainstZMomentum", true_momentum[2] * 0.001, ResolutionTemp[2] * 0.001 );
                FillHisto( "SpectrometerMomentumResolutionAgainstMomentum",  true_momentum.Mag() * 0.001, ResolutionTemp.Mag() * 0.001 );


                FillHisto( "KaonEndingPosition", true_position_end[2] * 0.001, true_position_end[0] );
                FillHisto( "TrueMomentumHist",   true_momentum.Mag() * 0.001 );
                FillHisto( "TruexMomentumHist",  true_momentum[0] * 0.001 );
                FillHisto( "TrueyMomentumHist",  true_momentum[1] * 0.001 );
                FillHisto( "TruezMomentumHist",  true_momentum[2] * 0.001 );


                FillHisto( "TrueTransverseMomentumHist",     true_momentum.Perp() );
                FillHisto( "TrueAzimuthalMomentumHist",      true_momentum.Phi() );
                FillHisto( "TruePolarMomentumHist",          true_momentum.Theta() );
                FillHisto( "TrueEnergyVsAzimuthal",          true_momentum.Phi(),   true_momentum.Mag()  * 0.001 );
                FillHisto( "TrueEnergyVsPolar",              true_momentum.Theta(), true_momentum.Mag()  * 0.001 );
                FillHisto( "TrueTranverseEnergyVsAzimuthal", true_momentum.Phi(),   true_momentum.Perp() * 0.001 );

                FillHisto( "CompareTrueMomentumHist",  true_momentum.Mag() * 0.001 );
                FillHisto( "CompareTruexMomentumHist", true_momentum[0] * 0.001 );
                FillHisto( "CompareTrueyMomentumHist", true_momentum[1] * 0.001 );
                FillHisto( "CompareTruezMomentumHist", true_momentum[2] * 0.001 );

                FillHisto( "CompareTrueTransverseMomentumHist", true_momentum.Perp() );
                FillHisto( "CompareTrueAzimuthalMomentumHist",  true_momentum.Phi() );
                FillHisto( "CompareTruePolarMomentumHist",      true_momentum.Theta() );

                FillHisto( "CompareTrueEnergyVsAzimuthal",          true_momentum.Phi(),   true_momentum.Mag() * 0.001 );
                FillHisto( "CompareTrueEnergyVsPolar",              true_momentum.Theta(), true_momentum.Mag() * 0.001 );
                FillHisto( "CompareTrueTranverseEnergyVsAzimuthal", true_momentum.Phi(),   true_momentum.Perp() * 0.001 );

                FillHisto( "TrueProductionPosition", true_position_start[2] * 0.001, true_position_start[0] );

                FillHisto( "TrueProductionPositionX", true_position_start[0] );
                FillHisto( "TrueProductionPositionY", true_position_start[1] );
                FillHisto( "TrueProductionPositionZ", true_position_start[2] * 0.001 );

                true_momentum.RotateY( -BeamAngleFromZAxis );
                momentum.RotateY( -BeamAngleFromZAxis );

                FillHisto( "TrueXZAngle",   true_xz_angle );
                FillHisto( "TrueYZAngle",   true_yz_angle );



                TLorentzVector true_muon_momentum = TrueCandidate->GetInitial4Momentum();
                //TLorentzVector TrueMuonMomentum = TrueCandidate->GetMomAtCheckPoint(2);
                TLorentzVector true_missing_mass = true_kaon_momentum - true_muon_momentum;
                true_muon_velocity_squared = true_muon_momentum.Vect() *  ( speed_of_light_2 / true_muon_momentum[3] );

                FillHisto("TrueMissingMass", true_missing_mass.Mag2() * 0.000001 );
                FillHisto( "Energyclusterenergy", true_muon_momentum[3] / cluster_energy[0] );
                FillHisto( "LkrEntryx", true_LKr_entry[0] );
                FillHisto( "LkrEntryy", true_LKr_entry[1] );
                FillHisto( "LkrEntryz", true_LKr_entry[2] );
            }
        }
        */
        //Reco stuff
        for ( int iLKr = 0; iLKr < accepted_number_of_LKr_candidates; iLKr++ )
        {
            FillHisto( "LKrtime", LKr_time[iLKr] );
            FillHisto( "LKrTimeDifference", LKr_time_difference[iLKr] );
        }
        FillHisto( "LKrcandidates", LKrEvent->GetNCandidates() );


        for ( int iMUV3 = 0; iMUV3 < accepted_number_of_MUV3_candidates; iMUV3++)
        {
            FillHisto( "closestxDistanceToMUV3", distance_to_MUV3[iMUV3][0] );
            FillHisto( "closestyDistanceToMUV3", distance_to_MUV3[iMUV3][1] );
            FillHisto( "closestzDistanceToMUV3", distance_to_MUV3[iMUV3][2] );
            FillHisto( "closestDistanceToMUV3", distance_to_MUV3[iMUV3].Mag() );
            FillHisto( "closestxPointToMUV3", closest_point_MUV3[iMUV3][0] );
            FillHisto( "closestyPointToMUV3", closest_point_MUV3[iMUV3][1] );
            FillHisto( "closestzPointToMUV3", closest_point_MUV3[iMUV3][2] );
            FillHisto( "MUV3Position", closest_point_MUV3[iMUV3][2] * 0.001 , closest_point_MUV3[iMUV3][0] );
            FillHisto( "MUV3time", MUV3_time[iMUV3] );
            FillHisto( "timedifference", MUV3_time_difference[iMUV3] );
        }
        FillHisto( "MUV3candidates", MUV3Event->GetNCandidates() );
        for ( int iCEDAR = 0; CEDAREvent->GetNCandidates(); iCEDAR++ )
        {
            FillHisto( "CEDARSectors", sectors[iCEDAR]);
        }
        FillHisto( "CEDARCandidates", CEDAREvent->GetNCandidates() );


        FillHisto( "XZAngle", reco_xz_angle );
        FillHisto( "YZAngle", reco_yz_angle );
        TVector3 Direction;
        Direction[0] = 1 * tan( reco_xz_angle );
        Direction[1] = 1 * tan( reco_yz_angle );
        Direction[2] = 1;
        Direction = Direction.Unit(); //Normalize
        TVector3 ReconstructedMomentum = Direction * momentum.Mag();

        ReconstructedMomentum.RotateY( BeamAngleFromZAxis );
        FillHisto( "ReconstructedMomentumHist", ReconstructedMomentum.Mag() * 0.001 );
        FillHisto( "ReconstructedxMomentumHist", ReconstructedMomentum[0] * 0.001 );
        FillHisto( "ReconstructedyMomentumHist", ReconstructedMomentum[1] * 0.001 );
        FillHisto( "ReconstructedzMomentumHist", ReconstructedMomentum[2] * 0.001 );
        ReconstructedMomentum.RotateY( -BeamAngleFromZAxis );

        FillHisto( "MissingMass", missing_mass.Mag2() / pow( 1000, 2) );

        momentum.RotateY(BeamAngleFromZAxis);   //Rotate the reference frame to be along the beam

        FillHisto( "MomentumHist",  momentum.Mag() * 0.001 );
        FillHisto( "xMomentumHist", momentum[0]    * 0.001 );
        FillHisto( "yMomentumHist", momentum[1]    * 0.001 );
        FillHisto( "zMomentumHist", momentum[2]    * 0.001 );



        FillHisto( "TransverseMomentumHist", momentum.Perp()  );
        FillHisto( "AzimuthalMomentumHist",  momentum.Phi()   );
        FillHisto( "PolarMomentumHist",      momentum.Theta() );

        FillHisto( "EnergyVsAzimuthal",          momentum.Phi(),   momentum.Mag()  * 0.001 );
        FillHisto( "EnergyVsPolar",              momentum.Theta(), momentum.Mag()  * 0.001 );
        FillHisto( "TranverseEnergyVsAzimuthal", momentum.Phi(),   momentum.Perp() * 0.001 );


        FillHisto( "spectrotime", Spectro_time );

        double Lkr_intersection_length;
        momentum.RotateY( -BeamAngleFromZAxis );  //Rotate reference frame back to along the detector
        /*
        if ( intersection( LKr1minx, Lkr1maxx, Lkr1miny, LKr1maxy, LKr1minz, Lkr1maxz, position_after, momentum_after ) == false )
        {
            Lkr_intersection_length = intersection_length( LKr1minx, Lkr1maxx, Lkr1miny, LKr1maxy, LKr1minz, Lkr1maxz, position_after, momentum_after );
        }
        else if ( intersection( LKr2minx, Lkr2maxx, Lkr2miny, LKr2maxy, LKr2minz, Lkr2maxz, position_after, momentum_after ) == false )
        {
            Lkr_intersection_length = intersection_length( LKr2minx, Lkr2maxx, Lkr2miny, LKr2maxy, LKr2minz, Lkr2maxz, position_after, momentum_after );
        }

        if  (   intersection( LKrhole1minx, Lkrhole1maxx, Lkrhole1miny, LKrhole1maxy, LKrhole1minz, Lkrhole1maxz, position_after, momentum_after ) == false
                && intersection( LKrhole2minx, Lkrhole2maxx, Lkrhole2miny, LKrhole2maxy, LKrhole2minz, Lkrhole2maxz, position_after, momentum_after ) == false
            )
        {
            FillHisto( "ClusterEnergyAndIntersectionLength", Lkr_intersection_length / cluster_energy[0] );
        }
        */
        if ( decay_area == 1 )    ///This is never true because the first if contains it to decay_area == 2.
        {
            FillHisto( "ClosestPointFromBeamAxis",   closest_point_from_baxis_before_fiducial.Mag() * 0.001 );
            FillHisto( "ClosestxPointFromBeamAxis",  closest_point_from_baxis_before_fiducial[0] );
            FillHisto( "ClosestyPointFromBeamAxis",  closest_point_from_baxis_before_fiducial[1] );
            FillHisto( "ClosestzPointFromBeamAxis",  closest_point_from_baxis_before_fiducial[2] * 0.001 );
            FillHisto( "MinimumDistanceToBeamAxis",  min_dist_to_baxis_before_fiducial );
            FillHisto( "ClosestDistanceToBeamAxis",  dist_to_baxis_before_fiducial.Mag() );
            FillHisto( "ClosestxDistanceToBeamAxis", dist_to_baxis_before_fiducial[0] );
            FillHisto( "ClosestyDistanceToBeamAxis", dist_to_baxis_before_fiducial[1] );
            FillHisto( "ClosestzDistanceToBeamAxis", dist_to_baxis_before_fiducial[2] );
            FillHisto( "DecayPoisition",             closest_point_from_baxis_before_fiducial[2] * 0.001 , closest_point_from_baxis_before_fiducial[0] );
            FillHisto("MissingMassVsZPosition",      closest_point_from_baxis_before_fiducial[2] * 0.001 , missing_mass.Mag2() * 0.000001 );

        }
        if ( decay_area == 2 )
        {
            FillHisto( "ClosestPointFromBeamAxis",   closest_point_from_baxis_after_fiducial.Mag() * 0.001 );
            FillHisto( "ClosestxPointFromBeamAxis",  closest_point_from_baxis_after_fiducial[0] );
            FillHisto( "ClosestyPointFromBeamAxis",  closest_point_from_baxis_after_fiducial[1] );
            FillHisto( "ClosestzPointFromBeamAxis",  closest_point_from_baxis_after_fiducial[2] * 0.001 );
            FillHisto( "MinimumDistanceToBeamAxis",  min_dist_to_baxis_after_fiducial );
            FillHisto( "ClosestDistanceToBeamAxis",  dist_to_baxis_after_fiducial.Mag() );
            FillHisto( "ClosestxDistanceToBeamAxis", dist_to_baxis_after_fiducial[0] );
            FillHisto( "ClosestyDistanceToBeamAxis", dist_to_baxis_after_fiducial[1] );
            FillHisto( "ClosestzDistanceToBeamAxis", dist_to_baxis_after_fiducial[2] );
            FillHisto( "DecayPoisition",             closest_point_from_baxis_after_fiducial[2] * 0.001 , closest_point_from_baxis_after_fiducial[0] );
            FillHisto( "MissingMassVsZPosition",     closest_point_from_baxis_after_fiducial[2] * 0.001 , missing_mass.Mag2() * 0.000001 );
        }
    }
}

void spectro2::PostProcess()
{

}

void spectro2::EndOfBurstUser()
{

}

void spectro2::EndOfRunUser()
{


    ///cout<<"start" << startime << "end:" << endtime << endl << "TotalTime:" << endtime - startime;

    TH1D* h58 = new TH1D("DetectorEfficiencyMomentum", "Detector Efficiency Momentum", NumberOfBins, 0, 0 );
    BookHisto(h58);
    fHisto.GetHisto("MomentumHist")->Copy(*h58);
    h58->Divide(fHisto.GetHisto("CompareTrueMomentumHist"));

    TH1D* h59 = new TH1D("DetectorxEfficiencyMomentum", "Detector Efficiency xMomentum", NumberOfBins, 0, 0 );
    BookHisto(h59);
    fHisto.GetHisto("xMomentumHist")->Copy(*h59);
    h59->Divide(fHisto.GetHisto("CompareTruexMomentumHist"));

    TH1D* h60 = new TH1D("DetectoryEfficiencyMomentum", "Detector Efficiency yMomentum", NumberOfBins, 0, 0 );
    BookHisto(h60);
    fHisto.GetHisto("yMomentumHist")->Copy(*h60);
    h60->Divide(fHisto.GetHisto("CompareTrueyMomentumHist"));

    TH1D* h61 = new TH1D("DetectorzEfficiencyMomentum", "Detector Efficiency zMomentum", NumberOfBins, 0, 0 );
    BookHisto(h61);
    fHisto.GetHisto("zMomentumHist")->Copy(*h61);
    h61->Divide(fHisto.GetHisto("CompareTruezMomentumHist"));



    TH1D* h62 = new TH1D("Reconstructed Momentum", "Reconstructed Momentum", NumberOfBins, 0, 0 );
    BookHisto(h62);
    fHisto.GetHisto("MomentumHist")->Copy(*h62);
    h62->Divide(fHisto.GetHisto("DetectorEfficiencyMomentum"));

    TH1D* h63 = new TH1D("Reconstructed x Momentum", "Reconstructed x Momentum", NumberOfBins, 0, 0 );
    BookHisto(h63);
    fHisto.GetHisto("xMomentumHist")->Copy(*h63);
    h63->Divide(fHisto.GetHisto("DetectorxEfficiencyMomentum"));

    TH1D* h64 = new TH1D("Reconstructed y Momentum", "Reconstructed y Momentum", NumberOfBins, 0, 0 );
    BookHisto(h64);
    fHisto.GetHisto("yMomentumHist")->Copy(*h64);
    h64->Divide(fHisto.GetHisto("DetectoryEfficiencyMomentum"));

    TH1D* h65 = new TH1D("Reconstructed z Momentum", "Reconstructed z Momentum", NumberOfBins, 0, 0 );
    BookHisto(h65);
    fHisto.GetHisto("zMomentumHist")->Copy(*h65);
    h65->Divide(fHisto.GetHisto("DetectorzEfficiencyMomentum"));

    TH1* h = fHisto.GetHisto("MissingMass");
    h->Fit("gaus");
    h->Draw();
    /*
    int n = 13;
    double x[n], y[n], ey[n], ry[n], ery[n];



    for( int j = 0; j < 6; j++ )
    {
        TString dim;
        TString num;
        for( int i = 0; i < 13; i++ )
        {
            string dimstr = "";
            switch ( j )
            {
                case 0:
                    dimstr = "";
                    x[i] = 5 * ( i + 2 );
                    break;
                case 1:
                    dimstr = "X";
                    x[i] = 0.026666666666666666 * ( i + 1 ) - 0.2;
                    break;
                case 2:
                    dimstr = "Y";
                    x[i] = 0.026666666666666666 * ( i + 1 ) - 0.2;
                    break;
                case 3:
                    dimstr = "Z";
                    x[i] = 5 * ( i + 1 );
                    break;
                case 4:
                    dimstr = "XZ";
                    x[i] = 0.0015333 * ( i + 1 ) - 0.008;
                    break;
                case 5:
                    dimstr = "YZ";
                    x[i] = ( 0.0016 ) * ( i + 1 ) - 0.012;
                    break;
                default:
                    dimstr = "";
                    x[i] = 5 * ( i + 1 );
                    break;
            }
            dim = TString( dimstr );
            string numstr = to_string( i );
            num = TString( numstr );
            TH1* res = fHisto.GetHisto( TString( "ResolutionTemp" ) + dim + num );
            double integral = res->Integral();
            double integral_temp = res->Integral( 1, 500 );

            double xmin = res->GetXaxis()->GetXmin();
            double xmax = res->GetXaxis()->GetXmax();
            double xbin = abs( xmin ) > abs( xmax ) ? abs( xmax ) : abs( xmin );

            int binx1 = res->GetXaxis()->FindBin( -xbin );
            int binx2 = res->GetXaxis()->FindBin( xbin );

            integral_temp = res->Integral( binx1, binx2 );
            int res_limits = 0;
            while( abs( integral_temp ) > abs( integral ) * 0.98 )
            {
                integral_temp = res->Integral( binx1 + res_limits, binx2 - res_limits );
                res_limits++;
            }
            double fitmin = res->GetXaxis()->GetBinLowEdge( binx1 + res_limits );
            double fitmax = res->GetXaxis()->GetBinUpEdge( binx2 - res_limits );

            cout << i << endl;
            res->Fit( "gaus", "Q", "", fitmin, fitmax );
            res->Draw();

            TF1* fit = res->GetFunction( "gaus" );
            y[i] = fit->GetParameter( 2 );
            ey[i] = fit->GetParError( 2 );
            ry[i] = y[i] / x[i] ;
            ery[i] = abs( ey[i] / x[i] ) ;
            cout << x[i] << " " << y[i] << endl;
        }



        TGraphErrors* graph = new TGraphErrors( n, x, y, 0, ey );
        BookHisto( dim + "MtmResolutionVs" + dim + "Mtm", graph );
        TGraphErrors* rgraph = new TGraphErrors( n, x, ry, 0, ery );
        BookHisto("Relative" + dim + "MtmResolutionVs" + dim + "Mtm", rgraph );
    }

    TF1* f = new TF1( "f", "sqrt([0]^2+([1]*x)^2)", 10, 75 );
    f->SetParameter( 0, 0.01 );
    f->SetParameter( 1, 0.001 );
    TGraph* g = fHisto.GetTGraph( "RelativeMtmResolutionVsMtm" );
    g->Fit( "f" );
    */
    SaveAllPlots();

}

void spectro2::DrawPlot()
{
    DrawAllPlots();
    //SaveAllPlotsPDF();
}
