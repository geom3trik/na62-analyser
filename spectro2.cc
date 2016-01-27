#include <stdlib.h>
#include <iostream>
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
{	//SOMETHING ELSE
    RequestTree( "GigaTracker", new TRecoGigaTrackerEvent );
    RequestTree( "Spectrometer", new TRecoSpectrometerEvent );
    RequestTree( "LKr", new TRecoLKrEvent );
    RequestTree( "MUV3", new TRecoMUV3Event );
    RequestTree( "CEDAR", new TRecoCedarEvent );

}

void spectro2::InitOutput()
{

}

void spectro2::CreateHist1D(TString name, TString title, int nbins, double low, double high)
{
    TH1D* h1 = new TH1D(name,title,nbins,low,high);
    BookHisto(h1);
}

void spectro2::SetHistAxisLabels(TString name, TString xlabel, TString ylabel)
{
    TH1* h = fHisto.GetHisto(name);
    if(h != NULL)
    {
        h -> GetXaxis() -> SetTitle(xlabel);
        h -> GetYaxis() -> SetTitle(ylabel);
    }

}

void spectro2::InitHist()
{
    // Reconstructed momentum histogram
    CreateHist1D("MomentumHist","Momentum",NumberOfBins,0,80);
    SetHistAxisLabels("MomentumHist","Momentum GeV","Number of Entries");

    // True momentum histogram plotted with same bins as reconstructed momentum
    CreateHist1D( "CompareTrueMomentumHist", " True Momentum", NumberOfBins, 0, 80 );
    SetHistAxisLabels("CompareTrueMomentumHist","Momentum GeV","Number of Entries");

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


    for(int i=1;i<=15;i++)
    {

        string intstr = to_string(i);
        TString num = TString(intstr);
        CreateHist1D(TString("ResolutionTemp") + num, "Title", 500,0,0);
        CreateHist1D(TString("ResolutionTempX") + num, "Title", 500,0,0);
        CreateHist1D(TString("ResolutionTempY") + num, "Title", 500,0,0);
        CreateHist1D(TString("ResolutionTempZ") + num, "Title", 500,0,0);
    }
    for(int i=1;i<=15;i++)
    {
        string intstr = to_string(i);
        TString num = TString(intstr);
        CreateHist1D(TString("xzResolutionTemp") + num, "Title", 500,0,0);
        CreateHist1D(TString("yzResolutionTemp") + num, "Title", 500,0,0);
    }



    TH2D* h56 = new TH2D( "MissingMassVsZPosition", "Missing Mass Vs Z position", NumberOfBins, 0, 0, NumberOfBins, 0, 0 );
    h56 -> GetXaxis() -> SetTitle( "GeV Squared" );
    h56 -> GetYaxis() -> SetTitle( "m" );
    BookHisto( h56 );

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
    for(ptr = canvases.begin(); ptr != canvases.end(); ptr++)
    {
        // Retrieve ROOT canvas and save as .pdf file
        ptr->second->GetCanvas()->SaveAs(TString(ptr->first + ".pdf"));
    }
}


TVector3 ClosestPointOnVectorToOtherVector( TVector3 Vector1Position, TVector3 Vector1Direction, TVector3 Vector2Position, TVector3 Vector2Direction) //14th November
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

TLorentzVector ClosestSpaceTimePointOnVectorToOtherVector( TVector3 Vector1Position, TVector3 Vector1Momentum, double Mass1, double Time1, TVector3 Vector2Position, TVector3 Vector2Momentum, double Mass2, double Time2) //14th November
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


void spectro2::Process( int iEvent )
{
	//if( fMCSimple.fStatus == MCSimple::kMissing ){printIncompleteMCWarning( iEvent );return;}
    //if( fMCSimple.fStatus == MCSimple::kEmpty ){printNoMCWarning();return;}

    // Get the events from each detector
	TRecoSpectrometerEvent *SpectrometerEvent = ( TRecoSpectrometerEvent* )GetEvent( "Spectrometer" );
    TRecoLKrEvent *LKrEvent = ( TRecoLKrEvent* )GetEvent( "LKr" );
    TRecoMUV3Event *MUV3Event = ( TRecoMUV3Event* )GetEvent( "MUV3" );
    TRecoCedarEvent *CEDAREvent = ( TRecoCedarEvent* )GetEvent( "CEDAR" );
    //Get truth event
    Event *MCTruthEvent = GetMCEvent();

	//Check to see if an event was detected in the spectrometer
	if( SpectrometerEvent->GetNCandidates() >= 1 )
	{
		//Loop through each detected particle in the spectrometer
		for ( int k = 0; k < SpectrometerEvent->GetNCandidates(); k++ )
		{
			//Get the candidate from the spectrometer event
			TRecoSpectrometerCandidate *SpectroCandidate = (TRecoSpectrometerCandidate*)SpectrometerEvent->GetCandidate(k);
			SpectroCandidate->SetEvent(SpectrometerEvent);

            //Get the three momentum of the candidate from the spectrometer before the magnet
            TVector3 momentum = SpectroCandidate->GetThreeMomentumBeforeMagnet();
            //Get the three position of the candidate from the spectrometer before the magnet
            TVector3 position_start = SpectroCandidate->GetPositionBeforeMagnet();

            //Lots of calculations
            double min_dist_to_baxis_before_fiducial = ((b.fiducial_entry - position_start).Dot(b.beam_axis.Cross(momentum))) / b.beam_axis.Cross(momentum).Mag();
            TVector3 closest_point_from_baxis_before_fiducial = ClosestPointOnVectorToOtherVector(position_start, momentum, b.fiducial_entry, b.beam_axis);
            TVector3 closest_point_of_beam_approach_before_fiducial = ClosestPointOnVectorToOtherVector(b.fiducial_entry, b.beam_axis, position_start, momentum);
            TVector3 dist_to_baxis_before_fiducial = - closest_point_of_beam_approach_before_fiducial + closest_point_from_baxis_before_fiducial;

            double min_dist_to_baxis_after_fiducial =  ((b.fiducial_entry - position_start).Dot(b.beam_axis_rotated.Cross(momentum))) / (b.beam_axis_rotated.Cross(momentum)).Mag();
            TVector3 closest_point_from_baxis_after_fiducial = ClosestPointOnVectorToOtherVector(position_start, momentum, b.fiducial_entry, b.beam_axis_rotated);
            TVector3 closest_point_of_beam_approach_after_fiducial = ClosestPointOnVectorToOtherVector(b.fiducial_entry, b.beam_axis_rotated, position_start, momentum);
            TVector3 dist_to_baxis_after_fiducial = -closest_point_of_beam_approach_after_fiducial + closest_point_from_baxis_after_fiducial;

            /////////////////////
            // Geometric Stuff //
            /////////////////////

            int decay_area = 0; //Initialise variable to determine where the particle decayed from. 0 is not in the beam, 1 is before the fiducial, 2 is in the fiducial.
            //Finds whether the muon comes from before or after the first magnet
            if(closest_point_from_baxis_before_fiducial[2] <= 104000 &&
              (closest_point_from_baxis_after_fiducial[2] < 104000 || abs(dist_to_baxis_before_fiducial.Mag()) < abs(dist_to_baxis_after_fiducial.Mag())))
            {
                decay_area = 1;
            }
            else if ( closest_point_from_baxis_after_fiducial[2] >= 104000 /*104000*/ && closest_point_from_baxis_after_fiducial[2] <= 166000 )
            {
                decay_area = 2;
            }
            double reco_xz_angle = SpectroCandidate->GetSlopeXBeforeMagnet();
            double reco_yz_angle = SpectroCandidate->GetSlopeYBeforeMagnet();

            ////////////////////////
            // Missing mass stuff //
            ////////////////////////

            //Calculate energy of kaon
            double kaon_mass = 493.667;
            double kaon_energy = sqrt((75e3*75e3) + (kaon_mass*kaon_mass));
            TLorentzVector kaon_momentum;
            kaon_momentum[2] = 75e3; //Set z component to beam energy
            kaon_momentum[3] = kaon_energy; //Set time component to kaon energy
            kaon_momentum.RotateY( -BeamAngleFromZAxis ); //Transform to beam frame

            //Calculate energy of muon
            double muon_mass = 105.6583715;
            double muon_energy = sqrt((momentum.Mag2()) + (muon_mass*muon_mass));
            TLorentzVector muon_momentum;
            muon_momentum.SetVect(momentum); //Set xyz components to muon 3-momentum
            muon_momentum[3] = muon_energy; //Set time component to muon energy

            //Calculate missing mass squared
            TLorentzVector missing_mass = kaon_momentum - muon_momentum;


            ///////////////////////////////////
            // Plot Histos with Restrictions //
            ///////////////////////////////////
            //Attempts to select only the muon
            if  ( SpectroCandidate->GetCharge() == 1 &&   //Positive Charge
                  SpectrometerEvent->GetNCandidates() == 1 && //Single Detection In Spectrometer
                  decay_area == 2  &&//Decay in Fiducial Region //I think this should be decay_area > 0
                  missing_mass.Mag2() / pow( 1000, 2 ) < 2000 &&//Have Missing Mass Correct for Decay
                  LKrEvent->GetNCandidates() >= 1 && // Be detected in LKr
                  MUV3Event->GetNCandidates() >= 1// Be detected in MUV3
                  //CEDAREvent->GetNCandidates() >=1 // Be detected in CEDAR
                  //CEDAREvent->GetNHits() >=10 //Have atleast 10 hits in the CEDAR
                )
            {

                //true stuff here

                if ( MCTruthEvent -> GetNKineParts() >= 1 )
                {

                    //Get the first decay product associated with the decay event (hopefully a muon)
                    KinePart *TrueCandidate = ( KinePart* )MCTruthEvent -> GetKineParts() -> At( 1 );
                    //Get the kaon associated with the decay event
                    KinePart *KaonCandidate = ( KinePart* )MCTruthEvent -> GetKineParts() -> At( 0 );

                    TLorentzVector true_kaon_momentum = KaonCandidate->GetFinal4Momentum();

                    TVector3 true_momentum = TrueCandidate->GetMomAtCheckPoint(2).Vect();
                    TVector3 true_position_start = TrueCandidate->GetProdPos().Vect();
                    TVector3 true_position_end = TrueCandidate->GetEndPos().Vect();

                    double true_xz_angle = atan(true_momentum[2] / true_momentum[0]);
                    double true_yz_angle = atan(true_momentum[1] / true_momentum[0]);

                    //Make sure the first true event decay product is a muon
                    if (true_momentum.Mag() != 0 && TrueCandidate -> GetPDGcode() == -13 && abs(true_momentum.Theta()) > 0 )
                    {
                        true_momentum.RotateY(BeamAngleFromZAxis);
                        momentum.RotateY(BeamAngleFromZAxis);
                        for(int i=1;i<=15;i++)
                        {
                            if(true_momentum.Mag() / 1000. < 5*i && true_momentum.Mag() / 1000. >= 5*(i-1))
                            {
                                string intstr = to_string(i);
                                TString num = TString(intstr);
                                FillHisto(TString("ResolutionTemp") + num, (true_momentum.Mag() - momentum.Mag()) / 1000.0);
                            }
                            if(true_momentum[0] / 1000. < 0.026666666666666666*i-0.2 && true_momentum[0] / 1000. >= 0.026666666666666666*(i-1) -0.2)
                            {
                                string intstr = to_string(i);
                                TString num = TString(intstr);
                                FillHisto(TString("ResolutionTempX") + num, (true_momentum[0] - momentum[0]) / 1000.0);
                            }
                            if(true_momentum[1] / 1000. < 0.026666666666666666*i -0.2 && true_momentum[1] / 1000. >= 0.026666666666666666*(i-1) -0.2)
                            {
                                string intstr = to_string(i);
                                TString num = TString(intstr);
                                FillHisto(TString("ResolutionTempY") + num, (true_momentum[1] - momentum[1]) / 1000.0);
                            }
                            if(true_momentum[2] / 1000. < 5*i && true_momentum[2] / 1000. >= 5*(i-1))
                            {
                                string intstr = to_string(i);
                                TString num = TString(intstr);
                                FillHisto(TString("ResolutionTempZ") + num, (true_momentum[2] - momentum[2]) / 1000.0);
                            }
                            //cout << "PARAMETER: " << param << endl;
                        }

                        true_momentum.RotateY(-BeamAngleFromZAxis);
                        momentum.RotateY(-BeamAngleFromZAxis);

                        for(int i=1;i<=15;i++)
                        {
                            if(true_xz_angle < pi/15*i - pi / 2 && true_xz_angle >= pi/15*(i-1) - pi / 2 )
                            {
                                string intstr = to_string(i);
                                TString num = TString(intstr);
                                FillHisto(TString("xzResolutionTemp") + num, ( true_xz_angle - reco_xz_angle ) );
                            }
                            if(true_yz_angle < pi/15*i - pi / 2 && true_yz_angle >= pi/15*(i-1) - pi / 2 )
                            {
                                string intstr = to_string(i);
                                TString num = TString(intstr);
                                FillHisto(TString("yzResolutionTemp") + num, ( true_yz_angle - reco_yz_angle ));
                            }
                        }

                        ////////////////////////////////////
                        // Calculate Momentum Resolutions //
                        ////////////////////////////////////
                        true_momentum.RotateY(BeamAngleFromZAxis);
                        momentum.RotateY(BeamAngleFromZAxis);
                        TVector3 ResolutionTemp = true_momentum - momentum;
                        FillHisto( "SpectrometerXMomentumResolution",  ResolutionTemp[0] );
                        FillHisto( "SpectrometerYMomentumResolution", ResolutionTemp[1] );
                        FillHisto( "SpectrometerZMomentumResolution", ResolutionTemp[2] / 1000. );
                        FillHisto( "SpectrometerMomentumResolution", ResolutionTemp.Mag() / 1000. );

                        FillHisto( "SpectrometerXMomentumResolutionAgainstXMomentum", true_momentum[0], ResolutionTemp[0] );
                        FillHisto( "SpectrometerYMomentumResolutionAgainstYMomentum", true_momentum[1], ResolutionTemp[1] );
                        FillHisto( "SpectrometerZMomentumResolutionAgainstZMomentum", true_momentum[2] / 1000, ResolutionTemp[2] / 1000. );
                        FillHisto( "SpectrometerMomentumResolutionAgainstMomentum",  true_momentum.Mag() / 1000, ResolutionTemp.Mag() / 1000. );


                        FillHisto( "KaonEndingPosition", true_position_end[2] / 1000., true_position_end[0] );
                        FillHisto( "TrueMomentumHist",   true_momentum.Mag() / 1000. );
                        FillHisto( "TruexMomentumHist",  true_momentum[0] / 1000. );
                        FillHisto( "TrueyMomentumHist",  true_momentum[1] / 1000. );
                        FillHisto( "TruezMomentumHist",  true_momentum[2] / 1000. );


                        FillHisto( "TrueTransverseMomentumHist",     true_momentum.Perp() );
                        FillHisto( "TrueAzimuthalMomentumHist",      true_momentum.Phi() );
                        FillHisto( "TruePolarMomentumHist",          true_momentum.Theta() );
                        FillHisto( "TrueEnergyVsAzimuthal",          true_momentum.Phi(),   true_momentum.Mag()  / 1000. );
                        FillHisto( "TrueEnergyVsPolar",              true_momentum.Theta(), true_momentum.Mag()  / 1000. );
                        FillHisto( "TrueTranverseEnergyVsAzimuthal", true_momentum.Phi(),   true_momentum.Perp() / 1000. );

                        FillHisto( "CompareTrueMomentumHist",  true_momentum.Mag() / 1000. );
                        FillHisto( "CompareTruexMomentumHist", true_momentum[0] / 1000. );
                        FillHisto( "CompareTrueyMomentumHist", true_momentum[1] / 1000. );
                        FillHisto( "CompareTruezMomentumHist", true_momentum[2] / 1000. );

                        FillHisto( "CompareTrueTransverseMomentumHist", true_momentum.Perp() );
                        FillHisto( "CompareTrueAzimuthalMomentumHist",  true_momentum.Phi() );
                        FillHisto( "CompareTruePolarMomentumHist",      true_momentum.Theta() );

                        FillHisto( "CompareTrueEnergyVsAzimuthal",          true_momentum.Phi(),   true_momentum.Mag()  / 1000. );
                        FillHisto( "CompareTrueEnergyVsPolar",              true_momentum.Theta(), true_momentum.Mag()  / 1000. );
                        FillHisto( "CompareTrueTranverseEnergyVsAzimuthal", true_momentum.Phi(),   true_momentum.Perp() / 1000. );

                        FillHisto( "TrueProductionPosition", true_position_start[2] / 1000., true_position_start[0] );

                        FillHisto( "TrueProductionPositionX", true_position_start[0]);
                        FillHisto( "TrueProductionPositionY", true_position_start[1]);
                        FillHisto( "TrueProductionPositionZ", true_position_start[2]/ 1000. );

                        true_momentum.RotateY(-BeamAngleFromZAxis);
                        momentum.RotateY(-BeamAngleFromZAxis);

                        FillHisto( "TrueXZAngle",   true_xz_angle );
                        FillHisto( "TrueYZAngle",   true_yz_angle );


                        TLorentzVector true_muon_momentum = TrueCandidate->GetInitial4Momentum();
                        //TLorentzVector TrueMuonMomentum = TrueCandidate->GetMomAtCheckPoint(2);
                        TLorentzVector true_missing_mass = true_kaon_momentum - true_muon_momentum;

                        FillHisto("TrueMissingMass", true_missing_mass.Mag2() / ( pow( 1000, 2 ) ) );
                    }
                }



                //Reco stuff

                FillHisto("XZAngle", reco_xz_angle);
                FillHisto("YZAngle", reco_yz_angle);
                TVector3 Direction;
                Direction[0] = 1 / tan( reco_xz_angle );
                Direction[1] = 1 / tan( reco_yz_angle );
                Direction[2] = 1;
                Direction = Direction.Unit(); //Normalize
                TVector3 ReconstructedMomentum = Direction * momentum.Mag();

                ReconstructedMomentum.RotateY( BeamAngleFromZAxis );
                FillHisto( "ReconstructedMomentumHist", ReconstructedMomentum.Mag() / 1000. );
                FillHisto( "ReconstructedxMomentumHist", ReconstructedMomentum[0] / 1000. );
                FillHisto( "ReconstructedyMomentumHist", ReconstructedMomentum[1] / 1000. );
                FillHisto( "ReconstructedzMomentumHist", ReconstructedMomentum[2] / 1000. );
                ReconstructedMomentum.RotateY( -BeamAngleFromZAxis );

                FillHisto( "MissingMass", missing_mass.Mag2() / pow( 1000, 2) );

                momentum.RotateY(BeamAngleFromZAxis);   //Rotate the reference frame to be along the beam

                FillHisto( "MomentumHist",  momentum.Mag() / 1000. );
                FillHisto( "xMomentumHist", momentum[0]    / 1000. );
                FillHisto( "yMomentumHist", momentum[1]    / 1000. );
                FillHisto( "zMomentumHist", momentum[2]    / 1000. );



                FillHisto( "TransverseMomentumHist", momentum.Perp()  );
                FillHisto( "AzimuthalMomentumHist",  momentum.Phi()   );
                FillHisto( "PolarMomentumHist",      momentum.Theta() );

                FillHisto( "EnergyVsAzimuthal",          momentum.Phi(),   momentum.Mag()  / 1000. );
                FillHisto( "EnergyVsPolar",              momentum.Theta(), momentum.Mag()  / 1000. );
                FillHisto( "TranverseEnergyVsAzimuthal", momentum.Phi(),   momentum.Perp() / 1000. );


                momentum.RotateY(-BeamAngleFromZAxis);  //Rotate reference frame back to along the detector
                if (decay_area == 1)    ///This is never true because the first if contains it to decay_area == 2.
                {
                    FillHisto( "ClosestPointFromBeamAxis",   closest_point_from_baxis_before_fiducial.Mag() / 1000. );
                    FillHisto( "ClosestxPointFromBeamAxis",  closest_point_from_baxis_before_fiducial[0] );
                    FillHisto( "ClosestyPointFromBeamAxis",  closest_point_from_baxis_before_fiducial[1] );
                    FillHisto( "ClosestzPointFromBeamAxis",  closest_point_from_baxis_before_fiducial[2] / 1000. );
                    FillHisto( "MinimumDistanceToBeamAxis",  min_dist_to_baxis_before_fiducial );
                    FillHisto( "ClosestDistanceToBeamAxis",  dist_to_baxis_before_fiducial.Mag() );
                    FillHisto( "ClosestxDistanceToBeamAxis", dist_to_baxis_before_fiducial[0] );
                    FillHisto( "ClosestyDistanceToBeamAxis", dist_to_baxis_before_fiducial[1] );
                    FillHisto( "ClosestzDistanceToBeamAxis", dist_to_baxis_before_fiducial[2] );
                    FillHisto( "DecayPoisition",             closest_point_from_baxis_before_fiducial[2] / 1000. , closest_point_from_baxis_before_fiducial[0] );
		    FillHisto("MissingMassVsZPosition",      closest_point_from_baxis_before_fiducial[2] / 1000. , missing_mass.Mag2()/ pow(1000,2));

                }
                if (decay_area == 2)
                {
                    FillHisto( "ClosestPointFromBeamAxis",   closest_point_from_baxis_after_fiducial.Mag() / 1000. );
                    FillHisto( "ClosestxPointFromBeamAxis",  closest_point_from_baxis_after_fiducial[0] );
                    FillHisto( "ClosestyPointFromBeamAxis",  closest_point_from_baxis_after_fiducial[1] );
                    FillHisto( "ClosestzPointFromBeamAxis",  closest_point_from_baxis_after_fiducial[2] / 1000. );
                    FillHisto( "MinimumDistanceToBeamAxis",  min_dist_to_baxis_after_fiducial );
                    FillHisto( "ClosestDistanceToBeamAxis",  dist_to_baxis_after_fiducial.Mag() );
                    FillHisto( "ClosestxDistanceToBeamAxis", dist_to_baxis_after_fiducial[0] );
                    FillHisto( "ClosestyDistanceToBeamAxis", dist_to_baxis_after_fiducial[1] );
                    FillHisto( "ClosestzDistanceToBeamAxis", dist_to_baxis_after_fiducial[2] );
                    FillHisto( "DecayPoisition",             closest_point_from_baxis_after_fiducial[2] / 1000. , closest_point_from_baxis_after_fiducial[0] );
		    FillHisto("MissingMassVsZPosition",	     closest_point_from_baxis_after_fiducial[2] / 1000. , missing_mass.Mag2()/ pow(1000,2));

                }
            }
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

    int n = 15;
    double x[n],xx[n],xy[n],xz[n], y[n], yx[n], yy[n], yz[n], ry[n], ryx[n], ryy[n], ryz[n],ey[n],eyx[n],eyy[n],eyz[n],ery[n],eryx[n],eryy[n],eryz[n],anglexz[n],angleyz[n];


    for(int i=1;i<=15;i++)
    {
        x[i-1] = 5*i;
        xx[i-1] = 0.026666666666666666*i-0.2;
        xy[i-1] = 0.026666666666666666*i-0.2;
        xz[i-1] = 5*i;
        Anglexz[i-1] = pi/15*i - pi / 2
        angleyz[i-1] = pi/15*i - pi / 2;

        string intstr = to_string(i);
        TString num = TString(intstr);
        TH1* res = fHisto.GetHisto(TString("ResolutionTemp") + num);
        TH1* resx = fHisto.GetHisto(TString("ResolutionTempX") + num);
        TH1* resy = fHisto.GetHisto(TString("ResolutionTempY") + num);
        TH1* resz = fHisto.GetHisto(TString("ResolutionTempZ") + num);
        TH1* angle_xzres = fHisto.GetHisto(TString("xzResolutionTemp") + num);
        TH1* angle_yzres = fHisto.GetHisto(TString("yzResolutionTemp") + num);

        double integral = res->Integral();
        double integral_temp = res->Integral(1,500);
        double xmin = res->GetXaxis()->GetXmin();
        double xmax = res->GetXaxis()->GetXmax();
        double xbin = abs(xmin) > abs(xmax) ? abs(xmax) : abs(xmin);
        int binx1 = res->GetXaxis()->FindBin(-xbin);
        int binx2 = res->GetXaxis()->FindBin(xbin);
        integral_temp = res->Integral(binx1,binx2);
        int res_limits = 0;
        while(abs(integral_temp) > abs(integral) * 0.98)
        {
            integral_temp = res->Integral(binx1+res_limits,binx2-res_limits);
            res_limits++;
        }
        double fitmin = res->GetXaxis()->GetBinLowEdge(binx1+res_limits);
        double fitmax = res->GetXaxis()->GetBinUpEdge(binx2-res_limits);

        double integralx = resx->Integral();
        double integralx_temp = resx->Integral(1,500);
        double xminx = resx->GetXaxis()->GetXmin();
        double xmaxx = resx->GetXaxis()->GetXmax();
        double xbinx = abs(xminx) > abs(xmaxx) ? abs(xmaxx) : abs(xminx);
        int binx1x = resx->GetXaxis()->FindBin(-xbinx);
        int binx2x = resx->GetXaxis()->FindBin(xbinx);
        integralx_temp = resx->Integral(binx1x,binx2x);
        int resx_limits = 0;
        while(abs(integralx_temp) > abs(integralx) * 0.98)
        {
            integralx_temp = resx->Integral(binx1x+resx_limits,binx2x-resx_limits);
            resx_limits++;
        }
        double fitminx = resx->GetXaxis()->GetBinLowEdge(binx1x+resx_limits);
        double fitmaxx = resx->GetXaxis()->GetBinUpEdge(binx2x-resx_limits);


        double integraly = resy->Integral();
        double integraly_temp = resy->Integral(1,500);
        double xminy = resy->GetXaxis()->GetXmin();
        double xmaxy = resy->GetXaxis()->GetXmax();
        double xbiny = abs(xminy) > abs(xmaxy) ? abs(xmaxy) : abs(xminy);
        int binx1y = resy->GetXaxis()->FindBin(-xbiny);
        int binx2y = resy->GetXaxis()->FindBin(xbiny);
        integraly_temp = resy->Integral(binx1y,binx2y);
        int resy_limits = 0;
        while(abs(integraly_temp) > abs(integraly) * 0.98)
        {
            integraly_temp = resy->Integral(binx1y+resy_limits,binx2y-resy_limits);
            resy_limits++;
        }
        double fitminy = resy->GetXaxis()->GetBinLowEdge(binx1y+resy_limits);
        double fitmaxy = resy->GetXaxis()->GetBinUpEdge(binx2y-resy_limits);


        double integralz = resz->Integral();
        double integralz_temp = resz->Integral(1,500);
        double xminz = resz->GetXaxis()->GetXmin();
        double xmaxz = resz->GetXaxis()->GetXmax();
        double xbinz = abs(xminz) > abs(xmaxz) ? abs(xmaxz) : abs(xminz);
        int binx1z = resz->GetXaxis()->FindBin(-xbinz);
        int binx2z = resz->GetXaxis()->FindBin(xbinz);
        integralz_temp = resz->Integral(binx1z,binx2z);
        int resz_limits = 0;
        while(abs(integralz_temp) > abs(integralz) * 0.98)
        {
            integralz_temp = resz->Integral(binx1z+resz_limits,binx2z-resz_limits);
            resz_limits++;
        }
        double fitminz = resz->GetXaxis()->GetBinLowEdge(binx1z+resz_limits);
        double fitmaxz = resz->GetXaxis()->GetBinUpEdge(binx2z-res_limits);

        double angle_xzintegral = angle_xzres->Integral();
        double angle_xzintegral_temp = angle_xzres->Integral(1,500);
        double angle_xzxmin = angle_xzres->GetXaxis()->GetXmin();
        double angle_xzxmax = angle_xzres->GetXaxis()->GetXmax();
        double angle_xzxbin = abs(angle_xzxmin) > abs(angle_xzxmax) ? abs(angle_xzxmax) : abs(angle_xzxmin);
        int angle_xzbinx1 = angle_xzres->GetXaxis()->FindBin(-angle_xzxbin);
        int angle_xzbinx2 = angle_xzres->GetXaxis()->FindBin(angle_xzxbin);
        angle_xzintegral_temp = angle_xzres->Integral(angle_xzbinx1,angle_xzbinx2);
        int angle_xzres_limits = 0;
        while(abs(angle_xzintegral_temp) > abs(angle_xzintegral) * 0.98)
        {
            angle_xzintegral_temp = angle_xzres->Integral(angle_xzbinx1+angle_xzres_limits,angle_xzbinx2-angle_xzres_limits);
            angle_xzres_limits++;
        }
        double angle_xzfitmin = angle_xzres->GetXaxis()->GetBinLowEdge(angle_xzbinx1+angle_xzres_limits);
        double angle_xzfitmax = angle_xzres->GetXaxis()->GetBinUpEdge(angle_xzbinx2-angle_xzres_limits);


        double angle_yzintegral = angle_yzres->Integral();
        double angle_yzintegral_temp = angle_yzres->Integral(1,500);
        double angle_yzxmin = angle_yzres->GetXaxis()->GetXmin();
        double angle_yzxmax = angle_yzres->GetXaxis()->GetXmax();
        double angle_yzxbin = abs(angle_yzxmin) > abs(angle_yzxmax) ? abs(angle_yzxmax) : abs(angle_yzxmin);
        int angle_yzbinx1 = angle_yzres->GetXaxis()->FindBin(-angle_yzxbin);
        int angle_yzbinx2 = angle_yzres->GetXaxis()->FindBin(angle_yzxbin);
        angle_yzintegral_temp = angle_yzres->Integral(angle_yzbinx1,angle_yzbinx2);
        int angle_yzres_limits = 0;
        while(abs(angle_yzintegral_temp) > abs(angle_yzintegral) * 0.98)
        {
            angle_yzintegral_temp = angle_yzres->Integral(angle_yzbinx1+angle_yzres_limits,angle_yzbinx2-angle_yzres_limits);
            angle_yzres_limits++;
        }
        double angle_yzfitmin = angle_yzres->GetXaxis()->GetBinLowEdge(angle_yzbinx1+angle_yzres_limits);
        double angle_yzfitmax = angle_yzres->GetXaxis()->GetBinUpEdge(angle_yzbinx2-angle_yzres_limits);



        res->Fit("gaus","Q","",fitmin,fitmax);
        resx->Fit("gaus","Q","",fitminx,fitmaxx);
        resy->Fit("gaus","Q","",fitminy,fitmaxy);
        resz->Fit("gaus","Q","",fitminz,fitmaxz);
        angle_xzres->Fit("gaus","Q","",angle_xzfitmin,angle_xzfitmax);
        angle_yzres->Fit("gaus","Q","",angle_yzfitmin,angle_yzfitmax);
        res->Draw();
        resx->Draw();
        resy->Draw();
        resz->Draw();
        angle_xzres->Draw();
        angle_yzres->Draw();
        TF1* fit = res->GetFunction("gaus");
        TF1* fitx = resx->GetFunction("gaus");
        TF1* fity = resy->GetFunction("gaus");
        TF1* fitz = resz->GetFunction("gaus");
        TF1* angle_xzfit = angle_xzres->GetFunction("gaus");
        TF1* angle_yzfitz = angle_yzres->GetFunction("gaus");

        y[i-1] = fit->GetParameter(2);
        ey[i-1] = fit->GetParError(2);
        yx[i-1] = fitx->GetParameter(2);
        eyx[i-1] = fitx->GetParError(2);
        yy[i-1] = fity->GetParameter(2);
        eyy[i-1] = fity->GetParError(2);
        yz[i-1] = fitz->GetParameter(2);
        eyz[i-1] = fitz->GetParError(2);
        angle_xzy[i-1] = angle_xzfit->GetParameter(2);
        eangle_xzy[i-1] = angle_xzfit->GetParError(2);
        angle_yzy[i-1] = angle_yzfitz->GetParameter(2);
        eangle_yzy[i-1] = angle_yzfitz->GetParError(2);


        ry[i-1] = y[i-1] / x[i-1] ;
        ery[i-1] = abs(ey[i-1] / x[i-1]) ;
        ryx[i-1] = yx[i-1] / xx[i-1];
        eryx[i-1] = abs(eyx[i-1] / xx[i-1]) ;
        ryy[i-1] = yy[i-1] / xy[i-1];
        eryy[i-1] = abs(eyy[i-1] / xy[i-1]) ;
        ryz[i-1] = yz[i-1] / xz[i-1];
        eryz[i-1] = abs(eyz[i-1] / xz[i-1]) ;
        rangle_xzy = angle_xzy[i-1] / anglexz[i-1];
        erangle_xzy = abs(eangle_xzy[i-1] / Anglexz[i-1]);
        rangle_yzy = angle_yzy[i-1] / angleyz[i-1];
        erangle_yzy = abs(eangle_yzy[i-1] / angleyz[i-1]);
    }

    TGraphErrors* graph = new TGraphErrors(n,x,y,0,ey);
    BookHisto("MtmResolutionVsMtm", graph);
    TGraphErrors* graphx = new TGraphErrors(n,xx,yx,0,eyx);
    BookHisto("XMtmResolutionVsXMtm", graphx);
    TGraphErrors* graphy = new TGraphErrors(n,xy,yy,0,eyy);
    BookHisto("YMtmResolutionVsYMtm", graphy);
    TGraphErrors* graphz = new TGraphErrors(n,xz,yz,0,eyz);
    BookHisto("ZMtmResolutionVsZMtm", graphz);
    TGraphErrors* angle_xzgraph = new TGraphErrors(n,anglexz,angle_xzy,0,eangle_xzy);
    BookHisto("angle_xzResolutionVsangle_xz", angle_xzgraph);
    TGraphErrors* angle_yzgraph = new TGraphErrors(n,angleyz,angle_yzy,0,eangle_yzy);
    BookHisto("angle_yzResolutionVsangle_yz", angle_yzgraph);



    TGraphErrors* rgraph = new TGraphErrors(n,x,ry,0,ery);
    BookHisto("RelativeMtmResolutionVsMtm", rgraph);
    //rgraph->Fit("[0]+[1]x");
    TGraphErrors* rgraphx = new TGraphErrors(n,xx,ryx,0,eryx);
    BookHisto("RelativeXMtmResolutionVsXMtm", rgraphx);
    //rgraphx->Fit("[0]+[1]x");
    TGraphErrors* rgraphy = new TGraphErrors(n,xy,ryy,0,eryy);
    BookHisto("RelativeYMtmResolutionVsYMtm", rgraphy);
    //rgraphy->Fit("[0]+[1]x");
    TGraphErrors* rgraphz = new TGraphErrors(n,xz,ryz,0,eryz);
    BookHisto("RelativeZMtmResolutionVsZMtm", rgraphz);
    //rgraphz->Fit("[0]+[1]x");
    TGraphErrors* rangle_xzgraph = new TGraphErrors(n,anglexz,rangle_xzy,0,erangle_xzy);
    BookHisto("Relativeangle_xzResolutionVsangle_xz", rangle_xzgraph);

    TGraphErrors* rangle_yzgraph = new TGraphErrors(n,anglexz,rangle_yzy,0,erangle_yzy);
    BookHisto("Relativeangle_yzResolutionVsangle_yz", rangle_yzgraph);


    SaveAllPlots();

}

void spectro2::DrawPlot()
{
    DrawAllPlots();
    //SaveAllPlotsPDF();
}
