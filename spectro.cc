#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include <math.h>
#include "spectro.hh"
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
spectro::spectro( Core::BaseAnalysis *ba ) : Analyzer( ba, "spectro" )
{
	RequestTree( "GigaTracker", new TRecoGigaTrackerEvent );
	RequestTree( "Spectrometer", new TRecoSpectrometerEvent );

}

void spectro::InitOutput()
{

}

void spectro::InitHist()
{
    // Reconstructed momentum histogram
    TH1D* h1 = new TH1D( "MomentumHist","Momentum", NumberOfBins, 0, 80 );
    h1 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h1 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h1 );

    // True momentum histogram plotted with same bins as reconstructed momentum
    TH1D* h101 = new TH1D( "CompareTrueMomentumHist", " True Momentum", NumberOfBins, 0, 80 );
    h101 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h101 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h101 );

    // Reconstructed x momentum histogram
    TH1D* h2 = new TH1D( "xMomentumHist", "Compare x Momentum", NumberOfBins, -0.3, 0.3 );
    h2 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h2 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h2 );

    // True x momentum histogram plotted with same bins as reconstructed x momentum
    TH1D* h102 = new TH1D( "CompareTruexMomentumHist", "Compare True x Momentum", NumberOfBins, -0.3, 0.3 );
    h102 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h102 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h102 );

    // Reconstructed y momentum histogram
    TH1D* h3 = new TH1D( "yMomentumHist", "y Momentum", NumberOfBins, -0.3, 0.3 );
    h3 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h3 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h3 );

    // True y momentum histogram plotted with same bins as reconstructed y momentum
    TH1D* h103 = new TH1D( "CompareTrueyMomentumHist", "Compare True y Momentum", NumberOfBins, -0.3, 0.3 );
    h103 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h103 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h103 );

    // Reconstructed z momentum histogram
    TH1D* h4 = new TH1D( "zMomentumHist", "z Momentum", NumberOfBins, 0, 80 );
    h4 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h4 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h4 );


    // True z momentum histogram plotted with same bins as reconstructed z momentum
    TH1D* h104 = new TH1D( "CompareTruezMomentumHist", "Compare True z Momentum", NumberOfBins, 0, 80 );
    h104 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h104 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h104 );

    // Reconstructed transverse momentum histogram
    TH1D* h13 = new TH1D( "TransverseMomentumHist", "Transverse Momentum", NumberOfBins, 0, 300 );
    h13 -> GetXaxis() -> SetTitle( "Momentum MeV" );
    h13 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h13 );

    // True transverse momentum histogram plotted with same bins as reconstructed transverse momentum
    TH1D* h113 = new TH1D( "CompareTrueTransverseMomentumHist", " Compare True Transverse Momentum", NumberOfBins, 0, 300 );
    h113 -> GetXaxis() -> SetTitle( "Momentum MeV" );
    h113 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h113 );

    // Reconstructed azimuthal angle histogram
    TH1D* h5 = new TH1D( "AzimuthalMomentumHist", "Azimuthal Angle", NumberOfBins, -pi, pi );
    h5 -> GetXaxis() -> SetTitle( "Azimuthal Angle Radians" );
    h5 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h5 );

    // True azimuthal angle histogram with same bins as reconstructed azimuthal angle histogram
    TH1D* h105 = new TH1D( "CompareTrueAzimuthalMomentumHist", "Compare True Azimuthal Angle", NumberOfBins, -pi, pi );
    h105 -> GetXaxis() -> SetTitle( "Azimuthal Angle Radians" );
    h105 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h105 );

    // Reconstructed polar angle histogram
    TH1D* h6 = new TH1D( "PolarMomentumHist", "Polar Angle", NumberOfBins, 0, 0.017 );
    h6 -> GetXaxis() -> SetTitle( "Polar Angle Radians" );
    h6 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h6 );

    // True polar angle histogram with same bins as reconstructed polar angle histogram
    TH1D* h106 = new TH1D( "CompareTruePolarMomentumHist", "Compare True Polar Angle", NumberOfBins, 0, 0.017 );
    h106 -> GetXaxis() -> SetTitle( "Polar Angle Radians" );
    h106 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h106 );


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











    TH1D* h18 = new TH1D( "TrueMomentumHist", " True Momentum", NumberOfBins, 0, 0 );
    h18 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h18 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h18 );

    TH1D* h10 = new TH1D( "TruexMomentumHist", "True x Momentum", NumberOfBins, 0, 0 );
    h10 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h10 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h10 );


    TH1D* h11 = new TH1D( "TrueyMomentumHist", "True y Momentum", NumberOfBins, 0, 0 );
    h11 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h11 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h11 );

    TH1D* h12 = new TH1D( "TruezMomentumHist", "True z Momentum", NumberOfBins, 0, 0 );
    h12 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h12 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h12 );

    TH1D* h19 = new TH1D( "TrueTransverseMomentumHist", " True Transverse Momentum", NumberOfBins, 0, 0 );
    h19 -> GetXaxis() -> SetTitle( "Momentum MeV" );
    h19 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h19 );

    TH1D* h20 = new TH1D( "TrueAzimuthalMomentumHist", " True Azimuthal Angle", NumberOfBins, 0, 0 );
    h20 -> GetXaxis() -> SetTitle( "Azimuthal Angle Radians" );
    h20 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h20 );

    TH1D* h21 = new TH1D( "TruePolarMomentumHist", " True Polar Angle", NumberOfBins, 0, 0 );
    h21 -> GetXaxis() -> SetTitle( "Polar Angle Radians" );
    h21 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h21 );

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

    TH1D* h25 = new TH1D( "ClosestPointFromBeamAxis", "Closest Point From Beam Axis", NumberOfBins, 0, 300 );
    h25 -> GetXaxis() -> SetTitle( "m" );
    h25 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h25 );

    TH1D* h26 = new TH1D( "ClosestxPointFromBeamAxis", " Closest x Point From Beam Axis", NumberOfBins, -200, 200 );
    h26 -> GetXaxis() -> SetTitle( "mm" );
    h26 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h26 );

    TH1D* h27 = new TH1D( "ClosestyPointFromBeamAxis", " Closest y Point From Beam Axis", NumberOfBins, -100, 100 );
    h27 -> GetXaxis() -> SetTitle( "mm" );
    h27 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h27 );

    TH1D* h28 = new TH1D( "ClosestzPointFromBeamAxis", " Closest z Point From Beam Axis", NumberOfBins, 0, 300 );
    h28 -> GetXaxis() -> SetTitle( "m" );
    h28 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h28 );

    TH1D* h29 = new TH1D( "MinimumDistanceToBeamAxis", " Minimum Distance From Beam Axis", NumberOfBins, -100, 100 );
    h29 -> GetXaxis() -> SetTitle( "mm" );
    h29 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h29 );

    TH1D* h30 = new TH1D( "ClosestDistanceToBeamAxis", " Minimum Distance(TWO) From Beam Axis", NumberOfBins, 0, 100 );
    h30 -> GetXaxis() -> SetTitle( "mm" );
    h30 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h30 );

    TH1D* h31 = new TH1D( "ClosestxDistanceToBeamAxis", " Minimum xDistance From Beam Axis", NumberOfBins, -100, 100 );
    h31 -> GetXaxis() -> SetTitle( "mm" );
    h31 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h31 );

    TH1D* h32 = new TH1D( "ClosestyDistanceToBeamAxis", " Minimum yDistance From Beam Axis", NumberOfBins, -100, 100 );
    h32 -> GetXaxis() -> SetTitle( "mm" );
    h32 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h32 );

    TH1D* h33 = new TH1D( "ClosestzDistanceToBeamAxis", " Minimum zDistance From Beam Axis", NumberOfBins, -0.1, 0.1 );
    h33 -> GetXaxis() -> SetTitle( "mm" );
    h33 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h33 );

    TH2D* h34 = new TH2D( "DecayPoisition", " Decay Poisition ", NumberOfBins, 0, 270, NumberOfBins, -180 , 180 );
    h34 -> GetXaxis() -> SetTitle( "z metres" );
    h34 -> GetYaxis() -> SetTitle( "x mm" );
    BookHisto( h34 );

    TH1D* h35 = new TH1D( "ClosestPointOfMuon", " ClosestPointOfMuon", NumberOfBins, 0, 0 );
    h35 -> GetXaxis() -> SetTitle( "m" );
    h35 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h35 );

    TH1D* h36 = new TH1D( "ClosestxPointOfMuon", " ClosestxPointOfMuon", NumberOfBins, 0, 0 );
    h36 -> GetXaxis() -> SetTitle( "mm" );
    h36 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h36 );

    TH1D* h37 = new TH1D( "ClosestyPointOfMuon", " ClosestyPointOfMuon", NumberOfBins, 0, 0 );
    h37 -> GetXaxis() -> SetTitle( "mm" );
    h37 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h37 );

    TH1D* h38 = new TH1D( "ClosestzPointOfMuon", " ClosestzPointOfMuon", NumberOfBins, 0, 0 );
    h38 -> GetXaxis() -> SetTitle( "m" );
    h38 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h38 );

    TH1D* h39 = new TH1D( "ClosestTimeOfMuon", " ClosestTimeOfMuon", NumberOfBins, 0, 0 );
    h39 -> GetXaxis() -> SetTitle( "s" );
    h39 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h39 );

    TH1D* h40 = new TH1D( "ClosestTimeOfKaon", " ClosestTimeOfKaon", NumberOfBins, 0, 0 );
    h40 -> GetXaxis() -> SetTitle( "s" );
    h40 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h40 );

    TH1D* h41 = new TH1D( "ClosestDistanceOfMuonToKaon", " ClosestDistanceOfMuonToKaon", NumberOfBins, 0, 0 );
    h41 -> GetXaxis() -> SetTitle( "mm" );
    h41 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h41 );

    TH1D* h42 = new TH1D( "ClosestxDistanceOfMuonToKaon", " ClosestxDistanceOfMuonToKaon", NumberOfBins, 0, 0 );
    h42 -> GetXaxis() -> SetTitle( "mm" );
    h42 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h42 );

    TH1D* h43 = new TH1D( "ClosestyDistanceOfMuonToKaon", " ClosestyDistanceOfMuonToKaon", NumberOfBins, 0, 0 );
    h43 -> GetXaxis() -> SetTitle( "mm" );
    h43 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h43 );

    TH1D* h44 = new TH1D( "ClosestzDistanceOfMuonToKaon", " ClosestzDistanceOfMuonToKaon", NumberOfBins, 0, 0 );
    h44 -> GetXaxis() -> SetTitle( "mm" );
    h44 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h44 );

    TH1D* h45 = new TH1D( "ClosestSpaceTimeInterval", " ClosestSpaceTimeInterval", NumberOfBins, 0, 0 );
    h45 -> GetXaxis() -> SetTitle( "m^2" );
    h45 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h45 );


    TH2D* h46 = new TH2D( "ParticleProductionPosition", "StartingPositionOfJetParticles", NumberOfBins, 0, 0, NumberOfBins, 0, 0 );
    h46 -> GetXaxis() -> SetTitle( "m" );
    h46 -> GetYaxis() -> SetTitle( "mm" );
    BookHisto( h46 );

    TH2D* h47 = new TH2D( "KaonEndingPosition", "EndingPositionOfKaons", NumberOfBins, 0, 0, NumberOfBins, 0, 0 );
    h47 -> GetXaxis() -> SetTitle( "m" );
    h47 -> GetYaxis() -> SetTitle( "mm" );
    BookHisto( h47 );

    TH1D* h48 = new TH1D( "TrueProuductionPositionx", " MuonProductionPointx", NumberOfBins, 0, 0 );
    h48 -> GetXaxis() -> SetTitle( "m" );
    h48 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h48 );

    TH1D* h49 = new TH1D( "TrueProuductionPositiony", " MuonProductionPointy", NumberOfBins, 0, 0 );
    h49 -> GetXaxis() -> SetTitle( "m" );
    h49 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h49 );

    TH1D* h50 = new TH1D( "TrueProuductionPositionz", " MuonProductionPointz", NumberOfBins, 0, 0 );
    h50 -> GetXaxis() -> SetTitle( "m" );
    h50 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h50 );

    TH2D* h51 = new TH2D( "TrueProuductionPosition", "MuonProductionPoint", NumberOfBins, 0, 0, NumberOfBins, 0, 0 );
    h51 -> GetXaxis() -> SetTitle( "m" );
    h51 -> GetYaxis() -> SetTitle( "mm" );
    BookHisto( h51 );

    TH1D* h52 = new TH1D("MissingMass", "Missing Mass Squared", NumberOfBins, 0, 0);
    h52 -> GetXaxis() -> SetTitle( "Missing Mass Squared" );
    h52 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h52 );

	TH1D* h53 = new TH1D("TrueMissingMass", "True Missing Mass Squared", 1000, 0, 0);
    h53 -> GetXaxis() -> SetTitle( "True Missing Mass Squared" );
    h53 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h53 );

    TH1D* h54 = new TH1D("SpectrometerxMomentumResolution", "Spectrometer x Momentum Resolution", NumberOfBins, 0, 0);
    h54 -> GetXaxis() -> SetTitle( "x Momentum MeV" );
    h54 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h54 );

    TH1D* h55 = new TH1D("SpectrometeryMomentumResolution", "Spectrometer y Momentum Resolution", NumberOfBins, 0, 0);
    h55 -> GetXaxis() -> SetTitle( "y Momentum MeV" );
    h55 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h55 );

    TH1D* h56 = new TH1D("SpectrometerzMomentumResolution", "Spectrometer z Momentum Resolution", NumberOfBins, 0, 0);
    h56 -> GetXaxis() -> SetTitle( "z Momentum GeV" );
    h56 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h56 );

    TH1D* h57 = new TH1D("SpectrometerMomentumResolution", "Spectrometer Momentum Resolution", NumberOfBins, 0, 0);
    h57 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h57 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h57 );

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

void spectro::DefineMCSimple()
{

}

void spectro::StartOfRunUser()
{

}

void spectro::StartOfBurstUser()
{

}

void spectro::SaveAllPlotsPDF()
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


void spectro::Process( int iEvent )
{
	//Add event to list of events
	event* reco_event = new event();
	reco_events.push_back(reco_event);

	event* true_event = new event();
	true_events.push_back(true_event);

	if( fMCSimple.fStatus == MCSimple::kMissing ){printIncompleteMCWarning( iEvent );return;}
	if( fMCSimple.fStatus == MCSimple::kEmpty ){printNoMCWarning();return;}

	TRecoSpectrometerEvent *SpectrometerEvent = ( TRecoSpectrometerEvent* )GetEvent( "Spectrometer" );
	TVector3    	ParticleTrueThreeMomentum, TrueParticleStartingThreePosition, TrueParticleEndingThreePosition,
               	 	ClosestPointFromBeamAxis, BeamPointFiducialEntry, DistanceToBeamAxis, ClosestPointOfBeamApproached, KaonThreeMomentum, KaonThreePositionGTK1,
                    ClosestPointFromBeamAxisBeforeFiducial, ClosestPointFromBeamAxisAfterFiducial, ClosestPointOfBeamApproachedBeforeFiducial, ClosestPointOfBeamApproachedAfterFiducial,
                    DistanceToBeamAxisBeforeFiducial, DistanceToBeamAxisAfterFiducial;
    TLorentzVector  TrueMuonFourMomentum, TrueFourMomentum, TrueParticleStartingFourPosition, TrueParticleEndingFourPosition, KaonFourMomentum, TrueKaonMomentum;
	double          MinimumDistanceToBeamAxis,MinimumDistanceToBeamAxisAfterFiducial,MinimumDistanceToBeamAxisBeforeFiducial, CheckIfEventCanBeMatchedToBeam = 0;

	//Check to see if an event was detected
	if( SpectrometerEvent -> GetNCandidates() >= 1 )
	{
		//Loop through each detected particle in the event
		for ( int k = 0; k < SpectrometerEvent -> GetNCandidates(); k++ )
		{
			//Get the candidate from the spectrometer event
			TRecoSpectrometerCandidate *SpectroCandidate = ( TRecoSpectrometerCandidate* )SpectrometerEvent->GetCandidate( k );
			SpectroCandidate->SetEvent( SpectrometerEvent );
			//Create a new particle and add it to the event
			particle* p = new particle();
			reco_event->add_particle(p);
			//Set the properties of the particle
			p->momentum = SpectroCandidate->GetThreeMomentumBeforeMagnet();
			p->position_start = SpectroCandidate->GetPositionBeforeMagnet();
			p->time_start = SpectroCandidate -> GetTime();
			p->charge = SpectroCandidate -> GetCharge();

            MinimumDistanceToBeamAxisBeforeFiducial = ( ( b.fiducial_entry - p->position_start ).Dot( b.beam_axis.Cross( p->momentum ) ) ) / ( b.beam_axis.Cross( p->momentum ) ).Mag() ;
            ClosestPointFromBeamAxisBeforeFiducial = ClosestPointOnVectorToOtherVector( p->position_start, p->momentum, b.fiducial_entry, b.beam_axis );
            ClosestPointOfBeamApproachedBeforeFiducial = ClosestPointOnVectorToOtherVector( b.fiducial_entry, b.beam_axis, p->position_start, p->momentum );
            DistanceToBeamAxisBeforeFiducial = -ClosestPointOfBeamApproachedBeforeFiducial + ClosestPointFromBeamAxisBeforeFiducial;
			MinimumDistanceToBeamAxisAfterFiducial = ( ( b.fiducial_entry - p->position_start ).Dot( b.beam_axis_rotated.Cross( p->momentum ) ) ) / ( b.beam_axis_rotated.Cross( p->momentum ) ).Mag() ;
			ClosestPointFromBeamAxisAfterFiducial = ClosestPointOnVectorToOtherVector( p->position_start, p->momentum, b.fiducial_entry, b.beam_axis_rotated );
			ClosestPointOfBeamApproachedAfterFiducial = ClosestPointOnVectorToOtherVector( b.fiducial_entry, b.beam_axis_rotated, p->position_start, p->momentum );
			DistanceToBeamAxisAfterFiducial = -ClosestPointFromBeamAxisAfterFiducial + ClosestPointOfBeamApproachedAfterFiducial;

            if ( ClosestPointFromBeamAxisBeforeFiducial( 2 ) <= 102000 && ( ClosestPointFromBeamAxisAfterFiducial( 2 ) < 102000 || abs( DistanceToBeamAxisBeforeFiducial.Mag() ) < abs ( DistanceToBeamAxisAfterFiducial.Mag() ) ) )
			{
				p->origin = ClosestPointFromBeamAxisBeforeFiducial;
				ClosestPointOfBeamApproached = ClosestPointOfBeamApproachedBeforeFiducial;
				p->beam_distance = DistanceToBeamAxisBeforeFiducial;
				p->minimum_beam_distance = MinimumDistanceToBeamAxisBeforeFiducial;
				CheckIfEventCanBeMatchedToBeam = 1;
            }

            else if ( ClosestPointFromBeamAxisAfterFiducial( 2 ) >= 102000 )
            {
                p->origin = ClosestPointFromBeamAxisAfterFiducial;
                ClosestPointOfBeamApproached = ClosestPointOfBeamApproachedAfterFiducial;
                p->beam_distance = DistanceToBeamAxisAfterFiducial;
                p->minimum_beam_distance = MinimumDistanceToBeamAxisAfterFiducial;
                CheckIfEventCanBeMatchedToBeam = 1;
            }

            /*Remove events that aren't a single positive particle being detected in the spectrometer (as this is not k->munu), Charge == 1 gets rid of 29 events, then && Candidates == 1 gets rid of another 12  */
			if ( p->charge == 1 && SpectrometerEvent->GetNCandidates() == 1 )
			{
				p->plot_momentum = true;

				double KaonMass = 493.667;
				double KaonEnergy = sqrt((75e3*75e3) + (KaonMass*KaonMass));
				TLorentzVector KaonMomentum;
				KaonMomentum[2] = 75e3;
				KaonMomentum[3] = KaonEnergy;
				KaonMomentum.RotateY( -BeamAngleFromZAxis );

				double MuonMass = 105.6583715;
				double MuonEnergy = sqrt((p->momentum.Mag2()) + (MuonMass*MuonMass));
				TLorentzVector MuonMomentum;
				MuonMomentum.SetVect(p->momentum);
				MuonMomentum[3] = MuonEnergy;

				TLorentzVector MissingMass = KaonMomentum - MuonMomentum;
				double MissingMass2 = MissingMass.Mag2();

                FillHisto( "MissingMass", MissingMass2 / ( pow( 1000, 2 ) ) );

				if ( CheckIfEventCanBeMatchedToBeam == 1 )
				{
                    p->plot_beam_distance = true;
					if ( abs( MissingMass2 ) / ( pow( 1000, 2 ) ) < 0.1 ) //This is very likey k->munu. One detected candidate in spec, correct charge, came from beam, kinematics consistent within resolution to this process.
					{
                        p->name = "Muon";
                        p->PDGcode = -13;
                        p->kmunu = true;
					}
				}
			}
		}
	}

	TRecoGigaTrackerEvent *GTKEvent = ( TRecoGigaTrackerEvent* )GetEvent( "GigaTracker" );
	//cout << "GTKEVENTNUMBER:" << GTKEvent -> GetNCandidates() << " " << "SPECTROMETEREVENTNUMBER:" << SpectrometerEvent -> GetNCandidates() << " ";
	//FOR SOME REASON  GTKEvent -> GetNCandidates() ALWAYS RETURNS 0 //
	if( GTKEvent -> GetNCandidates() >= 1 ) //Loop through every distinguishable detected event
	{
		for ( int k = 0; k < GTKEvent -> GetNHits(); k++)
		{
			// Create the Kaon and add it to the event
			particle* kaon = new particle();
			reco_event->add_particle(kaon);
			TRecoGigaTrackerCandidate *KaonCandidate = ( TRecoGigaTrackerCandidate* )GTKEvent->GetCandidate( 0 );
			KaonCandidate -> SetEvent( GTKEvent ); //THIS LINE CAUSES SEGMENTATION VIOLATION
			cout << "IS IT HERE????";
			//Set the properties of the particle
			KaonFourMomentum = KaonCandidate -> GetMomentum();
			kaon->momentum = KaonFourMomentum.Vect();
			kaon->position_start = KaonCandidate -> GetPosition( 0 );
			kaon->time_start = KaonCandidate -> GetTime1();
		}
	}

	/*
	TVector3 	MuonThreePosition, MuonMomentum, KaonThreePosition, KaonMomentum, ClosestPointOfMuon, ClosestPointOfKaon, ClosestDistanceFromMuonToKaon;
	TLorentzVector 	ClosestSpaceTimePointOfMuon, ClosestSpaceTimePointOfKaon;
	double 		MuonMass = 105.6583715, MuonDetectionTime, KaonMass = 493.667, KaonDetectionTime, SpaceTimeInterval, MinimumInterval, ClosestTimeOfMuon, ClosestTimeOfKaon, 				ClosestSpaceTimeInterval, ClosestTimeFromMuonToKaon;
	int 		K[10], KaonPosition[10], KaonDecayed;
	if ( iEvent == 100999 )
	{
		for( int i = 0; i < 10; i++ )
		{
			MinimumInterval = 99999999999999;
			for ( int k = 0; k < 101000; k++ )
			{
				MuonThreePosition.SetXYZ(EventParticleThreePositionBeforeMagnet[i][0],EventParticleThreePositionBeforeMagnet[i][1],EventParticleThreePositionBeforeMagnet[i][2]);
				MuonMomentum.SetXYZ(EventThreeMomentum[i][0], EventThreeMomentum[i][1], EventThreeMomentum[i][2] );
				MuonDetectionTime = EventParticleTime[i];
				KaonThreePosition.SetXYZ(EventKaonThreePositionGTK1[k][0],EventKaonThreePositionGTK1[k][1],EventKaonThreePositionGTK1[k][2]);
				KaonMomentum.SetXYZ(EventKaonThreeMomentum[k][0], EventKaonThreeMomentum[k][1], EventKaonThreeMomentum[k][2] );
				KaonDetectionTime = EventKaonTimeAtGTK1[k];
				ClosestSpaceTimePointOfMuon = ClosestSpaceTimePointOnVectorToOtherVector(MuonThreePosition, MuonMomentum, MuonMass, MuonDetectionTime, KaonThreePosition, KaonMomentum, KaonMass, KaonDetectionTime );
				SpaceTimeInterval = ClosestSpaceTimePointOfMuon.Mag();
				if ( SpaceTimeInterval < MinimumInterval )	//For each muon, find the kaon that has the smallest spacetime interval between them
				{
					MinimumInterval = SpaceTimeInterval;
					KaonPosition[i] = k; //Note which kaon corresponds to which muon.
				}
			}
		}
		for( int j = 0; j < 10; j++ )
		{
			KaonDecayed = KaonPosition[j];
			MuonThreePosition.SetXYZ(EventParticleThreePositionBeforeMagnet[j][0],EventParticleThreePositionBeforeMagnet[j][1],EventParticleThreePositionBeforeMagnet[j][2]);
			MuonMomentum.SetXYZ(EventThreeMomentum[j][0], EventThreeMomentum[j][1], EventThreeMomentum[j][2] );
			MuonDetectionTime = EventParticleTime[j];
			KaonThreePosition.SetXYZ(EventKaonThreePositionGTK1[KaonDecayed][0],EventKaonThreePositionGTK1[KaonDecayed][1],EventKaonThreePositionGTK1[KaonDecayed][2]);
			KaonMomentum.SetXYZ(EventKaonThreeMomentum[KaonDecayed][0], EventKaonThreeMomentum[KaonDecayed][1], EventKaonThreeMomentum[KaonDecayed][2] );
			KaonDetectionTime = EventKaonTimeAtGTK1[KaonDecayed];
			ClosestSpaceTimePointOfMuon = ClosestSpaceTimePointOnVectorToOtherVector(MuonThreePosition, MuonMomentum, MuonMass, MuonDetectionTime, KaonThreePosition, KaonMomentum, KaonMass, KaonDetectionTime );
			ClosestPointOfMuon = ClosestSpaceTimePointOfMuon.Vect();
			ClosestTimeOfMuon = ClosestSpaceTimePointOfMuon( 3 );
			ClosestSpaceTimePointOfKaon = ClosestSpaceTimePointOnVectorToOtherVector(KaonThreePosition, KaonMomentum, KaonMass, KaonDetectionTime, MuonThreePosition, MuonMomentum, MuonMass, MuonDetectionTime );
			ClosestPointOfKaon = ClosestSpaceTimePointOfKaon.Vect();
			ClosestTimeOfKaon = ClosestSpaceTimePointOfKaon( 3 );
			ClosestDistanceFromMuonToKaon = ClosestPointOfKaon - ClosestPointOfMuon;
			ClosestTimeFromMuonToKaon = ClosestTimeOfKaon - ClosestTimeOfMuon;
			ClosestSpaceTimeInterval =  ClosestDistanceFromMuonToKaon.Mag2() - pow( SpeedOfLight, 2 ) * pow( ClosestTimeFromMuonToKaon, 2 );
			FillHisto( "ClosestPointOfMuon",  ClosestPointOfMuon.Mag() );
			FillHisto( "ClosestxPointOfMuon",  ClosestPointOfMuon( 0 ) );
			FillHisto( "ClosestyPointOfMuon",  ClosestPointOfMuon( 1 ) );
			FillHisto( "ClosestzPointOfMuon",  ClosestPointOfMuon( 2 ) / 1000. );
			FillHisto( "ClosestTimeOfMuon",  ClosestTimeOfMuon );
			FillHisto( "ClosestTimeOfKaon",  ClosestTimeOfKaon );
			FillHisto( "ClosestDistanceOfMuonToKaon",  ClosestDistanceFromMuonToKaon.Mag() );
			FillHisto( "ClosestxDistanceOfMuonToKaon",  ClosestDistanceFromMuonToKaon( 0 ) );
			FillHisto( "ClosestyDistanceOfMuonToKaon",  ClosestDistanceFromMuonToKaon( 1 ) );
			FillHisto( "ClosestzDistanceOfMuonToKaon",  ClosestDistanceFromMuonToKaon( 2 ) );
			FillHisto( "ClosestSpaceTimeInterval",  ClosestSpaceTimeInterval );
		}
	}
*/
	Event *MCTruthEvent = GetMCEvent();
	if ( MCTruthEvent -> GetNKineParts() >= 3 )
	{
		for( int i = 0; i < MCTruthEvent -> GetNKineParts(); i++ )
		{
			particle* true_particle = new particle();
			true_event->add_particle(true_particle);
			KinePart *TrueCandidate = ( KinePart* )MCTruthEvent -> GetKineParts() -> At( i );

			KinePart *KaonCandidate = ( KinePart* )MCTruthEvent -> GetKineParts() -> At( 0 );
			TrueKaonMomentum = KaonCandidate->GetFinal4Momentum();

			true_particle->momentum = TrueCandidate -> GetMomAtCheckPoint( 2 ).Vect();
			true_particle->position_start = TrueCandidate -> GetProdPos().Vect();
			true_particle->position_end = TrueCandidate -> GetEndPos().Vect();
			if ( i == 0 )
				FillHisto( "KaonEndingPosition",true_particle->position_end[2] / 1000., true_particle->position_end[0] );
			FillHisto( "ParticleProductionPosition",true_particle->position_start[2] / 1000., true_particle->position_start[0] );
			if ( i == 1 && true_particle->momentum.Mag() != 0 && TrueCandidate -> GetPDGcode() == -13 && abs(true_particle->momentum.Theta()) > 0 )
			{
				TLorentzVector TrueMuonMomentum = TrueCandidate->GetInitial4Momentum();
				//TLorentzVector TrueMuonMomentum = TrueCandidate->GetMomAtCheckPoint(2);
				TLorentzVector TrueMissingMass = TrueKaonMomentum - TrueMuonMomentum;
				double TrueMissingMass2 = TrueMissingMass.Mag2();

				FillHisto("TrueMissingMass", TrueMissingMass2 / ( pow( 1000, 2 ) ) );

				cout    << true_particle->position_start[0] << " "
                        << true_particle->position_start[1] << " "
						<< true_particle->position_start[2] << " "
                        << true_particle->momentum.Mag() << " "
                        << true_particle->position_end[0] << " "
                        << true_particle->position_end[1] << " "
                        << true_particle->position_end[2] << " "
                        << endl
                        << MCTruthEvent -> GetNKineParts();

				for ( int k = 0; k < MCTruthEvent -> GetNKineParts(); k++ )
				{
					KinePart *CandidateN = ( KinePart* )MCTruthEvent -> GetKineParts() -> At( k );
					cout  << CandidateN -> GetParticleName() << CandidateN -> GetPDGcode() << " ";
				}
				if( i == 1 && true_particle->position_start[2] > 101000)
				{
					cout << endl;
					if ( TrueCandidate -> GetPDGcode() == -13) //If the kaon decays straight into a muon, do this shit.
					{
                        true_particle->plot_true_kmunu = true;
					}
				}
			}
		}
	}
}

void spectro::PostProcess()
{

}

void spectro::EndOfBurstUser()
{

}

void spectro::EndOfRunUser()
{


        for ( int i = 0; i < 100000; i++ )
        {
            int NumberDetected = reco_events[i]->particles.size();
            int TrueNumber = true_events[i]->particles.size();
            if  ( NumberDetected == 1 &&  reco_events[i]->particles[0]->plot_beam_distance == true )
            {
                FillHisto( "ClosestPointFromBeamAxis", reco_events[i]->particles[0]->origin.Mag() / 1000. );
                FillHisto( "ClosestxPointFromBeamAxis", reco_events[i]->particles[0]->origin[0] );
                FillHisto( "ClosestyPointFromBeamAxis", reco_events[i]->particles[0]->origin[1] );
                FillHisto( "ClosestzPointFromBeamAxis", reco_events[i]->particles[0]->origin[2] / 1000. );
                FillHisto( "MinimumDistanceToBeamAxis", reco_events[i]->particles[0]->minimum_beam_distance );
                FillHisto( "ClosestDistanceToBeamAxis", reco_events[i]->particles[0]->beam_distance.Mag() );
                FillHisto( "ClosestxDistanceToBeamAxis", reco_events[i]->particles[0]->beam_distance[0] );
                FillHisto( "ClosestyDistanceToBeamAxis", reco_events[i]->particles[0]->beam_distance[1] );
                FillHisto( "ClosestzDistanceToBeamAxis", reco_events[i]->particles[0]->beam_distance[2] );
                FillHisto( "DecayPoisition", reco_events[i]->particles[0]->origin[2] / 1000. , reco_events[i]->particles[0]->origin[0] );
            }
            if ( NumberDetected == 1 &&  reco_events[i]->particles[0]->plot_momentum == true )
            {
                reco_events[i]->particles[0]->momentum.RotateY(BeamAngleFromZAxis); //Switch to  reference system where beam is along z axis.
				FillHisto( "MomentumHist",  reco_events[i]->particles[0]->momentum.Mag() / 1000. );
				FillHisto( "xMomentumHist", reco_events[i]->particles[0]->momentum[0] / 1000. );
				FillHisto( "yMomentumHist", reco_events[i]->particles[0]->momentum[1] / 1000. );
				FillHisto( "zMomentumHist", reco_events[i]->particles[0]->momentum[2] / 1000. );
				FillHisto( "TransverseMomentumHist",  reco_events[i]->particles[0]->momentum.Perp() );
				FillHisto( "AzimuthalMomentumHist", reco_events[i]->particles[0]->momentum.Phi() );
				FillHisto( "PolarMomentumHist", reco_events[i]->particles[0]->momentum.Theta() );
				FillHisto( "EnergyVsAzimuthal", reco_events[i]->particles[0]->momentum.Phi(), reco_events[i]->particles[0]->momentum.Mag() / 1000. );
				FillHisto( "EnergyVsPolar", reco_events[i]->particles[0]->momentum.Theta(), reco_events[i]->particles[0]->momentum.Mag() / 1000. );
				FillHisto( "TranverseEnergyVsAzimuthal", reco_events[i]->particles[0]->momentum.Phi(), reco_events[i]->particles[0]->momentum.Perp() / 1000. );
                reco_events[i]->particles[0]->momentum.RotateY( -BeamAngleFromZAxis ); //Switch back to standard reference frame.
            }
            if ( TrueNumber >= 3 && true_events[i]->particles[1]->plot_true_kmunu == true )
            {
                true_events[i]->particles[1]->momentum.RotateY(BeamAngleFromZAxis);
                FillHisto( "TrueMomentumHist",  true_events[i]->particles[1]->momentum.Mag() / 1000. );
                FillHisto( "TruexMomentumHist", true_events[i]->particles[1]->momentum[0] / 1000. );
                FillHisto( "TrueyMomentumHist", true_events[i]->particles[1]->momentum[1] / 1000. );
                FillHisto( "TruezMomentumHist", true_events[i]->particles[1]->momentum[2] / 1000. );
                FillHisto( "TrueTransverseMomentumHist",  true_events[i]->particles[1]->momentum.Perp() );
                FillHisto( "TrueAzimuthalMomentumHist", true_events[i]->particles[1]->momentum.Phi() );
                FillHisto( "TruePolarMomentumHist", true_events[i]->particles[1]->momentum.Theta() );
                FillHisto( "TrueEnergyVsAzimuthal", true_events[i]->particles[1]->momentum.Phi(), true_events[i]->particles[1]->momentum.Mag() / 1000. );
                FillHisto( "TrueEnergyVsPolar", true_events[i]->particles[1]->momentum.Theta(), true_events[i]->particles[1]->momentum.Mag() / 1000. );
                FillHisto( "TrueTranverseEnergyVsAzimuthal", true_events[i]->particles[1]->momentum.Phi(), true_events[i]->particles[1]->momentum.Perp() / 1000. );
                FillHisto( "CompareTrueMomentumHist",  true_events[i]->particles[1]->momentum.Mag() / 1000. );
                FillHisto( "CompareTruexMomentumHist", true_events[i]->particles[1]->momentum[0] / 1000. );
                FillHisto( "CompareTrueyMomentumHist", true_events[i]->particles[1]->momentum[1] / 1000. );
                FillHisto( "CompareTruezMomentumHist", true_events[i]->particles[1]->momentum[2] / 1000. );
                FillHisto( "CompareTrueTransverseMomentumHist",  true_events[i]->particles[1]->momentum.Perp() );
                FillHisto( "CompareTrueAzimuthalMomentumHist", true_events[i]->particles[1]->momentum.Phi() );
                FillHisto( "CompareTruePolarMomentumHist", true_events[i]->particles[1]->momentum.Theta() );
                FillHisto( "CompareTrueEnergyVsAzimuthal", true_events[i]->particles[1]->momentum.Phi(), true_events[i]->particles[1]->momentum.Mag() / 1000. );
                FillHisto( "CompareTrueEnergyVsPolar", true_events[i]->particles[1]->momentum.Theta(), true_events[i]->particles[1]->momentum.Mag() / 1000. );
                FillHisto( "CompareTrueTranverseEnergyVsAzimuthal", true_events[i]->particles[1]->momentum.Phi(), true_events[i]->particles[1]->momentum.Perp() / 1000. );
                FillHisto( "TrueProuductionPositionx", true_events[i]->particles[1]->position_start[0]);
                FillHisto( "TrueProuductionPositiony", true_events[i]->particles[1]->position_start[1]);
                FillHisto( "TrueProuductionPositionz", true_events[i]->particles[1]->position_start[2]/ 1000. );
                true_events[i]->particles[1]->momentum.RotateY(-BeamAngleFromZAxis);

            }
            if  ( NumberDetected == 1  )
            {
                true_events[i]->particles[1]->momentum.RotateY(BeamAngleFromZAxis);
                reco_events[i]->particles[0]->momentum.RotateY(BeamAngleFromZAxis);
                TVector3 ResolutionTemp = true_events[i]->particles[1]->momentum - reco_events[i]->particles[0]->momentum;
                FillHisto( "SpectrometerxMomentumResolution",  ResolutionTemp[0] ) ;
                FillHisto( "SpectrometeryMomentumResolution", ResolutionTemp[1] );
                FillHisto( "SpectrometerzMomentumResolution", ResolutionTemp[2] / 1000. );
                FillHisto( "SpectrometerMomentumResolution", ResolutionTemp.Mag() / 1000. );
                true_events[i]->particles[1]->momentum.RotateY(-BeamAngleFromZAxis);
                reco_events[i]->particles[0]->momentum.RotateY(-BeamAngleFromZAxis);
            }

        }

        TH1D* h58 = h101;
        BookHisto(h58);
        h58->Divide(fHisto.GetHisto("MomentumHist"));
    SaveAllPlots();
}

void spectro::DrawPlot()
{
    DrawAllPlots();
    SaveAllPlotsPDF();
}
