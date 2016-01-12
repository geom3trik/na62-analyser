#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include <math.h>
#include <list>
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

void spectro::CreateHist1D(TString name, TString title, int nbins, double low, double high)
{
    TH1D* h1 = new TH1D(name,title,nbins,low,high);
    BookHisto(h1);
}

void spectro::SetHistAxisLabels(TString name, TString xlabel, TString ylabel)
{
    TH1* h = fHisto.GetHisto(name);
    if(h != NULL)
    {
        h -> GetXaxis() -> SetTitle(xlabel);
        h -> GetYaxis() -> SetTitle(ylabel);
    }

}

void spectro::InitHist()
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

    CreateHist1D( "TrueProuductionPositionx", " MuonProductionPointx", NumberOfBins, 0, 0 );
    SetHistAxisLabels("TrueProuductionPositionx","m","Number of Entries");

    CreateHist1D( "TrueProuductionPositiony", " MuonProductionPointy", NumberOfBins, 0, 0 );
    SetHistAxisLabels("TrueProuductionPositiony","m","Number of Entries");

    CreateHist1D( "TrueProuductionPositionz", " MuonProductionPointz", NumberOfBins, 0, 0 );
    SetHistAxisLabels("TrueProuductionPositionz","m","Number of Entries");

    TH2D* h51 = new TH2D( "TrueProuductionPosition", "MuonProductionPoint", NumberOfBins, 0, 0, NumberOfBins, 0, 0 );
    h51 -> GetXaxis() -> SetTitle( "m" );
    h51 -> GetYaxis() -> SetTitle( "mm" );
    BookHisto( h51 );

    CreateHist1D("MissingMass", "Missing Mass Squared", NumberOfBins, -0.2, 0.2);
    SetHistAxisLabels("MissingMass","Missing Mass Squared, GeV /c","Number of Entries");

	CreateHist1D("TrueMissingMass", "True Missing Mass Squared", 1000, 0, 0);
    SetHistAxisLabels("TrueMissingMass","True Missing Mass Squared","Number of Entries");

    CreateHist1D("SpectrometerxMomentumResolution", "Spectrometer x Momentum Resolution", NumberOfBins, 0, 0);
    SetHistAxisLabels("SpectrometerxMomentumResolution","x Momentum MeV","Number of Entries");

    CreateHist1D("SpectrometeryMomentumResolution", "Spectrometer y Momentum Resolution", NumberOfBins, 0, 0);
    SetHistAxisLabels("SpectrometeryMomentumResolution","y Momentum MeV","Number of Entries");

    CreateHist1D("SpectrometerzMomentumResolution", "Spectrometer z Momentum Resolution", NumberOfBins, 0, 0);
    SetHistAxisLabels("SpectrometerzMomentumResolution","z Momentum GeV","Number of Entries");

    CreateHist1D("SpectrometerMomentumResolution", "Spectrometer Momentum Resolution", NumberOfBins, 0, 0);
    SetHistAxisLabels("SpectrometerMomentumResolution","Momentum GeV","Number of Entries");

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

void spectro::SaveHistPDF()
{
    // Get iterator for CanvasOrganizer map
    NA62Analysis::NA62Map<TString,CanvasOrganizer*>::type::iterator ptr;
    // Get the CanvasOrganizer map
    NA62Analysis::NA62Map<TString,CanvasOrganizer*>::type canvases = GetCanvases();
    // Loop through all CanvasOrganizer map
    for(ptr = canvases.begin(); ptr != canvases.end(); ptr++)
    {
        // Retrieve ROOT canvas and save as .pdf file
        std::cout << ptr->second->GetName() << std::endl;
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


void spectro::Process( int iEvent )
{
	//Add event to list of events
	event* reco_event = new event();
	reco_events.push_back(reco_event);
    //Add true event to list of events
	event* true_event = new event();
	true_events.push_back(true_event);

	//if( fMCSimple.fStatus == MCSimple::kMissing ){printIncompleteMCWarning( iEvent );return;}
	//if( fMCSimple.fStatus == MCSimple::kEmpty ){printNoMCWarning();return;}

    // Get the spectrometer event
	TRecoSpectrometerEvent *SpectrometerEvent = ( TRecoSpectrometerEvent* )GetEvent( "Spectrometer" );
    /// Hopefully can replace these with class member variables in particle class?
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
			p->time_start = SpectroCandidate -> GetTime() * pow( 10, 1 ) / pow( 10 , 9 );
			p->charge = SpectroCandidate -> GetCharge();
            //Calculate stuff
            MinimumDistanceToBeamAxisBeforeFiducial = ( ( b.fiducial_entry - p->position_start ).Dot( b.beam_axis.Cross( p->momentum ) ) ) / ( b.beam_axis.Cross( p->momentum ) ).Mag() ;
            ClosestPointFromBeamAxisBeforeFiducial = ClosestPointOnVectorToOtherVector( p->position_start, p->momentum, b.fiducial_entry, b.beam_axis );
            ClosestPointOfBeamApproachedBeforeFiducial = ClosestPointOnVectorToOtherVector( b.fiducial_entry, b.beam_axis, p->position_start, p->momentum );
            DistanceToBeamAxisBeforeFiducial = -ClosestPointOfBeamApproachedBeforeFiducial + ClosestPointFromBeamAxisBeforeFiducial;
			MinimumDistanceToBeamAxisAfterFiducial = ( ( b.fiducial_entry - p->position_start ).Dot( b.beam_axis_rotated.Cross( p->momentum ) ) ) / ( b.beam_axis_rotated.Cross( p->momentum ) ).Mag() ;
			ClosestPointFromBeamAxisAfterFiducial = ClosestPointOnVectorToOtherVector( p->position_start, p->momentum, b.fiducial_entry, b.beam_axis_rotated );
			ClosestPointOfBeamApproachedAfterFiducial = ClosestPointOnVectorToOtherVector( b.fiducial_entry, b.beam_axis_rotated, p->position_start, p->momentum );
			DistanceToBeamAxisAfterFiducial = -ClosestPointFromBeamAxisAfterFiducial + ClosestPointOfBeamApproachedAfterFiducial;

            //I can't remember what this is doing
            if ( ClosestPointFromBeamAxisBeforeFiducial( 2 ) <= 104000 && ( ClosestPointFromBeamAxisAfterFiducial( 2 ) < 104000 || abs( DistanceToBeamAxisBeforeFiducial.Mag() ) < abs ( DistanceToBeamAxisAfterFiducial.Mag() ) ) )
			{
				p->origin = ClosestPointFromBeamAxisBeforeFiducial;
				ClosestPointOfBeamApproached = ClosestPointOfBeamApproachedBeforeFiducial;
				p->beam_distance = DistanceToBeamAxisBeforeFiducial;
				p->minimum_beam_distance = MinimumDistanceToBeamAxisBeforeFiducial;
				CheckIfEventCanBeMatchedToBeam = 0;
            }
            //Or this
            else if ( ClosestPointFromBeamAxisAfterFiducial( 2 ) >= 104000 && ClosestPointFromBeamAxisAfterFiducial( 2 ) <= 166000 )
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
			    //Did you add this Jack?
				p->plot_momentum = true;
                //Calculate energy of kaon
				double KaonMass = 493.667;
				double KaonEnergy = sqrt((75e3*75e3) + (KaonMass*KaonMass));
				TLorentzVector KaonMomentum;
				KaonMomentum[2] = 75e3; //Set z component to beam energy
				KaonMomentum[3] = KaonEnergy; //Set time component to kaon energy
				KaonMomentum.RotateY( -BeamAngleFromZAxis ); //Transform to beam frame

                //Calculate energy of muon
				double MuonMass = 105.6583715;
				double MuonEnergy = sqrt((p->momentum.Mag2()) + (MuonMass*MuonMass));
				TLorentzVector MuonMomentum;
				MuonMomentum.SetVect(p->momentum); //Set xyz components to muon 3-momentum
				MuonMomentum[3] = MuonEnergy; //Set time component to muon energy

                //Calculate missing mass squared
				TLorentzVector MissingMass = KaonMomentum - MuonMomentum;
				double MissingMass2 = MissingMass.Mag2();

                //What is this doing?
				if ( CheckIfEventCanBeMatchedToBeam == 1 )
				{

                    p->plot_beam_distance = true;
					if ( abs( MissingMass2 ) / ( pow( 1000, 2 ) ) < 1000000) //This is very likey k->munu. One detected candidate in spec, correct charge, came from beam, kinematics consistent within resolution to this process.
					{
                        FillHisto( "MissingMass", MissingMass2 / ( pow( 1000, 2 ) ) );
                        p->name = "Muon";
                        p->PDGcode = -13;
                        p->kmunu = true;
                        ///cout << "Muon:" << p->time_start <<endl;
					}
				}
			}
		}
	}

	//I wonder if the GigaTracker works yet

	/*
	TRecoGigaTrackerEvent *GTKEvent = ( TRecoGigaTrackerEvent* )GetEvent( "GigaTracker" );
	//cout << "GTKEVENTNUMBER:" << GTKEvent -> GetNCandidates() << " " << "SPECTROMETEREVENTNUMBER:" << SpectrometerEvent -> GetNCandidates() << " ";
	//FOR SOME REASON  GTKEvent -> GetNCandidates() ALWAYS RETURNS 0 //


	if( GTKEvent -> GetNCandidates() == 1 ) //Loop through every distinguishable detected event
	{

		for ( int k = 0; k < GTKEvent -> GetNHits(); k++)
		{
			// Create the Kaon and add it to the event
			particle* kaon = new particle();
			reco_event->add_particle(kaon);
			TRecoGigaTrackerCandidate *KaonCandidate = ( TRecoGigaTrackerCandidate* )GTKEvent->GetCandidate( 0 );
			KaonCandidate -> SetEvent( GTKEvent ); //THIS LINE CAUSES SEGMENTATION VIOLATION
			//Set the properties of the particle
			KaonFourMomentum = KaonCandidate -> GetMomentum();
			kaon->momentum = KaonFourMomentum.Vect();
			kaon->position_start = KaonCandidate -> GetPosition( 0 );
			kaon->time_start = KaonCandidate -> GetTime1() / pow( 10,18 ) ;
			kaon->detected = true;
			cout << "Kaon:" << kaon->time_start << endl;

		}
	}
	*/

    //Get truth event
	Event *MCTruthEvent = GetMCEvent();
	//Check if a event was detected
	if ( MCTruthEvent -> GetNKineParts() >= 1 )
	{
	    //Loop through all candidates in the event
		for( int i = 0; i < MCTruthEvent -> GetNKineParts(); i++ )
		{
		    //Create particle and add to the truth event
			particle* true_particle = new particle();
			true_event->add_particle(true_particle);
			//Get the current particle
			KinePart *TrueCandidate = ( KinePart* )MCTruthEvent -> GetKineParts() -> At( i );
            //Get the kaon associated with the decay event
			KinePart *KaonCandidate = ( KinePart* )MCTruthEvent -> GetKineParts() -> At( 0 );
			///This should really be replaced with particle class variable
			TrueKaonMomentum = KaonCandidate->GetFinal4Momentum();
            ///Like these
            //Which was checkpoint 2 again?
			true_particle->momentum = TrueCandidate -> GetMomAtCheckPoint( 2 ).Vect();
			true_particle->position_start = TrueCandidate -> GetProdPos().Vect();
			true_particle->position_end = TrueCandidate -> GetEndPos().Vect();

			if ( i == 0 ) //i.e. the kaon
				FillHisto( "KaonEndingPosition",true_particle->position_end[2] / 1000., true_particle->position_end[0] );

			FillHisto( "ParticleProductionPosition",true_particle->position_start[2] / 1000., true_particle->position_start[0] );

			if ( i == 1 && true_particle->momentum.Mag() != 0 && TrueCandidate -> GetPDGcode() == -13 && abs(true_particle->momentum.Theta()) > 0 )
			{
			    //Loop through each particle again?
                for ( int k = 0; k < MCTruthEvent -> GetNKineParts(); k++ )
				{
					KinePart *CandidateN = ( KinePart* )MCTruthEvent -> GetKineParts() -> At( k );
					///cout  << CandidateN -> GetParticleName() << CandidateN -> GetPDGcode() << " ";
				}

                ///More stuff that needs replacing with the way we did it for the real data
				TLorentzVector TrueMuonMomentum = TrueCandidate->GetInitial4Momentum();
				//TLorentzVector TrueMuonMomentum = TrueCandidate->GetMomAtCheckPoint(2);
				TLorentzVector TrueMissingMass = TrueKaonMomentum - TrueMuonMomentum;
				double TrueMissingMass2 = TrueMissingMass.Mag2();

				FillHisto("TrueMissingMass", TrueMissingMass2 / ( pow( 1000, 2 ) ) );

                /*
				cout    << true_particle->position_start[0] << " "
                        << true_particle->position_start[1] << " "
						<< true_particle->position_start[2] << " "
                        << true_particle->momentum.Mag() << " "
                        << true_particle->position_end[0] << " "
                        << true_particle->position_end[1] << " "
                        << true_particle->position_end[2] << " "
                        << endl
                        << MCTruthEvent -> GetNKineParts();
                */

                //What is this doing?
				true_particle->momentum.RotateY(BeamAngleFromZAxis);
				if( i == 1 && true_particle->position_start[2] > 101000 && abs(true_particle->momentum[0]) <= 235 && abs(true_particle->momentum[1]) <= 235 && true_particle->momentum.Theta() <= 0.020 )
				{
                    true_particle->momentum.RotateY(-BeamAngleFromZAxis);
					///cout << endl;
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
        ///What is the purpose of this?
        double startime = 999999999999999, endtime = 0, timeofevent;
        //Loop through all of the reconstructed events
        for ( int i = 0; i < reco_events.size(); i++ )
        {
            //Number of particles in this event
            int NumberDetected = reco_events[i]->num_of_particles();
            //Number of particles in this event
            int TrueNumber = true_events[i]->num_of_particles();
            //Loop through each particle in the reco_event
            for ( int j = 0; j < NumberDetected;j++ )
            {
                if  ( NumberDetected > 0
                        && reco_events[i]->particles[j]->plot_beam_distance == true
                        && reco_events[i]->particles[j]->kmunu == true
                        && abs(reco_events[i]->particles[j]->minimum_beam_distance) < 20 )
                {
                    FillHisto( "ClosestPointFromBeamAxis", reco_events[i]->particles[j]->origin.Mag() / 1000. );
                    FillHisto( "ClosestxPointFromBeamAxis", reco_events[i]->particles[j]->origin[0] );
                    FillHisto( "ClosestyPointFromBeamAxis", reco_events[i]->particles[j]->origin[1] );
                    FillHisto( "ClosestzPointFromBeamAxis", reco_events[i]->particles[j]->origin[2] / 1000. );
                    FillHisto( "MinimumDistanceToBeamAxis", reco_events[i]->particles[j]->minimum_beam_distance );
                    FillHisto( "ClosestDistanceToBeamAxis", reco_events[i]->particles[j]->beam_distance.Mag() );
                    FillHisto( "ClosestxDistanceToBeamAxis", reco_events[i]->particles[j]->beam_distance[0] );
                    FillHisto( "ClosestyDistanceToBeamAxis", reco_events[i]->particles[j]->beam_distance[1] );
                    FillHisto( "ClosestzDistanceToBeamAxis", reco_events[i]->particles[j]->beam_distance[2] );
                    FillHisto( "DecayPoisition", reco_events[i]->particles[j]->origin[2] / 1000. , reco_events[i]->particles[j]->origin[0] );

                    timeofevent = reco_events[i]->particles[j]->time_start;
                    if ( timeofevent > endtime)
                    {
                        endtime = timeofevent;
                    }
                    else if ( timeofevent < startime )
                    {
                        startime = timeofevent;
                    }
                }

                if ( NumberDetected > 0
                        && reco_events[i]->particles[j]->plot_momentum == true
                        && reco_events[i]->particles[j]->kmunu == true
                        && abs(reco_events[i]->particles[j]->minimum_beam_distance) < 20 )
                {
                    reco_events[i]->particles[0]->momentum.RotateY(BeamAngleFromZAxis); //Switch to  reference system where beam is along z axis.
                    FillHisto( "MomentumHist",  reco_events[i]->particles[j]->momentum.Mag() / 1000. );
                    FillHisto( "xMomentumHist", reco_events[i]->particles[j]->momentum[0] / 1000. );
                    FillHisto( "yMomentumHist", reco_events[i]->particles[j]->momentum[1] / 1000. );
                    FillHisto( "zMomentumHist", reco_events[i]->particles[j]->momentum[2] / 1000. );
                    FillHisto( "TransverseMomentumHist",  reco_events[i]->particles[j]->momentum.Perp() );
                    FillHisto( "AzimuthalMomentumHist", reco_events[i]->particles[j]->momentum.Phi() );
                    FillHisto( "PolarMomentumHist", reco_events[i]->particles[j]->momentum.Theta() );
                    FillHisto( "EnergyVsAzimuthal", reco_events[i]->particles[j]->momentum.Phi(), reco_events[i]->particles[j]->momentum.Mag() / 1000. );
                    FillHisto( "EnergyVsPolar", reco_events[i]->particles[j]->momentum.Theta(), reco_events[i]->particles[j]->momentum.Mag() / 1000. );
                    FillHisto( "TranverseEnergyVsAzimuthal", reco_events[i]->particles[j]->momentum.Phi(), reco_events[i]->particles[j]->momentum.Perp() / 1000. );
                    reco_events[i]->particles[0]->momentum.RotateY( -BeamAngleFromZAxis ); //Switch back to standard reference frame.
                }

				if ( TrueNumber >= 3
                        && true_events[i]->particles[1]->plot_true_kmunu == true
                        && abs(reco_events[i]->particles[j]->minimum_beam_distance) < 20  )
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

				if  ( NumberDetected == 1
                        && true_events[i]->particles[1]->plot_true_kmunu == true    ///Why is this here?
                        && reco_events[i]->particles[0]->plot_beam_distance == true
                        && abs(reco_events[i]->particles[j]->minimum_beam_distance) < 20  )
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
        }
        /*
        for ( int i = 0; i < reco_events.size(); i ++)
        {
            for ( int j = 0; j < reco_events.size(); j ++)
            {
                for ( int k1 = 0; k1 < reco_events[i]->particles.size(); k1++ )
                {
                    if ( reco_events[i]->particles[k1]-> kmunu == true )
                    {
                        for ( int k2 = 0; k2 < reco_events[j]->particles.size(); k2++ )
                        {
                            if ( reco_events[j]->particles[k2]-> detected == true )
                            {
                                reco_events[i]->particles[k1]->closest_spacetime_point = ClosestSpaceTimePointOnVectorToOtherVector( reco_events[i]->particles[k1]->position_start, reco_events[i]->particles[k1]->momentum, 105.6583715, reco_events[i]->particles[k1]->time_start, reco_events[j]->particles[k2]-> position_start, reco_events[j]->particles[k2]->momentum, 493.667, reco_events[j]->particles[k2]->time_start);
                                reco_events[j]->particles[k2]->closest_spacetime_point = ClosestSpaceTimePointOnVectorToOtherVector( reco_events[j]->particles[k2]->position_start, reco_events[j]->particles[k2]->momentum, 493.667, reco_events[j]->particles[k2]->time_start, reco_events[i]->particles[k1]-> position_start, reco_events[i]->particles[k1]->momentum, 105.6583715, reco_events[i]->particles[k1]->time_start);
                                reco_events[i]->particles[k1]->kaon_link = j;
                                if( abs( reco_events[i]->particles[k1]->minimum_spacetime_interval.Mag() ) > abs( (reco_events[j]->particles[k2]->closest_spacetime_point - reco_events[i]->particles[k1]->closest_spacetime_point).Mag() ) )
                                {
                                    reco_events[i]->particles[k1]->minimum_spacetime_interval = reco_events[j]->particles[k2]->closest_spacetime_point - reco_events[i]->particles[k1]->closest_spacetime_point;
                                }
                            }
                        }
                    }
                }
            }
        }

        for ( int i = 0; i < reco_events.size(); i ++)
        {
            for ( int k1 = 0; k1 < reco_events[i]->particles.size(); k1++ )
            {
                if ( reco_events[i]->particles[k1]-> kmunu == true )
                {
                    FillHisto( "ClosestPointOfMuon",  reco_events[i]->particles[k1]->closest_spacetime_point.Mag() );
                    FillHisto( "ClosestxPointOfMuon",  reco_events[i]->particles[k1]->closest_spacetime_point[0] );
                    FillHisto( "ClosestyPointOfMuon",  reco_events[i]->particles[k1]->closest_spacetime_point[1] );
                    FillHisto( "ClosestzPointOfMuon",  reco_events[i]->particles[k1]->closest_spacetime_point[2] / 1000. );
                    FillHisto( "ClosestTimeOfMuon",  reco_events[i]->particles[k1]->closest_spacetime_point[3] );
                    //FillHisto( "ClosestDistanceOfMuonToKaon",  reco_events[i]->particles[k1]->minimum_spacetime_interval.Vect.Mag() );
                    FillHisto( "ClosestxDistanceOfMuonToKaon",  reco_events[i]->particles[k1]->minimum_spacetime_interval[0] );
                    FillHisto( "ClosestyDistanceOfMuonToKaon",  reco_events[i]->particles[k1]->minimum_spacetime_interval[1] );
                    FillHisto( "ClosestzDistanceOfMuonToKaon",  reco_events[i]->particles[k1]->minimum_spacetime_interval[2] );
                    FillHisto( "ClosestSpaceTimeInterval",  reco_events[i]->particles[k1]->minimum_spacetime_interval.Mag() );
                }
                else if ( reco_events[i]->particles[k1]->detected == true )
                {
                    FillHisto( "ClosestTimeOfKaon",  reco_events[i]->particles[k1]->closest_spacetime_point[3] );
                }
            }
        }
		*/

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

    SaveAllPlots();
}

void spectro::DrawPlot()
{
    DrawAllPlots();
    //SaveAllPlotsPDF();
    SaveHistPDF();
}

