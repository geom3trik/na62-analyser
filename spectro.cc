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


    TH1D* h1 = new TH1D( "MomentumHist","Momentum", NumberOfBins, 0, 80 );
    h1 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h1 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h1 );

    TH1D* h101 = new TH1D( "CompareTrueMomentumHist", " True Momentum", NumberOfBins, 0, 80 );
    h101 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h101 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h101 );



    TH1D* h2 = new TH1D( "xMomentumHist", "Compare x Momentum", NumberOfBins, -0.3, 0.3 );
    h2 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h2 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h2 );

    TH1D* h102 = new TH1D( "CompareTruexMomentumHist", "Compare True x Momentum", NumberOfBins, -0.3, 0.3 );
    h102 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h102 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h102 );



    TH1D* h3 = new TH1D( "yMomentumHist", "y Momentum", NumberOfBins, -0.3, 0.3 );
    h3 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h3 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h3 );

    TH1D* h103 = new TH1D( "CompareTrueyMomentumHist", "Compare True y Momentum", NumberOfBins, -0.3, 0.3 );
    h103 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h103 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h103 );



    TH1D* h4 = new TH1D( "zMomentumHist", "z Momentum", NumberOfBins, 0, 80 );
    h4 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h4 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h4 );

    TH1D* h104 = new TH1D( "CompareTruezMomentumHist", "Compare True z Momentum", NumberOfBins, 0, 80 );
    h104 -> GetXaxis() -> SetTitle( "Momentum GeV" );
    h104 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h104 );



    TH1D* h13 = new TH1D( "TransverseMomentumHist", "Transverse Momentum", NumberOfBins, 0, 300 );
    h13 -> GetXaxis() -> SetTitle( "Momentum MeV" );
    h13 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h13 );

    TH1D* h113 = new TH1D( "CompareTrueTransverseMomentumHist", " Compare True Transverse Momentum", NumberOfBins, 0, 300 );
    h113 -> GetXaxis() -> SetTitle( "Momentum MeV" );
    h113 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h113 );



    TH1D* h5 = new TH1D( "AzimuthalMomentumHist", "Azimuthal Angle", NumberOfBins, -pi, pi );
    h5 -> GetXaxis() -> SetTitle( "Azimuthal Angle Radians" );
    h5 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h5 );

    TH1D* h105 = new TH1D( "CompareTrueAzimuthalMomentumHist", "Compare True Azimuthal Angle", NumberOfBins, -pi, pi );
    h105 -> GetXaxis() -> SetTitle( "Azimuthal Angle Radians" );
    h105 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h105 );



    TH1D* h6 = new TH1D( "PolarMomentumHist", "Polar Angle", NumberOfBins, 0, 0.017 );
    h6 -> GetXaxis() -> SetTitle( "Polar Angle Radians" );
    h6 -> GetYaxis() -> SetTitle( "Number of Entries" );
    BookHisto( h6 );

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

    NA62Analysis::NA62Map<TString,CanvasOrganizer*>::type::iterator ptr;

    NA62Analysis::NA62Map<TString,CanvasOrganizer*>::type canvases = GetCanvases();

    for(ptr = canvases.begin(); ptr != canvases.end(); ptr++)
    {
        ptr->second->GetCanvas()->SaveAs(TString(ptr->first + ".pdf"));
    }
    /*
    std::vector<TString>::iterator itOrder;
    NA62Analysis::NA62Map<TString,TH1*>::type::iterator ptr1;
    NA62Analysis::NA62Map<TString,TH2*>::type::iterator ptr2;
    NA62Analysis::NA62Map<TString,TGraph*>::type::iterator ptr3;
    NA62Analysis::NA62Map<TString,CanvasOrganizer*>::type::iterator ptr4;
    CanvasOrganizer *c;

    for(itOrder=fHisto.fHistoOrder.begin(); itOrder!=fHisto.fHistoOrder.end(); itOrder++)
    {
        if((ptr1=fHisto.fHisto.find(*itOrder))!=fHisto.end()){
            c = new CanvasOrganizer(TString("c_" + fAnalyzerName + "_") + *itOrder);
            c->AddHisto(ptr1->second);
            c->Draw();
            c->GetCanvas()->SaveAs(TString(ptr1->first + ".pdf"));
            fCanvas.insert(std::make_pair(c->GetName(), c));
        }
        else if((ptr2=fHisto.fHisto2.find(*itOrder))!=fHisto.fHisto2.end()){
            c = new CanvasOrganizer(TString("c_" + fAnalyzerName + "_") + *itOrder);
            c->AddHisto(ptr2->second);
            c->Draw();
            c->GetCanvas()->SaveAs(TString(ptr2->first + ".pdf"));
            fCanvas.insert(std::make_pair(c->GetName(), c));
        }
        else if((ptr3=fGraph.find(*itOrder))!=fGraph.end()){
            c = new CanvasOrganizer(TString("c_" + fAnalyzerName + "_") + *itOrder);
            c->AddHisto(ptr3->second);
            c->Draw();
            c->GetCanvas()->SaveAs(TString(ptr3->first + ".pdf"));
            fCanvas.insert(std::make_pair(c->GetName(), c));
        }
    }

    for(ptr4=fCanvas.begin(); ptr4!=fCanvas.end(); ptr4++)
    {
        ptr4->second->Draw();
    }
    */
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
    event* new_event = new event();

        events.push_back(new_event);
	if( fMCSimple.fStatus == MCSimple::kMissing ){printIncompleteMCWarning( iEvent );return;}
	if( fMCSimple.fStatus == MCSimple::kEmpty ){printNoMCWarning();return;}
 	//Get the spectrometer event
	TRecoSpectrometerEvent *SpectrometerEvent = ( TRecoSpectrometerEvent* )GetEvent( "Spectrometer" );
	//Define ParticleThreeMomentum and ParticleTrueThreeMomentum vectors
	TVector3    ParticleTrueThreeMomentum, TrueParticleStartingThreePosition, TrueParticleEndingThreePosition,
                ClosestPointFromBeamAxis, BeamPointFiducialEntry, DistanceToBeamAxis, ClosestPointOfBeamApproached, KaonThreeMomentum, KaonThreePositionGTK1,
                ClosestPointFromBeamAxisBeforeFiducial, ClosestPointFromBeamAxisAfterFiducial, ClosestPointOfBeamApproachedBeforeFiducial, ClosestPointOfBeamApproachedAfterFiducial, DistanceToBeamAxisBeforeFiducial, DistanceToBeamAxisAfterFiducial;
    //Define TrueFourMomentum Vector
    TLorentzVector  TrueMuonFourMomentum, TrueFourMomentum, TrueParticleStartingFourPosition, TrueParticleEndingFourPosition, KaonFourMomentum;

	double 		MinimumDistanceToBeamAxis,MinimumDistanceToBeamAxisAfterFiducial,MinimumDistanceToBeamAxisBeforeFiducial, CheckIfEventCanBeMatchedToBeam = 0;
                //EventParticleTime[101000], EventParticleThreePositionBeforeMagnet[101000][3], EventThreeMomentum[101000][3], KaonTimeAtGTK1, EventKaonThreeMomentum[101000][3], EventKaonThreePositionGTK1[101000][3], EventKaonTimeAtGTK1[101000];



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
            new_event->add_particle(p);
            //Set the properties of the particle
            p->momentum = SpectroCandidate->GetThreeMomentumBeforeMagnet();
            p->position_start = SpectroCandidate->GetPositionBeforeMagnet();
            p->time_start = SpectroCandidate -> GetTime();
            p->charge = SpectroCandidate -> GetCharge();


			//ParticleThreeMomentum = SpectroCandidate->GetThreeMomentumBeforeMagnet();
			//ParticleThreePositionBeforeMagnet = SpectroCandidate -> GetPositionBeforeMagnet();

			//BeamAxisDirection = { 0, 0, 1 };
			//BeamAxisDirection.RotateY( - BeamAngleFromZAxis );	//Rotate beam so it points along the dirction it should
			//BeamPointFiducialEntry = { 0, 0, 102000 };

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
                ClosestPointFromBeamAxis = ClosestPointFromBeamAxisBeforeFiducial;
                ClosestPointOfBeamApproached = ClosestPointOfBeamApproachedBeforeFiducial;
                    DistanceToBeamAxis = DistanceToBeamAxisBeforeFiducial;
                MinimumDistanceToBeamAxis = MinimumDistanceToBeamAxisBeforeFiducial;
                CheckIfEventCanBeMatchedToBeam = 1;
            }

            else if ( ClosestPointFromBeamAxisAfterFiducial( 2 ) >= 102000 )
            {
                ClosestPointFromBeamAxis = ClosestPointFromBeamAxisAfterFiducial;
                ClosestPointOfBeamApproached = ClosestPointOfBeamApproachedAfterFiducial;
                DistanceToBeamAxis = DistanceToBeamAxisAfterFiducial;
                MinimumDistanceToBeamAxis = MinimumDistanceToBeamAxisAfterFiducial;
                CheckIfEventCanBeMatchedToBeam = 1;
            }
            /*
            else if ( ClosestPointFromBeamAxisBeforeFiducial( 2 ) <= 102000  && abs( DistanceToBeamAxisBeforeFiducial.Mag() ) < abs ( DistanceToBeamAxisAfterFiducial.Mag() ) )
            {
                ClosestPointFromBeamAxis = ClosestPointFromBeamAxisBeforeFiducial;
                ClosestPointOfBeamApproached = ClosestPointOfBeamApproachedBeforeFiducial;
                    DistanceToBeamAxis = DistanceToBeamAxisBeforeFiducial;
                MinimumDistanceToBeamAxis = MinimumDistanceToBeamAxisBeforeFiducial;
                CheckIfEventCanBeMatchedToBeam = 1;
            }
            else if ( ClosestPointFromBeamAxisBeforeFiducial( 2 ) >= 102000 && ClosestPointFromBeamAxisAfterFiducial( 2 ) > 102000 )
            {
                ClosestPointFromBeamAxis = ClosestPointFromBeamAxisAfterFiducial;
                ClosestPointOfBeamApproached = ClosestPointOfBeamApproachedAfterFiducial;
                DistanceToBeamAxis = DistanceToBeamAxisAfterFiducial;
                MinimumDistanceToBeamAxis = MinimumDistanceToBeamAxisAfterFiducial;
                CheckIfEventCanBeMatchedToBeam = 1;
            }
            else if ( ClosestPointFromBeamAxisBeforeFiducial( 2 ) <= 102000 && abs( DistanceToBeamAxisBeforeFiducial.Mag() ) > abs ( DistanceToBeamAxisAfterFiducial.Mag() ) )
            {
                ClosestPointFromBeamAxis = ClosestPointFromBeamAxisAfterFiducial;
                ClosestPointOfBeamApproached = ClosestPointOfBeamApproachedAfterFiducial;
                DistanceToBeamAxis = DistanceToBeamAxisAfterFiducial;
                MinimumDistanceToBeamAxis = MinimumDistanceToBeamAxisAfterFiducial;
                CheckIfEventCanBeMatchedToBeam = 1;
            }
            */

			//Charge = SpectroCandidate -> GetCharge();
			//TimeAtBeforeMagnet = SpectroCandidate -> GetTime();

            /*Remove events that aren't a single positive particle being detected in the spectrometer (as this is not k->munu), Charge == 1 gets rid of 29 events, then && Candidates == 1 gets rid of another 12  */
			if ( p->charge == 1 && SpectrometerEvent->GetNCandidates() == 1 )
			{
				p->momentum.RotateY(BeamAngleFromZAxis);

				FillHisto( "MomentumHist",  p->momentum.Mag() / 1000. );
				FillHisto( "xMomentumHist", p->momentum[0] / 1000. );
				FillHisto( "yMomentumHist", p->momentum[1] / 1000. );
				FillHisto( "zMomentumHist", p->momentum[2] / 1000. );
				FillHisto( "TransverseMomentumHist",  p->momentum.Perp() );
				FillHisto( "AzimuthalMomentumHist", p->momentum.Phi() );
				FillHisto( "PolarMomentumHist", p->momentum.Theta() );
				FillHisto( "EnergyVsAzimuthal", p->momentum.Phi(), p->momentum.Mag() / 1000. );
				FillHisto( "EnergyVsPolar", p->momentum.Theta(), p->momentum.Mag() / 1000. );
				FillHisto( "TranverseEnergyVsAzimuthal", p->momentum.Phi(), p->momentum.Perp() / 1000. );

				if ( CheckIfEventCanBeMatchedToBeam == 1 )
				{
					FillHisto( "ClosestPointFromBeamAxis", ClosestPointFromBeamAxis.Mag() / 1000. );
					FillHisto( "ClosestxPointFromBeamAxis", ClosestPointFromBeamAxis( 0 ) );
					FillHisto( "ClosestyPointFromBeamAxis", ClosestPointFromBeamAxis( 1 ) );
					FillHisto( "ClosestzPointFromBeamAxis", ClosestPointFromBeamAxis( 2 ) / 1000. );
					FillHisto( "MinimumDistanceToBeamAxis", MinimumDistanceToBeamAxis );
					FillHisto( "ClosestDistanceToBeamAxis", DistanceToBeamAxis.Mag() );
					FillHisto( "ClosestxDistanceToBeamAxis", DistanceToBeamAxis( 0 ) );
					FillHisto( "ClosestyDistanceToBeamAxis", DistanceToBeamAxis( 1 ) );
					FillHisto( "ClosestzDistanceToBeamAxis", DistanceToBeamAxis( 2 ) );
					FillHisto( "DecayPoisition", ClosestPointFromBeamAxis( 2 ) / 1000. , ClosestPointFromBeamAxis ( 0 ) );
				}
			}
		}
	}

    particle* kaon = new particle();
    new_event->add_particle(kaon);
            //Set the properties of the particle

    TRecoGigaTrackerEvent *GTKEvent = ( TRecoGigaTrackerEvent* )GetEvent( "GigaTracker" );
    	//cout << "GTKEVENTNUMBER:" << GTKEvent -> GetNCandidates() << " " << "SPECTROMETEREVENTNUMBER:" << SpectrometerEvent -> GetNCandidates() << " ";
	//FOR SOME REASON  GTKEvent -> GetNCandidates() ALWAYS RETURNS 0 //
	if( GTKEvent -> GetNCandidates() >= 1 ) //Loop through every distinguishable detected event
    	{
		for ( int k = 0; k < GTKEvent -> GetNHits(); k++)
		{
			TRecoGigaTrackerCandidate *KaonCandidate = ( TRecoGigaTrackerCandidate* )GTKEvent->GetCandidate( 0 );
			KaonCandidate -> SetEvent( GTKEvent ); //THIS LINE CAUSES SEGMENTATION VIOLATION
			cout << "IS IT HERE????";
			KaonFourMomentum = KaonCandidate -> GetMomentum();
			kaon->momentum = KaonFourMomentum.Vect();
			kaon->position_start = KaonCandidate -> GetPosition( 0 );
			kaon->time_start = KaonCandidate -> GetTime1();
		}
	}

/*
	TVector3 MuonThreePosition, MuonMomentum, KaonThreePosition, KaonMomentum, ClosestPointOfMuon, ClosestPointOfKaon, ClosestDistanceFromMuonToKaon;
	TLorentzVector ClosestSpaceTimePointOfMuon, ClosestSpaceTimePointOfKaon;
	double MuonMass = 105.6583715, MuonDetectionTime, KaonMass = 493.667, KaonDetectionTime, SpaceTimeInterval, MinimumInterval, ClosestTimeOfMuon, ClosestTimeOfKaon, ClosestSpaceTimeInterval, ClosestTimeFromMuonToKaon;
	int K[10], KaonPosition[10], KaonDecayed;
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
	cout << "TESTING";
 	Event *MCTruthEvent = GetMCEvent();
	if ( MCTruthEvent -> GetNKineParts() >= 3 )
	{
		for( int i = 0; i < MCTruthEvent -> GetNKineParts(); i++ )
		{
			KinePart *TrueCandidate = ( KinePart* )MCTruthEvent -> GetKineParts() -> At( i );
			//TrueMuonFourMomentum = TrueCandidate -> GetMomSpectrometerEntry();
			TrueFourMomentum = TrueCandidate -> GetMomAtCheckPoint( 2 );
			ParticleTrueThreeMomentum = TrueFourMomentum.Vect();
			ParticleTrueThreeMomentum.RotateY( BeamAngleFromZAxis );
			TrueParticleStartingFourPosition = TrueCandidate -> GetProdPos();
			TrueParticleStartingThreePosition = TrueParticleStartingFourPosition.Vect();
			TrueParticleEndingFourPosition = TrueCandidate -> GetEndPos();
			TrueParticleEndingThreePosition = TrueParticleEndingFourPosition.Vect();
			if ( i == 0 )
			{
				FillHisto( "KaonEndingPosition",TrueParticleEndingThreePosition( 2 ) / 1000., TrueParticleEndingThreePosition( 0 ) );
			}
				FillHisto( "ParticleProductionPosition",TrueParticleStartingThreePosition( 2 ) / 1000., TrueParticleStartingThreePosition( 0 ) );
			if ( ParticleTrueThreeMomentum.Mag() != 0 )
			{
				cout << TrueParticleStartingThreePosition( 0 ) << " " << TrueParticleStartingThreePosition( 1 ) << " " << TrueParticleStartingThreePosition( 2 )
				<< " " << ParticleTrueThreeMomentum.Mag() << " " << TrueParticleEndingThreePosition( 0 ) << " " << TrueParticleEndingThreePosition( 1 ) << " " <<
				TrueParticleEndingThreePosition( 2 ) << " " << endl << MCTruthEvent -> GetNKineParts();
				for ( int k = 0; k < MCTruthEvent -> GetNKineParts(); k++ )
				{
					KinePart *CandidateN = ( KinePart* )MCTruthEvent -> GetKineParts() -> At( k );
					cout  << CandidateN -> GetParticleName() << CandidateN -> GetPDGcode() << " ";
				}


				if( i == 1 && TrueParticleStartingThreePosition( 2 ) > 101000 )
				{
					if ( TrueCandidate -> GetPDGcode() == -13) //If the kaon decays straight into a muon, do this shit.
					{
						cout << endl;
						FillHisto( "TrueMomentumHist",  ParticleTrueThreeMomentum.Mag() / 1000. );
						FillHisto( "TruexMomentumHist", ParticleTrueThreeMomentum(0) / 1000. );
						FillHisto( "TrueyMomentumHist", ParticleTrueThreeMomentum(1) / 1000. );
						FillHisto( "TruezMomentumHist", ParticleTrueThreeMomentum(2) / 1000. );
						FillHisto( "TrueTransverseMomentumHist",  ParticleTrueThreeMomentum.Perp() );
						FillHisto( "TrueAzimuthalMomentumHist", ParticleTrueThreeMomentum.Phi() );
						FillHisto( "TruePolarMomentumHist", ParticleTrueThreeMomentum.Theta() );
						FillHisto( "TrueEnergyVsAzimuthal", ParticleTrueThreeMomentum.Phi(), ParticleTrueThreeMomentum.Mag() / 1000. );
						FillHisto( "TrueEnergyVsPolar", ParticleTrueThreeMomentum.Theta(), ParticleTrueThreeMomentum.Mag() / 1000. );
						FillHisto( "TrueTranverseEnergyVsAzimuthal", ParticleTrueThreeMomentum.Phi(), ParticleTrueThreeMomentum.Perp() / 1000. );
						FillHisto( "CompareTrueMomentumHist",  ParticleTrueThreeMomentum.Mag() / 1000. );
						FillHisto( "CompareTruexMomentumHist", ParticleTrueThreeMomentum(0) / 1000. );
						FillHisto( "CompareTrueyMomentumHist", ParticleTrueThreeMomentum(1) / 1000. );
						FillHisto( "CompareTruezMomentumHist", ParticleTrueThreeMomentum(2) / 1000. );
						FillHisto( "CompareTrueTransverseMomentumHist",  ParticleTrueThreeMomentum.Perp() );
						FillHisto( "CompareTrueAzimuthalMomentumHist", ParticleTrueThreeMomentum.Phi() );
						FillHisto( "CompareTruePolarMomentumHist", ParticleTrueThreeMomentum.Theta() );
						FillHisto( "CompareTrueEnergyVsAzimuthal", ParticleTrueThreeMomentum.Phi(), ParticleTrueThreeMomentum.Mag() / 1000. );
						FillHisto( "CompareTrueEnergyVsPolar", ParticleTrueThreeMomentum.Theta(), ParticleTrueThreeMomentum.Mag() / 1000. );
						FillHisto( "CompareTrueTranverseEnergyVsAzimuthal", ParticleTrueThreeMomentum.Phi(), ParticleTrueThreeMomentum.Perp() / 1000. );
						FillHisto( "TrueProuductionPositionx", TrueParticleStartingThreePosition(0) );
						FillHisto( "TrueProuductionPositiony", TrueParticleStartingThreePosition(1) );
						FillHisto( "TrueProuductionPositionz", TrueParticleStartingThreePosition(2) / 1000. );
						FillHisto( "TrueProuductionPosition", TrueParticleStartingThreePosition(2) / 1000., TrueParticleStartingThreePosition(0) ); 						
					}
				}
			}
		}
	}
        //FillHisto( "TrueMuonxMomentumHist", TrueMuonFourMomentum(0) );
        //FillHisto( "TrueMuonyMomentumHist", TrueMuonFourMomentum(1) );
        //FillHisto( "TrueMuonzMomentumHist", TrueMuonFourMomentum(2) / 1000. );
	//FillHisto( "TrueMuonEnergyHist", TrueMuonFourMomentum(3) / 1000. );

}

void spectro::PostProcess()
{

}

void spectro::EndOfBurstUser()
{

}

void spectro::EndOfRunUser()
{
    SaveAllPlots();
}

void spectro::DrawPlot()
{
    DrawAllPlots();
    SaveAllPlotsPDF();
    //TCanvas* c = new TCanvas("canvas","canvas",800,450);
    //TH1D* h = (TH1D*)fHisto.GetHisto("MomentumHist");
    //iterator->second->Draw();
    //c->SaveAs("testy.pdf");
}
