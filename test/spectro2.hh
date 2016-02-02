#ifndef spectro_HH
#define spectro_HH

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include "MCSimple.hh"
#include "DetectorAcceptance.hh"
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <string>
#include <cstdlib>
#include <DetectorAcceptance.hh>


class TH1I;
class TH2F;
class TGraph;
class TTree;
class TF1;

const int       NumberOfBins = 2000, SpeedOfLight = 299792458;

const double    BeamAngleFromZAxis = -1.2E-3, pi = 3.141592653589793,speed_of_light_2 = 89875517873681764, kaon_mass = 493.667, kaon_mass_2 = 243707.106889, muon_mass = 105.6583715, muon_mass_2 = 11163.691468, speed_of_light_inverse = 1 / SpeedOfLight,
                refractive_index = 1.160183663671307728521380769209707, refractive_index_inverse = 1 / refractive_index, RICH_focal_length = 17, RICH_focal_length_inverse = 1 / RICH_focal_length,
                LKr1minx = 0, Lkr1maxx = 1000, Lkr1miny = -1260, LKr1maxy = 1260, LKr1minz = 240388, Lkr1maxz = 243222,
                LKr2minx = -0, Lkr2maxx = -1225 , Lkr2miny = -1260, LKr2maxy = 1260, LKr2minz = 240388, Lkr2maxz = 243222,
                LKrhole1minx = -118, Lkrhole1maxx = 118 , Lkrhole1miny = -118, LKrhole1maxy = 118, LKrhole1minz = 240388, Lkrhole1maxz = 240399,
                LKrhole2minx = -118, Lkrhole2maxx = 118 , Lkrhole2miny = -118, LKrhole2maxy = 118, LKrhole2minz = 243221, Lkrhole2maxz = 243221;

class beam
{
    public:
        beam()
          : fiducial_entry(0,0,104000), beam_axis(0,0,1), beam_axis_rotated(0,0,1)
        {
            beam_axis_rotated.RotateY(-BeamAngleFromZAxis);
        }

        TVector3 fiducial_entry, beam_axis, beam_axis_rotated;
};

class spectro2 : public NA62Analysis::Analyzer
{
	public:
		spectro2(NA62Analysis::Core::BaseAnalysis *ba);
		void InitHist();
		void InitOutput();
		void DefineMCSimple();
		void Process(int iEvent);

		void ProcessSpectrometer();
		void ProcessLKr();
		void ProcessMUV3();

		void StartOfBurstUser();
		void EndOfBurstUser();
		void StartOfRunUser();
		void EndOfRunUser();
		void PostProcess();
		void DrawPlot();
		void SaveAllPlotsPDF();
		void SaveHistPDF();
		void CreateHist1D(TString name, TString title = "Hist Title", int nbins = NumberOfBins, double low = 0, double high = 0);
		void SetHistAxisLabels(TString name, TString xlabel = "X Axis", TString ylabel = "Y Axis");
	protected:
        beam b;

};
#endif
