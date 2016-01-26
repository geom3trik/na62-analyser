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


class TH1I;
class TH2F;
class TGraph;
class TTree;
class TF1;

const int NumberOfBins = 100;
const double BeamAngleFromZAxis = -1.2E-3, pi = 3.141592653589793, SpeedOfLight = 299792458;

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
