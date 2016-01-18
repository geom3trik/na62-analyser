#ifndef spectro_HH
#define spectro_HH

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include "MCSimple.hh"
#include "DetectorAcceptance.hh"
#include <TCanvas.h>
#include <string>


class TH1I;
class TH2F;
class TGraph;
class TTree;

const int NumberOfBins = 100;
const double BeamAngleFromZAxis = -1.2E-3, pi = 3.141592653589793, SpeedOfLight = 299792458;

class particle
{
    public:

        particle()
          : name("unknown"), PDGcode(0), charge(0), time_start(0),kmunu(0),plot_beam_distance(0),minimum_beam_distance(0),plot_momentum(0),plot_true_kmunu(0),detected(0),kaon_link(0),beam_link(0),LKr_link(0),MUV2_link(0)//,minimum_spacetime_interval(9999999999,9999999999,9999999999,9999999999)
        {

        }

        std::string name;
        int PDGcode, kaon_link;
        TLorentzVector closest_spacetime_point,minimum_spacetime_interval;
        TVector3 position_start, position_end, momentum, origin,beam_distance;
        double charge, time_start,minimum_beam_distance;
        bool kmunu, plot_beam_distance,plot_momentum,plot_true_kmunu,detected, beam_link,LKr_link,MUV2_link;
};

//True particle class which inherits from particle class
//So we can put 4-vector stuff that only applies to true events here
//But still store the particle in the event particle list
class true_particle : public particle
{
    public:
        true_particle()
        {

        }
};

class event
{
    public:
        void add_particle()
        {
            particles.push_back(new particle());
        }
        void add_particle(particle* p)
        {
            particles.push_back(p);
        }
        particle* operator[](int i)
        {
            return particles[i];
        }
        //Returns the number of particles in this event
        int num_of_particles()
        {
            return particles.size();
        }
        std::vector<particle*> particles;
        TLorentzVector MissingMass;
};

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

//Need food. Going for dinner; 18:37.

class spectro : public NA62Analysis::Analyzer
{
	public:
		spectro(NA62Analysis::Core::BaseAnalysis *ba);
		void InitHist();
		void InitOutput();
		void DefineMCSimple();
		void Process(int iEvent);
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
        //Array of events
        std::vector<event*> reco_events;
        std::vector<event*> true_events;


};
#endif
