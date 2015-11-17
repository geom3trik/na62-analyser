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
          : name("unknown"), PDGcode(0), charge(0), time_start(0)
        {

        }

        std::string name;
        int PDGcode;
        TVector3 position_start, position_end, momentum;
        double charge, time_start;
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
        std::vector<particle*> particles;
};

class beam
{
    public:
        beam()
          : fiducial_entry(0,0,102000), beam_axis(0,0,1), beam_axis_rotated(0,0,1)
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
	protected:
        beam b;
        //Array of events
        std::vector<event> events;


};
#endif