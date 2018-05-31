#ifndef _PARA_HPP
#define _PARA_HPP

#include "macros.hpp"
#include <cassert>
#include <cstdlib>
#include <string>
#include <fstream>
#include "ioer.hpp"

namespace
{
    using std::fstream;
    using std::string;
    using ioer::keyval;
    using ioer::input_t;
    using ioer::output_t;
    using ioer::STDOUT;

    struct Para
    {
        int M;
        int Ne;
        double W;
        double gamma;
        double band_spread_over_gamma;
        double K;
        double mass;
        double omega;
        double deltaG;
        double nuclear_fric;
        double g;
        double a;
        double b;
        double kT;
        long Nstep;
        double dt;
        double dtq;
        int Anastep;
        int Ntraj;
        double x0avg;
        double v0avg;
        int surf0;
        double Er;
        double kT0;
        double thermal_tau;
        int bath_relax_Nstep;

        unsigned int random_seed;

        Para() = default;

        Para(const string& fname) {
            loadpara(fname);
        }

        void loadpara(const string& fname) {
            input_t inp(fname);

            inp.read_text();
            inp.extract_para("M", M);			    	    	
            inp.extract_para("W", W);	    		        	
            inp.extract_para("gamma", gamma);	    		    	
            inp.extract_para("band_spread_over_gamma", band_spread_over_gamma);
            inp.extract_para("K", K);
            inp.extract_para("mass", mass);	    		        	
            inp.extract_para("omega", omega);	    		    	
            inp.extract_para("deltaG", deltaG);	    		    	
            inp.extract_para("nuclear_fric", nuclear_fric);	    		
            inp.extract_para("g", g);	    		        	
            inp.extract_para("a", a);	
            inp.extract_para("b", b);	    		        	
            inp.extract_para("kT", kT);	    	            	
            inp.extract_para("Nstep", Nstep);	    		    	
            inp.extract_para("dt", dt);	           		    	
            inp.extract_para("dtq", dtq);	           		    	
            inp.extract_para("Anastep", Anastep);	    		    	
            inp.extract_para("Ntraj", Ntraj);	    		    	
            inp.extract_para("x0avg", x0avg);	    		    	
            inp.extract_para("v0avg", v0avg);	    		    	
            inp.extract_para("surf0", surf0);	    		    	
            inp.extract_para("kT0", kT0);  	    		    	
            inp.extract_para("thermal_tau", thermal_tau);   		    	
            inp.extract_para("bath_relax_Nstep", bath_relax_Nstep);
            inp.extract_para("random_seed", random_seed);   		    	
            inp.close();

            assert(static_cast<int>(dt / dtq) * dtq == dt);

            Ne = M / 2;
            Er = 0.5 * mass * omega * omega * g * g;
        }

        void showpara(output_t& outp=STDOUT)
        {
            outp.tabout("#", " M", M);
            outp.tabout("#", " W", W);
            outp.tabout("#", " gamma", gamma);
            outp.tabout("#", " band_spread_over_gamma", band_spread_over_gamma);
            outp.tabout("#", " K", K);
            outp.tabout("#", " mass", mass);
            outp.tabout("#", " omega", omega);
            outp.tabout("#", " deltaG", deltaG);
            outp.tabout("#", " nuclear_fric", nuclear_fric);
            outp.tabout("#", " g", g);
            outp.tabout("#", " a", a);
            outp.tabout("#", " b", b);
            outp.tabout("#", " kT", kT);
            outp.tabout("#", " Nstep", Nstep);
            outp.tabout("#", " dt", dt);
            outp.tabout("#", " dtq", dtq);
            outp.tabout("#", " Anastep", Anastep);
            outp.tabout("#", " Ntraj", Ntraj);
            outp.tabout("#", " x0avg", x0avg);
            outp.tabout("#", " v0avg", v0avg);
            outp.tabout("#", " surf0", surf0);
            outp.tabout("#", " kT0", kT0);
            outp.tabout("#", " Ne", Ne);
            outp.tabout("#", " Er", Er);
            outp.tabout("#", " thermal_tau", thermal_tau);
            outp.tabout("#", " bath_relax_Nstep", bath_relax_Nstep);
            outp.tabout("#", " random_seed", random_seed);
        }
    };
};

#endif // _PARA_HPP
