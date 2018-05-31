#ifndef _INTTABLE_MGR_HPP
#define _INTTABLE_MGR_HPP

#include <map>
#include <string>
#include <vector>
#include <string>
#include <cmath>
#include "inttable_mgr_base.hpp"
#include "ioer.hpp"

namespace {

    using std::string;
    using std::map;
    using std::vector;

    struct InttableMgr final : InttableMgr_Base {
        public:
            InttableMgr(const string& FNAME);
            ~InttableMgr() noexcept = default;

        private:
            void load_inttable_impl() override;
            double retrieve_impl(const string& key, double x) override;

        public:
            double kT, gamma;
            double emin, emax, de;
            uint32_t Ne;
    };

    /*
     * implementations
     */

    /*****************************************************************************/

    InttableMgr
        ::InttableMgr(const string& FNAME) :
            InttableMgr_Base(FNAME)
    {
    }

    /*****************************************************************************/

    void InttableMgr
        ::load_inttable_impl() 
        {
            // open file
            ioer::input_t dat(this->fname);
            // read header
            dat.read(   kT, gamma, 
                    this->xmin, this->xmax, this->dx, this->Nx,
                    emin, emax, de, Ne
                    );
            this->dx_inv = 1.0 / this->dx;
            // read data
            this->inttable_dict.insert(make_pair("n", vector<double>(Nx, 0.0)));
            this->inttable_dict.insert(make_pair("force", vector<double>(Nx, 0.0)));
            this->inttable_dict.insert(make_pair("fric", vector<double>(Nx, 0.0)));
            dat.read(
                    this->inttable_dict.at("n"),
                    this->inttable_dict.at("force"),
                    this->inttable_dict.at("fric")
                    );
            // close file
            dat.close();
        }

    /*****************************************************************************/

    double InttableMgr
        ::retrieve_impl(const string& key, double x) 
        {
            if (key == "n" or key == "force" or key == "fric") {
                return this->retrieve_raw(key, x);
            }
            else if(key == "fBCME") {

            }
        }

    /*****************************************************************************/

};

#endif // _INTTABLE_MGR_BASE_HPP
