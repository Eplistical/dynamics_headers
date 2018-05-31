#ifndef _INTTABLE_MGR_BASE_HPP
#define _INTTABLE_MGR_BASE_HPP

#include <map>
#include <string>
#include <vector>
#include <cassert>

namespace {

    using std::string;
    using std::map;
    using std::vector;

    struct InttableMgr_Base {
        public:
            InttableMgr_Base(const string& FNAME);
            virtual ~InttableMgr_Base() noexcept = default;

        public:
            void load_inttable();
            double retrieve(const string& key, double x);
            double retrieve_raw(const string& key, double x);

        public:
            map< string, vector<double> > inttable_dict;
            string fname;
            double xmin, xmax, dx, dx_inv;
            uint32_t Nx;

        private:
            virtual void load_inttable_impl() = 0;
            virtual double retrieve_impl(const string& key, double x) = 0;
    };

    /*
     * implementations
     */

    /*****************************************************************************/

    InttableMgr_Base
        ::InttableMgr_Base(const string& FNAME) :
            fname(FNAME)
        {
        }

    /*****************************************************************************/

    void InttableMgr_Base
        ::load_inttable() 
        {
            load_inttable_impl();
        }

    /*****************************************************************************/

    double InttableMgr_Base
        ::retrieve(const string& key, double x) 
        {
            retrieve_impl(key, x);
        }

    /*****************************************************************************/

    double InttableMgr_Base
        ::retrieve_raw(const string& key, double x)
        {
            const vector<double>& table(inttable_dict.at(key));
            const uint32_t ix(static_cast<uint32_t>((x - xmin) * dx_inv));
            assert(ix >= 0 and ix < Nx);

            if (ix == Nx - 1) {
                return table.at(ix);
            }
            else {
                const double x0(xmin + ix * dx);
                const double y0(table.at(ix)), y1(table.at(ix + 1));
                return y0 + (y1 - y0) * (x - x0) * dx_inv;
            }
        }

    /*****************************************************************************/

};

#endif // _INTTABLE_MGR_BASE_HPP
