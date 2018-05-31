#ifndef _BCME_PARTICLE_1D_HPP
#define _BCME_PARTICLE_1D_HPP

#include "cme_particle_1d.hpp"
#include "misc/randomer.hpp"

namespace {

    template <typename PotentialType, typename InttableMgrType>
        class BCME_Particle_1D final : public CME_Particle_1D<PotentialType>
    {
        public:
            using potential_t = PotentialType;
            using inttable_mgr_t = InttableMgrType;

        public:
            BCME_Particle_1D(double X, double V, double MASS, double KT,
                    int SURFACE,
                    double NUCLEAR_FRIC,
                    potential_t& POTENTIAL,
                    inttable_mgr_t& INTTABLE_MGR) noexcept :
                CME_Particle_1D<PotentialType>(X, V, MASS, KT, SURFACE, NUCLEAR_FRIC, POTENTIAL),
                inttable_mgr(INTTABLE_MGR)
                {
                }

            ~BCME_Particle_1D() noexcept = default;

        public:
            double get_N() {
                return this->surface;
            }

        private:
            void cal_force_fric(double& force, double& fric) {
                double fBCME;
                fBCME = inttable_mgr.retrieve("fBCME", this->x);

                force = this->potential.cal_force(this->x, this->surface) + fBCME;
                fric = this->nuclear_fric;
            }

            void do_evolve(double dt) override {
                // Velocity Verlet
                const double half_dt_mass_inv(0.5 * dt / this->mass);
                double force, fric, noise_force;

                cal_force_fric(force, fric);
                noise_force = randomer::normal(0.0, sqrt(2.0 * fric * this->kT / dt));
                this->v += (force - fric * this->v + noise_force) * half_dt_mass_inv;

                this->x += this->v * dt;

                cal_force_fric(force, fric);
                this->v += (force - fric * this->v + noise_force) * half_dt_mass_inv;

                // hopping
                this->hopper(dt);
            }

        private:
            inttable_mgr_t& inttable_mgr;
    };

};

#endif // _BCME_PARTICLE_1D_CPP
