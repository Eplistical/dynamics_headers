#ifndef _EF_PARTICLE_1D_HPP
#define _EF_PARTICLE_1D_HPP

#include "particle_1d_base.hpp"
#include "misc/randomer.hpp"

namespace {

    template <typename PotentialType, typename InttableMgrType>
        struct EF_Particle_1D final : public Particle_1D
    {
        public:
            using potential_t = PotentialType;
            using inttable_mgr_t = InttableMgrType;

        public:
            EF_Particle_1D(double X, double V, double MASS, double KT,
                    double NUCLEAR_FRIC,
                    potential_t& POTENTIAL,
                    inttable_mgr_t& INTTABLE_MGR) noexcept :
                Particle_1D(X, V, MASS, KT),
                nuclear_fric(NUCLEAR_FRIC), 
                potential(POTENTIAL), 
                inttable_mgr(INTTABLE_MGR)
            {
            }

            ~EF_Particle_1D() noexcept = default;

        public:
            // population on impurity
            double get_N() const noexcept {
                return inttable_mgr.retrieve("n", this->x);
            }

        private:
            void do_evolve(double dt) override {
                const double half_dt_mass_inv(0.5 * dt / this->mass);
                const double fric(nuclear_fric + inttable_mgr.retrieve("fric", this->x));
                const double noise_force(randomer::normal(0.0, sqrt(2.0 * fric * this->kT / dt)));

                this->v += (potential.cal_force(this->x, 0) + inttable_mgr.retrieve("force", this->x) - fric * this->v + noise_force) * half_dt_mass_inv;
                this->x += this->v * dt;
                this->v += (potential.cal_force(this->x, 0) + inttable_mgr.retrieve("force", this->x) - fric * this->v + noise_force) * half_dt_mass_inv;
            }

        public:
            double nuclear_fric;
            potential_t& potential;
            inttable_mgr_t& inttable_mgr;
    };

};

#endif // _EF_PARTICLE_CPP
