#ifndef _FP_PARTICLE_1D_HPP
#define _FP_PARTICLE_1D_HPP

#include "particle_1d_base.hpp"
#include "misc/randomer.hpp"
#include "misc/fermi.hpp"

namespace {

    template <typename PotentialType, typename InttableMgrType>
        struct FP_Particle_1D final : public Particle_1D
    {
        public:
            using potential_t = PotentialType;
            using inttable_mgr_t = InttableMgrType;

        public:
            FP_Particle_1D(double X, double V, double MASS, double KT,
                    double NUCLEAR_FRIC,
                    potential_t& POTENTIAL,
                    inttable_mgr_t& INTTABLE_MGR) noexcept :
                Particle_1D(X, V, MASS, KT),
                nuclear_fric(NUCLEAR_FRIC), 
                potential(POTENTIAL), 
                inttable_mgr(INTTABLE_MGR)
            {
            }

            ~FP_Particle_1D() noexcept = default;

        public:
            // population on impurity
            double get_N() const noexcept {
                return inttable_mgr.retrieve("n", this->x);
            }

        private:
            void cal_force_fric(double& force, double& fric) 
            {
                double h, dhdx, f, dfde;
                h = potential.cal_h(this->x);
                dhdx = potential.cal_dh_dx(this->x);
                f = misc::fermi(h * this->kT_inv);
                dfde = -this->kT_inv * f * (1.0 - f);

                force = potential.cal_force(this->x, 0) + dhdx * h;
                fric = nuclear_fric - 1.0 / potential.cal_gamma(this->x) * pow(dhdx, 2) * dfde;
            }

            void do_evolve(double dt) override {
                const double half_dt_mass_inv(0.5 * dt / this->mass);
                double force, fric, noise_force;

                cal_force_fric(force, fric);
                noise_force = randomer::normal(0.0, sqrt(2.0 * fric * this->kT / dt));
                this->v += (force - fric * this->v + noise_force) * half_dt_mass_inv;

                this->x += this->v * dt;
                
                cal_force_fric(force, fric);
                this->v += (force - fric * this->v + noise_force) * half_dt_mass_inv;
            }

        public:
            double nuclear_fric;
            potential_t& potential;
            inttable_mgr_t& inttable_mgr;
    };

};

#endif // _FP_PARTICLE_CPP
