#ifndef _CME_PARTICLE_1D_HPP
#define _CME_PARTICLE_1D_HPP

#include "particle_1d_base.hpp"
#include "misc/randomer.hpp"
#include "misc/fermi.hpp"

namespace {
    using std::mt19937;

    template <typename PotentialType>
        class CME_Particle_1D : public Particle_1D
    {
        public:
            using potential_t = PotentialType;

        public:
            CME_Particle_1D(double X, double V, 
                    double MASS, double KT,
                    int SURFACE,
                    double NUCLEAR_FRIC,
                    potential_t& POTENTIAL) noexcept :
                Particle_1D(X, V, MASS, KT),
                surface(SURFACE), nuclear_fric(NUCLEAR_FRIC),
                potential(POTENTIAL)
            {
            }

            virtual ~CME_Particle_1D() noexcept = default;

        public:
            // population on impurity
            double get_N() const noexcept {
                return surface;
            }

        private:
            void cal_force_fric(double& force, double& fric) {
                force = potential.cal_force(this->x, surface);
                fric = nuclear_fric;
            }

            virtual void do_evolve(double dt) override {
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
                hopper(dt);
            }

        protected:
            void hopper(double dt) {
                const double gamma(potential.cal_gamma(this->x));
                const double h(potential.cal_h(this->x));
                double prob;
                if (surface == 0) {
                    prob = gamma * misc::fermi(h * this->kT_inv) * dt;
                    if (randomer::rand() < prob) {
                        surface = 1;
                    }
                }
                else if (surface == 1) {
                    prob = gamma * (1.0 - misc::fermi(h * this->kT_inv)) * dt;
                    if (randomer::rand() < prob) {
                        surface = 0;
                    }
                }
            }
            
        public:
            int surface;
            double nuclear_fric;
            potential_t& potential;
    };

};

#endif // _CME_PARTICLE_1D_HPP
