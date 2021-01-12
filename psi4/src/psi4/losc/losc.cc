#include "losc.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/factory.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/v.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libdiis/diismanager.h"

#include <string>
#include <vector>

namespace psi {
namespace losc {
void LOSC::common_init() {
    // current implementation is only for c1 symmetry.
    if (molecule()->schoenflies_symbol() != "c1")
        throw PSIEXCEPTION("LOSC only supports C1 symmetry for now!");

    if (functional_->is_x_lrc()) {
        throw PSIEXCEPTION(
            "LOSC code does not support long range corrected functional yet!");
    }

    // set up reference at first!
    const string reference = to_lower_copy(options_.get_str("reference"));
    if (reference == "rks") {
        nspin_ = 1;
        same_a_b_dens_ = true;
        same_a_b_orbs_ = true;
    } else if (reference == "uks") {
        nspin_ = 2;
        same_a_b_dens_ = false;
        same_a_b_orbs_ = false;
    } else {
        throw PSIEXCEPTION(
            "unknown reference for LOSC calculation, uks or rks?");
    }

    // electron number
    nelec_[0] = Wfn::nalpha_;
    nelec_[1] = Wfn::nbeta_;

    // ==> allocate memory for matrices in base class <==
    // These allocations should be in base class, not in LOSC class.
    // But I don't want to touch psi4's source code.
    HF::Va_ = SharedMatrix(factory_->create_matrix("Va"));
    Wfn::Fa_ = SharedMatrix(factory_->create_matrix("Fa"));
    Wfn::Ca_ = SharedMatrix(factory_->create_matrix("Ca"));
    Wfn::Da_ = SharedMatrix(factory_->create_matrix("Da"));
    Wfn::epsilon_a_ = SharedVector(factory_->create_vector());
    Wfn::epsilon_a_->set_name("alpha orbital energies");
    if (nspin_ == 1) {
        HF::Vb_ = HF::Va_;
        Wfn::Fb_ = Wfn::Fa_;
        Wfn::Cb_ = Wfn::Ca_;
        Wfn::Db_ = Wfn::Da_;
        Wfn::epsilon_b_ = Wfn::epsilon_a_;
        Wfn::epsilon_b_->set_name("beta orbital energies");
    } else {
        HF::Vb_ = SharedMatrix(factory_->create_matrix("Vb"));
        Wfn::Fb_ = SharedMatrix(factory_->create_matrix("Fb"));
        Wfn::Cb_ = SharedMatrix(factory_->create_matrix("Cb"));
        Wfn::Db_ = SharedMatrix(factory_->create_matrix("Db"));
        Wfn::epsilon_b_ = SharedVector(factory_->create_vector());
        Wfn::epsilon_b_->set_name("beta orbital energies");
    }
}

LOSC::LOSC(SharedWavefunction dfa_wfn, shared_ptr<SuperFunctional> func)
    : psi::scf::HF(dfa_wfn, func, dfa_wfn->options(), PSIO::shared_object()),
      V_{HF::Va_, HF::Vb_},
      F_{Wfn::Fa_, Wfn::Fb_},
      C_{Wfn::Ca_, Wfn::Cb_},
      D_{Wfn::Da_, Wfn::Db_},
      eig_{Wfn::epsilon_a_, Wfn::epsilon_b_} {
    common_init();
}

LOSC::LOSC(SharedWavefunction dfa_wfn, shared_ptr<SuperFunctional> func,
           Options& options, shared_ptr<PSIO> psio)
    : psi::scf::HF(dfa_wfn, func, options, PSIO::shared_object()),
      V_{HF::Va_, HF::Vb_},
      F_{Wfn::Fa_, Wfn::Fb_},
      C_{Wfn::Ca_, Wfn::Cb_},
      D_{Wfn::Da_, Wfn::Db_},
      eig_{Wfn::epsilon_a_, Wfn::epsilon_b_} {
    common_init();
}

LOSC::~LOSC() {}

void LOSC::form_C() {
    for (size_t is = 0; is < nspin_; ++is) {
        diagonalize_F(F_[is], C_[is], eig_[is]);
    }
    find_occupation();
}

void LOSC::form_D() {
    // Fully occupied COs' coefficient matrix.
    // [nbasis x nocc].
    // TODO: take care of fractional occupation number.
    SharedMatrix C_occ[2];
    for (size_t is = 0; is < nspin_; ++is) {
        const size_t nocc = nelec_[is];
        // Copy the occupied block of the orbital matrix.
        C_occ[is] = std::make_shared<Matrix>("Occupied COs", nso_, nocc);
        auto t = C_[is];
        for (int p = 0; p < nso_; ++p) {
            for (int i = 0; i < nocc; ++i) {
                C_occ[is]->set(p, i, C_[is]->get(p, i));
            }
        }
        // D = C_occ * C_occ.T
        D_[is]->gemm(false, true, 1.0, C_occ[is], C_occ[is], 0.0);
    }
}

void LOSC::form_F() {
    // F_a = H + G_a + V_ext.
    for (size_t is = 0; is < nspin_; ++is) {
        F_[is]->copy(H_);
        F_[is]->add(G_[is]);
        for (const auto& Vext : external_potentials_) {
            F_[is]->add(Vext);
        }

        if (debug_) {
            F_[is]->print();
            if (is == 1) J_->print();
            K_[is]->print();
            if (functional_->needs_xc()) {
                V_[is]->print();
            }
            G_[is]->print();
        }
    }
}

/**
 * @note
 * Calling this function will allocate and compute the following matrices:
 * `J_`: Coulomb matrix.
 * `K_`: exchange matrix.
 * `V_losc_`: LOSC effective potential.
 * `G_`: G matrix.
 */
void LOSC::form_G() {
    // Compute J and K.
    // Push the C matrix on
    vector<SharedMatrix>& C = jk_->C_left();
    C.clear();
    C.push_back(Ca_subset("SO", "OCC"));
    C.push_back(Cb_subset("SO", "OCC"));
    // Run the JK object
    jk_->compute();
    // Pull the J and K matrices off
    const vector<SharedMatrix>& J = jk_->J();
    const vector<SharedMatrix>& K = jk_->K();
    const vector<SharedMatrix>& wK = jk_->wK();

    // save J and K.
    J_ = SharedMatrix(factory_->create_matrix("J"));
    if (nspin_ == 1) {
        J_->axpy(2, J[0]);
    } else {
        J_->copy(J[0]);
        J_->add(J[1]);
    }
    if (functional_->is_x_hybrid()) {
        for (size_t is = 0; is < nspin_; ++is) {
            K_[is] = K[is];
        }
    }

    // Compute LOSC-DFA potential: V = V_dfa + V_losc.
    if (functional_->needs_xc()) {
        form_V();
    }

    const double alpha = functional_->x_alpha();
    for (size_t is = 0; is < nspin_; ++is) {
        // G_a = J_a + J_b - alpha * K_a + V.
        G_[is] = SharedMatrix(factory_->create_matrix("G"));
        G_[is]->copy(J_);
        if (functional_->is_x_hybrid()) {
            G_[is]->axpy(-alpha, K[is]);
        }
        if (functional_->needs_xc()) {
            G_[is]->add(V_[is]);
        }
#ifdef USING_BrianQC
        if (brianEnable and brianEnableDFT) {
            // BrianQC multiplies with the exact exchange factors inside the
            // Fock building, so we must not do it here
            alpha = 1.0;
            beta = 1.0;
        }
        throw PSIEXCEPTION("losc form_G() does not support BrianQC yet.")
#endif
    }
}

/**
 * @note
 * Calling this function will allocate and compute `V_losc`, the LOSC effective
 * potential matrix.
 */
void LOSC::form_V_losc() {
    for (size_t is = 0; is < nspin_; ++is) {
        V_losc_[is] = SharedMatrix(factory_->create_matrix("V_losc"));
        // TODO: Add construnction of V_losc later.
    }
}

/**
 * @note
 * Calling this function will allocate and compute the following matrices:
 * `V_losc`: LOSC effective potential matrix.
 */
void LOSC::form_V() {
    // `potential_` is `HF::potential_`, whose `compute_V`
    // is not implemented and will throw error. No worries here,
    // because `HF::HF()` takes care to initialize it to `psi::RV`
    // or `psi::UV` according to the reference option.
    vector<SharedMatrix> D;
    vector<SharedMatrix> V;
    for (size_t is = 0; is < nspin_; ++is) {
        D.push_back(D_[is]);
        V.push_back(V_[is]);
    }
    potential_->set_D(D);
    potential_->compute_V(V);

    // V = V_dfa + V_losc.
    form_V_losc();
    for (size_t is = 0; is < nspin_; ++is) {
        V_[is]->add(V_losc_[is]);
    }
}

double LOSC::compute_E() {
    // General expression of electronic total energy:
    // E_elec  = \sum_s D_s * H + 1/2 * (J_a + J_b) * D_s - 1/2 * K_s * D_s
    // Summation over index a is for spin. Index a and b are for
    // alpha and beta spin respectively.

    double E_one_electron = 0.0;
    double E_coulomb = 0.0;
    double E_exchange = 0.0;
    const double alpha = functional_->x_alpha();
    const double factor = nspin_ == 1 ? 2.0 : 1.0;
    for (size_t is = 0; is < nspin_; ++is) {
        E_one_electron += factor * D_[is]->vector_dot(H_);
        E_coulomb += factor * D_[is]->vector_dot(J_);
        E_exchange -= factor * alpha * D_[is]->vector_dot(K_[is]);
    }

    double E_xc = 0.0;
    double E_vv10 = 0.0;
    if (functional_->needs_xc()) {
        E_xc = potential_->quadrature_values()["FUNCTIONAL"];
    }
    if (functional_->needs_vv10()) {
        E_vv10 = potential_->quadrature_values()["VV10"];
    }

    energies_["Nuclear"] = nuclearrep_;
    energies_["One-Electron"] = E_one_electron;
    energies_["Two-Electron"] = 0.5 * (E_coulomb + E_exchange);
    energies_["XC"] = E_xc;
    energies_["VV10"] = E_vv10;
    energies_["-D"] = scalar_variable("-D Energy");  // dispersion?
    double E_D = energies_["-D"];

    const double E_tot = nuclearrep_ + E_one_electron +
                         0.5 * (E_coulomb + E_exchange) + E_xc + E_D;
    return E_tot;
}

/**
 * @note
 * 1. Calling this function will release all the memory owned by the LOSC
 * object. This is okay. Because all the memember variables are
 * not exposed to users. So the user should not worry about the
 * availabilities of these internally used variables.
 * TODO:
 * Check if it really needs to manually release matrices here? They are actually
 * shared_ptr and will be released at the LOSC deconstructor.
 *
 * TODO:
 * 2. This function will set `Wavefunction::energy_` with
 * `HF::energies_["Total Energy"]`. However, it looks like
 * `HF::energies_["Total Energy"]` is never set in LOSC code.
 * Check where this quantity is set. This may give `Wavefunction::energy_`
 * with a unknown number.
 */
void LOSC::finalize() {
    // make sure all internally used variables are released.
    J_.reset();
    for (size_t is = 0; is < nspin_; ++is) {
        curvature_[is].reset();
        K_[is].reset();
        V_losc_[is].reset();
        G_[is].reset();
        D_old_[is].reset();
    }
    HF::finalize();
}

void LOSC::save_density_and_energy() {
    for (size_t is = 0; is < nspin_; ++is) {
        if (!D_old_[is]) {
            const string name = (is == 1 ? "D_old_a" : "D_old_b");
            D_old_[is] = SharedMatrix(factory_->create_matrix(name));
        }
        D_old_[is]->copy(D_[is]);
    }
}

void LOSC::damping_update(double factor) {
    // D = (1 - f) * D + f * D_old.
    for (size_t is = 0; is < nspin_; ++is) {
        D_[is]->scale(1.0 - factor);
        D_[is]->axpy(factor, D_old_[is]);
    }
}

bool LOSC::diis() {
    if (nspin_ == 1)
        return diis_manager_->extrapolate(1, F_[0].get());
    else if (nspin_ == 2)
        return diis_manager_->extrapolate(2, F_[0].get(), F_[1].get());
    else
        throw PSIEXCEPTION("LOSC diis: unknown spin.");
}

double LOSC::compute_orbital_gradient(bool save_fock, int max_diis_vectors) {
    // Conventional DIIS (X'[FDS - SDF]X, where X levels things out)
    SharedMatrix gradient[2];
    for (size_t is = 0; is < nspin_; ++is) {
        gradient[is] = form_FDSmSDF(F_[is], D_[is]);
    }

    if (save_fock) {
        // make sure the diis manager initialized.
        if (initialized_diis_manager_ == false) {
            if (scf_type_ == "DIRECT") {
                diis_manager_ = std::make_shared<DIISManager>(
                    max_diis_vectors, "LOSC DIIS vector",
                    DIISManager::LargestError, DIISManager::InCore);
            } else {
                diis_manager_ = std::make_shared<DIISManager>(
                    max_diis_vectors, "LOSC DIIS vector",
                    DIISManager::LargestError, DIISManager::OnDisk);
            }
            if (nspin_ == 1) {
                // RKS
                diis_manager_->set_error_vector_size(1, DIISEntry::Matrix,
                                                     gradient[0].get());
                diis_manager_->set_vector_size(1, DIISEntry::Matrix,
                                               F_[0].get());
            } else if (nspin_ == 2) {
                // UKS
                diis_manager_->set_error_vector_size(
                    2, DIISEntry::Matrix, gradient[0].get(), DIISEntry::Matrix,
                    gradient[1].get());
                diis_manager_->set_vector_size(2, DIISEntry::Matrix,
                                               F_[0].get(), DIISEntry::Matrix,
                                               F_[1].get());
            } else {
                throw PSIEXCEPTION(
                    "LOSC diis: unknown spin to calculate orbital gradients.");
            }
            initialized_diis_manager_ = true;
        }
        // save diis
        if (nspin_ == 1)
            diis_manager_->add_entry(2, gradient[0].get(), F_[0].get());
        else
            diis_manager_->add_entry(4, gradient[0].get(), F_[0].get(),
                                     gradient[1].get(), F_[1].get());
    }

    const double rms_v = gradient[0]->rms();
    const double max_v = gradient[0]->absmax();
    double rms[2] = {rms_v, rms_v};
    double max[2] = {max_v, max_v};
    if (nspin_ == 2) {
        rms[1] = gradient[1]->rms();
        max[1] = gradient[1]->absmax();
    }

    if (options_.get_bool("DIIS_RMS_ERROR")) {
        return std::sqrt(0.5 * (std::pow(rms[0], 2) + std::pow(rms[1], 2)));
    } else {
        return std::max(max[0], max[1]);
    }
}

int LOSC::soscf_update(double soscf_conv, int soscf_min_iter,
                       int soscf_max_iter, int soscf_print) {
    throw PSIEXCEPTION("Sorry, LOSC::soscf_update() is not implemented yet.");
}

bool LOSC::stability_analysis() {
    throw PSIEXCEPTION(
        "Sorry, LOSC::stability_analysis() is not implemented yet.");
}

vector<SharedMatrix> LOSC::onel_Hx(vector<SharedMatrix> x) {
    throw PSIEXCEPTION("Sorry, LOSC::onel_Hx() is not implemented yet.");
}

vector<SharedMatrix> LOSC::twoel_Hx(vector<SharedMatrix> x, bool combine,
                                    string return_basis) {
    throw PSIEXCEPTION("Sorry, LOSC::twoel_Hx() is not implemented yet.");
}

vector<SharedMatrix> LOSC::cphf_Hx(vector<SharedMatrix> x) {
    throw PSIEXCEPTION("Sorry, LOSC::cphf_Hx() is not implemented yet.");
}

vector<SharedMatrix> LOSC::cphf_solve(vector<SharedMatrix> x_vec,
                                      double conv_tol, int max_iter,
                                      int print_lvl) {
    throw PSIEXCEPTION("Sorry, LOSC::cphf_solve() is not implemented yet.");
}

SharedLOSC LOSC::c1_deep_copy(shared_ptr<BasisSet> basis) {
    throw PSIEXCEPTION("Sorry, LOSC::c1_deep_copy() is not implemented yet.");
}

}  // namespace losc
}  // namespace psi
