#include "losc.h"
#include "option_key.h"
#include "localization.h"
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
#include "psi4/libpsi4util/PsiOutStream.h"

#include <string>
#include <vector>
#include <memory>

namespace psi {
namespace losc {
void LOSC::common_init() {
    if (!functional_->needs_xc())
        throw PSIEXCEPTION(
            "You can only apply LOSC to a DFT functional, not HF!");

    if (functional_->is_x_lrc())
        throw PSIEXCEPTION(
            "LOSC code does not support long range corrected functionals yet!");

    // current implementation is only for c1 symmetry.
    if (molecule()->schoenflies_symbol() != "c1")
        throw PSIEXCEPTION("LOSC only supports C1 symmetry for now!");

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

    // initialize the size of vectors of data based on spin.
    K_.assign(nspin_, nullptr);
    V_losc_.assign(nspin_, nullptr);
    G_.assign(nspin_, nullptr);
    D_old_.assign(nspin_, nullptr);
    curvature_.assign(nspin_, nullptr);
    C_lo_.assign(nspin_, nullptr);
    local_occ_.assign(nspin_, nullptr);
}

/**
 * @note
 * The constructor of the base HF class takes a shallow copy of the input
 * wavefunction. We don't want the copy. We just need `dfa_wfn` as a reference
 * wavefunction that can be used for localization to obtain LOs.
 */
LOSC::LOSC(SharedHF dfa_wfn)
    : psi::scf::HF(
          std::make_shared<Wfn>(dfa_wfn->molecule(), dfa_wfn->basisset()),
          dfa_wfn->functional(), dfa_wfn->options(), PSIO::shared_object()),
      V_{HF::Va_, HF::Vb_},
      F_{Wfn::Fa_, Wfn::Fb_},
      C_{Wfn::Ca_, Wfn::Cb_},
      D_{Wfn::Da_, Wfn::Db_},
      eig_{Wfn::epsilon_a_, Wfn::epsilon_b_} {
    set_reference_wavefunction(dfa_wfn);
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
    for (size_t is = 0; is < nspin_; ++is) {
        const size_t nocc = nelec_[is];
        // Copy the occupied block of the orbital matrix.
        auto C_occ = std::make_shared<Matrix>("Occupied COs", nso_, nocc);
        for (int p = 0; p < nso_; ++p) {
            for (int i = 0; i < nocc; ++i) {
                C_occ->set(p, i, C_[is]->get(p, i));
            }
        }
        // D = C_occ * C_occ.T
        D_[is]->gemm(false, true, 1.0, C_occ, C_occ, 0.0);
    }
}

/**
 * @note
 * 1. In principle, this function calculates the total Fock matrix F = F_dfa +
 * F_losc for the LOSC-DFA. The construction of F_losc is based on F_dfa. Thus,
 * we need to build F_dfa first.
 * 2. During the construction of F, all the intermediate matrices, like G,
 * J, K, V_dfa and V_losc, are computed as well. There is no need to explicitly
 * expose functions to compute these intermediate matrices like what they did
 * in `scf::RHF` and `scf::UHF`. Here, we only override the `HF::form_F`
 * function.
 * 3. After calling this function, the total Fock matrices are updated in
 * `Wfn::Fa_` and `Wfn::Fb_`, and the total LOSC-DFA potential matrices V are
 * updated in `HF::Va_` and `HF::Vb_`.
 */
void LOSC::form_F() {
    // F_a = H + G_a + V_ext
    // G_a = J_a + J_b - alpha * K_a + V_dfa + V_losc

    // ==> build DFA Fock (without LOSC contribution) <==
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

    // Compute DFA potential V_dfa.
    if (functional_->needs_xc()) {
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
        // Note: `potential_` is just DFA potential, not LOSC-DFA potential.
        potential_->set_D(D);
        potential_->compute_V(V);
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
    }

    // Compute DFA Fock.
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

    // ==> Calculate LOSC-DFA Fock <==
    // Compute LOSC effective potential.
    build_V_losc();

    // Compute LOSC-DFA Fock and save LOSC-DFA V and G matrices.
    for (size_t is = 0; is < nspin_; ++is) {
        F_[is]->add(V_losc_[is]);
        G_[is]->add(V_losc_[is]);
        V_[is]->add(V_losc_[is]);
    }
}

/**
 * @note
 * 1. This is for frozen LO LOSC. Only build LO once.
 */
void LOSC::build_lo() {
    static size_t eval_num = 0;
    // Check if the LOs are constructed or not.
    bool do_lo = false;
    if (eval_num = 0 || C_lo_.size() == 0) do_lo = true;
    for (auto& lo : C_lo_) {
        if (!lo) do_lo = true;
    }
    if (!do_lo) return;

    int version = options_.get_int(option_localize_version);
    std::shared_ptr<LocalizerBase> localizer = nullptr;
    if (version == 2) {
        localizer = std::make_shared<LocalizerV2>(reference_wavefunction_);
    } else {
        throw PSIEXCEPTION(
            "Sorry, currently only support localization version 2.");
    }

    localizer->localize();
    C_lo_ = localizer->get_LO();
    eval_num = 1;
}

void LOSC::build_local_occupation() {
    local_occ_.clear();
    const size_t nso = basisset_->nbf();
    for (size_t is = 0; is < nspin_; ++is) {
        if (!C_lo_[is])
            throw PSIEXCEPTION(
                "LOSC implementation error: construct LOs first.");
        const size_t nlo = C_lo_[is]->ncol();
        auto SC = std::make_shared<Matrix>(nso, nlo);
        SC->gemm(false, false, 1.0, S_, C_lo_[is], 0.0);
        local_occ_.push_back(
            linalg::triplet(SC, D_[is], SC, true, false, false));
    }
}

void LOSC::build_V_losc() {
    build_lo();
    build_local_occupation();

    const size_t nso = basisset_->nbf();
    for (size_t is = 0; is < nspin_; ++is) {
        const size_t nlo = C_lo_[is]->ncol();
        // build A matrix.
        auto A = std::make_shared<Matrix>(nlo, nlo);
        for (size_t i = 0; i < nlo; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                // TODO: Add curvature contribution to V_losc later.
                const double K_ij = 0;
                const double L_ij = local_occ_[is]->get(i, j);
                if (i != j) {
                    A->set(i, j, -K_ij * L_ij);
                } else {
                    A->set(i, j, 0.5 * K_ij - K_ij * L_ij);
                    A->set(j, i, 0.5 * K_ij - K_ij * L_ij);
                }
            }
        }

        // calculate Losc correcting Hamiltonian matrix.
        auto CS = std::make_shared<Matrix>(nlo, nso);
        CS->gemm(true, false, 1.0, C_lo_[is], S_, 0.0);
        V_losc_[is] = linalg::triplet(CS, A, CS, true, false, false);
    }
}

// TODO: implement later.
double LOSC::compute_losc_energy() { return 0; }

/**
 * @note
 * We did not check if the internal matrices used to calculate the energy
 * is constructed or not. Just the evaluation.
 */
double LOSC::compute_E() {
    // General expression of electronic total energy for a DFA:
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

    // LOSC energy correction.
    const double E_corr = compute_losc_energy();

    energies_["Nuclear"] = nuclearrep_;
    energies_["One-Electron"] = E_one_electron;
    energies_["Two-Electron"] = 0.5 * (E_coulomb + E_exchange);
    energies_["XC"] = E_xc;
    energies_["VV10"] = E_vv10;
    energies_["-D"] = scalar_variable("-D Energy");
    energies_["LOSC"] = E_corr;
    double E_D = energies_["-D"];

    const double E_tot = nuclearrep_ + E_one_electron +
                         0.5 * (E_coulomb + E_exchange) + E_xc + E_D + E_corr;
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

SharedLOSC LOSC::c1_deep_copy(shared_ptr<BasisSet> basis) {
    throw PSIEXCEPTION("Sorry, LOSC::c1_deep_copy() is not implemented yet.");
}

}  // namespace losc
}  // namespace psi
