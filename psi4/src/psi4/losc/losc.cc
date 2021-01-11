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

#include <string>

namespace psi {
namespace losc {

using std::string;
using Wfn = Wavefunction;

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

LOSC::LOSC(SharedWavefunction dfa_wfn, std::shared_ptr<SuperFunctional> func)
    : psi::scf::HF(dfa_wfn, func, dfa_wfn->options(), PSIO::shared_object()) {
    common_init();
}

LOSC::LOSC(SharedWavefunction dfa_wfn, std::shared_ptr<SuperFunctional> func,
           Options& options, std::shared_ptr<PSIO> psio)
    : psi::scf::HF(dfa_wfn, func, options, PSIO::shared_object()) {}

LOSC::~LOSC() {}

/**
 * @note
 * Calling this function will allocate and compute following variables:
 * `C_`, `Wfn::Ca_`, `Wfn::Cb_`: the COs' coefficient matrix.
 * `eig_`: the orbital energies of COs.
 */
void LOSC::form_C() {
    for (size_t is = 0; is < nspin_; ++is) {
        const string C_name = (is == 1 ? "Ca" : "Cb");
        const string eig_name =
            (is == 1 ? "Alpha orbital energy" : "Beta orbital energy");
        C_[is] = SharedMatrix(factory_->create_matrix(C_name));
        eig_[is] = SharedVector(factory_->create_vector());
        eig_[is]->set_name(eig_name);

        diagonalize_F(F_[is], C_[is], eig_[is]);
    }

    // assign back to wavefunction.
    Wfn::Ca_ = C_[0];
    Wfn::Cb_ = (nspin_ == 1 ? C_[0] : C_[1]);

    find_occupation();
}
/**
 * @note
 * Calling this function will allocate and compute following variables:
 * `D_`, `Wfn::Da_`, `Wfn::Db_`: the density matrix.
 */
void LOSC::form_D() {
    // Fully occupied COs' coefficient matrix.
    // [nbasis x nocc].
    // TODO: take care of fractional occupation number.
    SharedMatrix C_occ[2];
    for (size_t is = 0; is < nspin_; ++is) {
        const size_t nocc = nelec_[is];
        // Copy the occupied block of the orbital matrix.
        C_occ[is] = std::make_shared<Matrix>("Occupied COs", nso_, nocc);
        for (int p = 0; p < nso_; ++p) {
            for (int i = 0; i < nocc; ++i) {
                C_occ[is]->set(p, i, C_[is]->get(p, i));
            }
        }
        // D = C_occ * C_occ.T
        const string name =
            (is == 1 ? "Density matrix: alpha" : "Density matrix: beta");
        D_[is] = SharedMatrix(factory_->create_matrix(name));
        D_[is]->gemm(false, true, 1.0, C_occ[is], C_occ[is], 0.0);
    }

    // assign back to wavefunction.
    Wfn::Da_ = D_[0];
    Wfn::Db_ = (nspin_ == 1 ? D_[0] : D_[1]);
}

/**
 * @note
 * Calling this function will allocate and compute following variables:
 * `F_`, `Wfn::Fa_`, `Wfn::Fb_`: the Fock matrix.
 */
void LOSC::form_F() {
    // F_a = H + G_a + V_ext.
    for (size_t is = 0; is < nspin_; ++is) {
        const string name =
            (is == 1 ? "Fock matrix: alpha" : "Fock matrix: beta");
        F_[is] = SharedMatrix(factory_->create_matrix(name));
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

    // assign back to the wavefunction
    Wfn::Fa_ = F_[0];
    Wfn::Fb_ = (nspin_ == 1 ? F_[0] : F_[1]);
}

/**
 * @note
 * Calling this function will allocate and compute the following matrices:
 * `J_`: Coulomb matrix.
 * `K_`: exchange matrix.
 * `V_losc_`: LOSC effective potential.
 * `V_`: LOSC-DFA effective potential.
 * `G_`: G matrix.
 */
void LOSC::form_G() {
    // Compute J and K.
    // Push the C matrix on
    std::vector<SharedMatrix>& C = jk_->C_left();
    C.clear();
    C.push_back(Ca_subset("SO", "OCC"));
    C.push_back(Cb_subset("SO", "OCC"));
    // Run the JK object
    jk_->compute();
    // Pull the J and K matrices off
    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();
    const std::vector<SharedMatrix>& wK = jk_->wK();

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
        // G_a = J_a + J_b - alpha * K_a.
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
void LOSC::form_V_losc() {}

/**
 * @note
 * Calling this function will allocate and compute the following matrices:
 * `V_losc`: LOSC effective potential matrix.
 * `V_` and `Wfn::V_`, the LOSC-DFA potential matrix.
 */
void LOSC::form_V() {
    HF::Va_ = SharedMatrix(factory_->create_matrix("Va"));
    HF::Vb_ = HF::Va_;
    V_[0] = Va_;
    if (nspin_ == 2) {
        HF::Vb_ = SharedMatrix(factory_->create_matrix("Vb"));
        V_[1] = Vb_;
    }

    // `potential_` is `HF::potential_`, whose `compute_V`
    // is not implemented and will throw error. No worries here,
    // because `HF::HF()` takes care to initialize it to `psi::RV`
    // or `psi::UV` according to the reference option.
    potential_->set_D({Da_, Db_});
    potential_->compute_V({Va_, Vb_});

    form_V_losc();

    // V = V_dfa + V_losc.
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
    for (size_t is = 0; is < nspin_; ++is) {
        E_one_electron += D_[is]->vector_dot(H_);
        E_coulomb += D_[is]->vector_dot(J_);
        E_exchange -= alpha * D_[is]->vector_dot(K_[is]);
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
    for (size_t is = 0; is < 2; ++is) {
        K_[is].reset();
        V_losc_[is].reset();
        G_[is].reset();
        V_[is].reset();
        F_[is].reset();
        C_[is].reset();
        D_[is].reset();
        eig_[is].reset();
    }
    HF::finalize();
}

void LOSC::damping_update(double) {
    throw PSIEXCEPTION("Sorry, LOSC::damping_update() is not implemented yet.");
}

SharedLOSC LOSC::c1_deep_copy(std::shared_ptr<BasisSet> basis) {
    throw PSIEXCEPTION("Sorry, LOSC::c1_deep_copy() is not implemented yet.");
}

}  // namespace losc
}  // namespace psi
