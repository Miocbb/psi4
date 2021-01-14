#include "losc.h"
#include "localization.h"
#include "option_key.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libqt/qt.h"
#include <memory>
#include <algorithm>
#include <random>

namespace psi {
namespace losc {

// ==> LOSC Localizer Base <== //
/**
 * @note
 * 1. COs in the wavefunction object is ALWAYS used as the LO basis.
 * 2. The window may be applied to the localization. Thus not all the COs
 * are used for localization.
 * 3. This function will allocate memory for `C_`, only when the localization is
 * done with window.
 */
void LocalizerBase::init_lo_basis() {
    // TODO: currently just use all COs and does not support LO window.
    C_ = {wfn_->Ca(), wfn_->Cb()};
    nlo_ = {C_[0]->ncol(), C_[1]->ncol()};
}

/**
 * @note
 * 1. This function allocate memory for `U_[spin]` matrix.
 * 2. The guess is set based on the options in the wavefunction `wfn_`.
 * 3. If the guess is for beta spin, try to use the U matrix from alpha
 * chanel first. This assumes alpha localization is done first and
 * alpha and beta LOs are similar to each other.
 */
void LocalizerBase::init_guess(size_t spin) {
    if (spin >= nspin_)
        throw PSIEXCEPTION("losc::LocalizerBase::guess(): invalid spin.");

    // make sure to allocate memory.
    if (!U_.at(spin))
        U_.at(spin) = std::make_shared<Matrix>(nlo_[spin], nlo_[spin]);
    // Beta localization for UKS. Try to use alpha localization U.
    if (spin == 1 && nlo_[0] == nlo_[1]) {
        U_.at(1) = U_.at(0)->clone();
        return;
    }
    // use guess setting from options.
    Options &opt = wfn_->options();
    string guess_str = to_lower_copy(opt.get_str(option_localize_guess));
    if (guess_str == "random")
        throw PSIEXCEPTION(
            "LOSC localization: sorry, does not support random initial guess.");
    else if (guess_str == "identity")
        U_.at(spin)->identity();
    else
        throw PSIEXCEPTION("LOSC localization: invalid localization guess.");
}

LocalizerBase::LocalizerBase(SharedHF &wfn) : wfn_{wfn} {
    // current implementation is only for c1 symmetry.
    if (wfn_->molecule()->schoenflies_symbol() != "c1")
        throw PSIEXCEPTION(
            "LOSC localization only supports C1 symmetry for now!");

    // RKS or UKS?
    nspin_ = (wfn_->same_a_b_orbs() && wfn_->same_a_b_dens() ? 1 : 2);

    // initialize the size of vectors of data based on spin.
    U_.assign(nspin_, nullptr);
    C_.assign(nspin_, nullptr);
    L_.assign(nspin_, nullptr);
    nlo_.assign(nspin_, 0);

    // initialize LO basis.
    init_lo_basis();
}

LocalizerBase::~LocalizerBase() {}

// ==> Localizer version 2 <== //

/**
 * @note
 * 1. Dipole matrix under AO will be computed during the construction.
 */
LocalizerV2::LocalizerV2(SharedHF &wfn) : LocalizerBase{wfn} {
    // initialize DFA Hamiltonian and dipole matrix under AO.
    H_ao_ = {wfn_->Fa(), wfn_->Fb()};
    D_ao_ = wfn_->mintshelper()->so_dipole();
    for (auto h : H_ao_) {
        if (!h)
            throw PSIEXCEPTION(
                "LOSC localization 2: wavefunction has empty Fock matrix.");
    }

    // initialize localization parameters.
    auto options = wfn_->options();
    para_c_ = options.get_double(option_localize_v2_c);
    para_gamma_ = options.get_double(option_localize_v2_gamma);
    maxiter_ = options.get_int(option_localize_max_iteration);
    convergence_ = options.get_int(option_localize_convergence);
    converged_ = false;
    js_random_permutation_ = true;
}

LocalizerV2::~LocalizerV2() {}

void LocalizerV2::print_header() {
    outfile->Printf("  ==> LOSC Localizer Version 2 <==\n\n");
    outfile->Printf("    Convergence = %11.3E\n", convergence_);
    outfile->Printf("    Maxiter     = %11d\n", maxiter_);
    outfile->Printf("\n");
}

/**
 * @brief Get the optimized rotation angle value for one pair of orbitals.
 */
void LocalizerV2::js_optimize_one_pair(const size_t i, const size_t j,
                                       const vector<SharedMatrix> &D_lo,
                                       const SharedMatrix &H_lo,
                                       double &theta_val, double &delta_val) {
    // Spacial part: constant x1 and x2.
    double x1 = 0.0;
    double x2 = 0.0;
    double Dii_Djj = 0.0;
    double Dij = 0.0;
    for (size_t xyz = 0; xyz < 3; ++xyz) {
        SharedMatrix DMat = D_lo[xyz];
        Dii_Djj = DMat->get(i, i) - DMat->get(j, j);
        Dij = DMat->get(i, j);
        x1 += -0.5 * Dii_Djj * Dii_Djj + 2 * Dij * Dij;
        x2 -= 2 * Dij * Dii_Djj;
    }

    // Energy part: constant y1 and y2.
    double y1 = 0.0;
    double y2 = 0.0;
    double Eii_Ejj = 0.0;
    double Eij = 0.0;
    SharedMatrix EMat = H_lo;
    Eii_Ejj = EMat->get(i, i) - EMat->get(j, j);
    Eij = EMat->get(i, j);
    y1 = -0.5 * Eii_Ejj * Eii_Ejj + 2 * Eij * Eij;
    y2 = -2 * Eij * Eii_Ejj;

    /* constant values */
    const double a1 = (1 - para_gamma_) * x1 + para_gamma_ * para_c_ * y1;
    const double a2 = (1 - para_gamma_) * x2 + para_gamma_ * para_c_ * y2;

    /* theta */
    double theta = 0.0;
    // atan2 returns a value in range [-pi, pi], which is of a
    // larger range than what atan gives ([-pi/2, pi/2]).
    // This is actually a desired behavior, as we may want the
    // derivative of U to be stable.
    // Then U should not jump from -pi/4 to pi/4 (or pi/4 to -pi/4).
    theta = 0.5 * atan2(-a1 - sqrt(a1 * a1 + a2 * a2), a2);
    // try not to get too large rotation angles.
    // however we don't want restrict it in [-pi/4, pi/4],
    // as that will cause derivative discontinuity for numerical dU/dP.
    if (theta > 3. * M_PI / 8.) theta -= M_PI / 2.;
    if (theta < -3. * M_PI / 8.) theta += M_PI / 2.;
    theta_val = theta;
    double delta = a1 * cos(4 * theta) + a2 * sin(4 * theta) - a1;
    delta_val = delta;  // unit in bohr^2.
}

/**
 * @brief Apply rotation to the pair of orbitals.
 */
void LocalizerV2::js_rotate_one_pair(const size_t i, const size_t j,
                                     const double theta, const SharedMatrix &U,
                                     const vector<SharedMatrix> &D_lo,
                                     const SharedMatrix &H_lo) {
    double costheta = cos(theta);
    double sintheta = sin(theta);
    const int nlo = U->ncol();

    auto U_Mat = U->pointer();
    C_DROT(nlo, &U_Mat[i][0], 1, &U_Mat[j][0], 1, costheta, sintheta);

    // rotate Dipole Matrix
    for (int xyz = 0; xyz < 3; ++xyz) {
        auto DMat = D_lo[xyz]->pointer();
        C_DROT(nlo, &DMat[i][0], 1, &DMat[j][0], 1, costheta, sintheta);
        C_DROT(nlo, &DMat[0][i], nlo, &DMat[0][j], nlo, costheta, sintheta);
    }

    // rotate localization Hamiltonian matrix.
    auto *HMat = H_lo->pointer();
    C_DROT(nlo, &HMat[i][0], 1, &HMat[j][0], 1, costheta, sintheta);
    C_DROT(nlo, &HMat[0][i], nlo, &HMat[0][j], nlo, costheta, sintheta);
}

void LocalizerV2::do_localization(int spin) {
    const size_t nlo = nlo_[spin];
    const size_t nso = wfn_->basisset()->nbf();
    SharedMatrix C = C_[spin];
    SharedMatrix U = U_[spin];
    SharedMatrix H_ao = H_ao_[spin];

    if (nlo < 1) return;

    SharedMatrix CU = std::make_shared<Matrix>(nso, nlo);
    CU->gemm(false, false, 1.0, C, U, 0.0);

    // calculate dipole on LO initial guess.
    // D_lo = U^T * C_lo_basis^T * D_ao * C_lo_basis * U
    vector<SharedMatrix> D_lo;
    for (size_t xyz = 0; xyz < 3; xyz++) {
        auto D_lo_init =
            linalg::triplet(CU, D_ao_[xyz], CU, true, false, false);
        D_lo.push_back(D_lo_init);
    }

    // calculate Hamiltonian matrix on LO initial guess.
    // H_lo = U^T * C_lo_basis^T * H_ao * C_lo_basis * U
    SharedMatrix H_lo = linalg::triplet(CU, H_ao, CU, true, false, false);

    // free CU, which is not needed any more.
    CU.reset();

    // Main loop of localization
    std::mt19937 g(0);
    size_t iter = 0;
    double cycle_delta = 10000.0;
    while (iter < maxiter_ && std::abs(cycle_delta) > convergence_) {
        cycle_delta = 0.0;
        vector<size_t> order;
        for (size_t i = 0; i < nlo; ++i) order.push_back(i);

        // random permutation
        if (js_random_permutation_) {
            std::shuffle(order.begin(), order.end(), g);
        }

        // Jacobi Sweep 2x2 rotation
        for (size_t ni = 0; ni < nlo; ++ni) {
            size_t i = order[ni];
            for (size_t nj = 0; nj < ni; ++nj) {
                size_t j = order[nj];
                double delta = 0.0;
                double theta = 0.0;
                // get the optimized rotation angle.
                js_optimize_one_pair(i, j, D_lo, H_lo, theta, delta);
                // apply rotation to related matrix.
                js_rotate_one_pair(i, j, theta, U, D_lo, H_lo);
                cycle_delta += delta;
            }
        }
        ++iter;
    }

    if (iter < maxiter_ || abs(cycle_delta) < convergence_) converged_ = true;

    // calculate the LO coefficient matrix.
    // L = C * U
    auto L = std::make_shared<Matrix>(nso, nlo);
    L->gemm(false, false, 1.0, C, U, 0.0);
    L_.push_back(L);
}

void LocalizerV2::localize() {
    L_.clear();
    for (size_t is = 0; is < nspin_; ++is) {
        init_guess(is);
        do_localization(is);
    }
}

}  // namespace losc
}  // namespace psi
