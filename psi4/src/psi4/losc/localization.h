#ifndef __PSI4_LOSC_LOCALIZATION_H_
#define __PSI4_LOSC_LOCALIZATION_H_

#include <cstddef>
#include "losc.h"
#include "psi4/libmints/typedefs.h"

namespace psi {
namespace losc {
class LocalizerBase {
   protected:
    SharedLOSC losc_wfn_;
    /**
     * Maximum iteration number of Jacobi-sweep algorithm for localization.
     */
    size_t js_max_iter_;
    /**
     * If use random permutation in Jacobi-sweep algorithm or not.
     */
    bool js_random_permutation_;
    /**
     * Convergence tolerance tolerance in Jacobi-sweep algorithm.
     */
    double js_convergence_;
    /**
     * @brief The unitary matrix that transfer LO basis coefficient matrix into
     * LO coefficient matrix.
     * @details Dimension: [nlo, nlo].
     * The relation between the LO \f$\psi\f$ and LO basis \f$\phi\f$ via the U
     * matrix is the following,
     * \f[
     * \psi_i = \sum_j U_{ij} \phi_j.
     * \f]
     */
    SharedMatrix U_[2];
    /**
     * @brief LO basis coefficient matrix under AO.
     * @details
     * Dimension: [nbasis, nlo]
     *
     * Relation between LO basis and AO via the coefficient matrix \f$C\f$ is
     * \f[
     * \psi_i = \sum_\mu C_{\mu i} \phi_{\mu},
     * \f]
     * where \f$C\f$ matrix is stored in `C_lo_basis_`.
     *
     * @note Usually, the COs are used as the LOs' basis.
     */
    SharedMatrix C_lo_basis_[2];
    /**
     * Set up the initial guess of the U matrix from a given U matrix.
     *
     * @note
     * The unitarity of the input U matrix is not verified.
     */
    SharedMatrix set_guess_u(int spin, SharedMatrix U_guess);
    /**
     * Set up the initial guess of the U matrix as a random unitary matrix.
     */
    SharedMatrix set_guess_u_random(int spin);
    /**
     * Set up the initial guess of the U matrix as an identity matrix.
     */
    SharedMatrix set_guess_u_idenity(int spin);

   public:
    /**
     * Constructor of localization base.
     */
    LocalizerBase(SharedLOSC losc_wfn);
    ~LocalizerBase();

    /**
     * Get the U matrix for the given spin.
     */
    SharedMatrix get_u(int spin);
    /**
     * Calculate the localize orbitals (LOs).
     *
     * @param spin[in]: the spin of the calculated curvature matrix.
     * 0 for alpha, 1 for beta.
     * @return the LOs' coefficients on atomic basis set.
     */
    virtual SharedMatrix compute(int spin) = 0;
};

/**
 * LOSC localization version 2.
 *
 * @see The first paper of LOSC https://doi.org/10.1093/nsr/nwx111.
 */
class LocalizerV2 : public LocalizerBase {
   protected:
    /**
     * KS Halmiltonian matrix of the DFA (NOT LOSC-DFA) on AO.
     */
    SharedMatrix H_ao_;
    /**
     * Dipole matrix under AO in order of x, y and z directions.
     */
    SharedMatrix Dipole_ao_[3];
    /**
     * Parameter \f$C\f$ in Losc2 localization cost function.
     */
    double para_c_;
    /**
     * Parameter \f$\gamma\f$ in Losc2 localization cost function.
     */
    double para_gamma_;

   public:
    /**
     * Constructor of localization version 2.
     */
    LocalizerV2(SharedLOSC losc_wfn);
    /**
     * Deconstructor of localization version 2.
     */
    ~LocalizerV2();
    /**
     * Calculate the localize orbitals (LOs) from localization version 2.
     *
     * @param spin[in]: the spin of the calculated curvature matrix.
     * 0 for alpha, 1 for beta.
     * @return the LOs' coefficients on atomic basis set.
     */
    SharedMatrix compute(int spin) override;
};

}  // namespace losc
}  // namespace psi

#endif
