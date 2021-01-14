#ifndef __PSI4_LOSC_LOCALIZATION_H_
#define __PSI4_LOSC_LOCALIZATION_H_

#include "losc.h"
#include "psi4/libmints/typedefs.h"
#include <vector>
#include <memory>
#include <cstddef>
#include <string>

namespace psi {
namespace losc {
using std::vector;

class LocalizerBase {
   protected:
    SharedWavefunction wfn_;

    // The unitary transformation matrix for CO to LO.
    // LO_i = CO_j Uji.
    vector<SharedMatrix> U_;

    // The delocalized orbitals which are usually the LO basis.
    // [nso, nlo]
    vector<SharedMatrix> C_;

    // LO coefficient matrix.
    // [nso, nlo]
    vector<SharedMatrix> L_;

    // number of LOs.
    vector<int> nlo_;

    // 1 for RKS and 2 for UKS.
    size_t nspin_;

    size_t maxiter_;
    double convergence_;
    bool converged_;

    /**
     * Initialize the LO basis `C_` and number of LOs `nlo_` from the input
     * wavefunction object.
     */
    void init_lo_basis();

    /**
     * Set the initial guess of U matrix for a spin (either alpha or beta).
     */
    void init_guess(size_t spin);

   public:
    // Note: the size of all the returned vectors of data associated with spin
    // is either 1 for RKS or 2 for UKS.

    /**
     * Constructor of localization base.
     */
    LocalizerBase(SharedWavefunction wfn);
    ~LocalizerBase();

    /**
     * Calculate the localize orbitals (LOs).
     */
    virtual void localize() = 0;

    vector<int> nlo() { return nlo_; }
    /**
     * Get the U matrix.
     */
    vector<SharedMatrix> get_U() { return U_; }

    /**
     * Get the LO coefficient matrix.
     */
    vector<SharedMatrix> get_LO() { return L_; }
};

/**
 * LOSC localization version 2.
 *
 * @see The first paper of LOSC https://doi.org/10.1093/nsr/nwx111.
 */
class LocalizerV2 : public LocalizerBase {
   protected:
    /**
     * KS Halmiltonian matrix under AO for the DFA (NOT LOSC-DFA).
     */
    vector<SharedMatrix> H_ao_;

    /**
     * Dipole matrix under AO in x, y and z directions.
     */
    vector<SharedMatrix> D_ao_;

    /**
     * Parameter C in Losc2 localization cost function.
     */
    double para_c_;

    /**
     * Parameter gamma in Losc2 localization cost function.
     */
    double para_gamma_;

    bool js_random_permutation_;

    void print_header();

    void js_optimize_one_pair(const size_t i, const size_t j,
                              const vector<SharedMatrix> &D_lo,
                              const SharedMatrix &H_lo, double &theta_val,
                              double &delta_val);
    void js_rotate_one_pair(const size_t i, const size_t j, const double theta,
                            const SharedMatrix &U,
                            const vector<SharedMatrix> &D_lo,
                            const SharedMatrix &H_lo);
    void do_localization(int spin);

   public:
    /**
     * Constructor of localization version 2.
     */
    LocalizerV2(SharedWavefunction wfn);

    /**
     * Deconstructor of localization version 2.
     */
    ~LocalizerV2();

    /**
     * Calculate the localize orbitals (LOs) from localization version 2.
     */
    void localize() override;
};

}  // namespace losc
}  // namespace psi

#endif
