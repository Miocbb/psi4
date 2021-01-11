/**
 * @brief
 * Decleration of LOSC class.
 */

#ifndef __PSI4_LOSC_LOSC_H__
#define __PSI4_LOSC_LOSC_H__

#include "psi4/libpsio/psio.hpp"
#include "psi4/libscf_solver/hf.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libmints/basisset.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.hpp"

#include <memory>

namespace psi {
namespace losc {

class LOSC;
using SharedLOSC = std::shared_ptr<LOSC>;

/**
 * Class of the localized orbital scaling correction (LOSC).
 *
 * @note
 * 1. The RKS and UKS LOSC-DFA calculation are combined into one decleration
 * to avoid code duplications.
 * 2. There should be at least a function named `c1_deep_copy()`
 * in `class psi::losc::LOSC`. This is the only one function exposed
 * to the py-side code, like in `psi::scf::UHF` and `psi::scf::RHF`.
 */
class LOSC : public psi::scf::HF {
   protected:
    int nelec_[2];  // Total electron number.
    size_t nspin_;  // 1 for restricted, 2 for unrestricted.

    // LOSC matrix pointers.
    SharedMatrix curvature_[2];  // LOSC curvature matrix.
    SharedMatrix J_;             // Coulomb matrix.
    SharedMatrix K_[2];          // exact exchange matrix.
    SharedMatrix V_losc_[2];     // LOSC effective potential matrix.
    SharedMatrix G_[2];          // G matrix: G = J + K + V_dfa + V_losc.

    // ==> convenient shared pointers to matrices <==
    // convenient shared pointers should be assigned back
    // the base class shared pointers afther the construction
    // of the corresponding matrices.

    // convenient pointers from psi::scf::HF.
    SharedMatrix V_[2];  // LOSC-DFA potential matrix: V = V_dfa + V_losc.

    // convenient pointers from psi::Wavefunction.
    SharedMatrix F_[2];    // LOSC-DFA Fock matrix.
    SharedMatrix C_[2];    // CO coefficient matrix.
    SharedMatrix D_[2];    // density matrix.
    SharedVector eig_[2];  // CO eigenvalues.

    void common_init();

   public:
    /**
     * Constructor of LOSC based on the wavefunction from a DFA.
     *
     * @param dfa_wfn[in]: the wavefunction of a DFA.
     * @param functinal[in]: functional object for the DFA.
     */
    LOSC(SharedWavefunction dfa_wfn,
         std::shared_ptr<SuperFunctional> functional);

    /**
     * Constructor of LOSC based on the wavefunction from a DFA.
     *
     * @param dfa_wfn[in]: the wavefunction of a DFA.
     * @param functinal[in]: functional object for the DFA.
     * @param options[in]: input options in psi4.
     * @param psio[in]: what is this?
     */
    LOSC(SharedWavefunction dfa_wfn,
         std::shared_ptr<SuperFunctional> functional, Options &options,
         std::shared_ptr<PSIO> psio);

    ~LOSC() override;

    /**
     * Calculate the LOSC effective potential matrix.
     */
    virtual void form_V_losc();

    /**
     * calculate the LOSC-DFA COs' coefficient matrix.
     */
    void form_C() override;

    /**
     * calculate the LOSC-DFA density matrix.
     */
    void form_D() override;

    /**
     * calculate the LOSC-DFA Fock (Hamiltonian) matrix.
     */
    void form_F() override;

    /**
     * calculate the LOSC-DFA G matrix.
     *
     * G matrix involves the effective contribution from LOSC, i.e.,
     * G = J + K + Vxc + V_losc.
     */
    void form_G() override;

    /**
     * Calculate the LOSC-DFA potential matrix.
     *
     * V = Vxc + V_losc
     */
    void form_V() override;

    /**
     * Compute LOSC-DFA total energy.
     *
     * E = E_dfa + E_losc
     */
    double compute_E() override;

    /**
     * TODO: what is this for?
     */
    void finalize() override;

    /**
     * TODO: damping for what? diis related?
     */
    void damping_update(double) override;

    /**
     * TODO: add details later.
     */
    SharedLOSC c1_deep_copy(std::shared_ptr<BasisSet> basis);
};

}  // namespace losc
}  // namespace psi

#endif
