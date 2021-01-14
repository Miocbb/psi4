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
#include <vector>
#include <string>
namespace psi {
namespace losc {

class LOSC;
class LocalizerBase;
class CurvatureBase;

using std::shared_ptr;
using std::string;
using std::vector;
using Wfn = Wavefunction;
using SharedLOSC = shared_ptr<LOSC>;
using SharedHF = shared_ptr<psi::scf::HF>;

/**
 * Class of the localized orbital scaling correction (LOSC).
 *
 * The LOSC class is for the post-SCF (post-LOSC-DFA) or frozen LO SCF
 * (SCF-LOSC-DFA) correction for an associated DFA.
 * For SCF-LOSC-DFA, see "J. Phys. Chem. Lett. 2020, 11, 23, 10269–10277".
 *
 * @note
 * 1. The RKS and UKS LOSC-DFA calculation are combined into one decleration
 * to avoid code duplications.
 */
class LOSC : public psi::scf::HF {
   private:
    /**
     * A convenient class that manages the spin related data
     * (matrix, vector ...) in the base class.
     */
    template <class T>
    class SpinData {
       private:
        vector<T*> data_;

       public:
        SpinData(T& a, T& b) : data_{&a, &b} {};
        T& operator[](size_t spin) { return *(data_.at(spin)); }
    };

    using SpinMatrix = SpinData<SharedMatrix>;
    using SpinVector = SpinData<SharedVector>;

    // ==> convenient shared pointers to matrices in the base class <==
    SpinMatrix V_;    // LOSC-DFA total Vxc. Linked to `HF::Va_` and `HF::Vb_`.
    SpinMatrix F_;    // LOSC-DFA Fock. Linked to `Wfn::Fa_` and `Wfn::Fb_`.
    SpinMatrix C_;    // COs' coefficient. Linked to `Wfn::Ca_` and `Wfn::Cb_`.
    SpinMatrix D_;    // Density matrix. Linked to `Wfn::Da_` and `Wfn::Db_`.
    SpinVector eig_;  // COs' energy. Linked to `Wfn::epsilon_a_` and
                      // `Wfn::epsilon_b_`.

    // Construct local occupation and save it in `local_occ_`.
    void build_local_occupation();
    // Construct LOSC effective potential and save it in `V_losc_`.
    void build_V_losc();
    // Compute and return the LOSC energy correction.
    double compute_losc_energy();

   protected:
    SharedHF dfa_wfn_;  // the parent DFA wavefunction
    size_t nspin_;      // 1 for restricted, 2 for unrestricted.
    int nelec_[2];      // Total electron number.

    // Matrix used and managed in LOSC class.
    SharedMatrix J_;                  // Coulomb matrix.
    vector<SharedMatrix> K_;          // exact exchange matrix.
    vector<SharedMatrix> V_losc_;     // LOSC effective potential matrix.
    vector<SharedMatrix> G_;          // G matrix: G = J + K + V_dfa + V_losc.
    vector<SharedMatrix> D_old_;      // Save density matrix for damping.
    vector<SharedMatrix> curvature_;  // LOSC curvature matrix.
    vector<SharedMatrix> C_lo_;       // LO coefficients.
    vector<SharedMatrix> local_occ_;  // Local occupation matrix.
    vector<int> nlo_;                 // number of LOs.

    void common_init();

   public:
    /**
     * Constructor of LOSC based on the DFA wavefunction.
     *
     * @param dfa_wfn[in]: the DFA wavefunction.
     */
    LOSC(SharedHF dfa_wfn);

    ~LOSC() override;

    // ==> Functions for building matrices <==
    /**
     * Diagonalize LOSC-DFA Fock matrix to update CO coefficient matrix.
     */
    void form_C() override;

    /**
     * Calculate the density matrix based on current COs.
     */
    void form_D() override;

    /**
     * Calculate the LOSC-DFA Fock (Hamiltonian) matrix based on current COs.
     */
    void form_F() override;

    // ==> Functions for SCF iteration <==
    /**
     * Compute LOSC-DFA total energy based on current COs. E = E_dfa + E_losc.
     */
    double compute_E() override;

    /**
     * TODO: what is this for?
     */
    void finalize() override;

    /**
     * Save density matrix and energy of current SCF step.
     */
    void save_density_and_energy() override;

    /**
     * Update density matrix with damping.
     *
     * @param factor: damping factor.
     */
    void damping_update(double factor) override;

    /**
     * Run diis algorithm for SCF.
     */
    bool diis() override;

    /**
     * Compute the gradient of orbitals. Needed for diis.
     */
    double compute_orbital_gradient(bool save_fock,
                                    int max_diis_vectors) override;

    /**
     * TODO: add details later.
     */
    SharedLOSC c1_deep_copy(shared_ptr<BasisSet> basis);
};

}  // namespace losc
}  // namespace psi

#endif
