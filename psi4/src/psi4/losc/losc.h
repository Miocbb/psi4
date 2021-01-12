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
using std::shared_ptr;
using std::string;
using std::vector;
using SharedLOSC = shared_ptr<LOSC>;
using Wfn = Wavefunction;

/**
 * Class of the localized orbital scaling correction (LOSC).
 *
 * @note
 * 1. The RKS and UKS LOSC-DFA calculation are combined into one decleration
 * to avoid code duplications.
 * 2. The public interface of LOSC should be the "same" as `psi::scf:RHF`
 * and `psi::scf::RHF`. This consistenty is to avoid modifications in
 * python side code.
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

   protected:
    int nelec_[2];  // Total electron number.
    size_t nspin_;  // 1 for restricted, 2 for unrestricted.

    // Matrix used and managed in LOSC class.
    SharedMatrix curvature_[2];  // LOSC curvature matrix.
    SharedMatrix J_;             // Coulomb matrix.
    SharedMatrix K_[2];          // exact exchange matrix.
    SharedMatrix V_losc_[2];     // LOSC effective potential matrix.
    SharedMatrix G_[2];          // G matrix: G = J + K + V_dfa + V_losc.
    SharedMatrix D_old_[2];      // The density matrix from last iteration step.
                                 // Used for update density matrix with damping.

    /**
     * Common initialization in LOSC constructors.
     */
    void common_init();

   public:
    /**
     * Constructor of LOSC based on the wavefunction from a DFA.
     *
     * @param dfa_wfn[in]: the wavefunction of a DFA.
     * @param functinal[in]: functional object for the DFA.
     */
    LOSC(SharedWavefunction dfa_wfn, shared_ptr<SuperFunctional> functional);

    /**
     * Constructor of LOSC based on the wavefunction from a DFA.
     *
     * @param dfa_wfn[in]: the wavefunction of a DFA.
     * @param functinal[in]: functional object for the DFA.
     * @param options[in]: input options in psi4.
     * @param psio[in]: what is this?
     */
    LOSC(SharedWavefunction dfa_wfn, shared_ptr<SuperFunctional> functional,
         Options& options, shared_ptr<PSIO> psio);

    ~LOSC() override;

    // ==> Functions for building matrices <==
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

    // ==> Functions for SCF iteration <==
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
     * TODO: What is this?
     */
    int soscf_update(double soscf_conv, int soscf_min_iter, int soscf_max_iter,
                     int soscf_print) override;

    /**
     * TODO: What is this?
     */
    bool stability_analysis() override;

    // ==> Hessian-vector computers and solvers <==
    // TODO: Not implemented yet!
    vector<SharedMatrix> onel_Hx(vector<SharedMatrix> x) override;
    vector<SharedMatrix> twoel_Hx(vector<SharedMatrix> x, bool combine = true,
                                  string return_basis = "MO") override;
    vector<SharedMatrix> cphf_Hx(vector<SharedMatrix> x) override;
    vector<SharedMatrix> cphf_solve(vector<SharedMatrix> x_vec,
                                    double conv_tol = 1.e-4, int max_iter = 10,
                                    int print_lvl = 1) override;

    /**
     * TODO: add details later.
     */
    SharedLOSC c1_deep_copy(shared_ptr<BasisSet> basis);
};

}  // namespace losc
}  // namespace psi

#endif
