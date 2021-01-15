#ifndef __PSI4_LOSC_CURVATURE_H_
#define __PSI4_LOSC_CURVATURE_H_

#include "losc.h"
#include "psi4/libmints/typedefs.h"

namespace psi {
namespace losc {

/**
 * Curvature base class for LOSC. The curvature is always based on a KS
 * wavefunction and a set of LOs.
 */
class CurvatureBase {
   protected:
    SharedHF wfn_;
    vector<SharedMatrix> C_lo_;
    vector<int> nlo_;
    SharedMatrix Qmn_;  // <Q|mn> full matrix, no symmetry used.
    SharedMatrix Vpq_;  // <P|Q>^-1 full matrix, no symmetry used.
    vector<SharedMatrix> Qll_; // <Q|ll> with l being LOs.

   public:
    /**
     * Explicit constructor of curvature base class.
     *
     * @param wfn: a HF or derived wavefunction.
     * @param C_lo: the localized orbitals associated with `wfn`. Its dimension
     * is [nso, nlo].
     */
    CurvatureBase(SharedHF& wfn, vector<SharedMatrix>& C_lo);
    ~CurvatureBase();

    void prepare_df();

    /**
     * Calculate the curvature matrix.
     */
    virtual vector<SharedMatrix> compute() = 0;
};

/**
 * LOSC Curvature version 1.
 *
 * @see The first paper of LOSC https://doi.org/10.1093/nsr/nwx111.
 */
class CurvatureV1 : public CurvatureBase {
   private:
    /**
     * Paramerter \f$C_x\f$ in curvature.
     * See Eq. (10) in the LOSC paper
     */
    double para_cx_;
    /**
     * Paramerter \f$\tau\f$ in curvature.
     * See Eq. (10) in the LOSC paper
     */
    double para_tau_;

    /**
     * Compute the Coulomb part of curvature for provided spin of electron.
     */
    vector<SharedMatrix> compute_kappa_J();
    /**
     * Compute the exchange-correlation part of curvature for provided spin of
     * electron.
     */
    vector<SharedMatrix> compute_kappa_xc();

   public:
    /**
     * Constructor of curvature v1 object from a LOSC wavefunction.
     */
    CurvatureV1(SharedHF& wfn, vector<SharedMatrix>& C_lo);
    ~CurvatureV1();
    /**
     * Compute the curvature v1 for the provided spin.
     */
    vector<SharedMatrix> compute() override;
};

/**
 * LOSC Curvature version 2 used in LOSC2.
 *
 * @note Curvature version 2 is based on version 1.
 * @see The paper of LOSC2: https://doi.org/10.1093/nsr/nwx111.
 */
class CurvatureV2 : public CurvatureBase {
   private:
    /**
     * Paramerter \f$\zeta\f$ in curvature v2.
     */
    double para_zeta_;
    /**
     * Paramerter \f$C_x\f$ in curvature v1.
     */
    double para_cx_;
    /**
     * Paramerter \f$\tau\f$ in curvature v1.
     */
    double para_tau_;

   public:
    /**
     * Constructor of curvature v2 object from a LOSC wavefunction.
     */
    CurvatureV2(SharedHF& wfn, vector<SharedMatrix>& C_lo);
    ~CurvatureV2();
    /**
     * Compute the curvature v2 for the provided spin.
     */
    vector<SharedMatrix> compute() override;
};

}  // namespace losc
}  // namespace psi

#endif
