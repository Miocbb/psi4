#ifndef __PSI4_LOSC_CURVATURE_H_
#define __PSI4_LOSC_CURVATURE_H_

#include "losc.h"
#include "psi4/libmints/typedefs.h"

namespace psi {
namespace losc {

class CurvatureBase {
   protected:
    SharedLOSC losc_wfn_;

   public:
    CurvatureBase(SharedLOSC losc_wfn);
    ~CurvatureBase();

    /**
     * Calculate the curvature matrix.
     *
     * @param spin[in]: the spin of the calculated curvature matrix.
     * 0 for alpha, 1 for beta.
     */
    virtual SharedMatrix compute(int spin) = 0;
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
    SharedMatrix compute_kappa_J(int spin);
    /**
     * Compute the exchange-correlation part of curvature for provided spin of
     * electron.
     */
    SharedMatrix compute_kappa_xc(int spin);

   public:
    /**
     * Constructor of curvature v1 object from a LOSC wavefunction.
     */
    CurvatureV1(SharedLOSC losc_wfn);
    ~CurvatureV1();
    /**
     * Compute the curvature v1 for the provided spin.
     */
    SharedMatrix compute(int spin) override;
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
    CurvatureV2(SharedLOSC losc_wfn);
    ~CurvatureV2();
    /**
     * Compute the curvature v2 for the provided spin.
     */
    SharedMatrix compute(int spin) override;
};

}  // namespace losc
}  // namespace psi

#endif
