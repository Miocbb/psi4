#include "curvature.h"
#include "losc.h"
#include "option_key.h"

#include "psi4/libpsi4util/exception.h"

namespace psi {
namespace losc {

// ==> CurvatureBase <== //
CurvatureBase::CurvatureBase(SharedLOSC losc_wfn) : losc_wfn_{losc_wfn} {}

CurvatureBase::~CurvatureBase() {}

SharedMatrix CurvatureBase::compute(int spin) {
    throw PSIEXCEPTION(
        "Sorry, CurvatureBase::compute function is not implemented.");
}

// ==> CurvatureV1 <== //
CurvatureV1::CurvatureV1(SharedLOSC losc_wfn) : CurvatureBase(losc_wfn) {
    Options &options = losc_wfn_->options();
    para_cx_ = options.get_double(option_curvature_v1_cx);
    para_tau_ = options.get_double(option_curvature_v1_tau);
}

CurvatureV1::~CurvatureV1() {}

SharedMatrix CurvatureV1::compute_kappa_J(int spin) {
    throw PSIEXCEPTION("Curvature V1 kappa_J: working on it.");
}

SharedMatrix CurvatureV1::compute_kappa_xc(int spin) {
    throw PSIEXCEPTION("Curvature V1 kappa_xc: Working on it.");
}

SharedMatrix CurvatureV1::compute(int spin) {
    throw PSIEXCEPTION("Curvature V1: Working on it.");
}

// ==> CurvatureV2 <== //
CurvatureV2::CurvatureV2(SharedLOSC losc_wfn) : CurvatureBase(losc_wfn) {
    Options &options = losc_wfn_->options();
    para_cx_ = options.get_double(option_curvature_v1_cx);
    para_tau_ = options.get_double(option_curvature_v1_tau);
    para_zeta_ = options.get_double(option_curvature_v2_zeta);
}

CurvatureV2::~CurvatureV2() {}

SharedMatrix CurvatureV2::compute(int spin) {
    throw PSIEXCEPTION("Curvature V2: Working on it.");
}

}  // namespace losc
}  // namespace psi
