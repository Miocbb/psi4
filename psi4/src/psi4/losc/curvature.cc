#include "curvature.h"
#include "losc.h"
#include "option_key.h"

#include "psi4/libpsi4util/exception.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"

namespace psi {
namespace losc {

// ==> CurvatureBase <== //
CurvatureBase::CurvatureBase(SharedHF& wfn, vector<SharedMatrix>& C_lo)
    : wfn_{wfn}, C_lo_{C_lo} {
    for (auto v : C_lo_) nlo_.push_back(v->ncol());
}

CurvatureBase::~CurvatureBase() {}

vector<SharedMatrix> CurvatureBase::compute() {
    throw PSIEXCEPTION(
        "Sorry, CurvatureBase::compute function is not implemented.");
}

// ==> CurvatureV1 <== //
CurvatureV1::CurvatureV1(SharedHF& wfn, vector<SharedMatrix>& C_lo)
    : CurvatureBase(wfn, C_lo) {
    Options& options = wfn_->options();
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

vector<SharedMatrix> CurvatureV1::compute() {
    // TODO: Currently it is a zero matrix. Implement later.
    vector<SharedMatrix> curvature;
    for (size_t is = 0; is < nlo_.size(); ++is) {
        curvature.push_back(std::make_shared<Matrix>(nlo_[is], nlo_[is]));
    }
    return curvature;
}

// ==> CurvatureV2 <== //
CurvatureV2::CurvatureV2(SharedHF& wfn, vector<SharedMatrix>& C_lo)
    : CurvatureBase(wfn, C_lo) {
    Options& options = wfn_->options();
    para_cx_ = options.get_double(option_curvature_v1_cx);
    para_tau_ = options.get_double(option_curvature_v1_tau);
    para_zeta_ = options.get_double(option_curvature_v2_zeta);
}

CurvatureV2::~CurvatureV2() {}

vector<SharedMatrix> CurvatureV2::compute() {
    // TODO: implement later.
    CurvatureV1 cv1_helper = CurvatureV1(wfn_, C_lo_);
    auto cv1 = cv1_helper.compute();
    return cv1;
}

}  // namespace losc
}  // namespace psi
