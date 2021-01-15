#include "curvature.h"
#include "losc.h"
#include "option_key.h"

#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/basisset.h"

namespace psi {
namespace losc {

// ==> CurvatureBase <== //
CurvatureBase::CurvatureBase(SharedHF& wfn, vector<SharedMatrix>& C_lo)
    : wfn_{wfn}, C_lo_{C_lo} {
    for (auto v : C_lo_) nlo_.push_back(v->ncol());
}

CurvatureBase::~CurvatureBase() {}

void CurvatureBase::prepare_df() {
    auto mints = MintsHelper(wfn_);
    auto zero_basis = BasisSet::zero_ao_basis_set();
    auto fitbasis = wfn_->get_basisset("DF_BASIS_SCF");
    auto ao_basis = wfn_->basisset();
    // <Q|mn>
    Qmn_ = mints.ao_eri(fitbasis, zero_basis, ao_basis, ao_basis);
    // <P|1/r|Q>^-1
    Vpq_ = mints.ao_eri(fitbasis, zero_basis, fitbasis, zero_basis);
    Vpq_->power(-1);

    const size_t nfit = fitbasis->nbf();
    const size_t nao = ao_basis->nbf();
    const size_t nao2 = nao * nao;

    // build <Q|ll>
    for (auto C : C_lo_) {
        const size_t nlo = C->ncol();
        auto Qll = std::make_shared<Matrix>(nfit, nlo);
        Matrix Qln(nlo, nao);
        for (size_t Q = 0; Q < nfit; ++Q) {
            // <Q|ln> = C_ml^T * <Q|mn>
            double* Qmn_ptr = &Qmn_->pointer()[0][Q * nao2];  // [nao x nao]
            double* Qln_ptr = Qln.pointer()[0];                // [nlo x nao]
            double* ml_ptr = C->pointer()[0];                 // [nao x nlo]
            C_DGEMM('T', 'N', nlo, nao, nao, 1.0, ml_ptr, nlo, Qmn_ptr, nao, 0.0,
                    Qln_ptr, nao);
            for (size_t l = 0; l < nlo; ++l) {
                // <Q|ll> = <Q|ln> dot C_nl
                double *lm = &Qln.pointer()[l][0];
                double *ml = &C->pointer()[0][l];
                const double v = C_DDOT(nao, lm, 1, ml, nlo);
                Qll->set(Q, l, v);
            }
        }
        Qll_.push_back(Qll);
    }
}

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

vector<SharedMatrix> CurvatureV1::compute_kappa_J() {
    prepare_df();
    vector<SharedMatrix> kappa_J;
    for (size_t is = 0; is < C_lo_.size(); ++is) {
        auto J = linalg::triplet(Qll_[is], Vpq_, Qll_[is], true, false, false);
        kappa_J.push_back(J);
    }
    return kappa_J;
}

vector<SharedMatrix> CurvatureV1::compute_kappa_xc() {
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
