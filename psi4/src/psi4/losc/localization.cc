#include "localization.h"
#include "losc.h"
#include "psi4/libmints/typedefs.h"

namespace psi {
namespace losc {

// ==> Localizer Base <== //
LocalizerBase::LocalizerBase(SharedLOSC losc_wfn) {}
LocalizerBase::~LocalizerBase() {}
SharedMatrix LocalizerBase::get_u(int spin) {}
SharedMatrix LocalizerBase::compute(int spin) {}

// ==> Localizer version 2 <== //
LocalizerV2::LocalizerV2(SharedLOSC losc_wfn) : LocalizerBase{losc_wfn} {}
LocalizerV2::~LocalizerV2() {}
SharedMatrix LocalizerV2::compute(int spin) {}

}  // namespace losc
}  // namespace psi
