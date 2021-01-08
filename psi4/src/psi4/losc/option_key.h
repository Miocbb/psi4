
#ifndef __PSI4_LOSC_OPTION_KEY_H__
#define __PSI4_LOSC_OPTION_KEY_H__

namespace psi {
namespace losc {
// general options
const char *option_print_level = "losc_print";

// localization options
const char *option_localize_convergence = "losc_localize_convergence";
const char *option_localize_max_iteration = "losc_localize_max_iter";
const char *option_localize_random_permute = "losc_localize_random_permute";
const char *option_localize_guess = "losc_localize_guess";
const char *option_localize_v2_c = "losc_localize_v2_c";
const char *option_localize_v2_gamma = "losc_localize_v2_gamma";

// curvature options
const char *option_curvature_version = "losc_curvature_version";
const char *option_curvature_v1_cx = "losc_curvature_v1_cx";
const char *option_curvature_v1_tau = "losc_curvature_v1_tau";
const char *option_curvature_v2_zeta = "losc_curvature_v2_zeta";

} // namespace losc
} // namespace psi
#endif
