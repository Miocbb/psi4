#ifndef __PSI4_LOSC_OPTION_KEY_H__
#define __PSI4_LOSC_OPTION_KEY_H__

namespace psi {
namespace losc {
// general options
static const char* option_print_level = "losc_print";

// localization options
static const char* option_localize_version = "losc_localizer_version";
static const char* option_localize_convergence = "losc_localizer_convergence";
static const char* option_localize_max_iteration = "losc_localizer_max_iter";
static const char* option_localize_random_permute =
    "losc_localize_random_permute";
static const char* option_localize_guess = "losc_localizer_guess";
static const char* option_localize_v2_c = "losc_localizer_v2_c";
static const char* option_localize_v2_gamma = "losc_localizer_v2_gamma";

// curvature options
static const char* option_curvature_version = "losc_curvature_version";
static const char* option_curvature_v1_cx = "losc_curvature_v1_cx";
static const char* option_curvature_v1_tau = "losc_curvature_v1_tau";
static const char* option_curvature_v2_zeta = "losc_curvature_v2_zeta";

}  // namespace losc
}  // namespace psi
#endif
