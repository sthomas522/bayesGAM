// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;


RcppExport SEXP _rcpp_module_boot_stan_fit4glm_continuous_mixed_with_qr_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4glm_continuous_with_qr_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4glm_discrete_mixed_with_qr_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4glm_discrete_with_qr_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4multresponse_semipar_array_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4multresponse_semipar_array_discrete_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4multresponse_semipar_array_mixed_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4multresponse_semipar_array_mixed_discrete_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4multresponse_semipar_array_mixed_randomint_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4multresponse_semipar_array_mixed_randomint_discrete_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_stan_fit4glm_continuous_mixed_with_qr_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4glm_continuous_mixed_with_qr_mod, 0},
    {"_rcpp_module_boot_stan_fit4glm_continuous_with_qr_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4glm_continuous_with_qr_mod, 0},
    {"_rcpp_module_boot_stan_fit4glm_discrete_mixed_with_qr_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4glm_discrete_mixed_with_qr_mod, 0},
    {"_rcpp_module_boot_stan_fit4glm_discrete_with_qr_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4glm_discrete_with_qr_mod, 0},
    {"_rcpp_module_boot_stan_fit4multresponse_semipar_array_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4multresponse_semipar_array_mod, 0},
    {"_rcpp_module_boot_stan_fit4multresponse_semipar_array_discrete_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4multresponse_semipar_array_discrete_mod, 0},
    {"_rcpp_module_boot_stan_fit4multresponse_semipar_array_mixed_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4multresponse_semipar_array_mixed_mod, 0},
    {"_rcpp_module_boot_stan_fit4multresponse_semipar_array_mixed_discrete_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4multresponse_semipar_array_mixed_discrete_mod, 0},
    {"_rcpp_module_boot_stan_fit4multresponse_semipar_array_mixed_randomint_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4multresponse_semipar_array_mixed_randomint_mod, 0},
    {"_rcpp_module_boot_stan_fit4multresponse_semipar_array_mixed_randomint_discrete_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4multresponse_semipar_array_mixed_randomint_discrete_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_bayesGAM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}