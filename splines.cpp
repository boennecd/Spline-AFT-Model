// [[Rcpp::depends(RcppArmadillo)]]
// #define DO_CHECKS
#define ARMA_NO_DEBUG

#include "splines.h"
#include <algorithm> // lower_bound
#include <cmath> // isnan
#include <stdexcept> // invalid_argument
#include <memory.h>

inline void check_splines
  (const arma::vec &boundary_knots, const arma::vec &interior_knots,
   const int order) {
#ifdef DO_CHECKS
  if(order<1, 0)
    throw std::invalid_argument("order<1");
  if(boundary_knots.size() != 2L)
    throw std::invalid_argument("boundary_knots should have length 2");
  if(interior_knots.size()>0 && boundary_knots(0)>min(interior_knots))
    throw std::invalid_argument("boundary_knots(0)>min(interior_knots)");
  if(interior_knots.size()>0 && boundary_knots(1)<max(interior_knots))
    throw std::invalid_argument("boundary_knots(1)<max(interior_knots)");
  // TODO: check if interior_knots are in ascending order?
#endif
}

inline void throw_invalid_out(
    std::string const &cl, unsigned const dim, unsigned const dim_ex){
#ifdef DO_CHECKS
  std::stringstream msg;
  msg << cl << ": invalid 'out' (dim is " << dim << "; expected "
      << dim_ex << ')';
  throw std::invalid_argument(msg.str());
#endif
}

namespace splines {

vec basisMixin::operator()(double const x, int const ders) const {
  vec out(get_n_basis());
  operator()(out, x, ders);
  return out;
}

mat basisMixin::basis(const vec &x, const int ders,
                      const double centre) const {
#ifdef DO_CHECKS
  if (ders < 0)
    throw std::invalid_argument("ders<0");
#endif
  uword const n_basis(get_n_basis()),
              n_x    (x.n_elem);
  rowvec centering =
    (std::isnan(centre) || ders > 0 ?
     zeros(n_basis) : operator()(centre, 0)).t();

  mat out(n_x, n_basis);
  vec wrk(n_basis);
  for (uword i = 0; i < n_x; i++){
    operator()(wrk, x[i], ders);
    out.row(i) = wrk.t() - centering;
  }

  return out;
}

SplineBasis::SplineBasis(const int order): order(order), knots() {
#ifdef DO_CHECKS
  if (order<1)
    throw std::invalid_argument("order<1");
#endif
}


SplineBasis::SplineBasis(const vec knots, const int order):
  order(order), knots(knots) {
#ifdef DO_CHECKS
  if (order<1)
    throw std::invalid_argument("order<1");
#endif
}

void SplineBasis::operator()(
    vec &out, double const x, const int ders) const {
  out.zeros();
#ifdef DO_CHECKS
  if(out.n_elem != SplineBasis::get_n_basis())
    throw_invalid_out(
      "splineBasis", out.n_elem, SplineBasis::get_n_basis());
#endif

  set_cursor(x);
  int io = curs - order;
  if (io < 0 || io > nknots) {
    /* Do nothing. x is already zero by default
    for (size_t j = 0; j < (size_t)order; j++) {
      out(j+io) = double(0); // R_NaN;
    }*/
  } if(ders > 0L){ // faster first order derivative
    derivs(wrk, x, ders);
    for (uword i = 0; i < wrk.n_elem; i++)
      out(i + io) = wrk(i);

  } else { /* fast method for value */
    basis_funcs(wrk, x);
    for (uword i = 0; i < wrk.n_elem; i++)
      out(i + io) = wrk(i);
  }
}

int SplineBasis::set_cursor(const double x) const {
  /* don't assume x's are sorted */
  curs = -1; /* Wall */
  boundary = 0;
  for (int i = 0; i < nknots; i++) {
    if (knots(i) >= x)
      curs = i;
    if (knots(i) > x)
      break;
  }
  if (curs > ncoef) {
    int const lastLegit = ncoef;
    if (x == knots(lastLegit)){
      boundary = 1;
      curs = lastLegit;
    }
  }
  return curs;
}

void SplineBasis::diff_table(const double x, const int ndiff) const {
  for (int i = 0; i < ndiff; i++) {
    rdel(i) = knots(curs + i) - x;
    ldel(i) = x - knots(curs - (i + 1));
  }
}

void SplineBasis::basis_funcs(vec &b, const double x) const {
  diff_table(x, ordm1);
  b(0) = 1;
  for (int j = 1; j <= ordm1; j++) {
    double saved(0);
    for (int r = 0; r < j; r++) { // do not divide by zero
      double const den = rdel(r) + ldel(j - 1 - r);
      if(den != 0) {
        double const term = b(r)/den;
        b(r) = saved + rdel(r) * term;
        saved = ldel(j - 1 - r) * term;
      } else {
        if(r != 0 || rdel(r) != 0)
          b(r) = saved;
        saved = 0.;
      }
    }
    b(j) = saved;
  }
}

void SplineBasis::derivs(vec &b, const double x, int const ders) const {
  diff_table(x, ordm1);
  b(0) = 1;
  for (int j = 1; j <= ordm1; j++) {
    bool const needs_derivs = j >= order - ders;
    double saved(0);
    for (int r = 0; r < j; r++) { // do not divide by zero
      double const den = rdel(r) + ldel(j - 1 - r);
      if(den != 0) {
        if(needs_derivs){
          // follow https://math.stackexchange.com/a/2119661/253239
          // and http://mat.fsv.cvut.cz/gcg/sbornik/prochazkova.pdf
          double const term = j * b(r) / den;
          b(r) = saved - term;
          saved = term;
          continue;
        }

        double const term = b(r)/den;
        b(r) = saved + rdel(r) * term;
        saved = ldel(j - 1 - r) * term;
      } else {
        if(r != 0 || rdel(r) != 0)
          b(r) = saved;
        saved = 0.;
      }
    }
    b(j) = saved;
  }
}

inline arma::vec get_SplineBasis_knots
(vec const &boundary_knots, vec const &interior_knots,
 int const order) {
  check_splines(boundary_knots, interior_knots, order);

  uword const nknots = interior_knots.size() + 2 * order;
  vec knots(nknots);
  for(uword i = 0; i < (uword)order; i++) {
    knots(i) = boundary_knots(0);
    knots(nknots - i - 1) = boundary_knots(1);
  }
  if (interior_knots.size() > 0)
    for(uword i = 0; i < interior_knots.size(); i++)
      knots(i + order) = interior_knots(i);

  return knots;
}

bs::bs(const vec &bk, const vec &ik, const bool inter, const int ord):
  SplineBasis(get_SplineBasis_knots(bk, ik, ord), ord),
  boundary_knots(bk), interior_knots(ik),
  intercept(inter),
  df((int)intercept + order - 1 + interior_knots.size()) {
  check_splines(boundary_knots, interior_knots, order);
}

void bs::operator()(vec &out, double const x, const int ders) const {
#ifdef DO_CHECKS
  if(out.n_elem != bs::get_n_basis())
    throw_invalid_out("bs", out.n_elem, bs::get_n_basis());
#endif
  if (x < boundary_knots(0) || x > boundary_knots(1)) {
    double const k_pivot =
        x < boundary_knots(0) ?
        0.75 * boundary_knots(0) + 0.25 * knots(order) :
        0.75 * boundary_knots(1) + 0.25 * knots(knots.n_elem - order - 2),
      delta = x - k_pivot;

    auto add_term = [&](int const d, double const f = 1){
      bs::operator()(wrks, k_pivot, d);
      out += f * wrks;
    };

    out.zeros();
    if (ders == 0) {
      add_term(0);
      add_term(1, delta);
      add_term(2, delta * delta/2.);
      add_term(3, delta * delta * delta /6.);

    } else if (ders == 1) {
      add_term(1);
      add_term(2, delta);
      add_term(3, delta * delta / 2.);

    } else if (ders == 2) {
      add_term(2);
      add_term(3, delta);

    } else if (ders == 3)
      add_term(3);

    return;
  }

  if(intercept)
    SplineBasis::operator()(out, x, ders);
  else {
    SplineBasis::operator()(wrk, x, ders);
    for(uword i = 1; i < wrk.n_elem; ++i)
      out[i - 1L] = wrk[i];
  }
}

ns::ns(const vec &boundary_knots, const vec &interior_knots,
       const bool intercept, const int order):
  bspline(boundary_knots, interior_knots, true, order),
  intercept(intercept),
  tl0(trans(bspline(boundary_knots(0), 0))),
  tl1(trans(bspline(boundary_knots(0), 1))),
  tr0(trans(bspline(boundary_knots(1), 0))),
  tr1(trans(bspline(boundary_knots(1), 1)))
  { }

void ns::operator()(vec &out, double const x, const int ders) const {
#ifdef DO_CHECKS
  if(out.n_elem != ns::get_n_basis())
    throw_invalid_out("ns", out.n_elem, ns::get_n_basis());
#endif
  if(x < bspline.boundary_knots(0)) {
    if (ders==0){
      out  = tl1;
      out *= x - bspline.boundary_knots(0);
      out += tl0;

    } else if (ders == 1)
      out = tl1;
    else
      out.zeros();

    return;

  } else if (x > bspline.boundary_knots(1)) {
    if (ders==0){
      out  = tr1;
      out *= x - bspline.boundary_knots(1);
      out += tr0;

    } else if (ders==1)
      out = tr1;
    else
      out.zeros();

    return;
  }

  out = trans(bspline(x, ders));
}

vec ns::trans(const vec &x) const {
  vec out = q_matrix * (intercept ? x : x(span(1, x.n_elem - 1)));
  return out(span(2, out.size() - 1));
}
} // namespace splines

// the R interface

// [[Rcpp::export(rng = false)]]
Rcpp::XPtr<splines::bs>
  get_bs_ptr(const arma::vec &knots, arma::vec const &boundary_knots,
             bool const intercept){
  if(boundary_knots.n_elem < 2)
    throw std::invalid_argument("get_ns_ptr: invalid boundary_knots size");

  return Rcpp::XPtr<splines::bs>(
    new splines::bs(boundary_knots, knots, intercept), true);
}

// [[Rcpp::export(rng = false)]]
Rcpp::XPtr<splines::ns>
  get_ns_ptr(const arma::vec &knots, arma::vec const &boundary_knots,
             bool const intercept){
  if(boundary_knots.n_elem < 2)
    throw std::invalid_argument("get_ns_ptr: invalid boundary_knots size");

  return Rcpp::XPtr<splines::ns>(
    new splines::ns(boundary_knots, knots, intercept), true);
}

// evaluates the spline at x.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix eval_spline_basis(arma::vec const &x, SEXP basis_ptr,
                                      int const ders = 0);

// evaluates the spline at x like eval_spline_basis but using pre-allocated
// memory.
// [[Rcpp::export(rng = false)]]
void eval_spline_basis_fill(arma::vec const &x, SEXP basis_ptr,
                            Rcpp::NumericMatrix out, int const ders = 0){
  Rcpp::XPtr<splines::basisMixin> basis(basis_ptr);
  size_t const n_x = x.size(),
             n_col = basis->get_n_basis();

  if(static_cast<unsigned>(out.nrow()) != n_x or
       static_cast<unsigned>(out.ncol()) != n_col)
    throw std::invalid_argument("Invalid out matrix is passed");
  arma::vec wrk(n_col);

  double *oi = &out[0];
  for(unsigned i = 0; i < n_x; ++i, ++oi){
    basis->operator()(wrk, x[i], ders);
    for(unsigned j = 0; j < n_col; ++j)
      oi[j * n_x] = wrk[j];
  }
}

Rcpp::NumericMatrix eval_spline_basis(arma::vec const &x, SEXP basis_ptr,
                                      int const ders){
  Rcpp::XPtr<splines::basisMixin> basis(basis_ptr);
  size_t const n_x = x.size(),
             n_col = basis->get_n_basis();
  Rcpp::NumericMatrix out(n_x, n_col);
  eval_spline_basis_fill(x, basis_ptr, out, ders);
  return out;
}

/*** R
# simple test that ns splines works
local({
  xs <- seq(0, 10, length.out = 100)
  knots <- c(3, 5, 7)
  b_knots <- c(1, 9)
  truth <- unclass(
    splines::ns(xs, knots = knots, Boundary.knots = b_knots, intercept = FALSE))
  f <- get_ns_ptr(knots = knots, boundary_knots = b_knots, intercept = FALSE)

  stopifnot(isTRUE(all.equal(eval_spline_basis(xs, f),
                             truth, check.attributes = FALSE)))

  out <- matrix(0., length(xs), NCOL(truth))
  eval_spline_basis_fill(xs, f, out)
  stopifnot(isTRUE(all.equal(out,
                             truth, check.attributes = FALSE)))

  # check derivs
  truth <- t(sapply(xs, function(x) numDeriv::jacobian(
    splines::ns, x, knots = knots, Boundary.knots = b_knots,
    intercept = FALSE)))

  stopifnot(isTRUE(all.equal(eval_spline_basis(xs, f, ders = 1L),
                             truth, check.attributes = FALSE)))

  eval_spline_basis_fill(xs, f, out, ders = 1L)
  stopifnot(isTRUE(all.equal(out,
                             truth, check.attributes = FALSE)))
})

# simple test that bs splines works
local({
  xs <- seq(0, 10, length.out = 100)
  knots <- c(3, 5, 7)
  b_knots <- c(1, 9)
  truth <- suppressWarnings(unclass(
    splines::bs(xs, knots = knots, Boundary.knots = b_knots,
                intercept = FALSE)))
  f <- get_bs_ptr(knots = knots, boundary_knots = b_knots, intercept = FALSE)

  stopifnot(isTRUE(all.equal(eval_spline_basis(xs, f), truth,
                             check.attributes = FALSE)))

  out <- matrix(0., length(xs), NCOL(truth))
  eval_spline_basis_fill(xs, f, out)
  stopifnot(isTRUE(all.equal(out,
                             truth, check.attributes = FALSE)))

  # check derivs
  truth <- t(sapply(xs, function(x) suppressWarnings(numDeriv::jacobian(
    splines::bs, x, knots = knots, Boundary.knots = b_knots,
    intercept = FALSE))))

  stopifnot(isTRUE(all.equal(eval_spline_basis(xs, f, ders = 1L),
                             truth, check.attributes = FALSE)))

  eval_spline_basis_fill(xs, f, out, ders = 1L)
  stopifnot(isTRUE(all.equal(out,
                             truth, check.attributes = FALSE)))
})
*/

// the class we use the C++ version saft_fit

class saft_fitter {
  arma::mat const X;
  arma::vec const y;
  arma::vec const event;

  double const n_events = arma::sum(event);

  std::unique_ptr<splines::basisMixin> base;

  unsigned const n_obs = X.n_cols,
                n_beta = X.n_rows,
                n_gamm = base->get_n_basis() + 1L;

  // working memory
  mutable arma::vec wrk_vec = arma::vec(base->get_n_basis());

  inline double dot_prod_beta
    (unsigned const idx, arma::vec const &beta) const noexcept {
    double out(0);
    double const * __restrict__ xi = X.colptr(idx),
                 * __restrict__ bi = beta.memptr();
    for(unsigned i = 0; i < n_beta; ++i)
      out += *xi++ * *bi++;
    return out;
  }

  inline double eval_basis
    (unsigned idx, double const point, arma::vec const &gamma,
     int const ders = 0L) const noexcept {
    double const * __restrict__ gi = gamma.begin();
    double out(*gi++ * (ders == 0)); // the intercept
    base->operator()(wrk_vec, point, ders);
    double const * wi = wrk_vec.begin();
    for(unsigned i = 1; i < n_gamm; ++i)
      out += *wi++ * *gi++;

    return out;
  }

public:
  enum supported_basis { ns, bs };

  saft_fitter(arma::mat const &X, arma::vec const &y, arma::vec const &event,
              arma::vec const &knots, arma::vec const &boundary_knots,
              supported_basis const which_base):
    X(X.t()), // we work with X transpose
    y(y), event(event),
    base(([&]() -> std::unique_ptr<splines::basisMixin> {
      switch(which_base){
      case ns:
        return std::unique_ptr<splines::basisMixin>(
          new splines::ns(boundary_knots, knots, false));
      case bs:
        return std::unique_ptr<splines::basisMixin>(
          new splines::bs(boundary_knots, knots, false));
      default:
        throw std::invalid_argument("which_base not supported");
      }

      return { nullptr };
    })()) {
    if(y.size() != n_obs)
      throw std::invalid_argument("invalid y");
    if(event.size() != n_obs)
      throw std::invalid_argument("invalid event");
  }

  double eval_log_likelihood(arma::vec const &beta, arma::vec const &gamma,
                             arma::vec nodes, arma::vec weights){
    if(beta.size() != n_beta)
      throw std::invalid_argument("invalid beta");
    if(gamma.size() != n_gamm)
      throw std::invalid_argument("invalid gamma");

    // rescale and relocate the nodes and weights
    nodes = (nodes + 1) * .5;
    weights *= .5;

    unsigned const n_nodes = nodes.size();

    // compute the log-likelihood observation by observation
    double out(0.);
    for(unsigned i = 0; i < n_obs; ++i){
      double const eta = dot_prod_beta(i, beta),
                     w = std::exp(eta),
                   w_y = w * y[i];

      if(event[i] > 0)
        // add the log hazard terms
        out -= eta + eval_basis(i, w_y, gamma);

      // add the terms from the log survival function
      for(unsigned j = 0; j < n_nodes; ++j){
        double const w_y_x = w_y * nodes[j];
        out += w_y * weights[j] * std::exp(eval_basis(i, w_y_x, gamma));
      }
    }

    return out;
  }

  arma::vec eval_log_likelihood_gr(
      arma::vec const &beta, arma::vec const &gamma, arma::vec nodes,
      arma::vec weights){
    if(beta.size() != n_beta)
      throw std::invalid_argument("invalid beta");
    if(gamma.size() != n_gamm)
      throw std::invalid_argument("invalid gamma");

    // rescale and relocate the nodes and weights
    nodes = (nodes + 1) * .5;
    weights *= .5;

    unsigned const n_nodes = nodes.size();
    arma::vec out(n_beta + n_gamm, arma::fill::zeros);
    double * const __restrict__ d_beta_begin  = out.memptr(),
           * const __restrict__ d_gamma_begin = d_beta_begin + n_beta;

    // hazard terms
    *d_gamma_begin = -n_events;

    // compute the log-likelihood observation by observation
    for(unsigned i = 0; i < n_obs; ++i){
      double const eta = dot_prod_beta(i, beta),
                     w = std::exp(eta),
                   w_y = w * y[i];

      if(event[i] > 0){
        // add the log hazard terms
        double const fac_beta = 1 + eval_basis(i, w_y, gamma, 1L) * w_y,
                     * xj = X.colptr(i);
        for(unsigned j = 0; j < n_beta; ++j)
          d_beta_begin[j] -= fac_beta * *xj++;

        base->operator()(wrk_vec, w_y);
        double const * wj = wrk_vec.begin();
        for(unsigned j = 1; j < n_gamm; ++j)
          d_gamma_begin[j] -= *wj++;
      }

      // add the terms from the log survival function
      for(unsigned j = 0; j < n_nodes; ++j){
        double const w_y_x = w_y * nodes[j],
                       fac = w_y * weights[j] *
                         std::exp(eval_basis(i, w_y_x, gamma));
        *d_gamma_begin += fac;
        double const * wj = wrk_vec.begin();
        for(unsigned j = 1; j < n_gamm; ++j)
          d_gamma_begin[j] += fac * *wj++;

        double const d_log_lambda0 = eval_basis(i, w_y_x, gamma, 1L),
                          fac_beta = fac * (1 + d_log_lambda0 * w_y_x);
        double const * xj = X.colptr(i);
        for(unsigned j = 0; j < n_beta; ++j)
          d_beta_begin[j] += fac_beta * *xj++;
      }
    }

    return out;
  }
};

// [[Rcpp::export(rng = false)]]
Rcpp::XPtr<saft_fitter> get_saft_fitter_ptr
  (arma::mat const &X, arma::vec const &y, arma::vec const &event,
   arma::vec const &knots, arma::vec const &boundary_knots,
   std::string const which_base){
  saft_fitter::supported_basis which_base_arg;
  if(which_base == "ns")
    which_base_arg = saft_fitter::supported_basis::ns;
  else if(which_base == "bs")
    which_base_arg = saft_fitter::supported_basis::bs;
  else
    throw std::invalid_argument("which_base is not supported");

  return Rcpp::XPtr<saft_fitter>(new saft_fitter(
    X, y, event, knots, boundary_knots, which_base_arg));
}

// [[Rcpp::export(rng = false)]]
double eval_log_likelihood_cpp
  (arma::vec const &beta, arma::vec const &gamma,
   arma::vec const &nodes, arma::vec const &weights,
   SEXP ptr){
  Rcpp::XPtr<saft_fitter> saft_fitter_ptr(ptr);

  return saft_fitter_ptr->eval_log_likelihood(beta, gamma, nodes, weights);
}

// [[Rcpp::export(rng = false)]]
arma::vec eval_log_likelihood_grad_cpp
  (arma::vec const &beta, arma::vec const &gamma,
   arma::vec const &nodes, arma::vec const &weights,
   SEXP ptr){
  Rcpp::XPtr<saft_fitter> saft_fitter_ptr(ptr);

  return saft_fitter_ptr->eval_log_likelihood_gr(beta, gamma, nodes, weights);
}
