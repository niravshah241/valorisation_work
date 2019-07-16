/** @defgroup grid Rectangular and triangular grid generation
 */

/** @defgroup rectgrid Rectangular 2-d grid
 * @ingroup grid
 */

/** @defgroup cubegrid Cubical n-dimensional grid
 * @ingroup grid
 */

/** @defgroup triagrid Triangular grid
 * @ingroup grid
 */

/** @defgroup vector Vector based functions
 */

/** @defgroup datafunc Data functions
 */

/** @defgroup visual Data visualization routines
 */

/** @defgroup discfunc Discretization methods
 */

/** @defgroup fv Finite volume discretizations
 * @ingroup discfunc
 */

/** @defgroup ldg Local discontinuous galerkin discretizations
 * @ingroup discfunc
 */

/** @defgroup fem Finite element discretizations
 * @ingroup discfunc
 */

/** @defgroup models Model implementations
 */

/** @defgroup model_examples Example of model implementations
 * @ingroup models
 */

/** @defgroup advection_output Linear dynamical system with advection and an output functional
 *
 * @ingroup model_examples
 */

/** @defgroup richards_fv Non-linear evolution equation with geometry transformation and an example of the richards' equation
 *
 * @par Included in the following presentations:
 *   - Algoritmy 2009
 *   - MoRePaS 2009
 *
 * @par Included in the following papers:
 *    - Diploma thesis
 *    - Algoritmy proceedings 2009
 *
 * @sa demo_richards_fv()
 *
 * @ingroup model_examples
 */

/** @defgroup convdiff Linear convection-diffusion examples with discretizations in matlab or dune-rb
 * @ingroup model_examples
 */

/** @defgroup lin_evol Linear evolution equations
 * @ingroup models
 */

/** @defgroup nonlin_evol Non-linear evolution equations
 * @ingroup models
 */

/** @defgroup lin_ds Dynamical systems for linear partial differential equations
 * @ingroup models
 */

/** @defgroup basisgen Reduced basis generation routines
 */

/** @defgroup ei Empirical interpolation of discrete operators
 *
 * @section Prerequisites
 *
 * The above functions can be used to empirically interpolate a discrete
 * operator @f$ L(\mu) : W_h \to W_h @f$ with an @f$H@f$-independent
 * Dof-dependence and acting on an arbitrary discrete function space @f$W_h@f$.
 * For details on the theory, we refer to
 * <a href="http://wwwmath.uni-muenster.de/num/publications/2010/DHO10/DHO10_nonlin_evol_preprint.pdf">[DHO10]</a>.
 *
 * @section Details
 *
 * The empirical interpolation is done in two steps:
 *
 *   -# Generate trajectories of solution vector by computing detailed
 *      simulations for parameters from a specific parameter set @f$M \in {\cal
 *      P} \in \mathbb{R}^p@f$ and compute corresponding operator evalutions on
 *      these solution vectors with a call like
 *      @code
 *        LU_fnames = ei_operator_collect_files(model, model_data, Mtrain, params);
 *      @endcode
 *   -# Compute the collateral reduced basis space and interpolation points for
 *      an operator @f$L_q@f$ interpolating evaluations on subtrajectories
 *      lying in the time slice 't'.
 *      @code
 *        params.index          = q;
 *        params.time_split_map = [ 1:model.nt+1; vec_of_slice_indices ];
 *        params.time_index     = t;
 *        detailed_data = ei_detailed(model, model_data, LU_fnames, params);
 *      @endcode
 */

/** @defgroup scripts Scripts
 */

/** @defgroup test Regression tests
 */

/** @defgroup demos Demos
 */

