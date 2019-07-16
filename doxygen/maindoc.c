/** @mainpage RBmatlab 0.11.04
 * @section intro Introduction
 *
 * RBmatlab is a MATLAB&copy; library providing routines for the solution of
 * numerical schemes based on partial differential equations. It includes
 *
 * -# low-level methods for
 *  - @ref grid "Generation of grids",
 *  - @ref vector "Vector based computations",
 *  - @ref datafunc "Useful data functions for initial and boundary conditions"
 *  and
 *  - @ref visual "Data visualization"
 * -# @ref discfunc "Discretization methods" for
 *   - @ref fv "Finite volume discretizations" and
 *   - @ref ldg "Local discontinuous Galerkin discretizations"
 * -# high-level Reduced basis routines for
 *   - @ref basisgen "Generation of a reduced basis and online matrices" and
 *   - @ref ei "Empirical interpolation of discrete operators"
 * -# @ref models "Implementations of models" for
 *   - @ref lin_evol "Linear evolution equations",
 *   - @ref nonlin_evol "Non-linear evolution equations" and
 *   - @ref lin_ds "Dynamical systems for linear partial differential equations"
 * -# @ref demos
 *
 * @par Stationary problems:
 * Although, the current development in RBmatlab was concentrated on evolution
 * schemes, the framework also allows stationary schemes.
 *
 * @par Usage with other numerical solvers
 * The reduced basis generation algorithms can also make use of detailed
 * simulations from other software packages, if those stick to the interface of
 * the RBmatlab framework. An example is the Dune module dune-rb for which a
 * model implementation exists
 * @link ./models/convdiff/convdiff_dune_model.m here @endlink.
 *
 * @section installation Installation
 *
 * After you have downloaded the <a
 * href="http://morepas.org/software/rbmatlab/0.11.04/rbmatlab-release-0.11.04.tar.gz">
 * Tarball</a> from the website and extracted it into a directory @c
 * /path/to/rbmatlab, installation is very simple:
 *
 * Make sure, that the following code is executed during matlab startup,
 * e.g. by putting it into the <tt>$HOME/matlab/startup.m</tt> file or by
 * executing it every time you start matlab by cut'n paste.
 * In particular you must set the path to the RBmatlab directory and choose a
 * temporary-data directory as summarized below
 *
 * @code
 * setenv('RBMATLABHOME','/path/to/rbmatlab');
 * setenv('RBMATLABTEMP','/tmp/matlab');
 * addpath(getenv('RBMATLABHOME'));
 * startup_rbmatlab
 * @endcode
 *
 * The environment variable @c RBMATLABTEMP must point to an existing
 * directory, where RBmatlab can create subdirectories for storing temporary
 * data. Large space, i.e. several gigabytes should be available.
 *
 * Alternatively you could also set the environment variables in the shell from
 * which you start MATLAB. Then you need to make sure, that the following code
 * is executed during MATLAB startup, e.g. by putting it into the @c
 * $HOME/matlab/startup.m file or by executing it every time you start MATLAB
 * by cut'n paste.
 *
 * @code
 * addpath(getenv('RBMATLABTEMP'));
 * startup_rbmatlab
 * @endcode
 *
 * @section learning How To learn RBmatlab
 *
 * Users new to RBmatlab should first get a copy of the partially available
 * RBmatlab-HowTo and look at the @ref demos "Demos" or the @ref scripts
 * "Scripts" section. Furthermore, there is a directory with @ref test
 * "regression tests" which could also be helpful.
 *
 * @section mtoc How To write doxygen documentation
 *
 * If you want to write your own documentation for your classes and functions,
 * you need to format the comment strings in your M-Files according to this
 * @ref doxygen.m "example" file. The documentation strings should be
 * compatible with ones expected by the MATLAB documentation tools, such that
 * @code help function @endcode should work as well at the MATLAB prompt.
 *
 * The documentation can be build automatically with the script
 * <tt>./make_docu.sh</tt> in the <tt>doxygen/</tt> directory. Currently, this
 * script needs a doxygen installation of version 1.6.2 or newer and uses a
 * filter program which is shipped for 'Linux x86_64' architectures only, but
 * can be compiled for other architectures as well.
 *
 * @section license License issues
 *
 * RBmatlab is published under the terms of the <a
 * href="http://www.opensource.org/licenses/afl-3.0.php">Academic Free License
 * 3.0</a>.
 *
 * Software using source files of this project or significant parts of it,
 * should include the following attribution notice:
 *
 * @code
 * % --------------------------------------------------------------------
 * % ATTRIBUTION NOTICE:
 * % This product includes software developed for the RBmatlab project at
 * % (C) Universities of Stuttgart and Münster, Germany.
 * %
 * % RBmatlab is a MATLAB software package for model reduction with an
 * % emphasis on Reduced Basis Methods. The project is maintained by
 * % M. Dihlmann, M. Drohmann, B. Haasdonk, M. Ohlberger and M. Schaefer.
 * % For Online Documentation and Download we refer to www.morepas.org.
 * % --------------------------------------------------------------------
 * @endcode
 */

