//===========================================================================
/*!
 *  \brief Minimizing functions of "N"-dimensional variables.
 *
 *
 *  \author  T. Trnka
 *  \date    2014
 *
 *  \par Copyright (c) 1998-2000:
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-25558<BR>
 *      Fax:   +49-234-32-14209<BR>
 *      eMail: shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
 *      www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
 *      <BR>
 *
 *
 *  <BR><HR>
 *  This file is part of Shark. This library is free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software
 *  Foundation; either version 3, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this library; if not, see <http://www.gnu.org/licenses/>.
 *
 *
 *
 */
//===========================================================================

#ifndef SHARK_ALGORITHMS_GRADIENTDESCENT_IMPL_LINMIN_INL
#define SHARK_ALGORITHMS_GRADIENTDESCENT_IMPL_LINMIN_INL

#include <shark/Core/Math.h>

namespace shark{ namespace detail{

//===========================================================================
/*!
 *  \brief Minimizes a function of "N" variables.
 *
 *  Performs a minimization of a function of \f$ N \f$ variables, i.e.
 *  given as input the vectors \f$ P \f$ and \f$ n \f$ and the function
 *  \f$ f \f$, the function finds the scalar \f$ \lambda \f$
 *  that minimizes \f$ f(P + \lambda n) \f$.
 *  \f$ P \f$ is replaced by \f$ P + \lambda n \f$ and \f$ n \f$
 *  by \f$ \lambda n \f$.
 *
 *  \param  point     N-dimensional initial starting point for the search, is set to the point where the function takes on a minimum.
 *  \param  searchDirection    N-dimensional search direction, is replaced by the actual vector displacement that \em p was moved.
 *  \param  value  The function value at the new point \em p.
 *  \param  func  The function that will be minimized.
 *  \param  ax    guess for the lower bracket
 *  \param  bx    guess for the upper bracket
 *  \return       none.
 *  \throw SharkException the type of the exception will be "size mismatch" and indicates that \em p is not one-dimensional or that \em p and \em searchDirection don't have the same length
 */
template<class VectorT,class VectorU,class Function>
void linmin
(
	VectorT& point,
	const VectorU& searchDirection,
	double& value,
	Function & func,
	double ax,
	double bx
)
{
	SIZE_CHECK(point.size()==searchDirection.size());

	const double   GOLD   = 1.618034;
	const double   GLIMIT = 100.;
	const double   TINY   = 1.0e-20;
	const unsigned ITMAX  = 100;
	const double   ZEPS   = 1.0e-10;
	const double   TOL    = 2.0e-4;
        const double   CGOLD  = 0.3819660;

	unsigned iter;
	double   fa, fb, fc, fp, cx;
	double   ulim, dum;
	double   a, b, d(0.), e, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm, etemp;

	unsigned        n = point.size();
	VectorT xt(n);

	//===================================================================

	fa = fp = func.eval(point);

	noalias(xt) = point + searchDirection;
	fb = func.eval(xt);

	//ensure that fb <= fa
	if (fb > fa) {
		dum = ax;
		ax = bx;
		bx = dum;
		dum = fb;
		fb = fa;
		fa = dum;
	}

	//evaluate function a golden selection
	cx = bx + GOLD * (bx - ax);
	noalias(xt) = point + cx*searchDirection;
	fc = func.eval(xt);

	while (fb > fc) {//find interval containing the minimum
		r = (bx - ax) * (fb - fc);
		q = (bx - cx) * (fb - fa);
		u = bx - ((bx - cx) * q - (bx - ax) * r) /
			(2. * copySign(std::max(std::abs(q - r), TINY), q - r));
		ulim = bx + GLIMIT * (cx - bx);
		if ((bx - u) *(u - cx) > 0.) {//u is in the interval (cx,bx)
			noalias(xt) = point + u*searchDirection;
			fu = func.eval(xt);
			if (fu < fc) {//continue in the interval [bx,u]
				ax = bx;
				bx = u;
				fa = fb;
				fb = fu;
				break;
			}
			else if (fu > fb) {
				cx = u;
				fc = fu;
				break;
			}
			//evaluate new golden selection point
			u = cx + GOLD * (cx - bx);
			noalias(xt) = point + u*searchDirection;
			fu = func.eval(xt);
		}
		else if ((cx - u) *(u - ulim) > 0.) {//u is in the interval [ulim,cx]
			noalias(xt) = point + u*searchDirection;
			fu = func.eval(xt);
			if (fu < fc) {//continue in the interval [ax,cx]
				bx = cx;
				cx = u;
				u  = cx + GOLD * (cx - bx);
				fb = fc;
				fc = fu;
				noalias(xt) = point + u*searchDirection;
				fu = func.eval(xt);
			}
		}
		else if ((u - ulim) *(ulim - cx) >= 0.) {//ulim is in the interval [cx,u]
			u = ulim;
			noalias(xt) = point + u*searchDirection;
			fu = func.eval(xt);
		}
		else {
			u = cx + GOLD * (cx - bx);
			noalias(xt) = point + u*searchDirection;
			fu = func.eval(xt);
		}
		//exchange interval limits and corresponding function values
		ax = bx;
		bx = cx;
		cx = u;
		fa = fb;
		fb = fc;
		fc = fu;
	}

	//=======================================================================

	e = 0.;
	if (ax < cx) {
		a = ax;
		b = cx;
	}
	else {
		a = cx;
		b = ax;
	}

	x = w = v = bx;
	noalias(xt) = point + x*searchDirection;
	fw = fv = fx = func.eval(xt);

	for (iter = 0; iter < ITMAX; iter++) {
		xm = 0.5 * (a + b);//interval centrum
		tol2 = 2. * (tol1 = TOL * std::abs(x) + ZEPS);
		if (std::abs(x - xm) <= (tol2 - 0.5 *(b - a))) {
			break;
		}
		// calculate delta x
		if (std::abs(e) > tol1) {
                        r = (x - w) * (fx - fv);
                        q = (x - v) * (fx - fw);
                        p = (x - v) * q - (x - w) * r;
                        q = 2.0 * (q - r);
                        if (q > 0) {
                                p = -p;
                        }
                        q = std::abs(q);
			etemp = e;
			e = d;
                        if (std::abs(p) >= std::abs(0.5 * q  *etemp) || p <= q * (a - x) || p >= q * (b - x)) {
                                e = (x >= xm) ? a - x : b - x;
                                d = CGOLD * e;
                        }
                        else {
                                d = p / q;
                                u = x + d;
                                if (u - a < tol2 || b - u < tol2) {
                                        d = copySign(tol1, xm - x);
                                }
                        }
		}
		else {
                        e = (x >= xm) ? a - x : b - x;
                        d = CGOLD * e;
		}
		//update u
		if (std::abs(d) >= tol1) {
			u = x + d;
		}
		else {
			u = x + copySign(tol1, d);
		}
                noalias(xt) = point + u*searchDirection;
                fu = func.eval(xt);
		//reduce interval length
		if (fu <= fx) {
			if (u >= x)
				a = x;
			else
				b = x;
			v  = w;
			w  = x;
			x  = u;
			fv = fw;
			fw = fx;
			fx = fu;
		}
		else {
			if (u < x)
				a = u;
			else
				b = u;
			if (fu <= fw || w == x) {
				v  = w;
				w  = u;
				fv = fw;
				fw = fu;
			}
			else if (fu <= fv || v == x || v == w) {
				v  = u;
				fv = fu;
			}
		}
	}

        value = fx;
        noalias(point) += searchDirection * x;
}

}}
#endif
