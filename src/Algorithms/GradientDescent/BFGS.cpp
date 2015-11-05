/*!
 * 
 *
 * \brief       BFGS
 * 
 * The Broyden, Fletcher, Goldfarb, Shannon (BFGS) algorithm is a
 * quasi-Newton method for unconstrained real-valued optimization.
 * 
 * 
 *
 * \author      O. Krause 
 * \date        2010
 *
 *
 * \par Copyright 1995-2015 Shark Development Team
 * 
 * <BR><HR>
 * This file is part of Shark.
 * <http://image.diku.dk/shark/>
 * 
 * Shark is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published 
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Shark is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with Shark.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
 #define SHARK_COMPILE_DLL
#include <shark/Algorithms/GradientDescent/BFGS.h>

using namespace shark;

void BFGS::initModel(){
	m_hessian.resize(m_dimension, m_dimension);
	m_hessian.clear();
	for (size_t i = 0; i < m_dimension; ++i)
	{
		m_hessian(i, i) = 1.;
	}
}

void BFGS::setHessian(const RealMatrix &hessian){
	m_hessian = hessian;
	axpy_prod(m_hessian,m_derivative,m_searchDirection);
	m_searchDirection *= -1;
	m_initialStepLength = 1.0;
}

void BFGS::computeSearchDirection(){
	RealVector gamma = m_derivative - m_lastDerivative;
	RealVector delta = m_best.point - m_lastPoint;
	double d = inner_prod(gamma,delta);
	
	RealVector Hg(m_dimension,0.0);
	axpy_prod(m_hessian,gamma,Hg);
	double gHg = inner_prod(gamma,Hg);

	//compute damped step
	//Method doi:10.1016/j.jcp.2013.08.044
	//Sigma parameters doi:10.1007/s10957-013-0448-8
	double sigma_2 = 0.6, sigma_3 = 3.0;
	double tau = d / gHg;

	double theta = 1.0;
	if (tau > 1 + sigma_3)
	{
		theta = sigma_3 / (tau - 1);
	}
	else if (tau < 1 - sigma_2)
	{
		theta = sigma_2 / (1 - tau);
	}

	delta = theta * delta + (1 - theta) * Hg;
	d = inner_prod(gamma,delta);

	//update hessian
	double scale = (gHg / d + 1) / d;

	m_hessian += scale * outer_prod(delta,delta)
		  - (outer_prod(Hg,delta)+outer_prod(delta,Hg))/d;

	//compute search direction
	axpy_prod(m_hessian,m_derivative,m_searchDirection);
	m_searchDirection *= -1;
}

//from ISerializable
void BFGS::read( InArchive & archive )
{
	AbstractLineSearchOptimizer::read(archive);
	archive>>m_hessian;
}

void BFGS::write( OutArchive & archive ) const
{
	AbstractLineSearchOptimizer::write(archive);
	archive<<m_hessian;
}