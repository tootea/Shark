//===========================================================================
/*!
 * 
 *
 * \brief      Contains the lapack bindings for the symmetric eigenvalue problem syev.
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
//===========================================================================
#ifndef SHARK_LINALG_BLAS_KERNELS_LAPACK_SYEV_HPP
#define SHARK_LINALG_BLAS_KERNELS_LAPACK_SYEV_HPP

#include "fortran.hpp"
#include "../traits.hpp"
#include "../../vector.hpp"

#define SHARK_LAPACK_SSYEV FORTRAN_ID(ssyev)
#define SHARK_LAPACK_DSYEV FORTRAN_ID(dsyev)
#define SHARK_LAPACK_SSYEVD FORTRAN_ID(ssyevd)
#define SHARK_LAPACK_DSYEVD FORTRAN_ID(dsyevd)
#define SHARK_LAPACK_SSYEVR FORTRAN_ID(ssyevr)
#define SHARK_LAPACK_DSYEVR FORTRAN_ID(dsyevr)

#define EIGENSYMM_USE_SYEVD

typedef int fortran_int_t;

extern "C" {
void SHARK_LAPACK_SSYEV( const char* jobz, const char* uplo, const fortran_int_t* n,
        float* a, const fortran_int_t* lda, float* w, float* work,
        const fortran_int_t* lwork, fortran_int_t* info );
void SHARK_LAPACK_DSYEV( const char* jobz, const char* uplo, const fortran_int_t* n,
        double* a, const fortran_int_t* lda, double* w, double* work,
        const fortran_int_t* lwork, fortran_int_t* info );

void SHARK_LAPACK_SSYEVD( const char* jobz, const char* uplo,
        const fortran_int_t* n, float* a, const fortran_int_t* lda, float* w,
        float* work, const fortran_int_t* lwork, fortran_int_t* iwork,
        const fortran_int_t* liwork, fortran_int_t* info );
void SHARK_LAPACK_DSYEVD( const char* jobz, const char* uplo,
        const fortran_int_t* n, double* a, const fortran_int_t* lda,
        double* w, double* work, const fortran_int_t* lwork,
        fortran_int_t* iwork, const fortran_int_t* liwork,
        fortran_int_t* info );

void SHARK_LAPACK_SSYEVR( const char* jobz, const char* range, const char* uplo,
        const fortran_int_t* n, float* a, const fortran_int_t* lda,
        const float* vl, const float* vu, const fortran_int_t* il,
        const fortran_int_t* iu, const float* abstol, fortran_int_t* m,
        float* w, float* z, const fortran_int_t* ldz, fortran_int_t* isuppz,
        float* work, const fortran_int_t* lwork, fortran_int_t* iwork,
        const fortran_int_t* liwork, fortran_int_t* info );
void SHARK_LAPACK_DSYEVR( const char* jobz, const char* range, const char* uplo,
        const fortran_int_t* n, double* a, const fortran_int_t* lda,
        const double* vl, const double* vu, const fortran_int_t* il,
        const fortran_int_t* iu, const double* abstol, fortran_int_t* m,
        double* w, double* z, const fortran_int_t* ldz, fortran_int_t* isuppz,
        double* work, const fortran_int_t* lwork, fortran_int_t* iwork,
        const fortran_int_t* liwork, fortran_int_t* info );
};

namespace shark { namespace blas { namespace bindings {

enum LAPACK_UPLO {
	Upper, Lower
};

enum LAPACK_EVJOB {
    ValuesOnly, Vectors
};

enum LAPACK_EVR_RANGE {
    All, ValueInterval, IndexInterval
};

static inline char lapack_evjob_char(LAPACK_EVJOB Job)
{
    switch (Job) {
        case ValuesOnly:
            return 'N';
        case Vectors:
            return 'V';
    }
    return '\0'; // undefined LAPACK_EVJOB value
}

static inline char lapack_evr_range_char(LAPACK_EVR_RANGE Range)
{
    switch (Range) {
        case All:
            return 'A';
        case ValueInterval:
            return 'V';
        case IndexInterval:
            return 'I';
    }
    return '\0'; // undefined LAPACK_EVR_RANGE value
}

static inline char lapack_uplo_char(LAPACK_UPLO Uplo)
{
    switch (Uplo) {
        case Upper:
            return 'U';
        case Lower:
            return 'L';
    }
    return '\0'; // undefined LAPACK_UPLO value
}

inline int syev( char Job, char Uplo, fortran_int_t N,
        float* A, fortran_int_t lda, float* w, float* work,
        fortran_int_t lwork )
{
    fortran_int_t info = 0;
    SHARK_LAPACK_SSYEV(&Job, &Uplo, &N, A, &lda, w, work, &lwork, &info);
    return info;
}

inline int syev( char Job, char Uplo, fortran_int_t N,
        double* A, fortran_int_t lda, double* w,
        double* work, fortran_int_t lwork )
{
    fortran_int_t info = 0;
    SHARK_LAPACK_DSYEV(&Job, &Uplo, &N, A, &lda, w, work, &lwork, &info);
    return info;
}

inline int syevd( char Job, char Uplo, fortran_int_t N,
        float* A, fortran_int_t lda, float* w,
        float* work, fortran_int_t lwork,
        fortran_int_t *iwork, fortran_int_t liwork)
{
    fortran_int_t info = 0;
    SHARK_LAPACK_SSYEVD(&Job, &Uplo, &N, A, &lda, w, work, &lwork, iwork, &liwork, &info);
    return info;
}

inline int syevd( char Job, char Uplo, fortran_int_t N,
        double* A, fortran_int_t lda, double* w,
        double* work, fortran_int_t lwork,
        fortran_int_t *iwork, fortran_int_t liwork)
{
    fortran_int_t info = 0;
    SHARK_LAPACK_DSYEVD(&Job, &Uplo, &N, A, &lda, w, work, &lwork, iwork, &liwork, &info);
    return info;
}

inline int syevr( char Job, char Range, char Uplo, fortran_int_t N,
        float* A, fortran_int_t lda, float vl, float vu, fortran_int_t il, fortran_int_t iu,
        float abstol, fortran_int_t* m, float* w,
        float* z, fortran_int_t ldz, fortran_int_t* isuppz,
        float* work, fortran_int_t lwork,
        fortran_int_t* iwork, fortran_int_t liwork)
{
    fortran_int_t info = 0;
    SHARK_LAPACK_SSYEVR(&Job, &Range, &Uplo, &N, A, &lda, &vl, &vu, &il, &iu, &abstol,
            m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);
    return info;
}

inline int syevr( char Job, char Range, char Uplo, fortran_int_t N,
        double* A, fortran_int_t lda, double vl, double vu, fortran_int_t il, fortran_int_t iu,
        double abstol, fortran_int_t* m, double* w,
        double* z, fortran_int_t ldz, fortran_int_t* isuppz,
        double* work, fortran_int_t lwork,
        fortran_int_t* iwork, fortran_int_t liwork)
{
    fortran_int_t info = 0;
    SHARK_LAPACK_DSYEVR(&Job, &Range, &Uplo, &N, A, &lda, &vl, &vu, &il, &iu, &abstol,
            m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);
    return info;
}


template<typename T>
inline int syev( LAPACK_EVJOB Job, LAPACK_UPLO Uplo, fortran_int_t N,
        T* A, fortran_int_t lda, T* w, T* work,
        fortran_int_t lwork )
{
    return syev(lapack_evjob_char(Job), lapack_uplo_char(Uplo), N, A, lda, w, work, lwork);
}

template<typename T>
inline int syevd( LAPACK_EVJOB Job, LAPACK_UPLO Uplo, fortran_int_t N,
        T* A, fortran_int_t lda, T* w,
        T* work, fortran_int_t lwork,
        fortran_int_t *iwork, fortran_int_t liwork)
{
    return syevd(lapack_evjob_char(Job), lapack_uplo_char(Uplo), N, A, lda, w, work, lwork, iwork, liwork);
}

template<typename T>
inline int syevr( LAPACK_EVJOB Job, LAPACK_EVR_RANGE Range, LAPACK_UPLO Uplo, fortran_int_t N,
        T* A, fortran_int_t lda, T vl, T vu, fortran_int_t il, fortran_int_t iu,
        T abstol, fortran_int_t* m, T* w,
        T* z, fortran_int_t ldz, fortran_int_t* isuppz,
        T* work, fortran_int_t lwork,
        fortran_int_t* iwork, fortran_int_t liwork)
{
    return syevr(lapack_evjob_char(Job), lapack_evr_range_char(Range), lapack_uplo_char(Uplo), N, A, lda,
                 vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork);
}


template <typename SymmA, typename VectorW>
inline int syev( LAPACK_EVJOB Job, LAPACK_UPLO Uplo,
          matrix_expression<SymmA> &A, vector_expression<VectorW> &w)
{
    typedef typename SymmA::value_type data_type;
    data_type opt_work_size_real;
    fortran_int_t N = A().size1();
    SIZE_CHECK(N == A().size2());

    syev(Job, Uplo, N, traits::storage(A()), traits::leading_dimension(A()),
         traits::storage(w()), &opt_work_size_real, -1);
    fortran_int_t opt_work_size = static_cast<fortran_int_t>(opt_work_size_real);

    shark::blas::vector<data_type> workspace(opt_work_size);

    return syev(Job, Uplo, N, traits::storage(A()), traits::leading_dimension(A()),
         traits::storage(w()), workspace.storage(), opt_work_size);
}

template <typename SymmA, typename VectorW>
inline int syevd( LAPACK_EVJOB Job, LAPACK_UPLO Uplo,
          matrix_expression<SymmA> &A, vector_expression<VectorW> &w)
{
    typedef typename SymmA::value_type data_type;
    data_type opt_work_size_real;
    fortran_int_t opt_iwork_size;
    fortran_int_t N = A().size1();
    SIZE_CHECK(N == A().size2());

    syevd(Job, Uplo, N, traits::storage(A()), traits::leading_dimension(A()),
         traits::storage(w()), &opt_work_size_real, -1, &opt_iwork_size, -1);
    fortran_int_t opt_work_size = static_cast<fortran_int_t>(opt_work_size_real);

    shark::blas::vector<data_type> workspace(opt_work_size);
    shark::blas::vector<fortran_int_t> iworkspace(opt_iwork_size);

    return syevd(Job, Uplo, N, traits::storage(A()), traits::leading_dimension(A()),
         traits::storage(w()), workspace.storage(), opt_work_size, iworkspace.storage(), opt_iwork_size);
}

template <typename SymmA, typename VectorW, typename MatrixZ, typename VectorS>
inline int syevr( LAPACK_EVJOB Job, LAPACK_EVR_RANGE Range, LAPACK_UPLO Uplo, matrix_expression<SymmA> &A,
                  typename SymmA::value_type vl, typename SymmA::value_type vu,
                  fortran_int_t il, fortran_int_t iu, typename SymmA::value_type abstol,
                  fortran_int_t &m, vector_expression<VectorW> &w,
                  matrix_expression<MatrixZ> &Z, vector_expression<VectorS> &s)
{
    typedef typename SymmA::value_type data_type;
    data_type opt_work_size_real;
    fortran_int_t opt_iwork_size;
    fortran_int_t N = A().size1();
    SIZE_CHECK(N == A().size2());

    syevd(Job, Uplo, N, traits::storage(A()), traits::leading_dimension(A()),
         traits::storage(w()), &opt_work_size_real, -1, &opt_iwork_size, -1);
    fortran_int_t opt_work_size = static_cast<fortran_int_t>(opt_work_size_real);

    shark::blas::vector<data_type> workspace(opt_work_size);
    shark::blas::vector<fortran_int_t> iworkspace(opt_iwork_size);

    return syevr(Job, Range, Uplo, N, traits::storage(A()), traits::leading_dimension(A()),
                 vl, vu, il, iu, abstol, &m, traits::storage(w()), traits::storage(Z()),
                 traits::leading_dimension(Z()), traits::storage(s()),
                 workspace.storage(), opt_work_size, iworkspace.storage(), opt_iwork_size);
}

template <typename SymmA, typename VectorW, typename MatrixZ, typename VectorS>
inline int syevr( LAPACK_EVJOB Job, LAPACK_UPLO Uplo, matrix_expression<SymmA> &A,
                  typename SymmA::value_type abstol,
                  vector_expression<VectorW> &w, matrix_expression<MatrixZ> &Z, vector_expression<VectorS> &s)
{
    fortran_int_t m_tmp;
    return syevr(Job, All, Uplo, A, 0, 0, 0, 0, abstol, m_tmp, w, Z, s);
}

template <typename SymmA, typename VectorW, typename MatrixZ>
inline int syevr( LAPACK_EVJOB Job, LAPACK_UPLO Uplo, matrix_expression<SymmA> &A,
                  typename SymmA::value_type abstol,
                  vector_expression<VectorW> &w, matrix_expression<MatrixZ> &Z)
{
    fortran_int_t m_tmp, N = A().size1();
    shark::blas::vector<fortran_int_t> s_tmp(N);

    return syevr(Job, All, Uplo, A, 0, 0, 0, 0, abstol, m_tmp, w, Z, s_tmp);
}


template <typename MatrA, typename VectorB>
void syev(
	matrix_expression<MatrA>& matA,
	vector_expression<VectorB>& eigenValues
) {
	SIZE_CHECK(matA().size1() == matA().size2());
	SIZE_CHECK(matA().size1() == eigenValues().size());
	
	std::size_t n = matA().size1();
	LAPACK_UPLO uplo;
	//lapack is column major storage.
	if(boost::is_same<typename MatrA::orientation, blas::row_major>::value){
		uplo = Upper;
	} else {
		uplo = Lower;
	}

	#ifdef EIGENSYMM_USE_SYEVD
		bindings::syevd(Vectors, uplo, matA(), eigenValues());
	#elif defined(EIGENSYMM_USE_SYEVR)
		matrix<value_type, column_major> tmp2(n, n);
		bindings::syevr(Vectors, uplo, matA(), 1e-6, eigenValues(), tmp2);
		matA() = tmp2;
	#else
		bindings::syev(Vectors, uplo, matA(), eigenValues());
	#endif

	matA() = trans(matA);
	
	//reverse eigenvectors and eigenvalues
	for (int i = 0; i < (int)n-i-1; i++)
	{
		int l =  n-i-1;
		std::swap(eigenValues()( l ),eigenValues()( i ));
	}
	for (int j = 0; j < (int)n; j++) {
		for (int i = 0; i < (int)n-i-1; i++)
		{
			int l =  n-i-1;
			std::swap(matA()( j , l ), matA()( j , i ));
		}
	}
}

}}}

#undef SHARK_LAPACK_SSYEV
#undef SHARK_LAPACK_DSYEV
#undef SHARK_LAPACK_SSYEVD
#undef SHARK_LAPACK_DSYEVD
#undef SHARK_LAPACK_SSYEVR
#undef SHARK_LAPACK_DSYEVR

#endif
