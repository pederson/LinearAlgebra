#include "Matrix.hpp"
#include "SparseVector.hpp"
#include "SparseMatrix.hpp"


using namespace std;

// g++ -std=c++11 Matrix.cpp LinearSolvers.hpp -o matrixtest

#include <typeinfo>
#include "LinearSolvers.hpp"

int main(int argc, char * argv[]){

	unsigned int nrows=20, ncols = 5;

	//**************** MATRIX TESTS *********************//
	Matrix A(nrows, ncols);
	for (auto i=0; i<nrows; i++)
	{
		for (auto j=0; j<ncols; j++)
		{
			A(i, j) = i*j;
		}
	}

	// assignment operator
	Matrix B = A;
	cout << "B = " << B << endl;

	A.print_summary();

	// submatrix extraction
	Matrix S = A(4, 8, 3, 4);
	cout << "S = " << S << endl;

	// submatrix assignment
	Matrix C(5, 2);
	C.fill(-10);
	B(0, 4, 2, 3) = C;
	B(15, 19, 0, 1) = S;
	cout << "B = " << B << endl;

	// submatrix assignment via matrix proxy
	// B(0, 1, 0, 1) = S(3, 4, 0, 1);
	// cout << "B = " << B << endl;

	// scalar multiplication, subtraction, addition, division
	Matrix D = (((C*0.5)/2.0) + 10) - 2.5;
	cout << "D = " << D << endl;
	cout << "C = " << C << endl;

	// vector proxy
	cout << "B(:, 3) = " << B.col(3) << endl;
	cout << "B(2, :) = " << B.row(2) << endl;

	// dot products
	cout << "B(:,1)' * B(:,2) = " << Vector_Proxy::dot(B.col(1), B.col(2)) << endl;

	// matrix-matrix mult
	Matrix E = B(0, 2, 0, 4);
	Matrix F = B(14, 18, 2, 4);
	Matrix G = E*F;
	cout << "E = " << E << endl;
	cout << "F = " << F << endl;
	cout << "G = " << G << endl;

	//***************************************************//




	//**************** VECTOR TESTS *********************//
	Vector v1(5);
	for (auto i=0; i<5; i++)
	{
		v1(i) = i;
	}

	// assignment operator
	Vector v2 = v1;
	cout << "v2 = " << v2 << endl;

	v2.print_summary();

	// // subvector extraction
	Vector v3 = v2(1, 2);
	cout << "v3 = " << v3 << endl;

	// // subvector assignment
	Vector v4(2);
	v4.fill(-10);
	v2(0, 1) = v4;
	cout << "v2 = " << v2 << endl;

	// // subvector assignment via vector proxy
	// // B(0, 1, 0, 1) = S(3, 4, 0, 1);
	// // cout << "B = " << B << endl;

	// scalar multiplication, subtraction, addition, division
	Vector v5 = (((v2*0.5)/2.0) + 10) - 2.5;
	cout << "v2 = " << v2 << endl;
	cout << "v5 = " << v5 << endl;

	// vector proxy
	cout << " v2(1:2) = " << v2(1,2) << endl;
	cout << "v2(2, :) = " << v2.row(2) << endl;

	// dot products
	cout << "v2' * v5= " << Vector::dot(v2, v5) << endl;

	// vector-matrix mult
	v5.transpose();
	Vector v6 = v5*F;
	cout << "v5 = " << v5 << endl;
	cout << "F = " << F << endl;
	cout << "v6 = " << v6 << endl;

	// matrix-vector mult
	v6.transpose();
	cout << "F = " << F << endl;
	cout << "v6 = " << v6 << endl;
	Vector v7 = F*v6;
	cout << "v7 = " << v7 << endl;

	// outer product
	Matrix m6 = v3*(v5);
	cout << "v3: " << v3 << endl;
	cout << "v5: " << v5 << endl;
	cout << "outer product: " << m6 << endl;


	//***************************************************//




	//************* DIRECT LINEAR SOLVER TESTS ***********//

	// generate random matrix
	Matrix m = randmat(5,3);
	cout << "m = " << m << endl;

	// generate random normal vector
	Vector v = randvecn(3);
	cout << "v = " << v << endl;

	// generate identity matrix
	Matrix I = eye(3,3);
	cout << "I = " << I*5.0 << endl;

	// orthogonalize matrix
	Matrix Q, Qt, R;
	qr_gram_schmidt(m,Q,R);
	cout << "************************** CLASSICAL GRAM-SCHMIDT:" << endl;
	cout << "Q: " << Q << endl;
	cout << "R: " << R << endl;
	cout << "Q'*Q : " << ~Q*Q << endl;
	cout << "Q*R : " << Q*R << endl;

	// check the solution to a problem
	Vector b = m*v;
	Vector y = unitary_solve(Q,b);
	cout << "y: " << y << endl;
	Vector x = upper_triangular_solve(R, y);
	cout << "************************** UPPER TRIANGULAR:" << endl;
	cout << "actual solution: " << v << endl;
	cout << "computed soution: " << x << endl;

	// check the lower-triangular solver
	Vector b2 = ~R*v;
	Vector x2 = lower_triangular_solve(~R, b2);
	cout << "************************** LOWER TRIANGULAR:" << endl;
	cout << "actual solution: " << v << endl;
	cout << "computed solution: " << x2 << endl;

	// check the diagonal solver
	Vector b3 = I*v;
	Vector x3 = diagonal_solve(I, b3);
	cout << "************************** DIAGONAL:" << endl;
	cout << "actual solution: " << v << endl;
	cout << "computed solution: " << x3 << endl;

	// check householder
	Matrix Uh, Rh, Qh;
	qr_householder(m, Uh, Rh, Qh);
	cout << "************************** HOUSEHOLDER:" << endl;
	cout << "U: " << Uh << endl;
	cout << "Q: " << Qh << endl;
	cout << "R: " << Rh << endl;
	cout << "Q'*Q : " << ~Qh*Qh << endl;
	cout << "Q*R : " << Qh*Rh << endl;

	// check lu decomp
	Matrix L, U;
	Matrix newmat = randmatn(5,5);
	lu(newmat, L, U);
	cout << "************************** LU DECOMP:" << endl;
	cout << "L: " << L << endl;
	cout << "U: " << U << endl;
	cout << "L*U : " << L*U << endl;
	cout << "newmat : " << newmat << endl;

	// check cholesky decomp
	newmat = hilb(5);
	cholesky(newmat, L);
	cout << "************************** CHOLESKY DECOMP:" << endl;
	cout << "A: " << newmat << endl;
	cout << "L: " << L << endl;
	cout << "L*L' : " << L*~L << endl;

	// check randomized basis method
	Matrix Qr;
	unsigned int targetrank = 4;
	rand_basis(newmat, Qr, targetrank);
	cout << "************************** RANDOMIZED BASIS:" << endl;
	cout << "A: " << newmat << endl;
	cout << "Q: " << Qr << endl;
	cout << "Q'*Q : " << ~Qr*Qr << endl;

	// check hessenberg form
	Matrix hess;
	hessenberg(newmat, hess);
	cout << "************************** HESSENBERG FORM:" << endl;
	cout << "hess: " << hess << endl;

	// test QR algorithm once
	Matrix Eg = hilb(4);
	cout << "hilb(4): " << Eg << endl;
	Matrix T;
	hessenberg(Eg, T);
	cout << "hess: " << T << endl;
	Matrix Tnew;
	qr_alg_tridiag_unshifted(T, Tnew);
	cout << "************************** QR ALGORITHM:" << endl;
	cout << "Tnew: " << Tnew << endl;


	// test QR double shift algorithm once
	Matrix Eg2 = randmatn(8,8);
	// Eg2(0,0) = 7; Eg2(0,1) = 3; Eg2(0,2) = 4; Eg2(0,3) = -11; Eg2(0,4) = -9; Eg2(0,5) = -2;
	// Eg2(1,0) = -6; Eg2(1,1) = 4; Eg2(1,2) = -5; Eg2(1,3) = 7; Eg2(1,4) = 1; Eg2(1,5) = 12;
	// Eg2(2,0) = -1; Eg2(2,1) = -9; Eg2(2,2) = 2; Eg2(2,3) = 2; Eg2(2,4) = 9; Eg2(2,5) = 1;
	// Eg2(3,0) = -8; Eg2(3,1) = 0; Eg2(3,2) = -1; Eg2(3,3) = 5; Eg2(3,4) = 0; Eg2(3,5) = 8;
	// Eg2(4,0) = -4; Eg2(4,1) = 3; Eg2(4,2) = -5; Eg2(4,3) = 7; Eg2(4,4) = 2; Eg2(4,5) = 10;
	// Eg2(5,0) = 6; Eg2(5,1) = 1; Eg2(5,2) = 4; Eg2(5,3) = -11; Eg2(5,4) = -7; Eg2(5,5) = -1;
	Matrix T2, Tnew2;
	hessenberg(Eg2, T2);
	qr_alg_double_shifted(T2, Tnew2);
	complex<double> eig1, eig2;
	eig2x2(Tnew2(Tnew2.rows()-2, Tnew2.rows()-1, Tnew2.cols()-2, Tnew2.cols()-1), eig1, eig2);
	cout << "************************** QR DOUBLE SHIFT ALGORITHM:" << endl;
	cout << "A: " << Eg2 << endl;
	cout << "Tnew: " << Tnew2 << endl;
	cout << "subTnew: " << Tnew2(0,1,0,1) << endl;
	cout << "eigs: " << eig1 << ", " << eig2 << endl;

	// check symmetric eigenvalue decomp
	Matrix eigs;
	eig_symm(Eg, eigs);
	cout << "************************** EIGENVALUES:" << endl;
	cout << "eigs: " << eigs << endl;


	// check eigenvalue decomp
	Matrix Teigsc;
	vector<complex<double>> eigsc;
	eig(Eg2, Teigsc, eigsc);
	cout << "************************** COMPLEX EIGENVALUES:" << endl;
	cout << "eigs: " << endl;
	for (auto i=0; i<eigsc.size(); i++) cout << eigsc[i] << ", " ;
	cout << '\n' << endl;


	
	// check singular values of a matrix
	Matrix Eg3 = hilb(6);
	Matrix Sing;
	singular_values_sq(Eg3, Sing);
	cout << "************************** SINGULAR VALUES SQUARED:" << endl;
	cout << "S: " << Sing << endl;
	cout << endl;



	// invert a square matrix
	Matrix Einv;
	invert(Eg2, Einv);
	cout << "************************** MATRIX INVERSE:" << endl;
	cout << "M*M_inv: " << Eg2*Einv << endl;
	cout << "M_inv*M: " << Einv*Eg2 << endl;
	cout << endl;

	//***************************************************//



	//************* ITERATIVE LINEAR SOLVER TESTS ***********//

	// jacobi
	cout << "************************** JACOBI ITERATION:" << endl;
	Matrix ddom = hilb(6);
	Vector dg(6); dg.fill(10);
	Matrix dgm = diag(dg);
	ddom = ddom + dgm;
	Vector rndx1 = randvecn(6);
	Vector ddomb = ddom*rndx1;
	Vector ddomx;
	jacobi(ddom, ddomb, ddomx, 100);
	cout << "A: " << ddom << endl;
	cout << "b: " << ddomb << endl;
	cout << "x_exact: " << rndx1 << endl;
	cout << "x_calc : " << ddomx << endl;
	cout << "error: " << (rndx1-ddomx).norm() << endl;

	// gauss-seidel
	cout << "************************** GAUSS-SEIDEL ITERATION:" << endl;
	Vector ddomx1;
	gauss_seidel(ddom, ddomb, ddomx1, 100);
	cout << "A: " << ddom << endl;
	cout << "b: " << ddomb << endl;
	cout << "x_exact: " << rndx1 << endl;
	cout << "x_calc : " << ddomx1 << endl;
	cout << "error: " << (rndx1-ddomx1).norm() << endl;

	// successive over-relaxation
	cout << "************************** SOR ITERATION:" << endl;
	Vector ddomx2;
	sor(ddom, ddomb, ddomx2, 1.5, 100);
	cout << "A: " << ddom << endl;
	cout << "b: " << ddomb << endl;
	cout << "x_exact: " << rndx1 << endl;
	cout << "x_calc : " << ddomx2 << endl;
	cout << "error: " << (rndx1-ddomx2).norm() << endl;

	// steepest descent
	cout << "************************** STEEPEST DESCENT:" << endl;
	Matrix spd = hilb(6);
	Vector rndx = randvecn(6);
	Vector solnb = spd*rndx;
	Vector solncalc(6);
	solncalc.fill(0);
	steepest_descent(spd, solnb, solncalc);
	cout << "A: " << spd << endl;
	cout << "b: " << solnb << endl;
	cout << "x_exact: " << rndx << endl;
	cout << "x_calc : " << solncalc << endl;
	cout << "error: " << (rndx-solncalc).norm() << endl;

	// conjugate gradient
	cout << "************************** CONJUGATE GRADIENT:" << endl;
	Vector solncalc2(6);
	solncalc2.fill(0);
	conjugate_gradient(spd, solnb, solncalc2);
	cout << "A: " << spd << endl;
	cout << "b: " << solnb << endl;
	cout << "x_exact: " << rndx << endl;
	cout << "x_calc : " << solncalc2 << endl;
	cout << "error: " << (rndx-solncalc2).norm() << endl;

	// conjugate resid
	cout << "************************** CONJUGATE RESIDUAL:" << endl;
	Vector solncalc7(6);
	solncalc7.fill(0);
	conjugate_residual(spd, solnb, solncalc7, 100);
	cout << "A: " << spd << endl;
	cout << "b: " << solnb << endl;
	cout << "x_exact: " << rndx << endl;
	cout << "x_calc : " << solncalc7 << endl;
	cout << "error: " << (rndx-solncalc7).norm() << endl;

	cout << "************************** BICG:" << endl;
	Vector solncalc4(6);
	solncalc4.fill(0);
	bicg(spd, solnb, solncalc4, 100);
	cout << "A: " << spd << endl;
	cout << "b: " << solnb << endl;
	cout << "x_exact: " << rndx << endl;
	cout << "x_calc : " << solncalc4 << endl;
	cout << "error: " << (rndx-solncalc4).norm() << endl;

	cout << "************************** BICR:" << endl;
	Vector solncalc8(6);
	solncalc8.fill(0);
	bicg(spd, solnb, solncalc8, 100);
	cout << "A: " << spd << endl;
	cout << "b: " << solnb << endl;
	cout << "x_exact: " << rndx << endl;
	cout << "x_calc : " << solncalc8 << endl;
	cout << "error: " << (rndx-solncalc8).norm() << endl;


	cout << "************************** BICGSTAB:" << endl;
	Vector solncalc3(6);
	solncalc3.fill(0);
	bicgstab(spd, solnb, solncalc3, 100);
	cout << "A: " << spd << endl;
	cout << "b: " << solnb << endl;
	cout << "x_exact: " << rndx << endl;
	cout << "x_calc : " << solncalc3 << endl;
	cout << "error: " << (rndx-solncalc3).norm() << endl;


	cout << "************************** GMRES(k):" << endl;
	Vector solncalc5(6);
	solncalc5.fill(0);
	gmres_k(spd, solnb, solncalc5, 10, 100);
	cout << "A: " << spd << endl;
	cout << "b: " << solnb << endl;
	cout << "x_exact: " << rndx << endl;
	cout << "x_calc : " << solncalc5 << endl;
	cout << "error: " << (rndx-solncalc5).norm() << endl;


	//***************************************************//



	//*************** SPARSE VECTOR TESTS ***************//
	cout << "************************** SPARSE VECTORS:" << endl;
	SparseVector sv(20, 0);
	sv(10) = 3.0;
	sv(19) = -10.0;
	sv(0) = 111.0;
	sv(1) = 1.0;
	cout << "\nsv: " << endl;
	cout << "length: " << sv.length() << endl;
	cout << "nnz: " << sv.support_size() << endl;
	cout << "norm: " << sv.norm() << endl;
	cout << sv << endl;

	SparseVector sv2(20, -1);
	sv2(0) = 111.0;
	sv2(1) = 3.0;
	sv2(11) =  3.0;
	sv2(18) =  -10.0;
	cout << "sv2: " << sv2 << endl;

	// sparse addition
	SparseVector sv3 = sv2 + sv;
	cout << "sv2+sv: " << sv3 << endl;
	
	// sparse subtraction
	SparseVector sv4 = sv - sv2;
	cout << "sv-sv2: " << sv4 << endl;

	// scalar multiplication, subtraction, addition, division
	SparseVector sv5 = (((sv2*0.5)/2.0) + 10) - 2.5;
	cout << "(((sv2*0.5)/2.0) + 10) - 2.5 : " << sv5 << endl;

	// dot product
	double sparse_dot = SparseVector::dot(sv,sv2);
	cout << "sv'*sv2 : " << sparse_dot << endl;

	// verify with dense dot product
	Vector dsv = sv.densify();
	Vector dsv2 = sv2.densify();
	double dense_dot = Vector::dot(dsv, dsv2);
	cout << "dense sv'*sv2 : " << dense_dot << endl;





	//*************** SPARSE MATRIX TESTS ***************//
	cout << "************************** SPARSE MATRICES:" << endl;
	SparseMatrix sm(5,5);
	sm.set(0,0, 3.0);
	sm.set(1,3, -10.0);
	sm.set(2,2, 111.0);
	sm.set(3,3, 1.0);
	sm.set(1,1, 2.0);
	cout << "\nsm: " << endl;
	cout << "nnz: " << sm.nnz() << endl;
	cout << sm << endl;

	
	// matrix transpose
	SparseMatrix sm2(5,5);
	sm2.set(1,1, 1.0);
	sm2.set(2,2, 1.0);
	sm2.set(3,3, 1.0);
	sm2.set(4,4, 1.0);
	sm2.set(0,2, 1.0);
	sm2.set(0,4, 1.0);
	cout << "sm2: " << sm2 << endl;

	// sparse matrix diag
	SparseMatrix sm1 = sprandmatn(6,6);
	cout << "sm1: " << sm1 << endl;
	Vector sm1d = diag(sm1);
	cout << "sm1 diag: " << sm1d << endl;

	// create SparseMatrix from diag
	SparseMatrix sm1diag = spdiag(sm1d);
	cout << "sm1d into sparse matrix: " << sm1diag << endl;

	// sparse matrix trace
	cout << "trace(sm1diag): " << trace(sm1diag) << "\n" << endl;

	/*
	// matrix transpose
	cout << "sm2: " << sm2 << endl;
	cout << "transposed: " << endl;
	sm2.transpose();
	cout << sm2.rows() << ", " << sm2.cols() << endl;
	cout << sm2.nnz() << endl;
	cout << sm2 << endl;
	*/

	// sparsematrix-vector product
	Vector dv2(5);
	dv2.fill(3);
	Vector dv3 = sm*dv2;
	cout << "SparseMatrix-Vector product: " << dv3 << endl;

	// sparse addition
	SparseMatrix sm3 = sm2 + sm;
	cout << "sm2+sm: " << sm3 << endl;
	
	// sparse subtraction
	cout << "sm2: " << sm2 << endl;
	SparseMatrix sm4 = sm2 - sm;
	cout << "sm2-sm: " << sm4 << endl;

	// scalar multiplication, subtraction, addition, division
	SparseMatrix sm5 = (((sm2*0.5)/2.0) + 10) - 2.5;
	cout << "(((sm2*0.5)/2.0) + 10) - 2.5 : " << sm5 << endl;

	// verify dense matrix conversion
	Matrix dm = sm.densify();
	cout << "dense sm : " << dm << endl;

	// solve upper triangular system
	SparseMatrix sput = speye(6,6) + strictly_upper(sm1);
	Vector spsoln = randvecn(6);
	Vector spb = sput*spsoln;
	Vector spsolncalc = upper_triangular_solve(sput, spb);
	cout << "sparse A: " << sput << endl;
	cout << "b: " << spb << endl;
	cout << "x_exact: " << spsoln << endl;
	cout << "x_calc : " << spsolncalc << endl;
	cout << "error: " << (spsoln-spsolncalc).norm() << endl;

	// solve lower triangular system
	SparseMatrix splt = speye(6,6) + strictly_lower(sm1);
	spb = splt*spsoln;
	Vector spsolncalclt = lower_triangular_solve(splt, spb);
	cout << "sparse A: " << splt << endl;
	cout << "b: " << spb << endl;
	cout << "x_exact: " << spsoln << endl;
	cout << "x_calc : " << spsolncalclt << endl;
	cout << "error: " << (spsoln-spsolncalclt).norm() << endl;

	// sparse jacobi iteration
	cout << "************************** SPARSE JACOBI ITERATION:" << endl;
	SparseMatrix spjac = sm1 + 10*speye(6,6);
	spb = spjac*spsoln;
	Vector spsolncalc_jac(6); spsolncalc_jac.fill(0);
	jacobi(spjac, spb, spsolncalc_jac, 100);
	cout << "sparse A: " << spjac << endl;
	cout << "b: " << spb << endl;
	cout << "x_exact: " << spsoln << endl;
	cout << "x_calc : " << spsolncalc_jac << endl;
	cout << "error: " << (spsoln-spsolncalc_jac).norm() << endl;

	// sparse gauss seidel iteration
	cout << "************************** SPARSE GAUSS-SEIDEL ITERATION:" << endl;
	SparseMatrix spgs = sm1 + 10*speye(6,6);
	spb = spgs*spsoln;
	Vector spsolncalc_gs(6); spsolncalc_gs.fill(0);
	gauss_seidel(spgs, spb, spsolncalc_gs, 100);
	cout << "sparse A: " << spgs << endl;
	cout << "b: " << spb << endl;
	cout << "x_exact: " << spsoln << endl;
	cout << "x_calc : " << spsolncalc_gs << endl;
	cout << "error: " << (spsoln-spsolncalc_gs).norm() << endl;

	// sparse SOR iteration
	cout << "************************** SPARSE SOR ITERATION:" << endl;
	SparseMatrix spsor = sm1 + 10*speye(6,6);
	spb = spsor*spsoln;
	Vector spsolncalc_sor(6); spsolncalc_sor.fill(0);
	sor(spsor, spb, spsolncalc_sor, 1.5, 100);
	cout << "sparse A: " << spsor << endl;
	cout << "b: " << spb << endl;
	cout << "x_exact: " << spsoln << endl;
	cout << "x_calc : " << spsolncalc_sor << endl;
	cout << "error: " << (spsoln-spsolncalc_sor).norm() << endl;

	// sparse bicgstab iteration
	cout << "************************** SPARSE BICGSTAB:" << endl;
	SparseMatrix spbicgstab = sm1 + 10*speye(6,6);
	spb = spsor*spsoln;
	Vector spsolncalc_bicgstab(6); spsolncalc_bicgstab.fill(0);
	bicgstab(spbicgstab, spb, spsolncalc_bicgstab, 100);
	cout << "sparse A: " << spsor << endl;
	cout << "b: " << spb << endl;
	cout << "x_exact: " << spsoln << endl;
	cout << "x_calc : " << spsolncalc_bicgstab << endl;
	cout << "error: " << (spsoln-spsolncalc_bicgstab).norm() << endl;


	return 0;
}