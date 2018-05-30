#include <libra.h>

using namespace std;
// using namespace libra;
// g++ -std=c++14 -I./ libra_test.cpp -o matrixtest

#include <typeinfo>

int main(int argc, char * argv[]){


	//**************** GENERALIZED VECTOR TESTS *********************//
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	cout << "******* GENERALIZED VECTOR TESTS ***********" << endl;
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	std::vector<double> ga1 = {0,1,2,3,4,5,6,7,8,9};
	std::vector<std::complex<double>> ga2 = {9,8,7,6,5,4,3,2,1,0};
	libra::write_vector(ga2);
	std::cout << "Generalized Inner Product: " << libra::inner_product(ga2, ga2) << std::endl;
	std::cout << "Generalized Norms: " << std::endl;
	std::cout << "                  infinity: " << libra::norm_inf(ga2) << std::endl;
	std::cout << "                       one: " << libra::norm_1(ga2) << std::endl;
	std::cout << "                       two: " << libra::norm_2(ga2) << std::endl;
	std::cout << "                     three: " << libra::norm_3(ga2) << std::endl;
	std::cout << "                      four: " << libra::norm_4(ga2) << std::endl;
	std::cout << "                      five: " << libra::norm_5(ga2) << std::endl;
	libra::fill(ga2, 1.0);
	libra::write_vector(ga2);
	std::cout << "Vector length: " << libra::length(ga2) << std::endl;
	auto subv = libra::vector_view(ga2, ga2.begin()+4, ga2.end());
	libra::fill(subv, 0.0);
	auto subsubv = libra::vector_view(subv, subv.begin()+1, subv.end()-1);
	libra::fill(subsubv, 2.0);
	std::cout <<"vector_view length: " << libra::length(subsubv) << std::endl;
	libra::write_vector(subv);
	std::cout << "vector_view norm_3: " << libra::norm_3(subsubv) << std::endl;
	libra::write_vector(ga2);
	auto subv1 = libra::vector_view(ga1, ga1.begin()+3, ga1.end()-2);
	libra::fill_randn(subv1);
	libra::write_vector(ga1);
	std::cout << "max value: " << libra::max(ga1) << " at position: " << libra::argmax(ga1) << std::endl;
	std::cout << "min value: " << libra::min(ga1) << " at position: " << libra::argmin(ga1) << std::endl;
	

	//**************** TENSOR TESTS *********************//
	cout << endl;
	cout << endl;
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	cout << "************ DENSE TENSOR TESTS ************" << endl;
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;

	libra::GenericTensor<int, libra::LinearMemory, 3, 2, libra::dynamic_size, libra::dynamic_size> t(2, 3);
	cout << "Tensor has rank: " << t.rank() << endl;
	cout << "Tensor has ndynamic: " << t.ndynamic() << endl;
	cout << "dimensions: " ;
	auto tdims = t.dims();
	for (auto i=tdims.begin(); i!=tdims.end(); i++) cout << *i << ", " ;
	cout << endl;

	// for (auto it=t.cbegin(); it!=t.cend(); it++) *it = 5;

	// for (auto d=0; d<t.rank(); d++){

	// }
	// for (auto i=tdims.begin(); i!=tdims.end(); i++) cout << *i << ", " ;

	std::cout << "size(3): " << t.size<3>() << std::endl;

	// t(1,2,2,2) = 12;
	cout << "t(1,2,2,2): " << t(0,0,0,0) << std::endl;

	// libra::Tensor<int, 3, 2, 2, 2> tens = {1,2,3,4,5,6,7,8};
	// cout << "<Tensor>" << endl;
	// for (auto it1 = tens.begin(); it1!=tens.end(); it1++){
	// 	for (auto it = it1->begin(); it!= it1->end(); it++){
	// 		libra::write_vector(*it);
	// 	}
	// }
	// cout << "</Tensor>" << endl;

	// throw -1;

	//**************** VECTOR TESTS *********************//
	cout << endl;
	cout << endl;
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	cout << "************ DENSE VECTOR TESTS ************" << endl;
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	libra::Vector<double, libra::dynamic_size> dmatx = {1.0,0.0,0.0};
	libra::Vector<double, 3> dresult = {1.0,2.0,3.0};
	dresult = dmatx;
	libra::write_vector(dresult);
	libra::Vector<double, 3> dresult_copy = dresult;
	libra::write_vector(dresult_copy);
	
	libra::Vector<int, libra::dynamic_size> dvec = {3,2,9,0,20, 9, 1};
	libra::write_vector(dvec);


	//**************** MATRIX TESTS *********************//
	cout << endl;
	cout << endl;
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	cout << "************ DENSE MATRIX TESTS ************" << endl;
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	libra::Matrix<int, 3, 3> dmat = {{1,2,3},{4,5,6},{7,8,9}}; dmat[0][0] = 9;
	// for (auto it = dmat.begin(); it!= dmat.end(); it++) libra::fill_randn(*it);
	cout << "<Matrix>" << endl;
	for (auto it = dmat.begin(); it!= dmat.end(); it++){
		libra::write_vector(*it);
	}
	cout << "</Matrix>" << endl;
	dmat.vmult(dmatx, dresult);
	libra::write_vector(dresult);
	
	// libra::



	unsigned int nrows=20, ncols = 5;

	//**************** MATRIX TESTS *********************//
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	cout << "************ DENSE MATRIX TESTS ************" << endl;
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
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

	// writing to file
	S.dlmwrite("M_dense.txt",":");

	//***************************************************//


	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	cout << "**** GENERALIZED LINEAR SOLVER TESTS *******" << endl;
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	// Vector lvec = libra::diag(A);
	// cout << "lvec: " << lvec << std::endl;

	typedef std::complex<double> SolveType;
	libra::Matrix<SolveType, 4, 4> lmat;// = {1.0, 0.1, 0.1, 0.2, 1.0, 0.2, 0.3, 0.3, 1.0};
	for (auto it=lmat.begin(); it!=lmat.end(); it++) libra::fill_rand(*it);
	// lmat(0,0) = lmat(1,1) = lmat(2,2) = 1.0;
	// libra::write_matrix(lmat);
	libra::Vector<SolveType, libra::dynamic_size> lvec(4, 0.0);//; lvec.resize(4);// = {0.2, 0.3, 0.2};
	libra::fill_rand(lvec);
	libra::write_vector(lvec);
	
	libra::Vector<SolveType, libra::dynamic_size> lresult; lresult.resize(4);// = {0,0,0};
	libra::fill(lresult, 0);
	// lmat.vmult(lvec, lresult);
	libra::bicgstab_l(2, lmat, lvec, lresult, 10);
	// solve the system with BiCGSTAB

	cout << "Solution: " << endl;
	// libra::write_vector(lresult);
	libra::Vector<SolveType, libra::dynamic_size> lresult2 = lresult;
	lmat.vmult(lresult, lresult2);
	std::cout << "Reconstructed RHS: " << std::endl;
	// libra::write_vector(lresult2);
	std::cout << "Exact RHS: " << std::endl;
	// libra::write_vector(lvec);
	std::cout << "resulting norm: " << norm_2(lresult2-lvec) << std::endl;


	typedef std::complex<double> Ctype;
	if (libra::is_complex<Ctype>::value){
		cout << "im complex" << endl;
		typedef typename libra::complex_detector<Ctype>::type 	MyType;
		MyType p;
		cout << "my underlying type is: " << typeid(p).name() << endl;

	} 
	throw -1;




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
	unsigned int niters;

	// jacobi
	cout << "************************** JACOBI ITERATION:" << endl;
	Matrix ddom = hilb(6);
	Vector dg(6); dg.fill(10);
	Matrix dgm = diag(dg);
	ddom = ddom + dgm;
	Vector rndx1 = randvecn(6);
	Vector ddomb = ddom*rndx1;
	Vector ddomx(6); ddomx.fill(0);
	niters = jacobi(ddom, ddomb, ddomx, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "A: " << ddom << endl;
	cout << "b: " << ddomb << endl;
	cout << "x_exact: " << rndx1 << endl;
	cout << "x_calc : " << ddomx << endl;
	cout << "error: " << (rndx1-ddomx).norm() << endl;

	// gauss-seidel
	cout << "************************** GAUSS-SEIDEL ITERATION:" << endl;
	Vector ddomx1(6); ddomx1.fill(0);
	niters = gauss_seidel(ddom, ddomb, ddomx1, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "A: " << ddom << endl;
	cout << "b: " << ddomb << endl;
	cout << "x_exact: " << rndx1 << endl;
	cout << "x_calc : " << ddomx1 << endl;
	cout << "error: " << (rndx1-ddomx1).norm() << endl;

	// successive over-relaxation
	cout << "************************** SOR ITERATION:" << endl;
	Vector ddomx2(6); ddomx2.fill(0);
	niters = sor(ddom, ddomb, ddomx2, 1.5, 100);
	cout << "iterated: " << niters << " times" << endl;
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
	niters = steepest_descent(spd, solnb, solncalc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "A: " << spd << endl;
	cout << "b: " << solnb << endl;
	cout << "x_exact: " << rndx << endl;
	cout << "x_calc : " << solncalc << endl;
	cout << "error: " << (rndx-solncalc).norm() << endl;

	// conjugate gradient
	cout << "************************** CONJUGATE GRADIENT:" << endl;
	Vector solncalc2(6);
	solncalc2.fill(0);
	niters = conjugate_gradient(spd, solnb, solncalc2, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "A: " << spd << endl;
	cout << "b: " << solnb << endl;
	cout << "x_exact: " << rndx << endl;
	cout << "x_calc : " << solncalc2 << endl;
	cout << "error: " << (rndx-solncalc2).norm() << endl;

	// conjugate resid
	cout << "************************** CONJUGATE RESIDUAL:" << endl;
	Vector solncalc7(6);
	solncalc7.fill(0);
	niters = conjugate_residual(spd, solnb, solncalc7, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "A: " << spd << endl;
	cout << "b: " << solnb << endl;
	cout << "x_exact: " << rndx << endl;
	cout << "x_calc : " << solncalc7 << endl;
	cout << "error: " << (rndx-solncalc7).norm() << endl;

	cout << "************************** BICG:" << endl;
	Vector solncalc4(6);
	solncalc4.fill(0);
	niters = bicg(spd, solnb, solncalc4, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "A: " << spd << endl;
	cout << "b: " << solnb << endl;
	cout << "x_exact: " << rndx << endl;
	cout << "x_calc : " << solncalc4 << endl;
	cout << "error: " << (rndx-solncalc4).norm() << endl;

	cout << "************************** BICR:" << endl;
	Vector solncalc8(6);
	solncalc8.fill(0);
	niters = bicg(spd, solnb, solncalc8, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "A: " << spd << endl;
	cout << "b: " << solnb << endl;
	cout << "x_exact: " << rndx << endl;
	cout << "x_calc : " << solncalc8 << endl;
	cout << "error: " << (rndx-solncalc8).norm() << endl;


	cout << "************************** BICGSTAB:" << endl;
	Vector solncalc3(6);
	solncalc3.fill(0);
	niters = bicgstab(spd, solnb, solncalc3, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "A: " << spd << endl;
	cout << "b: " << solnb << endl;
	cout << "x_exact: " << rndx << endl;
	cout << "x_calc : " << solncalc3 << endl;
	cout << "error: " << (rndx-solncalc3).norm() << endl;


	cout << "************************** GMRES(k):" << endl;
	Vector solncalc5(6);
	solncalc5.fill(0);
	niters = gmres_k(spd, solnb, solncalc5, 10, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "A: " << spd << endl;
	cout << "b: " << solnb << endl;
	cout << "x_exact: " << rndx << endl;
	cout << "x_calc : " << solncalc5 << endl;
	cout << "error: " << (solnb-spd*solncalc5).norm() << endl;
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

	// sparse matrix swap
	cout << "Swapping sparse matrices" << endl;
	cout << "sm1: " << sm1 << endl;
	cout << "sm2: " << sm2 << endl;
	swap(sm1, sm2);
	cout << "sm1: " << sm1 << endl;
	cout << "sm2: " << sm2 << endl;
	swap(sm1, sm2);

	// matrix transpose
	cout << "Sparse matrix transpose" << endl;
	cout << "sm2: " << sm2 << endl;
	cout << "transposed: " ;
	sm2.transpose();
	cout << sm2 << endl;
	sm2.transpose();

	// sparsematrix-vector product
	Vector dv2(5);
	dv2.fill(3);
	Vector dv3 = sm*dv2;
	cout << "SparseMatrix-Vector product: " << dv3 << endl;

	// sparsematrix-sparsematrix product
	SparseMatrix ss1 = sprandmatn(100,100);
	SparseMatrix ss2 = sprandmatn(100,100);
	SparseMatrix smsm2 = ss1*ss2;
	Vector ssp(smsm2.cols());
	Vector result1 = ss1*(ss2*ssp);
	Vector result2 = smsm2*ssp;
	cout << "sparsematrix-sparsematrix product: resid: " << norm_2(result2-result1) << endl;
	//cout << smsm2 << endl;

	// sparsematrix-sparsematrix tmult
	SparseMatrix smtsm2 = ss1.Tmult(ss2);
	result1 = ss1.Tmult(ss2*ssp);
	result2 = smtsm2*ssp;
	cout << "sparsematrix-sparsematrix tmult: resid: " << norm_2(result2-result1) << endl;
	//cout << smtsm2 << endl;

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

	// write to matrix-market format
	sm.mmwrite("M_sparse.txt");

	// read from matrix-market format
	SparseMatrix smread = mmread("M_sparse.txt");
	cout << "read sparse matrix: " << smread << endl;

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

	// incomplete cholesky decomp
	cout << "\nincomplete cholesky decomp" << endl;
	SparseMatrix spch = 10*speye(20,20) + sprandmatsymm(20,20, 1.0);
	//cout << spch << endl;
	SparseMatrix Lch;
	icholesky(spch, Lch);
	Vector chv = randvecn(20);
	cout << "resid: " << (spch*chv - Lch*Lch.Tmult(chv)).norm() << "\n" << endl;


	// incomplete LU decomp
	cout << "\nincomplete LU decomp" << endl;
	SparseMatrix splu = 10*speye(20,20) + sprandmatn(20,20, 1.0);
	cout << "fill: " << splu.nnz() << "/" << splu.rows()*splu.cols() << endl;
	//cout << spch << endl;
	SparseMatrix Llu, Ulu;
	ilu(splu, Llu, Ulu);
	cout << "L fill: " << Llu.nnz() << "/" << Llu.rows()*Llu.cols() << endl;
	cout << "U fill: " << Ulu.nnz() << "/" << Ulu.rows()*Ulu.cols() << endl;
	// cout << "L: " << Llu << endl;
	// cout << "U: " << Ulu << endl;
	Vector luv = randvecn(20);
	cout << "resid: " << (splu*luv - Llu*(Ulu*luv)).norm() << "\n" << endl;

	// throw -1;

	// sparse jacobi iteration
	cout << "************************** SPARSE JACOBI ITERATION:" << endl;
	SparseMatrix spjac = sm1 + 10*speye(6,6);
	spb = spjac*spsoln;
	Vector spsolncalc_jac(6); spsolncalc_jac.fill(0);
	niters = jacobi(spjac, spb, spsolncalc_jac, 100);
	cout << "iterated: " << niters << " times" << endl;
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
	niters = gauss_seidel(spgs, spb, spsolncalc_gs, 100);
	cout << "iterated: " << niters << " times" << endl;
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
	niters = sor(spsor, spb, spsolncalc_sor, 1.25, 100);
	cout << "iterated: " << niters << " times" << endl;
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
	niters = bicgstab(spbicgstab, spb, spsolncalc_bicgstab, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "sparse A: " << spsor << endl;
	cout << "b: " << spb << endl;
	cout << "x_exact: " << spsoln << endl;
	cout << "x_calc : " << spsolncalc_bicgstab << endl;
	cout << "error: " << (spsoln-spsolncalc_bicgstab).norm() << endl;

	// sparse gmres iteration
	cout << "************************** SPARSE GMRES:" << endl;
	SparseMatrix spgmres = sm1 + 10*speye(6,6);
	spb = spsor*spsoln;
	Vector spsolncalc_gmres(6); spsolncalc_gmres.fill(0);
	niters = gmres_k(spbicgstab, spb, spsolncalc_gmres, 10, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "sparse A: " << spsor << endl;
	cout << "b: " << spb << endl;
	cout << "x_exact: " << spsoln << endl;
	cout << "x_calc : " << spsolncalc_gmres << endl;
	cout << "error: " << (spsoln-spsolncalc_gmres).norm() << endl;




	cout << "****************** PRECONDITIONED SOLVERS ******************" << endl;
	unsigned int psolvesize = 500;
	unsigned int stencilsize = 5;
	double fill = double(stencilsize)/double(psolvesize);
	SparseMatrix spsymm = sprandmatnsymm(psolvesize,psolvesize, fill) + 10*speye(psolvesize, psolvesize);
	Vector ps_exact = randvecn(psolvesize);
	Vector ps_b = spsymm*ps_exact;
	Vector ps_calc(psolvesize);
	cout << "nnz: " << spsymm.nnz() << "/" << psolvesize*psolvesize << " = " << double(spsymm.nnz())/double(psolvesize*psolvesize) << endl;

	cout << "************************** CG - JACOBI PC:" << endl;
	JacobiPreconditioner jpc(spsymm);
	ps_calc.fill(0);
	niters = conjugate_gradient(spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "cg resid: " << (ps_b - spsymm*ps_calc).norm() << endl;
	ps_calc.fill(0);
	niters = conjugate_gradient(&jpc, spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "pc cg resid: " << (ps_b - spsymm*ps_calc).norm() << endl;


	cout << "************************** CG - GAUSS-SEIDEL PC:" << endl;
	GSPreconditioner gspc(spsymm);
	ps_calc.fill(0);
	niters = conjugate_gradient(spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "cg resid: " << (ps_b - spsymm*ps_calc).norm() << endl;
	ps_calc.fill(0);
	niters = conjugate_gradient(&gspc, spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "pc cg resid: " << (ps_b - spsymm*ps_calc).norm() << endl;


	cout << "************************** CG - SYMMETRIC GAUSS-SEIDEL PC:" << endl;
	SGSPreconditioner sgspc(spsymm);
	ps_calc.fill(0);
	niters = conjugate_gradient(spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "cg resid: " << (ps_b - spsymm*ps_calc).norm() << endl;
	ps_calc.fill(0);
	niters = conjugate_gradient(&sgspc, spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "pc cg resid: " << (ps_b - spsymm*ps_calc).norm() << endl;


	cout << "************************** CG - SOR PC:" << endl;
	SORPreconditioner sorpc(spsymm, 1.5);
	ps_calc.fill(0);
	niters = conjugate_gradient(spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "cg resid: " << (ps_b - spsymm*ps_calc).norm() << endl;
	ps_calc.fill(0);
	niters = conjugate_gradient(&sorpc, spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "pc cg resid: " << (ps_b - spsymm*ps_calc).norm() << endl;


	cout << "************************** CG - SSOR PC:" << endl;
	SSORPreconditioner ssorpc(spsymm, 0.8);
	ps_calc.fill(0);
	niters = conjugate_gradient(spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "cg resid: " << (ps_b - spsymm*ps_calc).norm() << endl;
	ps_calc.fill(0);
	niters = conjugate_gradient(&ssorpc, spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "pc cg resid: " << (ps_b - spsymm*ps_calc).norm() << endl;


	cout << "************************** CG - INCOMPLETE CHOLESKY PC:" << endl;
	ICPreconditioner icpc(spsymm);
	ps_calc.fill(0);
	niters = conjugate_gradient(spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "cg resid: " << (ps_b - spsymm*ps_calc).norm() << endl;
	ps_calc.fill(0);
	niters = conjugate_gradient(&icpc, spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "pc cg resid: " << (ps_b - spsymm*ps_calc).norm() << endl;

	cout << "************************** CG - INCOMPLETE LU PC:" << endl;
	ILUPreconditioner ilupc(spsymm);
	ps_calc.fill(0);
	niters = conjugate_gradient(spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "cg resid: " << (ps_b - spsymm*ps_calc).norm() << endl;
	ps_calc.fill(0);
	niters = conjugate_gradient(&ilupc, spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "pc cg resid: " << (ps_b - spsymm*ps_calc).norm() << endl;

	cout << "************************** CG - AMG PC:" << endl;
	AMGPreconditioner amgpc(spsymm);
	ps_calc.fill(0);
	niters = conjugate_gradient(spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "cg resid: " << (ps_b - spsymm*ps_calc).norm() << endl;
	ps_calc.fill(0);
	niters = conjugate_gradient(&amgpc, spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "pc cg resid: " << norm_2(ps_b - spsymm*ps_calc) << endl;

	/*
	cout << "************************** CR - INCOMPLETE LU PC:" << endl;
	ps_calc.fill(0);
	niters = conjugate_residual(spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "cg resid: " << (ps_b - spsymm*ps_calc).norm() << endl;
	ps_calc.fill(0);
	niters = conjugate_residual(&ilupc, spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "pc cg resid: " << (ps_b - spsymm*ps_calc).norm() << endl;
	*/


	cout << "************************** BICGSTAB - INCOMPLETE LU PC:" << endl;
	ps_calc.fill(0);
	niters = bicgstab(spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "bicgstab resid: " << (ps_b - spsymm*ps_calc).norm() << endl;
	ps_calc.fill(0);
	niters = bicgstab(&ilupc, spsymm, ps_b, ps_calc, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "pc bicgstab resid: " << norm_2(ps_b - spsymm*ps_calc) << endl;

	cout << "************************** GMRES - INCOMPLETE LU PC:" << endl;
	ps_calc.fill(0);
	niters = gmres_k(spsymm, ps_b, ps_calc, 20, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "gmres resid: " << (ps_b - spsymm*ps_calc).norm() << endl;
	ps_calc.fill(0);
	niters = gmres_k(&ilupc, spsymm, ps_b, ps_calc, 20, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "pc gmres resid: " << norm_2(ps_b - spsymm*ps_calc) << endl;

	cout << "************************** GMRES - AMG PC:" << endl;
	ps_calc.fill(0);
	niters = gmres_k(spsymm, ps_b, ps_calc, 20, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "gmres resid: " << (ps_b - spsymm*ps_calc).norm() << endl;
	ps_calc.fill(0);
	niters = gmres_k(&amgpc, spsymm, ps_b, ps_calc, 20, 100);
	cout << "iterated: " << niters << " times" << endl;
	cout << "pc gmres resid: " << norm_2(ps_b - spsymm*ps_calc) << endl;







	cout << "************************** ALGEBRAIC MULTIGRID: " << endl;
	Vector offd(1000 - 1); offd.fill(-1);
	SparseMatrix spamg = 2*speye(1000,1000) + spdiag(offd, 1) + spdiag(offd,-1);
	Vector ps_amg(1000); ps_amg.fill(0);
	// Vector psamg_b(1000); psamg_b.fill(0); psamg_b(50) = 1;
	Vector psamg_b = randvecn(1000);
	cout << "amg resid before: " << (psamg_b - spamg*ps_amg).norm() << endl;
	// niters = bicgstab(spamg, psamg_b, ps_amg, 50);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "amg resid after bicgstab: " << (psamg_b - spamg*ps_amg).norm() << endl;
	niters = amg(spamg, psamg_b, ps_amg, 100, 1.0e-12);
	cout << "iterated: " << niters << " times" << endl;
	cout << "amg resid: " << (psamg_b - spamg*ps_amg).norm() << endl;
	ps_amg.dlmwrite("amg_soln.txt");

	return 0;
}