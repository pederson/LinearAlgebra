#ifndef _COMPOSITEOPERATORS_H
#define _COMPOSITEOPERATORS_H



namespace libra{

// the TimeStepOperator implements the equivalent operation T(x)
// of acting with operator T on vector x. 
// The TimeStepper template paramter is an object that
// is required to have functions:
// 					- solution(), which returns the solution vector
//                  - time_step(double deltat), which advances the solution vector by delta_t in time
template <typename TimeStepper>
struct TimeStepOperator{
private:
	TimeStepper * mTOp;
	double mdt;


public:
	TimeStepOperator(TimeStepper & top, double deltat)
	: mTOp(&top), mdt(deltat) {};

	double & dt() {return mdt;};
	const double & dt() const {return mdt;};

	template <typename VectorType, typename VectorTypeOut>
	void vmult(const VectorType & v, VectorTypeOut && vout) const {
		mTOp->solution() = v;
		mTOp->time_step(mdt);
		vout = mTOp->solution();
	}
};

template <typename TimeStepper>
TimeStepOperator<TimeStepper> create_timestep_operator(TimeStepper & t, double dt){
	return TimeStepOperator<TimeStepper>(t, dt);
}



//********************************************************//
//********************************************************//



// the Fourier Operator implements a fourier operator
// (F-exp(-i*arg)) to solve periodic forcing problems 
// of the form (F-exp(-i*arg)) = s
template <typename Operator>
struct FourierOperator{
private:
	Operator * mTOp;
	double marg;
	std::complex<double> ceye = std::complex<double>(0, 1);

public:
	FourierOperator(Operator & top, double argg)
	: mTOp(&top), marg(argg) {};

	template <typename VectorType, typename VectorTypeOut>
	void vmult(const VectorType & v, VectorTypeOut && vout) const {
		mTOp->vmult(v, vout);
		vout = vout - exp(-ceye*marg)*v;
	}
};

template <typename Operator>
FourierOperator<Operator> create_fourier_operator(Operator & t, double f){
	return FourierOperator<Operator>(t, f);
}



//********************************************************//
//********************************************************//


// the RQIOperator implements a
// Rayleigh Quotient Iteration operator (A-lambda*I), which is
// used to solve the eigenvalue problem A*x = lambda*x
template <typename Operator, typename EigType>
struct RQIOperator{
private:
	Operator * mTOp;
	EigType mlambda;

public:
	RQIOperator(Operator & top, EigType lam)
	: mTOp(&top), mlambda(lam) {};

	EigType & lambda() {return mlambda;};
	const EigType & lambda() const {return mlambda;};


	template <typename VectorType, typename VectorTypeOut>
	void vmult(const VectorType & v, VectorTypeOut && vout) const {
		mTOp->vmult(v, vout);
		vout = vout - mlambda*v;
	}
};

template <typename Operator, typename EigType>
RQIOperator<Operator, EigType> create_rqi_operator(Operator & t, EigType f){
	return RQIOperator<Operator, EigType>(t, f);
}


//********************************************************//
//********************************************************//


// the GeneralizedRQIOperator implements a generalized
// Rayleigh Quotient Iteration operator (A-B*lambda), which is
// used to solve the generalized eigenvalue problem A*x = lambda*B*x
template <typename OperatorA, typename OperatorB, typename EigType>
struct GeneralizedRQIOperator{
private:
	OperatorA * mAop;
	OperatorB * mBop;
	EigType mlambda;


public:
	GeneralizedRQIOperator(OperatorA & a, OperatorB & b, EigType lam)
	: mAop(&a), mBop(&b), mlambda(lam) {};

	EigType & lambda() {return mlambda;};
	const EigType & lambda() const {return mlambda;};

	template <typename VectorType, typename VectorTypeOut>
	void vmult(const VectorType & v, VectorTypeOut && vout) const {

		libra::Vector<double, libra::dynamic_size> hold(v.size());
		mBop->vmult(v, hold);
		mAop->vmult(v, vout);
		vout = vout - mlambda*hold;
	}
};

template <typename OperatorA, typename OperatorB, typename EigType>
GeneralizedRQIOperator<OperatorA, OperatorB, EigType> create_grqi_operator(OperatorA & a, OperatorB & b, EigType & lam){
	return GeneralizedRQIOperator<OperatorA, OperatorB, EigType>(a, b, lam);
}


//********************************************************//
//********************************************************//


// the DeflatedOperator creates an explicit rank-1 
// modification to the operator A with known eigenvalue lambda
// and eigenvector u, with a general vector v as (A-lambda*u*v').
// If the vector v is chosen to be the right eigenvector, it 
// is the Hotelling deflation. It's also possible to choose
// for v the left eigenvector
template <typename Operator, typename EigVecType, typename EigValType>
struct DeflatedOperator{
private:
	Operator * mTOp;
	EigValType mlambda;
	EigVecType * meigvec;

public:
	DeflatedOperator(Operator & top, EigVecType & evec, EigValType lam)
	: mTOp(&top), meigvec(&evec), mlambda(lam) {};

	template <typename VectorType, typename VectorTypeOut>
	void vmult(const VectorType & v, VectorTypeOut && vout) const {
		mTOp->vmult(v, vout);
		vout = vout - mlambda*vector::inner_product(*meigvec, v)*v;
	}
};

template <typename Operator, typename EigVecType, typename EigValType>
DeflatedOperator<Operator, EigVecType, EigValType> create_hotelling_deflated_operator(Operator & t, EigVecType & v, EigValType f){
	return DeflatedOperator<Operator, EigVecType, EigValType>(t, v, f);
}


} // end namespace libra

#endif