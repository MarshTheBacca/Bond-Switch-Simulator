#include "opt.h"
#include <omp.h>
//##### STEEPEST DESCENT #####

//Default constructor
template <typename M>
SteepestDescentMultiDim<M>::SteepestDescentMultiDim() {
    itMax=1000;
    ls=1e-3;
    tol=1e-6;
}

//Constructor with algorithm search parameters
template <typename M>
SteepestDescentMultiDim<M>::SteepestDescentMultiDim(int iterationLimit, double lineSearchInc, double convergenceTolerance) {
    itMax=iterationLimit;
    ls=lineSearchInc;
    tol=convergenceTolerance;
}

//Optimisation Algorithm
template <typename M>
VecF<int> SteepestDescentMultiDim<M>::operator()(M& model, VecF<double> &x) {
    cout << "Steepest Descent" << endl;
    //Initialise optimisation variables
    int it=0; //number of iterations
    double f0, f1, f; //previous, current and generic function value
    double df; //difference between current and previous function evaluations
    VecF<double> g(x); //derivative of function
    VecF<double> d(x); //search direction
    VecF<int> status(2); //optimisation status: convergence status, iterations
    f0=numeric_limits<double>::infinity();

    //Evaluate gradient and check non-zero before commencing descent
    g=model.gradient(x);
    if(vAsum(g)<tol){
        status[0]=1;
        status[1]=0;
        return status;
    }

    //Steepest descent algorithm
    for(int i=0; i<itMax; ++i){
//        cout << "Iteration : " << i << endl;
//        cout << "Line search : " << ls << endl;
//        cout << "--------------------------------" << endl;
//        cout << "G coordinate : " << endl;
//        cout << g[0] << endl;
//        cout << "--------------------------------" << endl;
        //line search
        d=-g*ls;
        f1=model.function(x);
        //#pragma omp parallel for private(x,f) shared(d, f1)
        //{
            for (;;) {
                x += d;
//            cout << "----------------------                nodes[i].netCnxs.addValue(nodesA[i].netCnxs[j]);----------" << endl;
//            cout << "X coordinate : " << endl;
//            cout << x[0] << endl;
//            cout << "--------------------------------" << endl;
                f = model.function(x);
                if (f > f1) {
                    //passed through minimum, backtrack
                    x -= d;
                    ++it;
                    break;
                } else f1 = f;
            }
        //}

        //check convergence
        df=abs(f1-f0);
        if(abs(f1)<tol || df<tol) break;
        else f0=f1;
        //recalculate gradient
        g=model.gradient(x);
    }

    if(it==itMax-1) status[0]=2;
    else status[0]=0;
    status[1]=it;

    return status;
}


//##### STEEPEST DESCENT WITH ARMIJO BACKTRACKING LINE SEARCH #####

//Default constructor
template <typename M>
SteepestDescentArmijoMultiDim<M>::SteepestDescentArmijoMultiDim() {
    itMax=1000;
    tau=0.5;
    tol=1e-6;
}

//Constructor with algorithm search parameters
template <typename M>
SteepestDescentArmijoMultiDim<M>::SteepestDescentArmijoMultiDim(int iterationLimit, double lineSearchInc, double convergenceTolerance) {
    itMax=iterationLimit;
    tau=lineSearchInc;
    tol=convergenceTolerance;
}

//Optimisation Algorithm
template <typename M>
VecF<int> SteepestDescentArmijoMultiDim<M>::operator()(M& model, VecF<double> &x) {

    //Initialise optimisation variables
    int it=0; //number of iterations
    double f0, f1, f; //previous, current and generic function value
    double df; //difference between current and previous function evaluations
    double alpha; //overall line search step size
    double gSq; //modulus square of gradient
    VecF<double> g(x); //derivative of function
    VecF<double> d(x); //coordinates after search direction
    VecF<int> status(2); //optimisation status: convergence status, iterations
    f0=numeric_limits<double>::infinity();

    //Evaluate gradient and check non-zero before commencing descent
    g=model.gradient(x);
    if(vAsum(g)<tol){
        status[0]=1;
        status[1]=0;
        return status;
    }

    //Steepest descent algorithm
    for(int i=0; i<itMax; ++i){
//        cout << "Iteration : " << i << endl;

        //backtracking line search
        alpha=1.0;
        gSq=vNormSq(g);
        f=model.function(x);
        for(;;){
            //cout << "    step = " << alpha << endl;
            d=x-g*alpha;
            //cout << "Parameters : " << d << " " << x << " " << g << " " << alpha << endl;
            f1=model.function(d);
            f1+=0.1*alpha*gSq;
            //cout << "f0 = " << f << " --> " << " f1 = " << f1 << endl;
            if(f1<=f){
                x=d;
                ++it;
//                cout << "Energy Dropped" << endl;
                break;
            }
            else alpha*=tau;
        }

        //check convergence
        df=abs(f1-f0);
        if(abs(f1)<tol || df<tol) break;
        else f0=f1;
        //recalculate gradient
        g=model.gradient(x);
    }

    if(it==itMax-1) status[0]=2;
    else status[0]=0;
    status[1]=it;

    return status;
}

//##### NEWTON METHOD FOR ROOT FINDING #####

//Default constructor
template <typename M>
Newton<M>::Newton() {

    itMax=1000;
    tol=1e-6;
}

//Construct with optimisation parameters
template <typename M>
Newton<M>::Newton(int iterationLimit, long double convergenceTolerance) {

    itMax=iterationLimit;
    tol=convergenceTolerance;
}

//Root finding
template <typename M>
VecF<int> Newton<M>::operator()(M &model, long double &x) {

    VecF<int> status(2);
    status[0]=0;
    for(int i=1; i<=itMax; ++i){
        long double f0=model.function(x);
        long double f1=model.gradient(x);
        long double xx=x-f0/f1;
        if(fabs(x-xx)<tol){
            status[0]=1;
            status[1]=i;
            x=xx;
            break;
        }
        x=xx;
    }

    return status;
}

//##### HALLEY METHOD FOR ROOT FINDING #####

//Default constructor
template <typename M>
Halley<M>::Halley() {
    itMax=1000;
    tol=1e-6;
}

//Construct with optimisation parameters
template <typename M>
Halley<M>::Halley(int iterationLimit, long double convergenceTolerance) {
    itMax=iterationLimit;
    tol=convergenceTolerance;
}

//Root finding
template <typename M>
VecF<int> Halley<M>::operator()(M &model, long double &x) {

    VecF<int> status(2);
    status[0]=0;
    for(int i=1; i<=itMax; ++i){
        long double f0=model.function(x);
        long double f1=model.gradient(x);
        long double f2=model.hessian(x);
        long double xx=x-(2.0*f0*f1)/(2.0*f1*f1-f0*f2);
        if(fabs(x-xx)<tol){
            status[0]=1;
            status[1]=i;
            x=xx;
            break;
        }
        x=xx;
    }

    return status;
}
