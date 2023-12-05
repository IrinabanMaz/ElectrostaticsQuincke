#pragma once
#include "StokesSolver/ScalarSphericalHarmonicSeries.h"
#include "StokesSolver/GMRES.h"
/**
 * @brief Evaluates the Single Layer operator for points away from the surface, acting on ScalSphereData.
*/
class LaplaceSingleLayerDirectDiscrete : public SphericalScalarFunction
{
private:

    SphereCoordSystem* system;
    ScalSphereData* data;
    QuadConstants* consts;

    inline double kernel(const SphereCoord& x, const SphereCoord& y) const
    {
        RectCoord r = x - y;
        double R = norm(r);

        return 0.25 / (consts->PI * R);
    }
public:

    /**
     * @brief Constructs the functor, assigned to a SphereData.
     * @param d The density the opeator is acting on.
     * @param c The center of the sphere to integrate over.
     * @param na the name.
    */
    LaplaceSingleLayerDirectDiscrete(ScalSphereData* d, QuadConstants* cts, SphereCoordSystem* sys = &defaultsystem, std::string na = "Direct Single Layer Discrete")
    {
        data = d;
        name = na;
        system = sys;
        consts = cts;
    }

    /**
     * @brief Evaluates the operator.
     * @return
    */
    double  operator()(const SphereCoord& x) const
    {



        double total = 0.0;
        //loop over trapezoidal(outside) and legendre(inside) nodes.
        for (int p = 0; p < consts->NUMTRAPNODES; p++)
            for (int i = 0; i < consts->NUMGLNODES; i++)
            {
                //generate point on unit sphere, note theta is the first argument, phi the second.
                SurfaceCoord s(consts->GLnodes[i], 2.0 * consts->PI * (double)p / (double)consts->NUMTRAPNODES);
                //convert to spherical coordinate.
                SphereCoord y(s, system);
                //add contribution at y to integral.
                total = total + kernel(x,y) *(*data)[i][p] * consts->GLweights[i] * (2.0 * consts->PI / consts->NUMTRAPNODES);
            }

        return total;
    }
};

/**
 * @brief Functor for a term in the Single Layer integral operator evaluated using a VSHSeries representation.
*/
class LaplaceSingleLayerHarmonicTerm : public SphericalScalarFunction
{
private:

    Legendre* P;

    int m_;
    int n_;

    double rhohatYr;
    double rhohatYi;

    YReal yr;
    YImag yi;

    /**
 * @brief Coefficient of YReal or YImag in the SingleLayerHarmonic series.
 * @param r Radius from a SphereCoord.
 * @param n Degree of the YReal or YImag term.
 * @param rhohatY Normalized Fourier coefficient of the associated YReal or YImag function.
 * @return @f{equation}{
f_{n,m}^V(r) = \frac{1}{(2n + 1)}r^{-(n+1)} \hat{\rho}_{n,m}}
*/
    inline double f(double r, int n, double rhohatY) const
    {

        return rhohatY / (double)(2 * n + 1) * pow(r , -n-1);
    }
    /**
     * @brief Computes the derivative of fV() with respect to r. Used in helper functions for SingleLayerHarmonic.
     * @param r Radius from a SphereCoord.
     * @param n Degree of the VReal or VImag term.
     * @param rhohatV Normalized Fourier coefficient of the associated VReal or VImag function.
     * @param rhohatW Normalized Fourier coefficient of the associated WReal or WImag function.
     * @return @f{equation}{
    f_{n,m}^{V'}(r) = \frac{-n(n+2)}{(2n + 1)(2n+3)}r^{-(n+3)} \hat{\rho}_{n,m}^V + \frac{n+1}{4n + 2}\left(n r^{-(n+1) - (n+2)r^{-(n+3)}}\right)\hat{\rho}_{n,m}^W
     @f}
    */
    inline double fprime(double r, int n, double rhohatY) const
    {
        return ( - n - 1) * rhohatY / (double)(2 * n + 1) * pow(r, -n - 2);
    }

public:

    LaplaceSingleLayerHarmonicTerm(Legendre* poly = nullptr) :yr(poly), yi(poly) { P = poly; }
    /**
     * @brief Constructs the term in the series.
     * @param series VSHSeries corresponding to data to evaluate.
     * @param m The order of the term.
     * @param n The degree of the term.
    */
    LaplaceSingleLayerHarmonicTerm(const SSHSeries& series, int m, int n) :
        yr(series.P , m, n),
        yi(series.P , m, n)
    {
        m_ = m;
        n_ = n;

        P = series.P;

        rhohatYr = series.terms[n][m].rhohatYr;
        rhohatYi = series.terms[n][m].rhohatYi;

    }

    LaplaceSingleLayerHarmonicTerm(const SSHSeries& series, int m, int n , Legendre*poly) :
        yr(poly, m, n),
        yi(poly, m, n)
    {
        m_ = m;
        n_ = n;

        P = poly;

        rhohatYr = series.terms[n][m].rhohatYr;
        rhohatYi = series.terms[n][m].rhohatYi;

    }

    LaplaceSingleLayerHarmonicTerm(const LaplaceSingleLayerHarmonicTerm& term) :
        yr(term.yr),
        yi(term.yi)
    {

        int m = m_ = term.m_;
        int n = n_ = term.n_;

        rhohatYr = term.rhohatYr;
        rhohatYi = term.rhohatYi;

    }

    LaplaceSingleLayerHarmonicTerm(const LaplaceSingleLayerHarmonicTerm& term , Legendre* poly) :
        yr(term.yr , poly),
        yi(term.yi , poly)
    {

        int m = m_ = term.m_;
        int n = n_ = term.n_;

        rhohatYr = term.rhohatYr;
        rhohatYi = term.rhohatYi;

        P = poly;

    }


    LaplaceSingleLayerHarmonicTerm& operator =(const LaplaceSingleLayerHarmonicTerm& term)
    {
        P = term.P;
        int m = m_ = term.m_;
        int n = n_ = term.n_;
        yr = YReal(P,m, n);
        yi = YImag(P,m, n);

        rhohatYr = term.rhohatYr;
        rhohatYi = term.rhohatYi;

        return *this;

    }



    /**
     * @brief Evaluates the term.
     * @return @f{equation}{
     f_{n,m}(r)\boldsymbol{Y}_n^m(\theta , \phi)
    @f}

     see f , YReal , YImag 
    */
    inline double  operator()(const SphereCoord& s) const
    {
        return f(s.r, n_, rhohatYr) * yr(s) + f(s.r, n_,rhohatYi) * yi(s);

    }

    inline double dR(const SphereCoord& s) const
    {
        return fprime(s.r, n_, rhohatYr) * yr(s) + fprime(s.r, n_, rhohatYi) * yi(s);
    }

   

   
    inline double dTheta(const SphereCoord& s) const
    {
        return f(s.r, n_, rhohatYr) * yr.dTheta(s) + f(s.r, n_, rhohatYi) * yi.dTheta(s);
    }


    inline double dPhi(const SphereCoord& s) const
    {
        return f(s.r, n_, rhohatYr) * yr.dPhi(s) + f(s.r, n_, rhohatYi) * yi.dPhi(s);
    }



};

/**
 * @brief Functor for Single Layer operator. Sums SingleLayerHarmonicTerm objects.
*/
class LaplaceSingleLayerHarmonic : public SphereScalFunctionSum
{
private:

    Legendre* P;

    std::vector<std::vector<LaplaceSingleLayerHarmonicTerm>> terms;

public:

    SphereCoordSystem* system; /**< The Spherical Coordinate system of the boundary integral. */

    /**
     * @brief Constructs the Functor.
     * @param n number of the terms in the series.
     * @param c center of the SphereCoord in the integral.
    */
    LaplaceSingleLayerHarmonic(int n,Legendre*poly ,  SphereCoordSystem* sys) : P(poly)
    {

        system = sys;
        terms.resize(n + 1);

        for (size_t i = 0; i < terms.size(); i++)
        {
            terms[i].resize(i + 1);
            for (int j = 0; j <= i; j++)
                append(&terms[i][j]);
        }
    }

    /**
     * @brief Generates the Functor from a SSHSeries representation.
     * @param series the Source boundary data.
    */
    LaplaceSingleLayerHarmonic(const SSHSeries& series) : P(series.P)
    {
        terms.resize(series.N + 1);
        for (int i = 0; i <= series.N; i++)
            terms[i].resize(i + 1);
        for (int n = 0; n <= series.N; n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = LaplaceSingleLayerHarmonicTerm(series, m, n);
                append(&terms[n][m]);
            }

        system = series.system;

    }
    /**
     * @brief Copy Consructor
     * @param slh Single Layer to copy.
    */
    LaplaceSingleLayerHarmonic(const LaplaceSingleLayerHarmonic& slh) : P(slh.P)
    {
        terms = slh.terms;
        for (int n = 0; n < terms.size(); n++)
            for (int m = 0; m <= n; m++)
            {
                
                append(&terms[n][m]);
            }

        system = slh.system;
    }

    /**
     * @brief Copy Consructor
     * @param slh Single Layer to copy.
    */
    LaplaceSingleLayerHarmonic(const LaplaceSingleLayerHarmonic& slh , Legendre * poly) : P(poly)
    {
        terms = slh.terms;
        for (int n = 0; n < terms.size(); n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = LaplaceSingleLayerHarmonicTerm(slh.terms[n][m], poly);
                append(&terms[n][m]);
            }

        system = slh.system;
    }



    /**
     * @brief Assignment operator.
     * @param slh Single Layer to assign.
    */
    LaplaceSingleLayerHarmonic& operator =(const LaplaceSingleLayerHarmonic& slh)
    {
        terms = slh.terms;
        system = slh.system;
        clear();

        for (int n = 0; n < terms.size(); n++)
            for (int m = 0; m <= n; m++)
                append(&terms[n][m]);

        return *this;
    }





    /**
     * @brief Generates the Functor from a Fourier Series representation.
     * @param series the Source boundary data.
    */
    void solve(SSHSeries& series)
    {
        clear();
        terms.clear();
        terms.resize(series.N + 1);
        for (int i = 0; i <= series.N; i++)
            terms[i].resize(i + 1);
        for (int n = 0; n <= series.N; n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = LaplaceSingleLayerHarmonicTerm(series, m, n,P);
                append(&terms[n][m]);
            }

        system = series.system;
    }

    /**
     * @brief Evaluates the Single Layer operator at a point.
     * @return The sum of terms. See SingleLayerHarmonicTerm::operator()().
    */
    double operator()(const SphereCoord& s) const
    {

        SphereCoord temp = recast(s, system);
        int N = terms.size();
        P->populate(cos(temp.s.theta), N, N);
        return SphereScalFunctionSum::operator()(temp);

    }

    double dR(const SphereCoord& s) const
    {
        double temp = 0.0;

        SphereCoord sc = recast(s, system);

        for (int i = 0; i < terms.size(); i++)
            for (int j = 0; j < terms[i].size(); j++)
            {
                temp += terms[i][j].dR(sc);
            }

        return temp;
    }

    double dTheta(const SphereCoord& s) const
    {
        double temp = 0.0;

        SphereCoord sc = recast(s, system);

        for (int i = 0; i < terms.size(); i++)
            for (int j = 0; j < terms[i].size(); j++)
            {
                temp += terms[i][j].dTheta(sc);
            }

        return temp;
    }

    double dPhi(const SphereCoord& s) const
    {
        double temp = 0.0;

        SphereCoord sc = recast(s, system);

        for (int i = 0; i < terms.size(); i++)
            for (int j = 0; j < terms[i].size(); j++)
            {
                temp += terms[i][j].dPhi(sc);
            }

        return temp;
    }



};


/**
 * @brief Evaluates the Single Layer operator for points away from the surface, acting on ScalSphereData.
*/
class LaplaceTractionDirectDiscrete : public SphericalScalarFunction
{
private:

    SphereCoordSystem* system;
    ScalSphereData* data;
    QuadConstants* consts;

    inline double kernel(const SphereCoord& x, const SphereCoord& y) const
    {
        RectCoord r = x - y;
        double R = norm(r);

        return -0.25 * dot(r , e_r(x)) / (consts->PI * R * R * R);
    }
public:

    /**
     * @brief Constructs the functor, assigned to a SphereData.
     * @param d The density the opeator is acting on.
     * @param c The center of the sphere to integrate over.
     * @param na the name.
    */
    LaplaceTractionDirectDiscrete(ScalSphereData* d, QuadConstants* cts, SphereCoordSystem* sys, std::string na = "Direct Single Layer Discrete")
    {
        data = d;
        name = na;
        system = sys;
        consts = cts;
    }

    /**
     * @brief Evaluates the operator.
     * @return
    */
    double  operator()(const SphereCoord& x) const
    {



        double total = 0.0;
        //loop over trapezoidal(outside) and legendre(inside) nodes.
        for (int p = 0; p < consts->NUMTRAPNODES; p++)
            for (int i = 0; i < consts->NUMGLNODES; i++)
            {
                //generate point on unit sphere, note theta is the first argument, phi the second.
                SurfaceCoord s(consts->GLnodes[i], 2.0 * consts->PI * (double)p / (double)consts->NUMTRAPNODES);
                //convert to spherical coordinate.
                SphereCoord y(s, system);
                //add contribution at y to integral.
                total = total + kernel(x, y) * (*data)[i][p] * consts->GLweights[i] * (2.0 * consts->PI / consts->NUMTRAPNODES);
            }

        return total;
    }
};

/**
 * @brief Evaluates the Single Layer operator for points away from the surface, acting on ScalSphereData.
*/
class LaplaceTractionHarmonic : public SphericalScalarFunction
{
private:

    LaplaceSingleLayerHarmonic* slh;

    
public:

    /**
     * @brief Constructs the functor, assigned to a SphereData.
     * @param d The density the opeator is acting on.
     * @param c The center of the sphere to integrate over.
     * @param na the name.
    */
    LaplaceTractionHarmonic(LaplaceSingleLayerHarmonic* s)
    {
        slh = s;
    }

    /**
     * @brief Evaluates the operator.
     * @return
    */
    double  operator()(const SphereCoord& x) const
    {



        SphereCoord temp = recast(x, slh->system);

        double er = dot(e_r(x), e_r(temp));
        double etheta = dot(e_r(x), e_theta(temp));
        double ephi = dot(e_r(x), e_phi(temp));

        return slh->dR(temp) * er 
             + slh->dTheta(temp) * etheta / temp.r 
             + slh->dPhi(temp) * ephi / temp.r / std::sin(temp.s.theta);
    }
};

/**
 * @brief Evaluates the Single Layer operator for points away from the surface, acting on ScalSphereData.
*/
class ElectricField : public SphericalVectorField
{
private:

    LaplaceSingleLayerHarmonic* slh;
    
    
public:

    /**
     * @brief Constructs the functor, assigned to a SphereData.
     * @param d The density the opeator is acting on.
     * @param c The center of the sphere to integrate over.
     * @param na the name.
    */
    ElectricField(LaplaceSingleLayerHarmonic* s)
    {
        slh = s;
    }

    /**
     * @brief Evaluates the operator.
     * @return
    */
    RectCoord  operator()(const SphereCoord& x) const
    {



        SphereCoord temp = recast(x, slh->system);

        return -slh->dR(temp) * e_r(temp);
             - slh->dTheta(temp) / temp.r * e_theta(temp)
             - slh->dPhi(temp) / temp.r / std::sin(temp.s.theta) * e_phi(temp);
    }
};

