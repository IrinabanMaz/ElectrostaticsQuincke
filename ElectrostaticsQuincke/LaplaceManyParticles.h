#pragma once


#pragma once
#include "QuinckeRotator.h"
#include <ctime>

class LaplaceParticleSystem : public SphericalScalarFunction
{
public:

	std::vector<LaplaceParticle> particles;
	SphereScalFunctionSum traction;
	VectorFieldSum E;

	LaplaceParticleSystem(std::vector<SphereCoordSystem> cs = std::vector<SphereCoordSystem>(), int n = 6)
	{
		particles.resize(cs.size());
		for (size_t i = 0; i < particles.size(); i++)
		{
			particles[i] = LaplaceParticle(n, cs[i]);

			traction.append(&particles[i].traction);
			E.append(&particles[i].E);
		}
	}

	LaplaceParticleSystem(std::vector<ScalSphereData> ds, std::vector<SphereCoordSystem> cs, int nst , RectCoord E0)
	{

		if (ds.size() != cs.size())
			std::cout << "Error: Data size and centers size in paticle system are not equal.\n";

		particles.resize(cs.size());
		for (size_t i = 0; i < particles.size(); i++)
		{
			particles[i] = LaplaceParticle(ds[i], nst, cs[i]);

			traction.append(&particles[i].traction);
			E.append(&particles[i].E);
		}


	}


	LaplaceParticleSystem(const LaplaceParticleSystem& ps)
	{
		particles = ps.particles;

		traction.clear();
		E.clear();

		for (int i = 0; i < particles.size(); i++)
		{
			traction.append(&particles[i].traction);
			E.append(&particles[i].E);
		}
	}


	LaplaceParticleSystem(const int& n)
	{
		particles.resize(n);
		for (int i = 0; i < particles.size(); i++)
		{
			traction.append(&particles[i].traction);
			E.append(&particles[i].E);
		}
	}

	void push_back(const LaplaceParticle& p)
	{
		particles.push_back(p);
		traction.clear();
		E.clear();
		for (int i = 0; i < particles.size(); i++)
		{
			traction.append(&particles[i].traction);
			E.append(&particles[i].E);
		}
	}

	void setPermitivities(std::vector<double> perms)
	{
		for (int i = 0; i < particles.size(); i++)
			particles[i].setPermitivity(perms[i]);
	}
	LaplaceParticleSystem operator +(const LaplaceParticleSystem& ps) const
	{
		LaplaceParticleSystem temp(*this);

		if (particles.size() != ps.particles.size())
			std::cout << "Error:Summing different numbers of particles.\n";

		for (size_t i = 0; i < particles.size(); i++)
			temp.particles[i] += ps.particles[i];

		return temp;

	}

	LaplaceParticleSystem operator -(const LaplaceParticleSystem& ps) const
	{
		LaplaceParticleSystem temp(*this);

		if (particles.size() != ps.particles.size())
			std::cout << "Error:Summing different numbers of particles.\n";

		for (size_t i = 0; i < particles.size(); i++)
			temp.particles[i] = particles[i] - ps.particles[i];

		return temp;

	}

	LaplaceParticleSystem operator -() const
	{
		LaplaceParticleSystem temp(*this);


		for (size_t i = 0; i < particles.size(); i++)
			temp.particles[i] = -particles[i];

		return temp;

	}

	LaplaceParticleSystem& operator +=(const LaplaceParticleSystem& ps)
	{
		for (int i = 0; i < particles.size(); i++)
			particles[i] += ps.particles[i];

		return *this;
	}

	LaplaceParticleSystem& operator *=(const double& a)
	{
		for (int i = 0; i < particles.size(); i++)
			particles[i] *= a;

		return *this;
	}

	LaplaceParticleSystem operator *(const double& a) const
	{
		LaplaceParticleSystem temp(*this);


		for (size_t i = 0; i < particles.size(); i++)
			temp.particles[i] *= a;

		return temp;

	}

	LaplaceParticleSystem operator *(std::vector<double> scales)
	{
		LaplaceParticleSystem temp(*this);


		for (size_t i = 0; i < particles.size(); i++)
			temp.particles[i] *= scales[i];

		return temp;
	}
	LaplaceParticleSystem& operator =(const LaplaceParticleSystem& ps)
	{
		particles = ps.particles;

		traction.clear();
		E.clear();

		for (int i = 0; i < particles.size(); i++)
		{
			traction.append(&particles[i].traction);
			E.append(&particles[i].E);
		}

		return *this;

	}




	double norm() const
	{
		double temp = 0.0;
		for (int i = 0; i < particles.size(); i++)
		{
			double normj = particles[i].norm();
			temp += normj * normj;
		}
		return sqrt(temp);

	}

	double LInfNorm()
	{
		double max = 0.0;

		for (int i = 0; i < particles.size(); i++)
		{
			double temp = particles[i].lInfNorm();
			max = temp > max ? temp : max;
		}
		return max;
	}

	double dot(const LaplaceParticleSystem& ps) const
	{
		double temp = 0.0;
		for (size_t i = 0; i < particles.size(); i++)
		{
			temp += particles[i].dot(ps.particles[i]);
		}
		return temp;

	}

	static double dot(const LaplaceParticleSystem& particles, const LaplaceParticleSystem& qarticles)
	{
		return particles.dot(qarticles);
	}

	int getN() const
	{
		return particles.size();
	}

	void axpy(LaplaceParticleSystem* ps, double scale)
	{
		for (int i = 0; i < particles.size(); i++)
			particles[i].axpy(&(ps->particles[i]), scale);

	}


	double operator()(const SphereCoord& x) const
	{
		double temp = 0.0;

		for (int i = 0; i < particles.size(); i++)
		{

				temp += particles[i].potential(x);
		}

		return temp;
	}

};


class LPSQuinckeOperator
{
private:
	double permitivity;
public:

	LPSQuinckeOperator(double eps): permitivity(eps){}
	LaplaceParticleSystem operator *(const LaplaceParticleSystem& ps)
	{
		LaplaceQuinckeOperator L(permitivity);

	    LaplaceParticleSystem temp(ps.getN());

		for (size_t i = 0; i < ps.particles.size(); i++)
		{
			std::time_t ti = std::time(0);
			temp.particles[i] = L * ps.particles[i];
			
			SphereCoordSystem system = ps.particles[i].system;
			LaplaceParticle tr(ps.particles[i].numSeriesTerms , system);

			ScalSphereData t(2 * ps.particles[i].numSeriesTerms + 1, ps.particles[i].numSeriesTerms + 1);

			double lambda = (ps.particles[i].permitivity - permitivity) / (ps.particles[i].permitivity + permitivity);

			for (size_t j = 0; j < ps.particles.size(); j++)
			{
				if (i != j)
				{
						t = t + discretize(&ps.particles[j].traction, &system, ps.particles[i].consts);
				}
			}
			temp.particles[i] += LaplaceParticle(t, ps.particles[i].numSeriesTerms , system) * lambda;
			std::cout << "Operator on particle " << i + 1 << " computed. Time elapsed: " << std::time(0) - ti << "\n";
		}
		return temp;
	}
};

class LPSNeumannOperator
{

public:

	LaplaceParticleSystem operator *(const LaplaceParticleSystem& ps)
	{
		LaplaceNeumannOperator L;

	    LaplaceParticleSystem temp(ps.getN());

		for (size_t i = 0; i < ps.particles.size(); i++)
		{
			std::time_t ti = std::time(0);
			temp.particles[i] = L * ps.particles[i];
			SphereCoordSystem system = ps.particles[i].system;
			LaplaceParticle tr(ps.particles[i].numSeriesTerms , system);

			ScalSphereData t(2 * ps.particles[i].numSeriesTerms + 1, ps.particles[i].numSeriesTerms + 1);



			for (size_t j = 0; j < ps.particles.size(); j++)
			{
				if (i != j)
				{
						t = t + discretize(&ps.particles[j].traction, &system, ps.particles[i].consts);
				}
			}
			temp.particles[i] += LaplaceParticle(t, ps.particles[i].numSeriesTerms , system);
			std::cout << "Operator on particle " << i + 1 << " computed. Time elapsed: " << std::time(0) - ti << "\n";
		}

		return temp;
	}
};


class LPSIdentityPreconditioner
{
public:
	LaplaceParticleSystem solve(const LaplaceParticleSystem& ps)
	{
		return ps;
	}
};




LaplaceParticleSystem SolveQuinckeManySphereFromKnown(std::vector<SphereCoordSystem> cs,  int seriesSize,double l , int numIters = 30, double r = 1e-9)
{

	int size = cs.size();
	LaplaceParticleSystem solution(cs, seriesSize);


	LPSQuinckeOperator L(l);



	LaplaceParticleSystem rh(cs, seriesSize);

	for (int i = 0; i < rh.particles.size(); i++)
		rh.particles[i].generateRandomSmooth();

	LaplaceParticleSystem rhs = L * rh;

	LPSIdentityPreconditioner I;

	solution = -rhs;

	time_t t0= time(0);
	GMRES(&L, &solution, &rhs, &I, numIters, 1, r);
	std::cout << "Time to solve: " << time(0) - t0 << std::endl;
	//std::cout << Integrate(&BIEsolution, 1.0, BIEsolution.particles[0].center) * 0.25 / PI << std::endl;
	//std::cout << Integrate(&rh, 1.0, rh.particles[0].center) * 0.25 / PI << std::endl;

	SphereCoordSystem awaysystem = SphereCoordSystem(RectCoord(8.0, 0.0, 0.0));
	LaplaceParticle solutionEvaluatedAway(discretize(&rh, &awaysystem, rh.particles[0].consts),
		seriesSize,
		RectCoord(8.0, 0.0, 0.0));

	LaplaceParticle approxSolutionEvaluatedAway(discretize(&solution, &awaysystem, solution.particles[0].consts),
		seriesSize,
		RectCoord(8.0, 0.0, 0.0));

    	
	std::cout << "Error in BIE method(L 2 on particle system): " << (rh - solution).norm() << std::endl;
	std::cout << "Error in BIE method(L infinity on particle 1): " << LInfDifference(rh.particles[0].data , solution.particles[0].data , rh.particles[0].consts) << std::endl;
 	std::cout << "Error in BIE method(L infinity on particle 2): " << LInfDifference(rh.particles[1].data , solution.particles[1].data , rh.particles[1].consts) << std::endl;
 
	return solution;


}

class NormalField : public SphericalScalarFunction
{
private:
	RectCoord E;

public:
	NormalField(const RectCoord& e)
	{
		E = e;
	}

	double operator()(const SphereCoord& x) const
	{
		return dot(E, e_r(x));
	}
};

LaplaceParticleSystem SolveQuinckeRandom(std::vector<SphereCoordSystem> cs, int seriesSize, double epsplus , double epsminus, RectCoord E0,int numIters = 30, double r = 1e-9)
{

	int size = cs.size();
	LaplaceParticleSystem solution(cs, seriesSize);

	double lambda = (epsminus - epsplus) / (epsplus + epsminus);
	LPSQuinckeOperator L(lambda);



	static LaplaceParticleSystem maxq(cs, seriesSize);
	static bool called = false;

    LaplaceParticleSystem q(cs, seriesSize);

	if (!called)
	{
		called = true;
		for (int i = 0; i < maxq.particles.size(); i++)
			maxq.particles[i].generateRandomSmooth();
	}
	
	
		for (int i = 0; i < q.particles.size(); i++)
		{
			ScalSphereData data = discretize(&maxq.particles[i].fourierdata,
				&maxq.particles[i].system,
				q.particles[i].consts);
			q.particles[i] = LaplaceParticle(data,seriesSize, q.particles[i].system);
	   }
	
	LaplaceParticleSystem appliedFieldTerm(cs, seriesSize);

	for (int i = 0; i < appliedFieldTerm.particles.size(); i++)
	{
		NormalField field(E0);
		ScalSphereData data = discretize(&field,
			&appliedFieldTerm.particles[i].system,
			appliedFieldTerm.particles[i].consts);
		appliedFieldTerm.particles[i] = LaplaceParticle(data, seriesSize, appliedFieldTerm.particles[i].system);
	}
	LaplaceParticleSystem rhs =  q * (1 / (epsplus + epsminus)) + appliedFieldTerm * lambda;

	LPSIdentityPreconditioner I;

	time_t t0 = time(0);
	GMRES(&L, &solution, &rhs, &I, numIters, 1, r);
	std::cout << "Time to solve: " << time(0) - t0 << std::endl;

	LaplaceParticleSystem Enplus;
	Enplus = (q - solution * epsminus) * (1.0 / (epsplus - epsminus));

	std::vector<SphereData> EData;
	for (int i = 0; i < cs.size(); i++)
	{
		EData.push_back(discretize(&solution.particles[i].E,
			&solution.particles[i].system,
			solution.particles[i].consts));
	}

	//std::cout << Integrate(&BIEsolution, 1.0, BIEsolution.particles[0].center) * 0.25 / PI << std::endl;
	//std::cout << Integrate(&rh, 1.0, rh.particles[0].center) * 0.25 / PI << std::endl;

	return solution;


}

LaplaceParticleSystem SolveNeumannManySphereFromKnown(std::vector<SphereCoordSystem> cs,  int seriesSize,double l , int numIters = 30, double r = 1e-9)
{

	int size = cs.size();
	LaplaceParticleSystem solution(cs, seriesSize);


	LPSNeumannOperator L;



	LaplaceParticleSystem rh(cs, seriesSize);

	for (int i = 0; i < rh.particles.size(); i++)
		rh.particles[i].generateRandomSmooth();


	LaplaceParticleSystem rhs = L * rh;

	LPSIdentityPreconditioner I;

	solution = rh * 1.5;

	time_t t = time(0);
	GMRES(&L, &solution, &rhs, &I, numIters, 1, r);

	std::cout << "Time to solve: " << time(0) - t << "\n";

	//std::cout << Integrate(&BIEsolution, 1.0, BIEsolution.particles[0].center) * 0.25 / PI << std::endl;
	//std::cout << Integrate(&rh, 1.0, rh.particles[0].center) * 0.25 / PI << std::endl;

	SphereCoordSystem awaysystem = SphereCoordSystem(RectCoord(8.0, 0.0, 0.0));
	LaplaceParticle solutionEvaluatedAway(discretize(&rh, &awaysystem, rh.particles[0].consts),
		seriesSize,
		RectCoord(8.0, 0.0, 0.0));

	LaplaceParticle approxSolutionEvaluatedAway(discretize(&solution, &awaysystem, solution.particles[0].consts),
		seriesSize,
		RectCoord(8.0, 0.0, 0.0));

    	
	std::cout << "Error in BIE method(L 2 on particle system): " << (rh - solution).norm() << std::endl;

	for(int i = 0; i < solution.particles.size(); i++)
	std::cout << "Error in BIE method(L infinity on particle " << i + 1 <<"): " << LInfDifference(rh.particles[i].data , solution.particles[i].data , rh.particles[i].consts) << std::endl;

 
	return solution;

}
