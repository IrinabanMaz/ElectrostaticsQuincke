#pragma once

#include "LaplaceOperators.h"

class LaplaceParticle
{
public:
	SphereCoordSystem system;
    QuadConstants consts;
	Legendre poly;
	int numSeriesTerms;

	double permitivity;

	ScalSphereData data;
	LaplaceSingleLayerDirectDiscrete quadpotential;
	LaplaceTractionDirectDiscrete quadtraction;

	SSHSeries fourierdata;
	LaplaceSingleLayerHarmonic potential;
	LaplaceTractionHarmonic traction;

	ElectricField E;

	LaplaceParticle(int p = 0, SphereCoordSystem sys = defaultsystem ) :
		consts(2 * p + 1,  p + 1),
		data(2 * p + 1,  p + 1),
		fourierdata(p, &consts, &poly, &system),
		quadpotential(&data,&consts, &system),
		potential(fourierdata, &poly),
		quadtraction(&data, &consts, &system),
		traction(&potential),
		E(&potential )
	{
		fourierdata.approximate(data, p);
		potential.solve(fourierdata);
		numSeriesTerms = p;
		system = sys;
	}
	LaplaceParticle(const ScalSphereData & d , int p , SphereCoordSystem sys = defaultsystem):
		consts(2 * p + 1, p + 1),
		data(d),
		fourierdata(p, &consts, &poly, &system),
		quadpotential(&data, &consts, &system),
		potential(fourierdata, &poly),
		quadtraction(&data, &consts, &system),
		traction(&potential),
		E(&potential )
	{
		numSeriesTerms = p;
		fourierdata.approximate(data, p);
		potential.solve(fourierdata);
		system = sys;
	}

	LaplaceParticle(const LaplaceParticle & p):
		consts(p.consts),
		data(p.data),
		fourierdata(p.fourierdata, &consts, &poly , &system),
		quadpotential(&data, &consts, &system),
		potential(fourierdata, &poly),
		quadtraction(&data, &consts, &system),
		traction(&potential),
		E(&potential)
	{
		permitivity = p.permitivity;
		
		system = p.system;
		numSeriesTerms = p.numSeriesTerms;
	}


	LaplaceParticle &  operator = (const LaplaceParticle& particle)
	{
		system = particle.system;
		consts = particle.consts;
		numSeriesTerms = particle.numSeriesTerms;
		data = particle.data;
		fourierdata = SSHSeries(particle.fourierdata, &consts, &poly , &system);
		potential.solve(fourierdata);

		permitivity = particle.permitivity;

		return *this;
	}
	void generateRandomData()
	{
		std::default_random_engine gen;
		std::uniform_real_distribution<double> r(1.0 , 10.0);

		size_t numtheta = data.size();
		size_t numphi = data[0].size();

		for (size_t t = 0; t < numtheta; t++)
			for (size_t p = 0; p < numphi; p++)
				data[t][p] = r(gen);

		fourierdata.approximate(data, numSeriesTerms);
		double check = validateData();
		potential.solve(fourierdata);

	}

	void generateRandomSmooth()
	{

		fourierdata.generateRandomSmooth(numSeriesTerms);
		data = discretize(&fourierdata, &system, consts);

		double check = validateData();
		potential.solve(fourierdata);

	}

	void setPermitivity(double eps)
	{
		permitivity = eps;
	}


	double validateData()
	{
		return LInfDifference(data, &fourierdata, consts , &system);
	}

	LaplaceParticle operator +(const LaplaceParticle& p) const
	{
		LaplaceParticle temp(*this);

	//	double check1 = temp.validateData();
	//	double check2 = p.validateData();

	//	if (check1 > 1e-4 || check2 > 1e-4)
	//		std::cout << " Fourier Data invalid at time of + operator.\n";

		for(int i = 0; i < consts.NUMGLNODES; i++)
			for(int j = 0; j < consts.NUMTRAPNODES; j++)
				temp.data[i][j] += p.data[i][j];

		temp.fourierdata += p.fourierdata;

		temp.potential.solve(temp.fourierdata);

	//	if (temp.validateData() > 1.0001*( check1 + check2))
	//		std::cout << " + operator is invalid.\n";

		return temp;
	}

	LaplaceParticle operator -(const LaplaceParticle& p) const
	{
		LaplaceParticle temp(*this);

		//double check1 = temp.validateData();
		//double check2 = p.validateData();

		//if (check1 > 1e-4 || check2 > 1e-4)
		//	std::cout << " Fourier Data invalid at time of - operator.\n";


		for (int i = 0; i < consts.NUMGLNODES; i++)
			for (int j = 0; j < consts.NUMTRAPNODES; j++)
				temp.data[i][j] -= p.data[i][j];

		temp.fourierdata -= p.fourierdata;

		temp.potential.solve(temp.fourierdata);

		//if (temp.validateData() > 1.0001 * (check1 + check2))
		//	std::cout << "binary - operator is invalid.\n";

		return temp;
	}

	LaplaceParticle operator -() const
	{
		LaplaceParticle temp(*this);

		double check = temp.validateData();
		for (int t = 0; t < consts.NUMGLNODES; t++)
			for (int p = 0; p < consts.NUMTRAPNODES; p++)
				temp.data[t][p] = -data[t][p];

		temp.fourierdata = fourierdata.uminus(&temp.poly);

		temp.numSeriesTerms = numSeriesTerms;

		temp.potential.solve(temp.fourierdata);

		//if (temp.validateData() > 1.0001 * check)
		//	std::cout << "unary - operator is invalid.\n";
		return temp;
	}

	LaplaceParticle& operator +=(const LaplaceParticle& particle)
	{

		//double check1 = validateData();
		//double check2 = particle.validateData();

		//if (check1 > 1e-4 || check2 > 1e-4)
		//	std::cout << " Fourier Data invalid at time of += operator.\n";

		for (int t = 0; t < (int)consts.NUMGLNODES; t++)
			for (int p = 0; p < (int)consts.NUMTRAPNODES; p++)
				data[t][p] += particle.data[t][p];


		fourierdata += particle.fourierdata;



		potential.solve(fourierdata);

		//if (validateData() > 1.0001 * (check1 + check2))
		//	std::cout << " += operator is invalid.\n";

		return *this;
	}

	LaplaceParticle& operator *=(const double& a)
	{
		//double check = validateData();

		//if (check > 1e-4)
		//	std::cout << " Fourier Data invalid at time of *= operator.\n";


		for (int t = 0; t < (int)consts.NUMGLNODES; t++)
			for (int p = 0; p < (int)consts.NUMTRAPNODES; p++)
				data[t][p] *= a;


		fourierdata *= a;


		//if (validateData() > 1.0001 * std::abs(a) * check)
		//	std::cout << " *= operator is invalid.\n";

		potential.solve(fourierdata);
		return *this;
	}

	LaplaceParticle operator *(const double& a) const
	{
		LaplaceParticle temp(*this);

		//double check = temp.validateData();

		//if (check > 1e-4)
		//	std::cout << " Fourier Data invalid at time of * operator.\n";


		for (int t = 0; t < (int)consts.NUMGLNODES; t++)
			for (int p = 0; p < (int)consts.NUMTRAPNODES; p++)
				temp.data[t][p] *= a;

		temp.fourierdata *= a;

		temp.numSeriesTerms = numSeriesTerms;

		temp.potential.solve(temp.fourierdata);

		temp.system = system;


		//if (temp.validateData() > 1.0001 *std::abs(a) * check)
		//	std::cout << " * operator is invalid.\n";

		return temp;
	}


	double norm() const
	{
		return sqrt(dot(*this));
	}

	double lInfNorm()
	{
		SphericalScalarFunction zero([](const SphereCoord& x) {return 0.0; });
		return LInfDifference(data, &zero, consts , &system);
	}

	double dot(const LaplaceParticle& particle) const
	{
		return L2InnerProductDiscrete(data, particle.data , 1.0 , consts);
	}

	static double dot(const LaplaceParticle& particle, const LaplaceParticle& qarticle)
	{
		return particle.dot(qarticle);
	}

	int getN() const
	{
		return numSeriesTerms;
	}

	void axpy(LaplaceParticle* particle, double scale)
	{

		//double check1 = validateData();
		//double check2 = particle->validateData();

		for (int t = 0; t < (int)consts.NUMGLNODES; t++)
			for (int p = 0; p < (int)consts.NUMTRAPNODES; p++)
				data[t][p] += particle->data[t][p] * scale;


		fourierdata.axpy(particle->fourierdata, scale);



		potential.solve(fourierdata);

		//if (validateData() > 1.0001 * (check1 + std::abs(scale) * check2))
			//std::cout << " axpy operator is invalid.";

	}



	void refreshData()
	{
		fourierdata.approximate(data, numSeriesTerms);


		potential.solve(fourierdata);
	}


};

class LaplaceNeumannOperator
{
public:
	LaplaceParticle operator *(const LaplaceParticle& particle) const
	{
		LaplaceParticle temp(particle.numSeriesTerms , particle.system);

		for (int t = 0; t < particle.consts.NUMGLNODES; t++)
			for (int p = 0; p < particle.consts.NUMTRAPNODES; p++)
			{
				SurfaceCoord s(particle.consts.GLnodes[t], 2.0 * particle.consts.PI * (double)p / (double)particle.consts.NUMTRAPNODES);
				SphereCoord x(s, &particle.system);

				temp.data[t][p] = -particle.data[t][p] + particle.traction(x);
			}

		temp.refreshData();

		return temp;
	}


};


class LaplaceQuinckeOperator
{

private:
	double permitivity;
public:

	LaplaceQuinckeOperator(const double & eps): permitivity(eps){}

	LaplaceParticle operator *(const LaplaceParticle& particle) const
	{
		LaplaceParticle temp(particle.numSeriesTerms , particle.system);

		temp.setPermitivity(particle.permitivity);
		double lambda = (particle.permitivity - permitivity) / (particle.permitivity + permitivity);
		for (int t = 0; t < particle.consts.NUMGLNODES; t++)
			for (int p = 0; p < particle.consts.NUMTRAPNODES; p++)
			{
				SurfaceCoord s(particle.consts.GLnodes[t], 2.0 * particle.consts.PI * (double)p / (double)particle.consts.NUMTRAPNODES);
				SphereCoord x(s, &particle.system);
				
				temp.data[t][p] = 0.5 * (1.0 + lambda) * particle.data[t][p] + lambda * particle.traction(x);
			}

		temp.refreshData();

		return temp;
	}


};
class LaplaceIdentityPreconditioner
{
public:
	LaplaceParticle solve(const LaplaceParticle& particle)
	{
		return particle;
	}
};

class PointChargePotential :public SphericalScalarFunction
{
private:
	RectCoord center;
	double charge;
	const double PI = 4.0 * atan(1.0);

public:

	PointChargePotential(double q, RectCoord x0) : charge(q), center(x0) {}

	double operator()(const SphereCoord& x) const
	{
		RectCoord r = x - center;
		double R = norm(r);

		return charge * 0.25 / PI  / R;
	}
};


class PointChargeNormalDerivative : public SphericalScalarFunction
{
private:
	RectCoord center;
	double charge;
	const double PI = 4.0 * atan(1.0);

public:

	PointChargeNormalDerivative(double q , RectCoord x0 ) : charge(q) , center(x0){}

	double operator()(const SphereCoord& x) const
	{
		RectCoord r = x - center;
		RectCoord er = e_r(x);
		double R = norm(r);

		return -charge* 0.25 / PI * dot(r, er) / pow(R, 3);
	}
};


LaplaceParticle SolveNeumannSingleSphere(double q , RectCoord x0, int seriesSize)
{

	LaplaceNeumannOperator L;

	const double PI = 4.0 * atan(1.0);

	QuadConstants consts(2 * seriesSize + 1, seriesSize + 1);


	PointChargePotential trueSolution(q, x0);
	PointChargeNormalDerivative r(q, x0);

	LaplaceParticle rh(discretize(&r , &defaultsystem , consts), seriesSize, defaultsystem);
	

	LaplaceIdentityPreconditioner I;
	LaplaceParticle solution = -rh;

	GMRES(&L, &solution, &rh, &I, 30, 1, 1e-10);


	SphereCoordSystem awaysystem = SphereCoordSystem(RectCoord(3.0, 0.0, 0.0));
	LaplaceParticle solutionEvaluatedAway(discretize(&trueSolution, &awaysystem, consts),
		                                  seriesSize,
		                                  RectCoord(3.0, 0.0, 0.0));

	LaplaceParticle approxSolutionEvaluatedAway(discretize(&solution.potential, &awaysystem, consts),
		                                        seriesSize,
		                                        RectCoord(3.0, 0.0, 0.0));

	std::cout << "Norm of solution: " << solutionEvaluatedAway.norm() << "\n";
	std::cout << "Norm of approximation: " << approxSolutionEvaluatedAway.norm() << "\n";
	std::cout << "Error in BIE method: " << (solutionEvaluatedAway - approxSolutionEvaluatedAway).norm() << std::endl;
	return solution;

}

LaplaceParticle SolveNeumannKnown(int seriesSize)
{

	LaplaceNeumannOperator L;

	const double PI = 4.0 * atan(1.0);

	QuadConstants consts(2 * seriesSize + 1,  seriesSize + 1);


	LaplaceParticle trueSolution(seriesSize);
	
	trueSolution.generateRandomSmooth();

	LaplaceParticle rh = L * trueSolution;


	LaplaceIdentityPreconditioner I;
	LaplaceParticle solution = -rh;

	GMRES(&L, &solution, &rh, &I, 30, 1, 1e-10);


	SphereCoordSystem awaysystem = SphereCoordSystem(RectCoord(3.0, 0.0, 0.0));

	LaplaceParticle solutionEvaluatedAway(discretize(&trueSolution.potential, &awaysystem, consts),
		seriesSize,
		RectCoord(3.0, 0.0, 0.0));

	LaplaceParticle approxSolutionEvaluatedAway(discretize(&solution.potential,&awaysystem, consts),
		seriesSize,
		RectCoord(3.0, 0.0, 0.0));


	std::cout << "Error in BIE method: " << LInfDifference(solution.data , trueSolution.data,solution.consts) << std::endl;
	return solution;

}


LaplaceParticle SolveQuinckeKnown(int seriesSize , double l)
{

	LaplaceQuinckeOperator L(l);

	const double PI = 4.0 * atan(1.0);

	QuadConstants consts(2 * seriesSize + 1, seriesSize + 1);


	LaplaceParticle trueSolution(seriesSize);
	
	trueSolution.generateRandomSmooth();

	LaplaceParticle rh = L * trueSolution;


	LaplaceIdentityPreconditioner I;
	LaplaceParticle solution = -rh;

	GMRES(&L, &solution, &rh, &I, 30, 1, 1e-10);


	SphereCoordSystem awaysystem = SphereCoordSystem(RectCoord(3.0, 0.0, 0.0));

	LaplaceParticle solutionEvaluatedAway(discretize(&trueSolution.potential, &awaysystem, consts),
		seriesSize,
		RectCoord(3.0, 0.0, 0.0));

	LaplaceParticle approxSolutionEvaluatedAway(discretize(&solution.potential, &awaysystem, consts),
		seriesSize,
		RectCoord(3.0, 0.0, 0.0));

		std::cout << "Error in BIE method (L2 ): " << (solution-trueSolution).norm() << std::endl;

	std::cout << "Error in BIE method (L Infinity): " << LInfDifference(solution.data , trueSolution.data,solution.consts) << std::endl;
	return solution;

}

