#include "LaplaceManyParticles.h"

void testSingleLayer(int p = 8);
void testTraction(int p = 8);

int main()
{

	std::vector<RectCoord> centers;

	

	RectCoord E(10.0, 10.0, 0.0);

	std::array<int, 6> ps = {4 , 6 , 8 , 10 , 12 , 14};
	std::array<double, 8> rs = { 2.001 , 2.005, 2.01 , 2.05 , 2.1 , 2.5, 3.0 , 4.0 };
	std::array<double, 6> errs;
/*
	for (int p : ps)
	{

		std::cout << "Using " << p << "th order Spherical Harmonics and " << 2 * p + 1 << " quadrature nodes in each direction.\n";
		SolveQuinckeManySphereFromKnown(centers, p , 0.5);
		std::cout << "\n";
	}
	*/
/*
	for (int p : ps)
	{
     	std::cout << "Using " << p << "th order Spherical Harmonics and " << 2 * p + 1 << " quadrature nodes in each direction.\n";
		SolveQuinckeKnown( p , 0.5);
		std::cout << "\n";

	} */
/* for (int p : ps)
	{
     	std::cout << "Using " << p << "th order Spherical Harmonics and " << 2 * p + 1 << " quadrature nodes in each direction.\n";
		SolveNeumannSingleSphere( 5.0 ,RectCoord(0.0 , 0.0 , 3.0), p);
		std::cout << "\n";

	}

	*/
	/*
 for (int p : ps)
	{
     	std::cout << "Using " << p << "th order Spherical Harmonics and " << 2 * p + 1 << " quadrature nodes in each direction.\n";
		SolveNeumannKnown(p);
		std::cout << "\n";

	}
	*/
	for (double r : rs)
	{
		centers.push_back(RectCoord(r, 0.0, 0.0));
		centers.push_back(RectCoord(-3.0, 0.0, 1.5));
		centers.push_back(RectCoord(0.0, 4.0, 2.0));
	    centers.push_back(RectCoord(0.0, -5.0, -1.25));
		centers.push_back(RectCoord(0.0, 0.0, 0.0));
	int maxp = 16;
	LaplaceParticleSystem comparison = SolveQuinckeRandom(centers, maxp, 5.0, 1.0, E);

	double compnorm = comparison.LInfNorm();
	int i = 0;
	
		for (int p : ps)
		{
		
			std::cout << "Using " << p << "th order Spherical Harmonics and " << 2 * p + 1 << " quadrature nodes in each direction.\n";
			LaplaceParticleSystem temp = SolveQuinckeRandom(centers, p, 5.0, 1.0, E);

			LaplaceParticleSystem currentsolution(centers, maxp);
			for (int i = 0; i < temp.particles.size(); i++)
			{
				ScalSphereData data = discretize(&temp.particles[i].fourierdata,
					&temp.particles[i].system,
					currentsolution.particles[i].consts);
				currentsolution.particles[i] = LaplaceParticle(data, maxp, temp.particles[i].system);
			}
			errs[i] = (currentsolution - comparison).LInfNorm() / compnorm;
			std::cout << "Error for p = " << p << ": " << errs[i++];
			std::cout << "\n";
		}

		std::cout << "Norm of solution: " << compnorm << std::endl;
		std::cout << "errors : ";
		for (int j = 0; j < 6; j++)
			std::cout << errs[j] << ", ";

		std::cout << std::endl;
		centers.clear();
	}
	return 0;
}

SphericalScalarFunction zero([](const SphereCoord& x) {return 0.0; });


void testSingleLayer(int n)
{
	LaplaceParticle p(n);

	p.generateRandomSmooth();

	std::array<double, 6> rs = { 0.1 , 0.5 , 1.0 , 2.0 , 4.0, 8.0 };

	for (double r : rs)
	{
		SphereCoordSystem system(RectCoord(2.0 + r, 0.0, 0.0));
		std::cout << "r = " << r << ".Relative max error in Single Layer between methods: "
			<< LInfDifference(&p.potential
				, &p.quadpotential
				, p.consts
			    , &system)/
			LInfDifference(&p.potential
				, &zero
				, p.consts
				,&system) << "\n";
	}
}

void testTraction(int n)
{
	LaplaceParticle p(n);

	p.generateRandomData();

	std::array<double, 6> rs = { 0.1 , 0.5 , 1.0 , 2.0 , 4.0, 8.0 };

	for (double r : rs)
	{
		SphereCoordSystem system(RectCoord(2.0 + r, 0.0, 0.0));
		std::cout << "r = " << r << ".Relative max error in Traction between methods: "
			<< LInfDifference(&p.traction
				, &p.quadtraction
				, p.consts
				, &system) /
			LInfDifference(&p.traction
				, &zero
				, p.consts
				, &system) << "\n";
	}
}