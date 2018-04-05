#ifndef TRACE_HPP
#define TRACE_HPP

#include <cmath>
#include <vector>
#include <cstdlib>
#include <stdexcept>

#define ONEDAY (24*3600)
#define ONEYEAR (365.25*ONEDAY)
#define ONEHOUR 3600


class Trace
{
public:
	Trace(double mu, double downtime, double length = ONEYEAR) : 
		mu(mu), i(0)
	{
		p = 0.0;
		nF = 0;
		while(p < length)
		{
			double u = (double)rand() / (double)RAND_MAX;
			double X = -log(u)*mu;
			errors.push_back(p+downtime/ONEYEAR+X);
			p += downtime/ONEYEAR+X;
			nF++;
		}
		horizon = length;
	}

	double next(double t)
	{
		while(t >= errors[i] && i < errors.size())
			i++;
		if(i >= errors.size())
		{
		/*	while (p < t)
			{
				double u = (double)rand() / (double)RAND_MAX;
				double X = -log(u)*mu;
				errors.push_back(p+X);
				p += X;
				nF++;
				i++;
			}
			//if (p > horizon)
			//	throw std::runtime_error("ERROR: application too long, increasing trace length");*/
			return horizon+0.1;
		}
		return errors[i];
	}

	double mu;
	unsigned int nF;
	std::vector<double> errors;

private:
	unsigned int i;
	double horizon;
	double p;
};

#endif
