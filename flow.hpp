#ifndef INCLUDE_FLOW
#define INCLUDE_FLOW

#include "math.hpp"
#include <vector>
#include <algorithm>
#include <random>

static float PI = 3.1415926535897932384626433832795f;

// ============================================================
enum EIntegrator {
	Integrator_Euler, Integrator_RK2, Integrator_RK4
};

// ============================================================
// base class for fields (every field has a sample function)
template<typename vin, typename vout>
class field
{
public:
	virtual vout sample(const vin& pos) const = 0;
};

// ============================================================
// vector field is a map from vec -> vec
template<typename vec>
class vfield : public field<vec, vec>
{
public:
	// takes an Euler step
	vec stepEuler(const vec& pos, float dt) const {
		// TODO: implement!
		return pos + sample(pos) * dt;
	}

	// takes a 2nd-order Runge-Kutta step
	vec stepRK2(const vec& pos, float dt) const {
		// TODO: implement!
		vec v1(sample(pos));
		vec v2(sample(pos + v1 * 0.5 * dt));
		return pos + v2 * dt;
	}

	// takes a 4th-order Runge-Kutta step
	vec stepRK4(const vec& pos, float dt) const {
		// TODO: implement!
		vec v1(sample(pos));
		vec v2(sample(pos + v1 * 0.5 * dt));
		vec v3(sample(pos + v2 * 0.5 * dt));
		vec v4(sample(pos + v3 * dt));
		return pos + (v1 / 6 + v2 / 3 + v3 / 3 + v4 / 6) * dt;
	}

	// takes a step with the selected integrator
	vec step(const vec& pos, float dt, EIntegrator integrator) const {
		switch (integrator)
		{
		case Integrator_Euler:
			return stepEuler(pos, dt);
		case Integrator_RK2:
			return stepRK2(pos, dt);
		case Integrator_RK4:
		default:
			return stepRK4(pos, dt);
		}
	}
	
	// traces a tangent curve and returns all particle positions
	void tangentCurve(vec pos, float dt, float tau, EIntegrator integrator, std::vector<vec>& curve) const
	{
		curve.clear();
		curve.push_back(pos);
		tau = std::abs(tau);
		while (tau > std::abs(dt))
		{
			pos = step(pos, dt, integrator);
			curve.push_back(pos);
			tau -= std::abs(dt);
		}
		if (dt < 0)
			tau = -tau;
		pos = step(pos, tau, integrator);
		curve.push_back(pos);
		// TODO: repeatedly call step() and store the intermediate locations in the std::vector 'curve'.
	}
};

// ============================================================
// vector field in 2D
class vfield2d : public vfield<vec2f>
{
};

#endif // !INCLUDE_FLOW
