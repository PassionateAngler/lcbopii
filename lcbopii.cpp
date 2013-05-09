#include <cmath>
#include <gmtl/gmtl.h>
#include "atom.h"
#include "lcbopii.h"

namespace simul
{
	//S_up / S_down parameters from TABLE I.
	float LCBOPII::S_down_sr(float q)
	{
		return S_down(x_from_q(q, 1.7, 2.2), 3.0);
	}
	float LCBOPII::S_down_lr(float q)
	{
		return S_down(x_from_q(q, 5.5, 6.0), 0.0);
	}
	float LCBOPII::S_down_db(float q)
	{
		return S_down(x_from_q(q, 0.0, 1.0), 0.0);
	}

	float LCBOPII::S_up_mr(float q)
	{
		return S_up(x_from_q(q, 1.7, 2.2), -2.0);
	}
	float LCBOPII::S_up_M(float q)
	{
		return S_up(x_from_q(q, 2.0, 3.0), 0.0);
	}
	float LCBOPII::S_up_gamma_0(float q)
	{
		return S_up(x_from_q(q, 0.34, 0.93), 0.0);
	}

	float LCBOPII::S_down_N(float q)
	{
		return S_down(x_from_q(q, 1.7, 2.2), -3.0);
	}
	float LCBOPII::S_down_sat(float q)
	{
		return S_down(x_from_q(q, 3.0, 4.0), 0.0);
	}
	float LCBOPII::S_up_gamma_2(float q)
	{
		return S_up(x_from_q(q, 0.30, 0.93), 0.0);
	}

	/**
	 * B. Short range potential
	 */
	// (5)
	float LCBOPII::V_sr(Atom *i, Atom *j)
	{
		Atom::position_type r_ij = j->r - i->r;
		float r = gmtl::length(r_ij);

		return V_sr_R(r) - B(i, j)*V_sr_A(r);
	}
	// (6)
	float LCBOPII::V_sr_R(float r)
	{
		static const float A_sr = 53026.92614;
		static const float alpha = 6.74750993;

		return A_sr*std::exp(-alpha*r);
	}
	// (7)
	float LCBOPII::V_sr_A(float r)
	{
		static const float B_sr_1 = 27618.35706;
		static const float beta_1 = 6.34503890;
		static const float B_sr_2 = 34.07142502;
		static const float beta_2 = 1.19712839;

		return B_sr_1*std::exp(-beta_1*r) + B_sr_2*std::exp(-beta_2*r);
	}

	// (8)
	float LCBOPII::B(Atom *i, Atom *j)
	{
		return 0.5*(b(i,j) + b(j,i)) + F_conj(i,j) + A(i,j) + T(i,j);
	}
	// (9)
	float LCBOPII::b(Atom *i, Atom *j)
	{
		float sum = 0.0;

		Atom * k;
		Atom::bond_type bonds_i = i->get_bonds();
		Atom::position_type r_ik;

		//Loop over all neighbours of i != j
		for(Atom::bond_type::iterator it_k = bonds_i.begin();
				it_k != bonds_i.end(); it_k++ )
		{
			k = *it_k;
			if(k->get_id() != j->get_id())
			{
				r_ik = k->r - i->r;
				sum += S_down_N(gmtl::length(r_ik));
			}
		}

		return 1.0/std::sqrt(1 + sum);
	}
	// (22)
	float LCBOPII::F_conj(Atom *i, Atom *j)
	{
		return 0.0;
	}
	// (32)
	float LCBOPII::A(Atom *i, Atom *j)
	{
		return 0.0;
	}
	// (35)
	float LCBOPII::T(Atom *i, Atom *j)
	{
		return 0.0;
	}
}
