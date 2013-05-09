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
		for(Atom::bond_type::iterator it = bonds_i.begin();
				it != bonds_i.end(); it++ )
		{
			k = *it;
			if(k->get_id() != j->get_id())
			{
				r_ik = k->r - i->r;
				sum += S_down_N(gmtl::length(r_ik));
			}
		}

		return 1.0/std::sqrt(1 + sum);
	}
	// (10)
	float LCBOPII::N_ijk(Atom *i, Atom *j, Atom *k)
	{
		Atom::position_type r_ij, r_ik;
		r_ij = j->r - i->r;
		r_ik = k->r - i->r;

		return N(i) - S_down_N(gmtl::length(r_ij)) - S_down_N(gmtl::length(r_ik));
	}
	// (11)
	float LCBOPII::N(Atom *i)
	{
		float sum = 0.0;
		Atom * j;
		Atom::bond_type bonds_i = i->get_bonds();
		Atom::position_type r_ij;
		for(Atom::bond_type::iterator it = bonds_i.begin();
				it != bonds_i.end(); it++ )
		{
			j = *it;
			r_ij = j->r - i->r;
			sum += S_down_N(gmtl::length(r_ij));
		}

		return sum;
	}
	// (12)
	float LCBOPII::G(Atom *i, Atom *j, Atom *k)
	{
		Atom::position_type r_ij = j->r - i->r;
		Atom::position_type r_ik = k->r - i->r;

		//y = cos_theta_ijk
		float y = gmtl::dot(r_ij, r_ik)/(gmtl::length(r_ij)*gmtl::length(r_ik));

		//z = N_ijk
		float z = N_ijk(i,j,k);

		float y_zero = y_0(z);

		return THETA(y_zero - y)*G_1(y) + THETA(y - y_zero)*G_2(y, z, y_zero);
	}
	// (13)
	float LCBOPII::y_0(float z)
	{
		static const float A_y0 = -0.4;
		static const float B_y0 = 0.01875;

		return A_y0 + B_y0*z*(1.0 + z);
	}
	// (14)
	float LCBOPII::G_1(float y)
	{
		static const float one_third = 1.0/3.0;

		static const float g_1_n[] = {
										 0.7233666272,
										 1.7334665088,
										 1.8701997632
									    };

		static const float g_2_n[] = {
								   	     0.73994527795,
								   	     -1.999211817,
								   	     -17.43251545,
								   	   	 -33.96127110,
								   	   	 -44.65392079
								 	 	};

		static const float g_3_n[] = {
										 -15.19,
										 -25.6168552398,
										 -21.51728397,
										 0.9899080993,
										 13.66416160
									    };

		float ret = 0.0;

		if(y >= -1 and y < 0.5 )
		{

		}
		else if(y >= 0.5 and y < one_third )
		{

		}
		else
		{

		}

		return ret;
	}
	// (15)
	float LCBOPII::G_2(float y, float z, float y0)
	{

	}
	// (16)
	float LCBOPII::g_z_max(float y0)
	{

	}
	// (17)
	float LCBOPII::g_z_2(float z)
	{

	}
	// (18)
	float LCBOPII::g_z_1(float y0, float g_z2)
	{

	}
	// (19)
	float LCBOPII::g_z_0(float y0, float g_zMax,
							float g_z1, float g_z2)
	{

	}
	// (20)
	float LCBOPII::H(Atom *i, Atom *j, Atom *k)
	{

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
