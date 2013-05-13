#include <cmath>
#include <gmtl/gmtl.h>
#include "atom.h"
#include "lcbopii.h"

namespace simul
{
	//S_up / S_down parameters from TABLE I.
	double LCBOPII::S_down_sr(double q)
	{
		return S_down(x_from_q(q, 1.7, 2.2), 3.0);
	}
	double LCBOPII::S_down_lr(double q)
	{
		return S_down(x_from_q(q, 5.5, 6.0), 0.0);
	}
	double LCBOPII::S_down_db(double q)
	{
		return S_down(x_from_q(q, 0.0, 1.0), 0.0);
	}

	double LCBOPII::S_up_mr(double q)
	{
		return S_up(x_from_q(q, 1.7, 2.2), -2.0);
	}
	double LCBOPII::S_up_M(double q)
	{
		return S_up(x_from_q(q, 2.0, 3.0), 0.0);
	}
	double LCBOPII::S_up_gamma_0(double q)
	{
		return S_up(x_from_q(q, 0.34, 0.93), 0.0);
	}

	double LCBOPII::S_down_N(double q)
	{
		return S_down(x_from_q(q, 1.7, 2.2), -3.0);
	}
	double LCBOPII::S_down_sat(double q)
	{
		return S_down(x_from_q(q, 3.0, 4.0), 0.0);
	}
	double LCBOPII::S_up_gamma_2(double q)
	{
		return S_up(x_from_q(q, 0.30, 0.93), 0.0);
	}

	// TABLE I. parameters of LCBOPII (p. 4)
	// Short range potential
	const double LCBOPII::A_sr = 53026.92614;
	const double LCBOPII::alpha = 6.74750993;

	const double LCBOPII::B_sr_1 = 27618.35706;
	const double LCBOPII::beta_1 = 6.34503890;
	const double LCBOPII::B_sr_2 = 34.07142502;
	const double LCBOPII::beta_2 = 1.19712839;

	const double LCBOPII::A_y0 = -0.4;
	const double LCBOPII::B_y0 = 0.01875;

	const double LCBOPII::one_third = 1.0/3.0;

	const double LCBOPII::g_min = 0.0020588719;
	//const double LCBOPII::g_gr = 0.0831047003;
	/**
	 * Wartosc g_gr dofitiowana do csplines
	 */
	const double LCBOPII::g_gr = 0.0831049;
	const double LCBOPII::g_max = 16.0;

	const double LCBOPII::g_1_n[] = {
									 0.7233666272,
									 1.7334665088,
									 1.8701997632
								    };

	/*const double LCBOPII::g_2_n[] = {
							   	     0.73994527795,
							   	     -1.999211817,
							   	     -17.43251545,
							   	   	 -33.96127110,
							   	   	 -44.65392079
							 	 	};*/

	/**
	 *  Parametry dofitowane dla obszaru -0.5 <= y < -1/3
	 *  przy pomocoy csplines dla starych parametrow funkcja
	 *  byla nieciagla na ww. obszarze
	 */
	const double LCBOPII::g_2_n[] = {
							   	     1.2313,
							   	     3.1061,
							   	     2.26709,
							   	     -0.527691,
							   	     -0.337354
							 	 	};

	const double LCBOPII::g_3_n[] = {
									 -15.19,
									 -25.6168552398,
									 -21.51728397,
									 0.9899080993,
									 13.66416160
								    };

	const double LCBOPII::A_g = 5.6304664723;
	const double LCBOPII::B_g = 0.1516943990;
	const double LCBOPII::C_g = 0.009832975891;

	const double LCBOPII::D_g = -0.189175977654;
	const double LCBOPII::E_g = 0.050977653631;

	/**
	 * B. Short range potential
	 */
	// (5)
	double LCBOPII::V_sr(Atom *i, Atom *j)
	{
		Atom::position_type r_ij = j->r - i->r;
		double r = gmtl::length(r_ij);

		return V_sr_R(r) - B(i, j)*V_sr_A(r);
	}
	// (6)
	double LCBOPII::V_sr_R(double r)
	{
		return A_sr*std::exp(-alpha*r);
	}
	// (7)
	double LCBOPII::V_sr_A(double r)
	{
		return B_sr_1*std::exp(-beta_1*r) + B_sr_2*std::exp(-beta_2*r);
	}

	// (8)
	double LCBOPII::B(Atom *i, Atom *j)
	{
		return 0.5*(b(i,j) + b(j,i)) + F_conj(i,j) + A(i,j) + T(i,j);
	}
	// (9)
	double LCBOPII::b(Atom *i, Atom *j)
	{
		double sum = 0.0;

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
	double LCBOPII::N_ijk(Atom *i, Atom *j, Atom *k)
	{
		Atom::position_type r_ij, r_ik;
		r_ij = j->r - i->r;
		r_ik = k->r - i->r;

		return N(i) - S_down_N(gmtl::length(r_ij)) - S_down_N(gmtl::length(r_ik));
	}
	// (11)
	double LCBOPII::N(Atom *i)
	{
		double sum = 0.0;
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
	double LCBOPII::G(Atom *i, Atom *j, Atom *k)
	{
		Atom::position_type r_ij = j->r - i->r;
		Atom::position_type r_ik = k->r - i->r;

		//y = cos_theta_ijk
		double y = gmtl::dot(r_ij, r_ik)/(gmtl::length(r_ij)*gmtl::length(r_ik));

		//z = N_ijk
		double z = N_ijk(i,j,k);

		return G(y,z);
	}
	// (12)
	double LCBOPII::G(double y, double z)
	{
		double y_zero = y_0(z);
		std::cout << "#Y0 " << y_zero << std::endl;
		std::cout << "#THETA(y_zero - y) G_1=> " << THETA(y_zero - y) << std::endl;
		std::cout << "#THETA(y - y_zero) G_2=>" << THETA(y - y_zero) << std::endl;
		return THETA(y_zero - y)*G_1(y) + THETA(y - y_zero)*G_2(y, z, y_zero);
	}
	// (13)
	double LCBOPII::y_0(double z)
	{
		return A_y0 + B_y0*z*(1.0 + z);
	}
	// (14)
	double LCBOPII::G_1(double y)
	{
		double ret = 0.0;

		if(y >= -1 and y < -0.5 )
		{
			for(int i=0; i<=2; i++)
			{
				ret += g_1_n[i]*std::pow(y,i);
			}
			ret *= std::pow(y + 1, 2);
			ret += g_min;
		}
		else if(y >= -0.5 and y < -one_third )
		{
			for(int i=0; i<=4; i++)
			{
				ret += (y + 0.5)*g_2_n[i]*std::pow(y,i);
			}
			ret += g_gr;
		}
		else
		{
			for(int i=0; i<=4; i++)
			{
				ret += g_3_n[i]*std::pow(y,i);
			}
			ret *= std::pow(y - 1, 2);
			ret += g_max;
		}

		return ret;
	}

	double LCBOPII::G_1_prim(double y0)
	{
		double G_1prim = 0.0;
		if(y0 >= -1 and y0 < -0.5 )
		{
			G_1prim = 2*(y0 + 1)*g_1_n[0];
			for(int i=1; i <= 2; i++)
			{
				G_1prim = (
								2*(y0 + 1)*std::pow(y0, i) +
								i*std::pow(y0 + 1, 2)*std::pow(y0, i - 1)
						   )*g_1_n[i];
			}
		}
		else if(y0 >= -0.5 and y0 < -one_third )
		{
			G_1prim = g_2_n[0];
			for(int i=1; i <= 4; i++)
			{
				G_1prim += ((y0 + 0.5)*i*std::pow(y0, i-1) + std::pow(y0, i))*g_2_n[i];
			}
		}
		else
		{
			G_1prim = 2*(y0 - 1)*g_3_n[0];
			for(int i=1; i <= 4; i++)
			{
				G_1prim += (
								2*(y0 - 1)*std::pow(y0, i) +
								std::pow(y0 - 1, 2)*i*std::pow(y0, i - 1)
							)*g_3_n[i];
			}
		}

		return G_1prim;
	}
	// (15)
	double LCBOPII::G_2(double y, double z, double y0)
	{
		double g_zMax = g_z_max(z, y0);
		std::cout << "#g_z_max " << g_zMax << std::endl;

		double g_z2 = g_z_2(z);
		//double g_z2 = 0.0;
		std::cout << "#g_z_2 " << g_z2 << std::endl;

		double g_z1 = g_z_1(y0, g_zMax, g_z2);
		std::cout << "#g_z_1 " << g_z1 << std::endl;

		double g_z0 = g_z_0(y0, g_zMax, g_z1, g_z2);
		std::cout << "#g_z_0 " << g_z0 << std::endl;

		return  g_zMax + std::pow((1-y), 2)*(
					g_z0 +
					g_z1*y +
					g_z2*std::pow(y , 2)
				);
	}
	double LCBOPII::G_2_prim(double y, double z, double y0)
	{
		double g_zMax = g_z_max(z, y0);

		double g_z2 = g_z_2(z);

		double g_z1 = g_z_1(y0, g_zMax, g_z2);

		double g_z0 = g_z_0(y0, g_zMax, g_z1, g_z2);

		double ret = -2*(1 - y)*g_z0;
		for(int i=1; i<=2; i++)
		{
			switch(i)
			{
			case 1:
				ret += (-2*(1 - y)*y + std::pow(1-y, 2))*g_z1;
				break;
			case 2:
				ret += (-2*(1 - y)*y*y + 2*std::pow(1-y, 2)*y)*g_z2;
				break;
			}
		}

		return ret;
	}
	// (16)
	double LCBOPII::g_z_max(double z, double y0)
	{
		return g_max
				- (A_g + B_g*z + C_g*std::pow(z, 4))*std::pow((1 - y0), 3);
	}
	// (17)
	double LCBOPII::g_z_2(double z)
	{
		double z4 = std::pow(z, 4);
		return (D_g*z4)/(1 + E_g*z4);
	}
	// (18)
	double LCBOPII::g_z_1(double y0, double g_zMax, double g_z2)
	{
		double G_1prim = G_1_prim(y0);

		std::cout << "#G_1_prim " << G_1prim << std::endl;

		// Funkcja poprawiona przez zemnie w drugim wyrazie
		// w liczniku powinno byc [G_1 - g_z_max] !!!!
		// orginalnie jest G_1
		return G_1prim/std::pow(y0 - 1, 2) -
				2 * (G_1(y0) - g_zMax)/std::pow(y0 - 1, 3) -
				2 * g_z2 * y0;
	}
	// (19)
	double LCBOPII::g_z_0(double y0, double g_zMax,
							double g_z1, double g_z2)
	{
		std::cout << "# G_1(y0)" << G_1(y0) << std::endl;
		return (G_1(y0) - g_zMax)/std::pow(y0 - 1, 2) -
				g_z1*y0 -
				g_z2*std::pow(y0, 2);
	}

	// (20)
	double LCBOPII::H(Atom *i, Atom *j, Atom *k)
	{

	}

	// (22)
	double LCBOPII::F_conj(Atom *i, Atom *j)
	{
		return 0.0;
	}
	// (32)
	double LCBOPII::A(Atom *i, Atom *j)
	{
		return 0.0;
	}
	// (35)
	double LCBOPII::T(Atom *i, Atom *j)
	{
		return 0.0;
	}
}
