#include <cmath>
#include <cstdint>
#include <cfloat>
#include <cstdio>
#include <Eigen/Dense>
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
	const double LCBOPII::g_gr = 0.0831047003;
	/**
	 * Wartosc g_gr dofitiowana do csplines
	 */
	//const double LCBOPII::g_gr = 0.0831049;

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
	 *	parameters for H(x) function
	 */
	const double LCBOPII::d = 0.14;
	const double LCBOPII::C_1 = 3.335;
	const double LCBOPII::C_4 = 220.0;
	const double LCBOPII::C_6 = -(std::pow(C_1, 2) + 12*C_4*std::pow(d, 2))
									/(30*std::pow(d, 4));

	const double LCBOPII::R_0 = LCBOPII::H_2(d);

	const double LCBOPII::R_1 = C_1 + std::pow(C_1, 2)*d
									+ 4*C_4*std::pow(d,3) + 6*C_6*std::pow(d,5);

	const double LCBOPII::L = R_0 - 2*C_1*d;
	const double LCBOPII::kappa = (2*C_1 - R_1)/L;

	/**
	 * parameters matrices for
	 * F(N_ij, N_ji, N_ij_conj) function
	 */
	const double LCBOPII::F_ij_0[][4] = {
			{0.0, 0.0207, -0.0046, -0.1278},
			{0.0207, 0.0, -0.0365, -0.1043},
			{-0.0046, -0.0365, 0.0, -0.0273},
			{-0.1278, -0.1043, -0.0273, 0.0}
	};

	const double LCBOPII::F_ij_1[][4] = {
			{0.0, 0.0584, 0.0416, -0.1278},
			{0.0584, 0.1379, 0.0062, -0.1243},
			{0.0416, 0.0062, 0.0936, -0.0393},
			{-0.1278, -0.1243, -0.0393, 0.0}
	};

	/**
	 * parameters for A_ij
	 */
	const double LCBOPII::alpha_0 = 0.95;

	/**
	 * parameters for T_ij
	 */
	const double LCBOPII::A_t  = -13.152909887;
	const double LCBOPII::B_t1 = -0.0486839616;
	const double LCBOPII::B_t2 = 3.8;
	const double LCBOPII::B_t3 = 0.62;
	const double LCBOPII::B_t4 = 0.005;

	/**
	 * Long range potential V_lr parameters
	 */
	const double LCBOPII::r_0 = 3.715735;
	const double LCBOPII::eps_2 = 0.002827918;
	const double LCBOPII::lambda_1 = 1.338162;
	const double LCBOPII::lambda_2 = 2.260479;
	const double LCBOPII::eps_1 = eps_2*std::pow(lambda_2/lambda_1, 2);
	const double LCBOPII::v_1 = eps_1 - eps_2;

	/**
	 * Middle - range potential V_mr parameters
	 */
	const double LCBOPII::r_1 = 4.0;
	const double LCBOPII::r_2 = 2.9;
	const double LCBOPII::A_mr_0 = -0.2345;
	const double LCBOPII::A_mr_1 = -0.67;
	const double LCBOPII::A_mr_2 = -4.94;

	/**
	 * Pair bonding energy Eb_ij (1)
	 */
	double LCBOPII::Eb_ij(Atom * i, Atom * j)
	{
		Atom::position_type r_ij = j->r - i->r;
		double r = r_ij.norm();
		double ret = 0.0;

		double S_d_sr = S_down_sr(r);
		double Bij =  B(i, j);

		// r < 1.7
		if(S_d_sr == 1.0){
			ret = V_sr(i, j, Bij);
			std::cout << std::endl;
			return ret;
		}

		// 1.7 < r  < 2.2
		if(S_d_sr > 0)
		{
			ret += S_down_sr(r)*V_sr(i,j, Bij);
		}

		// r > 1.7
		if(S_d_sr < 1) // S_up_sr > 0
		{
			std::cout << "LR" << std::endl;
			ret += (1 - S_d_sr)*V_lr(r);
			ret += S_up_mr(r)*V_mr(i,j, Bij)/Z_mr_i(i,j);
		}

		std::cout << std::endl;
		return ret;
	}

	/**
	 * B. Short range potential
	 */
	// (5)
	double LCBOPII::V_sr(Atom *i, Atom *j, double Bij)
	{
		Atom::position_type r_ij = j->r - i->r;
		double r = r_ij.norm();

		/*
		std::cout << "V_sr Atom i = " << i->get_id() << std::endl;
		std::cout << "V_sr_R(r) " << V_sr_R(r) << std::endl;
		std::cout << "first B(i,j) " << B(i,j) << std::endl;
		std::cout << "V_sr_A(r) " << V_sr_A(r) << std::endl;
		*/

		return V_sr_R(r) - Bij*V_sr_A(r);
		//return r;
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
		return 0.5*(b(i,j) + b(j,i)) + F_A_T(i, j);
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
				std::cout << "H("
								  << i->get_id() << ", "
						  	  	  << j->get_id() << ", "
						  	  	  << k->get_id() <<") "
						  	  	  << H(i, j, k) << std::endl;
				std::cout << "G("
								  << i->get_id() << ", "
						  	  	  << j->get_id() << ", "
						  	  	  << k->get_id() <<") "
						  	  	  << G(i, j, k) << std::endl;
				sum += S_down_N(r_ik.norm())*H(i, j, k)*G(i, j, k);
			}
		}

		sum = 1.0/std::sqrt(1 + sum);
		if(sum > 0){
			std::cout << "b(" << i->get_id() << ", " <<  j->get_id() << ") = " << sum << std::endl;
		}
		return sum;
	}

	double LCBOPII::F_A_T(Atom *i, Atom*j,
			bool enable_F, bool enable_A, bool enable_T)
	{
		double part;
		double ret = 0.0;
		//double weight = 1.0;
		double weight_i, weight_j, weight;
		weight_i = weight_j = 1.0;

		/**
		 *	K, L number of neighbours of "i" without "j" and vice versa
		 */
		const uint32_t K = i->get_bonds().size() - 1;
		const uint32_t L = j->get_bonds().size() - 1;

		const uint32_t MAX_CONFIGURATION = 0x01 << (K + L);
		uint32_t sigmas_k, sigmas_l;

		double Nij_el, Nji_el, Nij_max_el, Nji_max_el, Nij_min_el, Nji_min_el;
		unsigned int Nij_sigma, Nji_sigma;
		double Mij_sigma,  Mji_sigma;
		double Nij_sigma_k_l_conj;

		//std::cout << "K = "<< K << " L = " <<  L << std::endl;
		/**
		 * Higer L bits of "configuration" are neighbours of "j" sigmas
		 * Lower K bits of "configuration" are neighbours of "i" sigmas
		 */
		for(u_int32_t configuration = 0L; configuration < MAX_CONFIGURATION; configuration++)
		{
			//Current sigmas for "i" nn.
			//std::cout << "conf. " << configuration << std::endl;
			sigmas_k = configuration & (~(UINT32_MAX<<K));
			weight_i = W_ij(i, j, sigmas_k);
			//if(weight == 0.0) continue;

			Nij_sigma = N_ij_sigma_k(sigmas_k);
			Mij_sigma = M_ij_sigma_k(i, j, sigmas_k);
			Nij_el = N_ij_el(Nij_sigma, Mij_sigma);
			Nij_min_el = N_ij_min_el(Nij_sigma);
			Nij_max_el = N_ij_max_el(Nij_sigma);
			std::cout << "Nij_sigma " << Nij_sigma << std::endl;
			/*
			std::cout << "Mij_sigma" << Mij_sigma << std::endl;
			std::cout << "Nij_el"
			std::cout << "Nij_min_el"
			std::cout << "Nij_max_el"
			*/

			//Current sigmas for "j" nn.
			sigmas_l = configuration >> K;
			weight_j = W_ij(j, i, sigmas_l);
			std::cout << "sigmas k =  " << sigmas_k << " sigmas l = " << sigmas_l << std::endl;
			std::cout << "weight k " << weight_i << " weight l " << weight_j << std::endl;

			weight = weight_i*weight_j;
			if(weight == 0.0) continue;

			Nji_sigma =  N_ij_sigma_k(sigmas_l);
			Mji_sigma =  M_ij_sigma_k(j, i, sigmas_l);
			Nji_el =     N_ij_el(Nji_sigma, Mji_sigma);
			Nji_min_el = N_ij_min_el(Nji_sigma);
			Nji_max_el = N_ij_max_el(Nji_sigma);
			std::cout << "Nji_sigma " << Nji_sigma << std::endl;
			std::cout << std::endl;

			std::cout << "Mij_sigma = " << Mij_sigma << std::endl;
			std::cout << "Mji_sigma = " << Mji_sigma << std::endl;

            std::cout << "Nij_el = " << Nij_el << " Nji_el = " << Nji_el
            		  << " Nij_min_el = " << Nij_min_el
            		  << " Nji_min_el = " <<  Nji_min_el << std::endl;

            std::cout << "Nij_el + Nji_el - Nij_min_el - Nji_min_el = " <<
            		Nij_el + Nji_el - Nij_min_el - Nji_min_el << std::endl;
			// eq. 26
			Nij_sigma_k_l_conj = (Nij_el + Nji_el - Nij_min_el - Nji_min_el)/
				(Nij_max_el + Nji_max_el - Nij_min_el - Nji_min_el + DBL_MIN);


			// F_conj eq. 22 page 6
			if(enable_F)
			{
				part = weight*F_conj(Nij_sigma, Nji_sigma, Nij_sigma_k_l_conj);
				//part = weight*F_conj(Nij_sigma, Nji_sigma, 1);
				std::cout << "F_con("
						<< Nij_sigma << ", "
						<< Nji_sigma << ", "
						<< Nij_sigma_k_l_conj <<
						") = " << part << std::endl;
				ret += part;
			}

			// A(ntibonding) term eq. 33 page 6
			if(enable_A)
			{
				if(((Nij_sigma == Nji_sigma) and (Nij_sigma == 1 or Nij_sigma == 2))
					or (Nij_sigma == 1 and Nji_sigma == 2)
					or (Nij_sigma == 2 and Nji_sigma == 1))
				{
					part = weight*A(Nij_el - Nji_el);
					std::cout << "A = " << part << std::endl;
					ret += part;
				}
			}

			// T(orsion) term eq. 35 page 7
			if(enable_T)
			{
				if((Nij_sigma == Nji_sigma) and  Nij_sigma == 2)
				{
					part = weight*t_ij(i, sigmas_k, j, sigmas_l,
							           Nij_sigma_k_l_conj, Nij_el - Nji_el);
					std::cout << "T = " << part << std::endl;
					ret += part;
				}
			}
		}

		std::cout << "F_A_T for i = "<< i->get_id() << " and j = " <<  j->get_id() << " = " << ret << std::endl;
		return ret;
	}
	// (10)
	double LCBOPII::N_ijk(Atom *i, Atom *j, Atom *k)
	{
		Atom::position_type r_ij, r_ik;
		r_ij = j->r - i->r;
		r_ik = k->r - i->r;

		return N(i) - S_down_N(r_ij.norm()) - S_down_N(r_ik.norm());
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
			sum += S_down_N(r_ij.norm());
		}

		return sum;
	}
	// (12)
	double LCBOPII::G(Atom *i, Atom *j, Atom *k)
	{
		Atom::position_type r_ij = j->r - i->r;
		Atom::position_type r_ik = k->r - i->r;

		//y = cos_theta_ijk
		double y = r_ij.dot(r_ik)/(r_ij.norm()*r_ik.norm());

		std::cout << "cos(i,j,k)= " << y << std::endl;
		//z = N_ijk
		double z = N_ijk(i,j,k);

		return G(y,z);
	}
	// (12)
	double LCBOPII::G(double y, double z)
	{
		double y_zero = y_0(z);

		std::cout << "y_zero= " << z << std::endl;
		return THETA(y_zero - y)*G_1(y) + THETA(y - y_zero)*G_2(y, z, y_zero);
	}
	// (13)
	double LCBOPII::y_0(double z)
	{
		//std::cout << "y_0= " <<  A_y0 + B_y0*z*(1.0 + z)<< std::endl;

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

		//std::cout << "G_1= " << ret << std::endl;
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
		//std::cout << "#g_z_max " << g_zMax << std::endl;

		double g_z2 = g_z_2(z);
		//double g_z2 = 0.0;
		//std::cout << "#g_z_2 " << g_z2 << std::endl;

		double g_z1 = g_z_1(y0, g_zMax, g_z2);
		//std::cout << "#g_z_1 " << g_z1 << std::endl;

		double g_z0 = g_z_0(y0, g_zMax, g_z1, g_z2);
		//std::cout << "#g_z_0 " << g_z0 << std::endl;

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

		//std::cout << "#G_1_prim " << G_1prim << std::endl;

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
		//std::cout << "# G_1(y0)" << G_1(y0) << std::endl;
		return (G_1(y0) - g_zMax)/std::pow(y0 - 1, 2) -
				g_z1*y0 -
				g_z2*std::pow(y0, 2);
	}

	// (20)
	double LCBOPII::H(Atom *i, Atom *j, Atom *k)
	{
		Atom::position_type r_ij = j->r - i->r;
		Atom::position_type r_ik = k->r - i->r;

		double delta_r_ijk = r_ij.norm() - r_ik.norm();  // (page 3)

		return H(delta_r_ijk);
	}

	double LCBOPII::H(double x)
	{
		if(x < -d)
		{
			return H_1(x);
		}
		else if(x >= -d and x <= d)
		{
			return H_2(x);
		}
		else
		{
			return H_3(x);
		}
	}

	double LCBOPII::H_1(double x)
	{
		return L * (
			1 + kappa*(x + d)*std::pow(1/(1 + std::pow(kappa*(x + d),4)), 0.25)
		);
	}

	double LCBOPII::H_2(double x)
	{
		return 1 + C_1*x + 0.5*std::pow(C_1*x, 2)
					+ C_4*std::pow(x, 4) + C_6*std::pow(x, 6);
	}

	double LCBOPII::H_3(double x)
	{
		return R_0 + R_1*(x - d);
	}

	// (23)
	double LCBOPII::W_ij(Atom *i, Atom *j, uint32_t sigmas)
	{
		double mul = 1.0;
		uint8_t t;
		Atom * k;
		int k_idx = 0;
		uint8_t sigma_k;
		Atom::bond_type bonds_i = i->get_bonds();
		Atom::position_type r_ik;
		double len_ik, S_Nik;

		for(Atom::bond_type::iterator it = bonds_i.begin();
				it != bonds_i.end(); it++)
		{
			k = *it;
			if(k->get_id() == j->get_id()) continue;
			r_ik = k->r - i->r;
			len_ik = r_ik.norm();

			// sigma of k-th neigbour of "i"
			sigma_k = (sigmas >> k_idx) & 1L;
			S_Nik = S_down_N(len_ik);

			if(sigma_k)
			{
				mul *= S_Nik;
			}
			else
			{
				if(S_Nik == 1.0)
				{
					/**
					 * for sigma_k=0 and S_Nik=1 eq. 23 is always 0
					 */
					mul = 0.0;
					break;
				}
				else
				{
					mul *= (1 - S_Nik);
				}
			}
			k_idx++;
		}
		if(k_idx)
			return mul;
		else
			return 0.0;
	}

	double LCBOPII::M_ij_sigma_k(Atom *i, Atom *j, uint32_t sigmas)
	{
		uint8_t t;
		Atom * k;
		int k_idx = 0;
		uint8_t sigma_k;
		Atom::bond_type bonds_i = i->get_bonds();
		Atom::position_type r_ik;
		double len_ik, Nki;

		double ret = 0.0;
		for(Atom::bond_type::iterator it = bonds_i.begin();
				it != bonds_i.end(); it++)
		{
			k = *it;
			if(k->get_id() == j->get_id()) continue;
			r_ik = k->r - i->r;
			len_ik = r_ik.norm();

			// sigma of k-th neigbour of "i"
			sigma_k = (sigmas >> k_idx) & 1L;
			k_idx++;

			if(sigma_k)
			{
				ret += S_up_M(N(k) - S_down_N(len_ik));
				if(ret >= 3.0)
				{
					ret = 3.0;
					break;
				}
			}
			else
			{
				continue;
			}
		}
		return ret;
	}

	// (35)
	double LCBOPII::t_ij(Atom *i, uint32_t sigmas_k,
				     	    Atom *j, uint32_t sigmas_l,
				     	    double z, double delta_el)
	{
		double ret = 0.0;
		double y;
		Atom * nn;
		Atom::position_type r_ij, r_ik1, r_ik2, w_ijk_p, w_ijk_m, t_ijk,
							r_ji, r_jl1, r_jl2, w_jil_p, w_jil_m, t_jil, dummy_v;

		r_ij = j->r - i->r;
		r_ji = i->r - j->r;
		r_ij.normalize();
		r_ji.normalize();

		_t_ij_search_nn(i, j, &r_ik1, &r_ik2, sigmas_k);
		_t_ij_search_nn(j, i, &r_jl1, &r_jl2, sigmas_l);

		w_ijk_m = r_ik1 - r_ik2;
		w_ijk_p = r_ik1 + r_ik2;
		w_ijk_m.normalize();
		w_ijk_p.normalize();

		w_jil_m = r_jl1 - r_jl2;
		w_jil_p = r_jl1 + r_jl2;
		w_jil_m.normalize();
		w_jil_p.normalize();

		/**
		 * eq. 40 for "i"
		 */
		t_ijk = r_ij.cross(w_ijk_m);
		dummy_v = r_ij.cross(w_ijk_p);
		dummy_v *= r_ij.dot(w_ijk_m);
		t_ijk += dummy_v;

		/**
		 * eq. 40 for "j"
		 */
		t_jil = r_ji.cross(w_jil_m);
		dummy_v = r_ji.cross(w_jil_p);
		dummy_v *= r_ji.dot(w_jil_m);
		t_jil += dummy_v;

		t_ijk.normalize();
		t_jil.normalize();

		/**
		 * eq. 39
		 */
		y = t_ijk.dot(t_jil);

		/*
		std::cout << "#y = " << y << std::endl;
		std::cout << "#t_ijk = " << t_ijk.x() << " " << t_ijk.y() << " " << t_ijk.z() << std::endl;
		std::cout << "#t_jil = " << t_jil.x() << " " << t_jil.y()  << " " << t_jil.z() << std::endl;
		std::cout << "#z = " << z << std::endl;
		std::cout << "#delat_el = " << delta_el << std::endl;
		*/

		return t_ij(y, z, std::pow(delta_el, 2));
	}

	double LCBOPII::t_ij(double y, double z, double delta_el)
	{
		double y_q = std::pow(y, 2);
		double delta_el_q = std::pow(delta_el, 2);

		if(z <= 1.0/8.0)
		{
			return tau_1(z)*std::pow(y_q*(1 - y_q), 2);
		}
		else
		{
			return tau_2(z, delta_el_q)*(1 - y_q)*std::pow(2 - y_q, 2);
			//return -0.01837987358053528*(1 - y_q)*std::pow(2 - y_q, 2);
		}
	}

	void LCBOPII::_t_ij_search_nn(Atom *center, Atom *pass,
							  Atom::position_type * r_n1,
							  Atom::position_type * r_n2,
							  uint32_t sigmas_nn)
	{
		Atom * nn;
		uint8_t sigma_nn;
		unsigned int nn_idx = 0;
		unsigned int nn_cnt = 0;
		// searching for "i" neighbours "k1" and "k2"
		for(Atom::bond_type::iterator it = center->get_bonds().begin();
				it != center->get_bonds().end(); it++)
		{
			nn = *it;
			if(nn->get_id() == pass->get_id()) continue;

			// sigma of k-th neigbour of "i"
			sigma_nn = (sigmas_nn >> nn_idx) & 1L;
			nn_idx++;

			if(sigma_nn)
			{
				if(nn_cnt == 0)
				{
					*r_n1 = nn->r;
					*r_n1 -= center->r;
					(*r_n1).normalize();
				}
				else if(nn_cnt == 1)
				{
					*r_n2 = nn->r;
					*r_n2 -= center->r;
					(*r_n2).normalize();
				}
				else
				{
					return;
				}
				nn_cnt++;
			}
		}
	}

	/**
	 * C. Long-range potential
	 */
	double LCBOPII::V_lr(double r)										// (42)
	{
		double ret = S_down_lr(r);

		if(ret == 0.0) return 0.0;

		double r_r0 = r - r_0;

		if(THETA(-r_r0) > 0)
		{
			return ret*V_lr_1(r_r0);
		}
		else
		{
			return ret*V_lr_2(r_r0);
		}
	}

	/**
	 * D. Middle-range potential
	 */
	double LCBOPII::V_mr(Atom * i, Atom * j, double Bij)
	{
		double Ndb_ij = N_db_ij(i, j);
		double x_db = x_db_ij(Ndb_ij);

		if(Ndb_ij <= 1)
		{
			return 0.0;
		}
		else if(Ndb_ij > 1 and Ndb_ij <= 2)
		{

		}
		else if(Ndb_ij > 2 and Ndb_ij <=3)
		{

		}
		else
		{
			return 0.0;
		}

	};

	double LCBOPII::N_db_ij(Atom *i, Atom *j)
	{
		double ret = 4.0;
		Atom * k;
		Eigen::Vector3d r_ki;

		// searching for "i" neighbours "k1" and "k2"
		for(Atom::bond_type::iterator it = i->get_bonds().begin();
				it != i->get_bonds().end(); it++)
		{
			k = *it;
			if(j->get_id() == i->get_id()) continue;

			r_ki = i->r - k->r;

			ret -= S_down_N(r_ki.norm())*N_el_ki(k, i, r_ki.norm());
		}

		return ret;
	}

	double LCBOPII::N_el_ki(Atom * k, Atom * i, double r_ki)
	{
		double Nki = N(k) - S_down_N(r_ki);
		double Mki = M_ki(k, i);
		double SdMki = S_down_sat(Nki)*Mki;

		double ret = (4 - SdMki)/(Nki + 1 - SdMki);

		return ret;
	}

	double LCBOPII::M_ki(Atom * k, Atom * i)
	{
		double ret = 0.0;
		double Nmk;
		Atom * m;
		Eigen::Vector3d r_mk;

		for(Atom::bond_type::iterator it = k->get_bonds().begin();
				it != k->get_bonds().end(); it++)
		{
			m = *it;
			if(k->get_id() == i->get_id()) continue;

			r_mk = m->r - k->r;
			Nmk = N(m) - S_down_N(r_mk.norm());
			ret += S_down_N(r_mk.norm())*S_up_M(Nmk);
		}

		return ret;
	}

	double LCBOPII::gamma_ij(Atom * i, Atom *j, double Bij)
	{
		if(Bij == 0.0) return 1.0;

		double sum = 0.0;
		Atom * k;
		Eigen::Vector3d r_ij, r_ik;

		r_ij = j->r - i->r;

		double Nij = N(i) - S_down_N(r_ij.norm());

		//if(Nij == 0.0) return 1.0;

		for(Atom::bond_type::iterator it = i->get_bonds().begin();
				it != i->get_bonds().end(); it++)
		{
			k = *it;
			if(j->get_id() == i->get_id()) continue;
			r_ik = k->r - i->r;
			sum += S_down_N(r_ik.norm())*
						std::pow(1 + (r_ij.dot(r_ik)/(r_ij.norm()*r_ik.norm())), 4);
		}
		return 1.0/(1 + Bij*sum/Nij);
	}

	double LCBOPII::Z_mr_i(Atom * i, Atom *j)
	{
		/*TODO*/
		return 1.0;
	};
}
