#ifndef __LCBOPII_H
#define __LCBOPII_H
#include <cstdint>
#include <cmath>
#include <Eigen/Dense>
#include "atom.h"

/**
 * DIY implementation of LCBOPII potential
 * 	Phys. Rev. B 72, 214102
 * 	"Improved long-range reactive bond-order potential for carbon. I. Construction"
 * 	by. Jan. H. Los. et al.
 */
namespace simul
{
	class LCBOPII
	{
	public:
		// TABLE I.
		const static double A_sr ;
		const static double alpha;

		const static double B_sr_1;
		const static double beta_1;
		const static double B_sr_2;
		const static double beta_2;

		const static double A_y0; 
		const static double B_y0; 

		const static double one_third;

		// function G(x,y)
		const static double g_min;
		const static double g_gr;
		const static double g_max;

		const static double g_1_n[], g_2_n[], g_3_n[];

		const static double A_g;
		const static double B_g;
		const static double C_g;

		const static double D_g;
		const static double E_g;

		// function H(x)
		const static double d;
		const static double C_1;
		const static double C_4;

		const static double L;
		const static double R_0;
		const static double kappa;
		const static double R_1;
		const static double C_6;

		// function F(N_ij, N_ji, N_ij_conj)
		const static double F_ij_0[][4];
		const static double F_ij_1[][4];

		// function A_ij
		const static double alpha_0;

		// function T_ij
		const static double A_t;// = -13.152909887;
		const static double B_t1;
		const static double B_t2;
		const static double B_t3;
		const static double B_t4;

	private:

		template<typename T> static double THETA(T x)		// Heavyside step function
		{
			if(x)
				return (x >= 0 ? 1.0:0.0);
			else
				return 0.5;
		}

		double S_down(double x, double p);						// (2)
		double S_up(double x, double p);							// (3)
		double x_from_q(double q, double q_min, double q_max);	// (4)

		void _t_ij_search_nn(Atom *center, Atom *pass,
							  Atom::position_type * r_n1,
							  Atom::position_type * r_n2,
							  uint32_t sigmas_nn);

	public:
		/**
		 * Switch functions from TABLE I. page 4
		 */
		double S_down_sr(double q);
		double S_down_lr(double q);
		double S_down_db(double q);

		double S_up_mr(double q);
		double S_up_M(double q);
		double S_up_gamma_0(double q);

		double S_down_N(double q);
		double S_down_sat(double q);
		double S_up_gamma_2(double q);

		/**
		 * B. Short range potential
		 */
		double V_sr(Atom *i, Atom *j);							// (5)
		double V_sr_R(double r);								// (6)
		double V_sr_A(double r);								// (7)
		double B(Atom *i, Atom *j);								// (8)

		/**
		 * function returns F_ij_conj + A_ij + T_ij
		 */
		double F_A_T(Atom *i, Atom *j,
						bool enable_F = true,
						bool enable_A = true,
						bool enable_T = true);

		/**
		 * Term b_ij (p. 3)
		 */
		double b(Atom *i, Atom *j);								// (9)
		double N_ijk(Atom *i, Atom *j, Atom *k);				// (10)
		double N(Atom *i);										// (11)
		double G(Atom *i, Atom *j, Atom *k);					// (12)
		double G(double y, double z);							// (12 y,z)
		double y_0(double z);									// (13)
		double G_1(double y);									// (14)
		double G_1_prim(double y0);

		double G_2(double y, double z, double y0);				// (15)
		double G_2_prim(double y, double z, double y0);			// (15)

		double g_z_max(double z, double y0);						// (16)
		double g_z_2(double z);									// (17)
		double g_z_1(double y0, double g_zMax, double g_z2);		// (18)
		double g_z_0(double y0, double g_zMax,
					double g_z1, double g_z2);						// (19)

		double H(Atom *i, Atom *j, Atom *k);						// (20)
		double H(double x);										// (20)
		double H_1(double x);
		static double H_2(double x);
		double H_3(double x);

		double F_conj(unsigned int Nij_sigma,						// (22)
					    unsigned int Nji_sigma,
					    double Nij_sigma_k_l_conj);

		double W_ij(Atom *i, Atom *j, uint32_t sigmas);		  		// (23)
		unsigned int N_ij_sigma_k(uint32_t sigmas);			  		// (24), (25)
		double N_ij_el(unsigned int Nij_sigma, double Mij_sigma);  	// (27)
		double M_ij_sigma_k(Atom *i, Atom *j, uint32_t sigmas); 		// (28)
		double N_ij_min_el(unsigned int Nij_sigma);                 // (30)
		double N_ij_max_el(unsigned int Nij_sigma);                 // (30)

		double A(double delta_el);									    // (32)
		double t_ij(Atom *i, uint32_t sigmas_k,						// (35)
				     Atom *j, uint32_t sigmas_l,
				     double z, double delta_el);
	};


	inline double LCBOPII::S_down(double x, double p)
	{
		return THETA(-x) + THETA(x)*THETA(1-x)*(1 + 2*x + p*x*x)*(1-x)*(1-x);
	}

	inline double LCBOPII::S_up(double x, double p)
	{
		return 1.0 - S_down(x, p);
	}

	inline double LCBOPII::x_from_q(double q, double q_min, double q_max)
	{
		return (q - q_min)/(q_max - q_min);
	}

	inline unsigned int LCBOPII::N_ij_sigma_k(uint32_t sigmas)
	{
		return (sigmas & 1L) + ((sigmas>>1) & 1L) + ((sigmas>>2) & 1L);
	}

	inline double LCBOPII::N_ij_el(unsigned int  Nij_sigma, double Mij_sigma)
	{
		return (4.0 - Mij_sigma)/(Nij_sigma + 1.0 - Mij_sigma);
	}

	inline double LCBOPII::N_ij_min_el(unsigned int Nij_sigma)
	{
		return 4.0/(Nij_sigma + 1.0);
	}

	inline double LCBOPII::N_ij_max_el(unsigned int Nij_sigma)
	{
		return 4.0 - Nij_sigma;
	}

	inline 	double LCBOPII::F_conj(unsigned int Nij_sigma, unsigned int Nji_sigma,
											double Nij_sigma_k_l_conj)
	{
		return (1.0 - Nij_sigma_k_l_conj)*F_ij_0[Nij_sigma][Nji_sigma]
		          + Nij_sigma_k_l_conj*F_ij_1[Nij_sigma][Nji_sigma] ;
	}

	inline double LCBOPII::A(double delta_el)
	{
		return (alpha_0 * std::pow(delta_el, 2))/(1 + 10*std::abs(delta_el));
	}
}
#endif
