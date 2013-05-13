#ifndef __LCBOPII_H
#define __LCBOPII_H
#include <gmtl/gmtl.h>
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

		const static double g_min;
		const static double g_gr;
		const static double g_max;

		const static double g_1_n[], g_2_n[], g_3_n[];

		const static double A_g;
		const static double B_g;
		const static double C_g;

		const static double D_g;
		const static double E_g;


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
		 * Term b_ij (p. 3)
		 */
		double b(Atom *i, Atom *j);								// (9)
		double N_ijk(Atom *i, Atom *j, Atom *k);				// (10)
		double N(Atom *i);										// (11)
		double G(Atom *i, Atom *j, Atom *k);						// (12)
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

		double F_conj(Atom *i, Atom *j);						// (22)
		double A(Atom *i, Atom *j);								// (32)
		double T(Atom *i, Atom *j);								// (35)
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
}
#endif
