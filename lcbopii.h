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

	private:

		template<typename T> static float THETA(T x)		// Heavyside step function
		{
			if(x)
				return (x >= 0 ? 1.0:0.0);
			else
				return 0.5;
		}

		float S_down(float x, float p);						// (2)
		float S_up(float x, float p);							// (3)
		float x_from_q(float q, float q_min, float q_max);	// (4)

	public:
		/**
		 * Switch functions from TABLE I. page 4
		 */
		float S_down_sr(float q);
		float S_down_lr(float q);
		float S_down_db(float q);

		float S_up_mr(float q);
		float S_up_M(float q);
		float S_up_gamma_0(float q);

		float S_down_N(float q);
		float S_down_sat(float q);
		float S_up_gamma_2(float q);

		/**
		 * B. Short range potential
		 */
		float V_sr(Atom *i, Atom *j);							// (5)
		float V_sr_R(float r);								// (6)
		float V_sr_A(float r);								// (7)
		float B(Atom *i, Atom *j);								// (8)
		/**
		 * Term b_ij (p. 3)
		 */
		float b(Atom *i, Atom *j);								// (9)
		float N_ijk(Atom *i, Atom *j, Atom *k);				// (10)
		float N(Atom *i);										// (11)
		float H(Atom *i, Atom *j, Atom *k);						// (20)

		float F_conj(Atom *i, Atom *j);						// (22)
		float A(Atom *i, Atom *j);								// (32)
		float T(Atom *i, Atom *j);								// (35)
	};



	inline float LCBOPII::S_down(float x, float p)
	{
		return THETA(-x) + THETA(x)*THETA(1-x)*(1 + 2*x + p*x*x)*(1-x)*(1-x);
	}

	inline float LCBOPII::S_up(float x, float p)
	{
		return 1.0 - S_down(x, p);
	}

	inline float LCBOPII::x_from_q(float q, float q_min, float q_max)
	{
		return (q - q_min)/(q_max - q_min);
	}
}
#endif
