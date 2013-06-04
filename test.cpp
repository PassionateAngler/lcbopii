#include <iostream>
#include <Eigen/Dense>
#include <deque>
#include <cmath>
#include "atom.h"
#include "lcbopii.h"
#include "printers.h"

using namespace  simul;

void test_G(LCBOPII lcbopii)
{
	double y = -1.0;
	double step =2/50.0; // 1 - (-1)
	double z = 0.0;
	double y0;

	while(z <= 8)
	{
		y = -1.0;
		std::cout << "# z=" << z << std::endl;

		while(y <= 1.0)
		{
			std::cout << y << " " << lcbopii.G(y,z) << std::endl;
			//std::cout << y << " " << lcbopii.G_1(y) << std::endl;
			//std::cout << std::pow((1 - y), 2) << std::endl;
			y += step;
		}
		y0 = lcbopii.y_0(z);
		//std::cout << "G_1(y0) = " << lcbopii.G_1(y0) << " G_2(y0) = " << lcbopii.G_2(y0, z, y0) << std::endl;
		//std::cout << "G_1_prim(y0) = " << lcbopii.G_1_prim(y0)
		//	     << " G_2_prim(y0) = " << lcbopii.G_2_prim(y0, z, y0) << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		z += 1.0;
	}
}

void test_H(LCBOPII lcbopii)
{
	double x = -1.25;
	double step = (2.5)/400;

	while(x < 1.25)
	{

		//std::cout << "# H_2(d) =" << lcbopii.H_2(lcbopii.d) << std::endl;
		//std::cout << "# H_2_prim(d) =" << lcbopii.H_2_prim(lcbopii.d) << std::endl;
		//std::cout << "# H_3(d) =" << lcbopii.H_3(lcbopii.d) << std::endl;
		std::cout << x << " " << lcbopii.H(x) << std::endl;
		x += step;
	}
}

/*class Test
{
	virtual void dupa();
};
*/
/*
class TestAnalizer: public Atom::StructAnalizer
{
	virtual void bound_found(Atom *u, Atom *v)
	{
		std::cout << "Atom: " << u->get_id() << " --- " << v->get_id() << std::endl;
	}
	virtual void atom_found(Atom *u)
	{
		std::cout << "Odwiedzilem Atom: " << u->get_id() << std::endl;
	}

};
*/

int main()
{
	double r = 1.42;
	double omega = 2.0*M_PI/3.0;

	Atom * i = new Atom(0.0, 0, 0.0);
	Atom * k1 = new Atom(r*std::cos(omega), r*std::sin(omega), 0.0);
	Atom * k2 = new Atom(r*std::cos(2*omega), r*std::sin(2*omega), 0.0);


	Atom * j = new Atom(r, 0.0, 0.0);
	Atom * l1 = new Atom(r*std::cos(omega + M_PI) + r, r*std::sin(omega + M_PI), 0.0);
	Atom * l2 = new Atom(r*std::cos(2*omega + M_PI) + r, r*std::sin(2*omega + M_PI), 0.0);

	i->addBond(j);
	i->addBond(k1);
	i->addBond(k2);

	j->addBond(l1);
	j->addBond(l2);

	Eigen::Vector3d ri = i->r - j->r;

//	for(double theta = 0.0; theta<= M_PI/2.0; theta += M_PI/12.0)

	const double DELTA_THETA = M_PI/180.0;

	Eigen::Vector3d k1_org, k2_org;

	k1_org = k1->r;
	k2_org = k2->r;
	LCBOPII lcbopii;// = new LCBOPII();
	for(double theta = 0.0; theta<= M_PI/2.0; theta += DELTA_THETA )
	{

	k1->r = k1_org;
	k2->r = k2_org;
	Eigen::AngleAxisd rot_ri(theta, ri.normalized());
	Eigen::Vector3d r_k1_rot = rot_ri * (k1->r - i->r);
	k1->r += (r_k1_rot - k1->r);

	Eigen::Vector3d r_k2_rot = rot_ri * (k2->r - i->r);
	k2->r += (r_k2_rot - k2->r);

	/*
	SdfPrinter * aw = new SdfPrinter();
	std::cout << aw->get_sdf(k1);
	delete aw;
	*/

	//std::cout << std::endl;
	//aw->walk_BFS(k2);

	//std::cout << std::endl;
	//Atom::walk_BFS(i, new TestAnalizer());

	//test_H(lcbopii);
	/*
	std::cout << "#Fij= "<< lcbopii.F_A_T(i, j, true, false, false) << std::endl;
	std::cout << "#Aij= "<< lcbopii.F_A_T(i, j, false, true, false) << std::endl;
	std::cout << "#Tij= "<< lcbopii.F_A_T(i, j, false, false, true) << std::endl;
	*/

	//std::cout << theta/M_PI << " " << lcbopii.F_A_T(i, j)*lcbopii.V_sr(i,j) << std::endl;
	//std::cout << std::endl;
	/*
   std::cout << "H_2(d)" << lcbopii.H_2(lcbopii.d) << std::endl;
   std::cout << "H_3(d)" << lcbopii.H_3(lcbopii.d) << std::endl;*/
	/*
   std::cout << "#V_sr " << lcbopii.V_sr(i, j) << std::endl;
   std::cout << "#V_sr " << lcbopii.V_sr(i, k1) << std::endl;
   std::cout << "#V_sr " << lcbopii.V_sr(i, k2) << std::endl;
   std::cout << "#V_sr " << lcbopii.V_sr(j, i) << std::endl;
   std::cout << "#V_sr " << lcbopii.V_sr(j, l1) << std::endl;
   std::cout << "#V_sr " << lcbopii.V_sr(j, l2) << std::endl;

   */
   //double sum = lcbopii.V_sr(i, j) + lcbopii.V_sr(i, k1) + lcbopii.V_sr(i, k2)
//		   + lcbopii.V_sr(j, i) + lcbopii.V_sr(j, l1) + lcbopii.V_sr(j, l2);

   //std::cout << theta/M_PI << " " << sum/6 << std::endl;

   //std::cout << "V_sr " << lcbopii.V_sr(j, i) << std::endl;
	}

	//std::cout << lcbopii.tau_2(1.0, 0);
	return 1;
}
